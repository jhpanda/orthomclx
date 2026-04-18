from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Union

from orthomcl.compat import dataclass
from orthomcl.evalue import parse_evalue_cutoff

from orthomcl.io import format_score, round_score, write_mcl_input
from orthomcl.models import EdgeRecord


@dataclass(slots=True)
class IndexedPairsSummary:
    ortholog_count: int
    inparalog_count: int
    coortholog_count: int
    rbh_count: int
    mcl_count: int
    out_dir: Path


@dataclass(slots=True)
class RawOrtholog:
    seq_a: str
    seq_b: str
    taxon_a: str
    taxon_b: str
    unnormalized_score: float


@dataclass(slots=True)
class RawInparalog:
    seq_a: str
    seq_b: str
    taxon_id: str
    unnormalized_score: float


def ensure_binary(root: Path, binary_name: str, make_target: str) -> Path:
    binary = root / "build" / binary_name
    source = root / "src" / "c" / f"{binary_name}.c"
    if binary.exists() and source.exists() and binary.stat().st_mtime >= source.stat().st_mtime:
        return binary
    subprocess.run(["make", make_target], cwd=root, check=True)
    return binary


def run_indexed_pairs(
    compiled_dir: Union[str, Path],
    out_dir: Union[str, Path],
    percent_match_cutoff: float,
    evalue_cutoff: float,
    jobs: int = 1,
) -> IndexedPairsSummary:
    compiled_dir = Path(compiled_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    root = Path(__file__).resolve().parents[2]
    orth_bin = ensure_binary(root, "indexed_orthologs", "build-indexed-orthologs")
    inpara_bin = ensure_binary(root, "indexed_inparalogs", "build-indexed-inparalogs")
    co_bin = ensure_binary(root, "indexed_coorthologs", "build-indexed-coorthologs")
    rbh_bin = ensure_binary(root, "indexed_rbh", "build-indexed-rbh")

    stage_dir = out_dir / "indexed_stages"
    stage_dir.mkdir(parents=True, exist_ok=True)
    cutoff_mant, cutoff_exp = parse_evalue_cutoff(evalue_cutoff)

    print(f"[orthomclx] building orthologs (jobs={jobs}, pcut={percent_match_cutoff}, ecut={evalue_cutoff})")

    orthologs_raw = run_parallel_stage(
        binary=orth_bin,
        base_args=[
            str(compiled_dir / "similarities.bin"),
            str(compiled_dir / "proteins.tsv"),
            str(compiled_dir / "taxa.tsv"),
        ],
        extra_before_out=[],
        out_base_dir=stage_dir / "orthologs",
        output_name="orthologs.indexed.raw.tsv",
        merged_output=stage_dir / "orthologs.merged.raw.tsv",
        percent_match_cutoff=percent_match_cutoff,
        cutoff_mant=cutoff_mant,
        cutoff_exp=cutoff_exp,
        jobs=jobs,
    )
    ortholog_edges = normalize_raw_orthologs(read_raw_orthologs(orthologs_raw), "ortholog")
    orthologs_file = out_dir / "pairs" / "orthologs.txt"
    write_edge_records(orthologs_file, ortholog_edges)
    print(f"[orthomclx] orthologs complete: {len(ortholog_edges)} edges")

    print(f"[orthomclx] building inparalogs (jobs={jobs})")
    inparalogs_raw = run_parallel_stage(
        binary=inpara_bin,
        base_args=[
            str(compiled_dir / "similarities.bin"),
            str(compiled_dir / "proteins.tsv"),
            str(compiled_dir / "taxa.tsv"),
        ],
        extra_before_out=[str(orthologs_file)],
        out_base_dir=stage_dir / "inparalogs",
        output_name="inparalogs.indexed.raw.tsv",
        merged_output=stage_dir / "inparalogs.merged.raw.tsv",
        percent_match_cutoff=percent_match_cutoff,
        cutoff_mant=cutoff_mant,
        cutoff_exp=cutoff_exp,
        jobs=jobs,
    )
    inparalog_edges = normalize_raw_inparalogs(
        read_raw_inparalogs(inparalogs_raw),
        {edge.seq_a for edge in ortholog_edges} | {edge.seq_b for edge in ortholog_edges},
    )
    inparalogs_file = out_dir / "pairs" / "inparalogs.txt"
    write_edge_records(inparalogs_file, inparalog_edges)
    print(f"[orthomclx] inparalogs complete: {len(inparalog_edges)} edges")

    print(f"[orthomclx] building coorthologs (jobs={jobs})")
    coorthologs_raw = run_parallel_stage(
        binary=co_bin,
        base_args=[
            str(compiled_dir / "similarities.bin"),
            str(compiled_dir / "proteins.tsv"),
            str(compiled_dir / "taxa.tsv"),
        ],
        extra_before_out=[str(orthologs_file), str(inparalogs_file)],
        out_base_dir=stage_dir / "coorthologs",
        output_name="coorthologs.indexed.raw.tsv",
        merged_output=stage_dir / "coorthologs.merged.raw.tsv",
        percent_match_cutoff=percent_match_cutoff,
        cutoff_mant=cutoff_mant,
        cutoff_exp=cutoff_exp,
        jobs=jobs,
    )
    coortholog_edges = normalize_raw_orthologs(read_raw_orthologs(coorthologs_raw), "coortholog")
    coorthologs_file = out_dir / "pairs" / "coorthologs.txt"
    write_edge_records(coorthologs_file, coortholog_edges)
    print(f"[orthomclx] coorthologs complete: {len(coortholog_edges)} edges")

    rbh_dir = stage_dir / "rbh"
    rbh_dir.mkdir(parents=True, exist_ok=True)
    print(f"[orthomclx] building strict 1:1 RBH (jobs={jobs})")
    rbh_merged = run_parallel_rbh_stage(
        binary=rbh_bin,
        base_args=[
            str(compiled_dir / "similarities.bin"),
            str(compiled_dir / "proteins.tsv"),
            str(compiled_dir / "taxa.tsv"),
        ],
        out_base_dir=rbh_dir,
        merged_output=stage_dir / "rbh_1to1.merged.txt",
        percent_match_cutoff=percent_match_cutoff,
        cutoff_mant=cutoff_mant,
        cutoff_exp=cutoff_exp,
        jobs=jobs,
    )
    rbh_file = out_dir / "pairs" / "rbh_1to1.txt"
    rbh_file.parent.mkdir(parents=True, exist_ok=True)
    rbh_file.write_text(rbh_merged.read_text())
    with rbh_file.open() as handle:
        rbh_count = sum(1 for _ in handle)
    print(f"[orthomclx] rbh_1to1 complete: {rbh_count} pairs")

    all_edges = sorted(
        ortholog_edges + inparalog_edges + coortholog_edges,
        key=lambda edge: (edge.seq_a, edge.seq_b, edge.edge_type),
    )
    write_mcl_input(out_dir / "mclInput", all_edges)
    print(f"[orthomclx] mclInput written: {len(all_edges)} weighted edges")
    return IndexedPairsSummary(
        ortholog_count=len(ortholog_edges),
        inparalog_count=len(inparalog_edges),
        coortholog_count=len(coortholog_edges),
        rbh_count=rbh_count,
        mcl_count=len(all_edges),
        out_dir=out_dir,
    )


def run_parallel_rbh_stage(
    *,
    binary: Path,
    base_args: List[str],
    out_base_dir: Path,
    merged_output: Path,
    percent_match_cutoff: float,
    cutoff_mant: float,
    cutoff_exp: int,
    jobs: int,
) -> Path:
    jobs = max(1, jobs)
    out_base_dir.mkdir(parents=True, exist_ok=True)
    shard_dirs = [out_base_dir / f"shard_{index}" for index in range(jobs)]
    for shard_dir in shard_dirs:
        shard_dir.mkdir(parents=True, exist_ok=True)

    def run_one(shard_index: int) -> Path:
        shard_dir = shard_dirs[shard_index]
        cmd = (
            [str(binary)]
            + base_args
            + [str(shard_dir), str(percent_match_cutoff), str(cutoff_mant), str(cutoff_exp), str(shard_index), str(jobs)]
        )
        subprocess.run(cmd, check=True)
        return shard_dir / "rbh_1to1.txt"

    if jobs == 1:
        shard_files = [run_one(0)]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            shard_files = list(executor.map(run_one, range(jobs)))

    merge_plain_text_files(shard_files, merged_output)
    return merged_output


def run_parallel_stage(
    *,
    binary: Path,
    base_args: List[str],
    extra_before_out: List[str],
    out_base_dir: Path,
    output_name: str,
    merged_output: Path,
    percent_match_cutoff: float,
    cutoff_mant: float,
    cutoff_exp: int,
    jobs: int,
) -> Path:
    jobs = max(1, jobs)
    out_base_dir.mkdir(parents=True, exist_ok=True)
    shard_dirs = [out_base_dir / f"shard_{index}" for index in range(jobs)]
    for shard_dir in shard_dirs:
        shard_dir.mkdir(parents=True, exist_ok=True)

    def run_one(shard_index: int) -> Path:
        shard_dir = shard_dirs[shard_index]
        cmd = (
            [str(binary)]
            + base_args
            + extra_before_out
            + [str(shard_dir), str(percent_match_cutoff), str(cutoff_mant), str(cutoff_exp), str(shard_index), str(jobs)]
        )
        subprocess.run(cmd, check=True)
        return shard_dir / output_name

    if jobs == 1:
        shard_files = [run_one(0)]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            shard_files = list(executor.map(run_one, range(jobs)))

    merged_output.parent.mkdir(parents=True, exist_ok=True)
    merge_edge_files(shard_files, merged_output)
    return merged_output


def merge_edge_files(shard_files: List[Path], output_path: Path) -> None:
    merged: dict[tuple[str, ...], str] = {}
    for shard_file in shard_files:
        if not shard_file.exists():
            continue
        with shard_file.open() as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                key = tuple(fields[:-1])
                value = fields[-1]
                existing = merged.get(key)
                if existing is not None and existing != value:
                    raise ValueError(f"Conflicting scores for {key}: {existing} vs {value}")
                merged[key] = value
    with output_path.open("w") as handle:
        for key in sorted(merged):
            handle.write("\t".join(key) + f"\t{merged[key]}\n")


def merge_plain_text_files(shard_files: List[Path], output_path: Path) -> None:
    merged: set[str] = set()
    for shard_file in shard_files:
        if not shard_file.exists():
            continue
        with shard_file.open() as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line:
                    merged.add(line)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for line in sorted(merged):
            handle.write(f"{line}\n")


def read_edge_records(path: Path, edge_type: str) -> list[EdgeRecord]:
    edges: list[EdgeRecord] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            seq_a, seq_b, score = line.split("\t")
            edges.append(
                EdgeRecord(
                    seq_a=seq_a,
                    seq_b=seq_b,
                    taxon_a="",
                    taxon_b="",
                    unnormalized_score=float(score),
                    normalized_score=float(score),
                    edge_type=edge_type,
                )
            )
    return edges


def read_raw_orthologs(path: Path) -> list[RawOrtholog]:
    rows: list[RawOrtholog] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            seq_a, seq_b, taxon_a, taxon_b, score = line.split("\t")
            rows.append(
                RawOrtholog(
                    seq_a=seq_a,
                    seq_b=seq_b,
                    taxon_a=taxon_a,
                    taxon_b=taxon_b,
                    unnormalized_score=float(score),
                )
            )
    return rows


def read_raw_inparalogs(path: Path) -> list[RawInparalog]:
    rows: list[RawInparalog] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            seq_a, seq_b, taxon_id, score = line.split("\t")
            rows.append(
                RawInparalog(
                    seq_a=seq_a,
                    seq_b=seq_b,
                    taxon_id=taxon_id,
                    unnormalized_score=float(score),
                )
            )
    return rows


def normalize_raw_orthologs(rows: list[RawOrtholog], edge_type: str) -> list[EdgeRecord]:
    sums: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        key = tuple(sorted((row.taxon_a, row.taxon_b)))
        sums.setdefault(key, []).append(row.unnormalized_score)
    avgs = {key: sum(values) / len(values) for key, values in sums.items()}
    edges = [
        EdgeRecord(
            seq_a=row.seq_a,
            seq_b=row.seq_b,
            taxon_a=row.taxon_a,
            taxon_b=row.taxon_b,
            unnormalized_score=row.unnormalized_score,
            normalized_score=row.unnormalized_score / avgs[tuple(sorted((row.taxon_a, row.taxon_b)))],
            edge_type=edge_type,
        )
        for row in rows
    ]
    return sorted(edges, key=lambda edge: (edge.taxon_a, edge.taxon_b, edge.seq_a, edge.seq_b))


def normalize_raw_inparalogs(rows: list[RawInparalog], ortholog_ids: set[str]) -> list[EdgeRecord]:
    all_scores: dict[str, list[float]] = {}
    orth_scores: dict[str, list[float]] = {}
    for row in rows:
        all_scores.setdefault(row.taxon_id, []).append(row.unnormalized_score)
        if row.seq_a in ortholog_ids or row.seq_b in ortholog_ids:
            orth_scores.setdefault(row.taxon_id, []).append(row.unnormalized_score)
    avgs = {
        taxon: (sum(orth_scores[taxon]) / len(orth_scores[taxon])) if taxon in orth_scores else (sum(scores) / len(scores))
        for taxon, scores in all_scores.items()
    }
    edges = [
        EdgeRecord(
            seq_a=row.seq_a,
            seq_b=row.seq_b,
            taxon_a=row.taxon_id,
            taxon_b=row.taxon_id,
            unnormalized_score=row.unnormalized_score,
            normalized_score=row.unnormalized_score / avgs[row.taxon_id],
            edge_type="inparalog",
        )
        for row in rows
    ]
    return sorted(edges, key=lambda edge: (edge.taxon_a, edge.seq_a, edge.seq_b))


def write_edge_records(path: Path, edges: list[EdgeRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for edge in edges:
            handle.write(f"{edge.seq_a}\t{edge.seq_b}\t{format_score(edge.normalized_score)}\n")
