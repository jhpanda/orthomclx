from __future__ import annotations

import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Union

from orthomcl.io import (
    ensure_directory,
    read_edge_file,
    read_similarity_records,
    write_edge_file,
    write_mcl_input,
)
from orthomcl.models import EdgeRecord, SimilarityRecord


def evalue_key(record: SimilarityRecord) -> tuple[int, float]:
    return record.evalue_exp, record.evalue_mant


def score_from_records(left: SimilarityRecord, right: SimilarityRecord, zero_cutoff: float) -> float:
    if left.evalue_mant < zero_cutoff or right.evalue_mant < zero_cutoff:
        return (left.evalue_exp + right.evalue_exp) / -2.0
    return (math.log10(left.evalue_mant * right.evalue_mant) + left.evalue_exp + right.evalue_exp) / -2.0


def canonical_pair(seq_a: str, seq_b: str) -> tuple[str, str]:
    return (seq_a, seq_b) if seq_a < seq_b else (seq_b, seq_a)


def passes_thresholds(
    record: SimilarityRecord,
    percent_match_cutoff: float,
    evalue_exp_cutoff: int,
) -> bool:
    return record.percent_match >= percent_match_cutoff and record.evalue_exp <= evalue_exp_cutoff


def average(values: list[float]) -> float:
    if not values:
        raise ValueError("Cannot average an empty list")
    return sum(values) / len(values)


def build_pairs(
    similar_sequences_path: Union[str, Path],
    out_dir: Union[str, Path],
    percent_match_cutoff: float,
    evalue_exp_cutoff: int,
    engine: str = "python",
) -> dict[str, list[EdgeRecord]]:
    if engine == "auto":
        engine = "c" if shutil.which("cc") else "python"
    if engine == "c":
        return build_pairs_c(
            similar_sequences_path,
            out_dir,
            percent_match_cutoff,
            evalue_exp_cutoff,
        )
    if engine != "python":
        raise ValueError(f"Unsupported engine: {engine}")

    return build_pairs_python(
        similar_sequences_path,
        out_dir,
        percent_match_cutoff,
        evalue_exp_cutoff,
    )


def build_pairs_python(
    similar_sequences_path: Union[str, Path],
    out_dir: Union[str, Path],
    percent_match_cutoff: float,
    evalue_exp_cutoff: int,
) -> dict[str, list[EdgeRecord]]:
    similarities = read_similarity_records(similar_sequences_path)
    by_pair = {(record.query_id, record.subject_id): record for record in similarities}

    orthologs = build_orthologs(similarities, by_pair, percent_match_cutoff, evalue_exp_cutoff)
    inparalogs = build_inparalogs(
        similarities,
        by_pair,
        orthologs,
        percent_match_cutoff,
        evalue_exp_cutoff,
    )
    coorthologs = build_coorthologs(
        similarities,
        by_pair,
        orthologs,
        inparalogs,
        percent_match_cutoff,
        evalue_exp_cutoff,
    )

    output_dir = ensure_directory(out_dir)
    pairs_dir = ensure_directory(output_dir / "pairs")
    write_edge_file(pairs_dir / "orthologs.txt", orthologs)
    write_edge_file(pairs_dir / "inparalogs.txt", inparalogs)
    write_edge_file(pairs_dir / "coorthologs.txt", coorthologs)
    all_edges = sorted(
        orthologs + inparalogs + coorthologs,
        key=lambda edge: (edge.seq_a, edge.seq_b, edge.edge_type),
    )
    write_mcl_input(output_dir / "mclInput", all_edges)
    return {
        "orthologs": orthologs,
        "inparalogs": inparalogs,
        "coorthologs": coorthologs,
    }


def ensure_c_engine_built() -> Path:
    root = Path(__file__).resolve().parents[2]
    binary = root / "build" / "pairs_engine"
    source = root / "src" / "c" / "pairs_engine.c"
    if binary.exists() and binary.stat().st_mtime >= source.stat().st_mtime:
        return binary
    subprocess.run(["make", "build-c-engine"], cwd=root, check=True)
    return binary


def build_pairs_c(
    similar_sequences_path: Union[str, Path],
    out_dir: Union[str, Path],
    percent_match_cutoff: float,
    evalue_exp_cutoff: int,
) -> dict[str, list[EdgeRecord]]:
    binary = ensure_c_engine_built()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(binary),
            str(similar_sequences_path),
            str(out_dir),
            str(percent_match_cutoff),
            str(evalue_exp_cutoff),
        ],
        check=True,
    )
    return {
        "orthologs": read_edge_file(out_dir / "pairs" / "orthologs.txt", "ortholog"),
        "inparalogs": read_edge_file(out_dir / "pairs" / "inparalogs.txt", "inparalog", same_taxon=True),
        "coorthologs": read_edge_file(out_dir / "pairs" / "coorthologs.txt", "coortholog"),
    }


def build_orthologs(
    similarities: list[SimilarityRecord],
    by_pair: dict[tuple[str, str], SimilarityRecord],
    percent_match_cutoff: float,
    evalue_exp_cutoff: int,
) -> list[EdgeRecord]:
    best_by_query_taxon: dict[tuple[str, str], tuple[int, float]] = {}
    for record in similarities:
        if record.query_taxon == record.subject_taxon:
            continue
        key = (record.query_id, record.subject_taxon)
        record_key = evalue_key(record)
        current = best_by_query_taxon.get(key)
        if current is None or record_key < current:
            best_by_query_taxon[key] = record_key

    best_hits: dict[tuple[str, str], SimilarityRecord] = {}
    for record in similarities:
        if record.query_taxon == record.subject_taxon:
            continue
        if not passes_thresholds(record, percent_match_cutoff, evalue_exp_cutoff):
            continue
        cutoff = best_by_query_taxon.get((record.query_id, record.subject_taxon))
        if cutoff is None:
            continue
        if record.evalue_mant < 0.01 or evalue_key(record) == cutoff:
            best_hits[(record.query_id, record.subject_id)] = record

    temp_edges: list[tuple[str, str, str, str, float]] = []
    for (query_id, subject_id), left in best_hits.items():
        if query_id >= subject_id:
            continue
        right = best_hits.get((subject_id, query_id))
        if right is None:
            continue
        temp_edges.append(
            (
                query_id,
                subject_id,
                left.query_taxon,
                left.subject_taxon,
                score_from_records(left, right, zero_cutoff=0.01),
            )
        )

    avg_by_taxon_pair: dict[tuple[str, str], float] = {}
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for _, _, taxon_a, taxon_b, score in temp_edges:
        grouped[tuple(sorted((taxon_a, taxon_b)))].append(score)
    for key, scores in grouped.items():
        avg_by_taxon_pair[key] = average(scores)

    return sorted(
        [
            EdgeRecord(
                seq_a=seq_a,
                seq_b=seq_b,
                taxon_a=taxon_a,
                taxon_b=taxon_b,
                unnormalized_score=score,
                normalized_score=score / avg_by_taxon_pair[tuple(sorted((taxon_a, taxon_b)))],
                edge_type="ortholog",
            )
            for seq_a, seq_b, taxon_a, taxon_b, score in temp_edges
        ],
        key=lambda edge: (edge.taxon_a, edge.taxon_b, edge.seq_a, edge.seq_b),
    )


def build_inparalogs(
    similarities: list[SimilarityRecord],
    by_pair: dict[tuple[str, str], SimilarityRecord],
    orthologs: list[EdgeRecord],
    percent_match_cutoff: float,
    evalue_exp_cutoff: int,
) -> list[EdgeRecord]:
    best_inter_taxon: dict[str, tuple[int, float]] = {}
    for record in similarities:
        if record.query_taxon == record.subject_taxon:
            continue
        key = record.query_id
        record_key = evalue_key(record)
        current = best_inter_taxon.get(key)
        if current is None or record_key < current:
            best_inter_taxon[key] = record_key

    intra_candidates: dict[tuple[str, str], SimilarityRecord] = {}
    queries_with_any_similarity = {record.query_id for record in similarities}
    for record in similarities:
        if record.query_taxon != record.subject_taxon:
            continue
        if record.query_id == record.subject_id:
            continue
        if not passes_thresholds(record, percent_match_cutoff, evalue_exp_cutoff):
            continue
        threshold = best_inter_taxon.get(record.query_id)
        if threshold is None:
            if record.query_id in queries_with_any_similarity:
                intra_candidates[(record.query_id, record.subject_id)] = record
            continue
        if (
            record.evalue_mant < 0.001
            or record.evalue_exp < threshold[0]
            or (record.evalue_exp == threshold[0] and record.evalue_mant <= threshold[1])
        ):
            intra_candidates[(record.query_id, record.subject_id)] = record

    temp_edges: list[tuple[str, str, str, float]] = []
    for (query_id, subject_id), left in intra_candidates.items():
        if query_id >= subject_id:
            continue
        right = intra_candidates.get((subject_id, query_id))
        if right is None:
            continue
        temp_edges.append(
            (
                query_id,
                subject_id,
                left.query_taxon,
                score_from_records(left, right, zero_cutoff=0.01),
            )
        )

    ortholog_ids = {edge.seq_a for edge in orthologs} | {edge.seq_b for edge in orthologs}
    taxon_all_scores: dict[str, list[float]] = defaultdict(list)
    taxon_orth_scores: dict[str, list[float]] = defaultdict(list)
    for seq_a, seq_b, taxon, score in temp_edges:
        taxon_all_scores[taxon].append(score)
        if seq_a in ortholog_ids or seq_b in ortholog_ids:
            taxon_orth_scores[taxon].append(score)

    avg_by_taxon: dict[str, float] = {}
    for taxon, scores in taxon_all_scores.items():
        avg_by_taxon[taxon] = average(taxon_orth_scores[taxon]) if taxon_orth_scores[taxon] else average(scores)

    return sorted(
        [
            EdgeRecord(
                seq_a=seq_a,
                seq_b=seq_b,
                taxon_a=taxon,
                taxon_b=taxon,
                unnormalized_score=score,
                normalized_score=score / avg_by_taxon[taxon],
                edge_type="inparalog",
            )
            for seq_a, seq_b, taxon, score in temp_edges
        ],
        key=lambda edge: (edge.taxon_a, edge.seq_a, edge.seq_b),
    )


def build_coorthologs(
    similarities: list[SimilarityRecord],
    by_pair: dict[tuple[str, str], SimilarityRecord],
    orthologs: list[EdgeRecord],
    inparalogs: list[EdgeRecord],
    percent_match_cutoff: float,
    evalue_exp_cutoff: int,
) -> list[EdgeRecord]:
    inparalog_2way = {(edge.seq_a, edge.seq_b) for edge in inparalogs} | {
        (edge.seq_b, edge.seq_a) for edge in inparalogs
    }
    ortholog_2way = {(edge.seq_a, edge.seq_b) for edge in orthologs} | {
        (edge.seq_b, edge.seq_a) for edge in orthologs
    }
    ortholog_set = {canonical_pair(edge.seq_a, edge.seq_b) for edge in orthologs}

    candidate_pairs: set[tuple[str, str]] = set()
    for left_a, left_b in inparalog_2way:
        for ortho_a, ortho_b in ortholog_2way:
            if left_b != ortho_a:
                continue
            candidate_pairs.add(canonical_pair(left_a, ortho_b))

    for ip1_a, ip1_b in inparalog_2way:
        for ortho_a, ortho_b in ortholog_2way:
            if ip1_b != ortho_a:
                continue
            for ip2_a, ip2_b in inparalog_2way:
                if ortho_b != ip2_a:
                    continue
                candidate_pairs.add(canonical_pair(ip1_a, ip2_b))

    candidate_pairs = {pair for pair in candidate_pairs if pair not in ortholog_set and pair[0] != pair[1]}

    temp_edges: list[tuple[str, str, str, str, float]] = []
    for seq_a, seq_b in sorted(candidate_pairs):
        left = by_pair.get((seq_a, seq_b))
        right = by_pair.get((seq_b, seq_a))
        if left is None or right is None:
            continue
        if not passes_thresholds(left, percent_match_cutoff, evalue_exp_cutoff):
            continue
        if not passes_thresholds(right, percent_match_cutoff, evalue_exp_cutoff):
            continue
        temp_edges.append(
            (
                seq_a,
                seq_b,
                left.query_taxon,
                left.subject_taxon,
                score_from_records(left, right, zero_cutoff=0.00001),
            )
        )

    avg_by_taxon_pair: dict[tuple[str, str], float] = {}
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for _, _, taxon_a, taxon_b, score in temp_edges:
        grouped[tuple(sorted((taxon_a, taxon_b)))].append(score)
    for key, scores in grouped.items():
        avg_by_taxon_pair[key] = average(scores)

    return sorted(
        [
            EdgeRecord(
                seq_a=seq_a,
                seq_b=seq_b,
                taxon_a=taxon_a,
                taxon_b=taxon_b,
                unnormalized_score=score,
                normalized_score=score / avg_by_taxon_pair[tuple(sorted((taxon_a, taxon_b)))],
                edge_type="coortholog",
            )
            for seq_a, seq_b, taxon_a, taxon_b, score in temp_edges
        ],
        key=lambda edge: (edge.taxon_a, edge.taxon_b, edge.seq_a, edge.seq_b),
    )
