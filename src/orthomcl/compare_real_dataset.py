from __future__ import annotations

import argparse
import filecmp
import re
from dataclasses import dataclass
from pathlib import Path

from orthomcl.blast_parser import parse_blast_m8
from orthomcl.pairs import build_pairs


ROOT = Path(__file__).resolve().parents[2]


@dataclass(slots=True)
class FileSummary:
    path: Path
    size_bytes: int
    line_count: int


def count_lines(path: Path) -> int:
    count = 0
    with path.open("rb") as handle:
        for _ in handle:
            count += 1
    return count


def summarize_file(path: Path) -> FileSummary:
    return FileSummary(path=path, size_bytes=path.stat().st_size, line_count=count_lines(path))


def choose_blast_source(dataset_dir: Path, explicit: str | None) -> Path:
    if explicit:
        candidate = dataset_dir / explicit
        if not candidate.exists():
            raise FileNotFoundError(f"Requested blast source not found: {candidate}")
        return candidate

    script_path = dataset_dir / "orthomcl.sh"
    if script_path.exists():
        match = re.search(r'blastp_out="[^"]*/([^"/]+)"', script_path.read_text())
        if match:
            candidate = dataset_dir / match.group(1)
            if candidate.exists():
                return candidate

    for fallback in ["blast_out2", "blast_out"]:
        candidate = dataset_dir / fallback
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not determine BLAST source file")


def compare_files(left: Path, right: Path, sample_limit: int = 5) -> tuple[bool, list[str]]:
    if filecmp.cmp(left, right, shallow=False):
        return True, []

    mismatches: list[str] = []
    with left.open() as lf, right.open() as rf:
        for line_number, (l_line, r_line) in enumerate(zip(lf, rf), start=1):
            if l_line != r_line:
                mismatches.append(
                    f"line {line_number}: left={l_line.rstrip()!r} right={r_line.rstrip()!r}"
                )
                if len(mismatches) >= sample_limit:
                    break

    if len(mismatches) < sample_limit:
        left_count = count_lines(left)
        right_count = count_lines(right)
        if left_count != right_count:
            mismatches.append(f"line counts differ: left={left_count} right={right_count}")
    return False, mismatches


def print_summary(dataset_dir: Path, blast_source: Path) -> None:
    print(f"Dataset: {dataset_dir}")
    print(f"Chosen BLAST source: {blast_source.name}")

    important = [
        dataset_dir / "similarSequences.txt",
        dataset_dir / "similarSequence_orthomcl.txt",
        dataset_dir / "pairs" / "orthologs.txt",
        dataset_dir / "pairs" / "inparalogs.txt",
        dataset_dir / "pairs" / "coorthologs.txt",
        dataset_dir / "mclInput",
        dataset_dir / "groups.txt",
    ]
    for path in important:
        summary = summarize_file(path)
        print(f"{path.name}: {summary.line_count} lines, {summary.size_bytes} bytes")

    identical, _ = compare_files(
        dataset_dir / "similarSequences.txt",
        dataset_dir / "similarSequence_orthomcl.txt",
    )
    print(
        "saved similarSequences.txt vs similarSequence_orthomcl.txt: "
        + ("identical" if identical else "different")
    )


def run_parser_compare(dataset_dir: Path, blast_source: Path, out_dir: Path) -> bool:
    out_dir.mkdir(parents=True, exist_ok=True)
    generated = out_dir / "similarSequences.generated.txt"
    with generated.open("w") as handle:
        parse_blast_m8(blast_source, dataset_dir / "compliantFasta", handle)

    target = dataset_dir / "similarSequence_orthomcl.txt"
    identical, mismatches = compare_files(generated, target)
    print(f"parser compare against {target.name}: {'MATCH' if identical else 'DIFFER'}")
    if not identical:
        for mismatch in mismatches:
            print(f"  {mismatch}")
    return identical


def run_pairs_compare(dataset_dir: Path, parsed_input: Path, out_dir: Path, engine: str) -> bool:
    build_pairs(
        parsed_input,
        out_dir,
        percent_match_cutoff=30.0,
        evalue_exp_cutoff=-3,
        engine=engine,
    )

    checks = [
        (out_dir / "pairs" / "orthologs.txt", dataset_dir / "pairs" / "orthologs.txt"),
        (out_dir / "pairs" / "inparalogs.txt", dataset_dir / "pairs" / "inparalogs.txt"),
        (out_dir / "pairs" / "coorthologs.txt", dataset_dir / "pairs" / "coorthologs.txt"),
        (out_dir / "mclInput", dataset_dir / "mclInput"),
    ]
    all_ok = True
    for generated, target in checks:
        identical, mismatches = compare_files(generated, target)
        print(f"{generated.name} vs {target.name}: {'MATCH' if identical else 'DIFFER'}")
        if not identical:
            all_ok = False
            for mismatch in mismatches:
                print(f"  {mismatch}")
    return all_ok


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Summary-first comparison harness for the saved real OrthoMCL dataset."
    )
    parser.add_argument(
        "--dataset-dir",
        default=str(ROOT / "tests" / "orthomcl_final"),
        help="Path to the real dataset directory",
    )
    parser.add_argument(
        "--blast-source",
        choices=["blast_out", "blast_out2"],
        help="Force which BLAST file to use; otherwise auto-detect from orthomcl.sh",
    )
    parser.add_argument(
        "--work-dir",
        default=str(ROOT / "tmp" / "real_dataset_compare"),
        help="Directory for generated comparison outputs",
    )
    parser.add_argument(
        "--run-parser",
        action="store_true",
        help="Regenerate similarSequences.txt from the chosen BLAST file and compare it",
    )
    parser.add_argument(
        "--run-pairs",
        action="store_true",
        help="Regenerate pairs and mclInput from the parsed similarities and compare them",
    )
    parser.add_argument(
        "--parsed-input",
        help="Use this parsed similarity file for --run-pairs; defaults to generated parser output if available, otherwise similarSequence_orthomcl.txt",
    )
    parser.add_argument(
        "--engine",
        choices=["python", "c", "auto"],
        default="c",
        help="Execution engine for the pair-building stage",
    )
    args = parser.parse_args(argv)

    dataset_dir = Path(args.dataset_dir)
    work_dir = Path(args.work_dir)
    blast_source = choose_blast_source(dataset_dir, args.blast_source)

    print_summary(dataset_dir, blast_source)

    parser_ok = True
    generated_parsed = work_dir / "parser" / "similarSequences.generated.txt"
    if args.run_parser:
        parser_ok = run_parser_compare(dataset_dir, blast_source, work_dir / "parser")

    pairs_ok = True
    if args.run_pairs:
        if args.parsed_input:
            parsed_input = Path(args.parsed_input)
        elif args.run_parser:
            parsed_input = generated_parsed
        else:
            parsed_input = dataset_dir / "similarSequence_orthomcl.txt"
        pairs_ok = run_pairs_compare(dataset_dir, parsed_input, work_dir / "pairs", args.engine)

    if parser_ok and pairs_ok:
        print("Comparison completed without detected mismatches.")
        return 0
    print("Comparison completed with mismatches.")
    return 1

