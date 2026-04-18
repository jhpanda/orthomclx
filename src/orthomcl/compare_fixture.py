from __future__ import annotations

import argparse
import difflib
import filecmp
from pathlib import Path
from typing import List, Optional

from orthomcl.groups import mcl_to_groups_file
from orthomcl.pairs import build_pairs


ROOT = Path(__file__).resolve().parents[2]


def read_text(path: Path) -> str:
    return path.read_text() if path.exists() else ""


def assert_file_equal(actual: Path, expected: Path) -> None:
    if filecmp.cmp(actual, expected, shallow=False):
        return
    diff = difflib.unified_diff(
        read_text(expected).splitlines(keepends=True),
        read_text(actual).splitlines(keepends=True),
        fromfile=str(expected),
        tofile=str(actual),
    )
    raise AssertionError("".join(diff))


def compare_python_outputs(fixture_dir: Path, work_dir: Path) -> None:
    build_pairs(
        fixture_dir / "similarSequences.txt",
        work_dir / "python",
        percent_match_cutoff=50.0,
        evalue_exp_cutoff=-5,
    )
    mcl_to_groups_file(
        fixture_dir / "mclOutput",
        work_dir / "python" / "groups.txt",
        "OG000",
        1000,
    )
    expected_dir = fixture_dir / "expected"
    assert_file_equal(work_dir / "python" / "pairs" / "orthologs.txt", expected_dir / "orthologs.txt")
    assert_file_equal(work_dir / "python" / "pairs" / "inparalogs.txt", expected_dir / "inparalogs.txt")
    assert_file_equal(work_dir / "python" / "pairs" / "coorthologs.txt", expected_dir / "coorthologs.txt")
    assert_file_equal(work_dir / "python" / "mclInput", expected_dir / "mclInput")
    assert_file_equal(work_dir / "python" / "groups.txt", expected_dir / "groups.txt")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Compare Python rewrite outputs against the fixture outputs.")
    parser.add_argument(
        "--fixture",
        default=str(ROOT / "tests" / "fixtures" / "toy_pairs"),
        help="Fixture directory containing similarSequences.txt and expected/ outputs",
    )
    parser.add_argument(
        "--work-dir",
        default=str(ROOT / "tmp" / "compat"),
        help="Directory where generated outputs should be written",
    )
    args = parser.parse_args(argv)

    fixture_dir = Path(args.fixture)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    compare_python_outputs(fixture_dir, work_dir)
    print("Compatibility check passed.")
    return 0
