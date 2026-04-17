from __future__ import annotations

import sys
import tempfile
import unittest
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.pairs import build_pairs


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


class PairsTest(unittest.TestCase):
    def test_build_pairs_writes_expected_toy_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture_dir = ROOT / "tests" / "fixtures" / "toy_pairs"
            similar_sequences = tmp_path / "similarSequences.txt"
            write_file(similar_sequences, (fixture_dir / "similarSequences.txt").read_text())

            results = build_pairs(
                similar_sequences,
                tmp_path / "run",
                percent_match_cutoff=50.0,
                evalue_exp_cutoff=-5,
            )

            self.assertEqual(
                [(edge.seq_a, edge.seq_b) for edge in results["orthologs"]],
                [("a1", "b1"), ("a2", "b2")],
            )
            self.assertEqual(
                [(edge.seq_a, edge.seq_b) for edge in results["inparalogs"]],
                [("a1", "a2"), ("b1", "b2")],
            )
            self.assertEqual(
                [(edge.seq_a, edge.seq_b) for edge in results["coorthologs"]],
                [("a1", "b2"), ("a2", "b1")],
            )

            self.assertEqual(
                (tmp_path / "run" / "pairs" / "orthologs.txt").read_text(),
                (fixture_dir / "expected" / "orthologs.txt").read_text(),
            )
            self.assertEqual(
                (tmp_path / "run" / "pairs" / "inparalogs.txt").read_text(),
                (fixture_dir / "expected" / "inparalogs.txt").read_text(),
            )
            self.assertEqual(
                (tmp_path / "run" / "pairs" / "coorthologs.txt").read_text(),
                (fixture_dir / "expected" / "coorthologs.txt").read_text(),
            )
            self.assertEqual(
                (tmp_path / "run" / "mclInput").read_text(),
                (fixture_dir / "expected" / "mclInput").read_text(),
            )

    @unittest.skipUnless(shutil.which("cc") and shutil.which("make"), "C toolchain not available")
    def test_c_engine_matches_python_outputs_on_toy_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture_dir = ROOT / "tests" / "fixtures" / "toy_pairs"

            build_pairs(
                fixture_dir / "similarSequences.txt",
                tmp_path / "python",
                percent_match_cutoff=50.0,
                evalue_exp_cutoff=-5,
                engine="python",
            )
            build_pairs(
                fixture_dir / "similarSequences.txt",
                tmp_path / "c",
                percent_match_cutoff=50.0,
                evalue_exp_cutoff=-5,
                engine="c",
            )

            self.assertEqual(
                (tmp_path / "python" / "pairs" / "orthologs.txt").read_text(),
                (tmp_path / "c" / "pairs" / "orthologs.txt").read_text(),
            )
            self.assertEqual(
                (tmp_path / "python" / "pairs" / "inparalogs.txt").read_text(),
                (tmp_path / "c" / "pairs" / "inparalogs.txt").read_text(),
            )
            self.assertEqual(
                (tmp_path / "python" / "pairs" / "coorthologs.txt").read_text(),
                (tmp_path / "c" / "pairs" / "coorthologs.txt").read_text(),
            )
            self.assertEqual(
                (tmp_path / "python" / "mclInput").read_text(),
                (tmp_path / "c" / "mclInput").read_text(),
            )


if __name__ == "__main__":
    unittest.main()
