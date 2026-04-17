from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.compile_similarities import compile_similarities
from orthomcl.indexed_pairs import run_indexed_pairs


@unittest.skipUnless(shutil.which("cc") and shutil.which("make"), "C toolchain not available")
class IndexedPairsTest(unittest.TestCase):
    def test_indexed_pairs_single_worker_matches_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture_dir = ROOT / "tests" / "fixtures" / "toy_pairs"
            compiled_dir = tmp_path / "compiled"
            out_dir = tmp_path / "indexed_pairs"
            compile_similarities(fixture_dir / "similarSequences.txt", compiled_dir)

            summary = run_indexed_pairs(compiled_dir, out_dir, 50.0, -5, jobs=1)

            self.assertEqual(summary.ortholog_count, 2)
            self.assertEqual(summary.inparalog_count, 2)
            self.assertEqual(summary.coortholog_count, 2)
            self.assertEqual(summary.mcl_count, 6)
            self.assertEqual(
                (out_dir / "pairs" / "orthologs.txt").read_text(),
                (fixture_dir / "expected" / "orthologs.txt").read_text(),
            )
            self.assertEqual(
                (out_dir / "pairs" / "inparalogs.txt").read_text(),
                (fixture_dir / "expected" / "inparalogs.txt").read_text(),
            )
            self.assertEqual(
                (out_dir / "pairs" / "coorthologs.txt").read_text(),
                (fixture_dir / "expected" / "coorthologs.txt").read_text(),
            )
            self.assertEqual(
                (out_dir / "mclInput").read_text(),
                (fixture_dir / "expected" / "mclInput").read_text(),
            )

    def test_indexed_pairs_multi_worker_matches_single_worker(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture_dir = ROOT / "tests" / "fixtures" / "toy_pairs"
            compiled_dir = tmp_path / "compiled"
            single_out = tmp_path / "single"
            multi_out = tmp_path / "multi"
            compile_similarities(fixture_dir / "similarSequences.txt", compiled_dir)

            run_indexed_pairs(compiled_dir, single_out, 50.0, -5, jobs=1)
            run_indexed_pairs(compiled_dir, multi_out, 50.0, -5, jobs=2)

            self.assertEqual(
                (single_out / "pairs" / "orthologs.txt").read_text(),
                (multi_out / "pairs" / "orthologs.txt").read_text(),
            )
            self.assertEqual(
                (single_out / "pairs" / "inparalogs.txt").read_text(),
                (multi_out / "pairs" / "inparalogs.txt").read_text(),
            )
            self.assertEqual(
                (single_out / "pairs" / "coorthologs.txt").read_text(),
                (multi_out / "pairs" / "coorthologs.txt").read_text(),
            )
            self.assertEqual(
                (single_out / "mclInput").read_text(),
                (multi_out / "mclInput").read_text(),
            )


if __name__ == "__main__":
    unittest.main()
