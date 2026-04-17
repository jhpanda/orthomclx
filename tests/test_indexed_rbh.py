from __future__ import annotations

import shutil
import subprocess
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
class IndexedRbhTest(unittest.TestCase):
    def test_indexed_rbh_stage_emits_strict_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture = ROOT / "tests" / "fixtures" / "toy_pairs" / "similarSequences.txt"
            compiled_dir = tmp_path / "compiled"
            out_dir = tmp_path / "rbh"
            compile_similarities(fixture, compiled_dir)

            subprocess.run(["make", "build-indexed-rbh"], cwd=ROOT, check=True)
            subprocess.run(
                [
                    str(ROOT / "build" / "indexed_rbh"),
                    str(compiled_dir / "similarities.bin"),
                    str(compiled_dir / "proteins.tsv"),
                    str(compiled_dir / "taxa.tsv"),
                    str(out_dir),
                    "50",
                    "-5",
                ],
                check=True,
            )

            self.assertEqual(
                (out_dir / "rbh_1to1.txt").read_text(),
                "a1\tb1\n"
                "a2\tb2\n",
            )

    def test_indexed_pairs_includes_rbh_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture = ROOT / "tests" / "fixtures" / "toy_pairs" / "similarSequences.txt"
            compiled_dir = tmp_path / "compiled"
            out_dir = tmp_path / "pairs"
            compile_similarities(fixture, compiled_dir)

            summary = run_indexed_pairs(compiled_dir, out_dir, 50.0, -5, jobs=1)

            self.assertEqual(summary.rbh_count, 2)
            self.assertEqual(
                (out_dir / "pairs" / "rbh_1to1.txt").read_text(),
                "a1\tb1\n"
                "a2\tb2\n",
            )

    def test_indexed_pairs_parallel_rbh_matches_single_worker_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture = ROOT / "tests" / "fixtures" / "toy_pairs" / "similarSequences.txt"
            compiled_dir = tmp_path / "compiled"
            single_dir = tmp_path / "pairs_single"
            parallel_dir = tmp_path / "pairs_parallel"
            compile_similarities(fixture, compiled_dir)

            single = run_indexed_pairs(compiled_dir, single_dir, 50.0, -5, jobs=1)
            parallel = run_indexed_pairs(compiled_dir, parallel_dir, 50.0, -5, jobs=2)

            self.assertEqual(single.rbh_count, parallel.rbh_count)
            self.assertEqual(
                (single_dir / "pairs" / "rbh_1to1.txt").read_text(),
                (parallel_dir / "pairs" / "rbh_1to1.txt").read_text(),
            )


if __name__ == "__main__":
    unittest.main()
