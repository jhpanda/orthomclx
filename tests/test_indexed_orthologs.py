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


@unittest.skipUnless(shutil.which("cc") and shutil.which("make"), "C toolchain not available")
class IndexedOrthologsTest(unittest.TestCase):
    def test_indexed_ortholog_stage_emits_summary_and_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture = ROOT / "tests" / "fixtures" / "toy_pairs" / "similarSequences.txt"
            compiled_dir = tmp_path / "compiled"
            out_dir = tmp_path / "indexed"
            compile_similarities(fixture, compiled_dir)

            subprocess.run(["make", "build-indexed-orthologs"], cwd=ROOT, check=True)
            subprocess.run(
                [
                    str(ROOT / "build" / "indexed_orthologs"),
                    str(compiled_dir / "similarities.bin"),
                    str(compiled_dir / "proteins.tsv"),
                    str(compiled_dir / "taxa.tsv"),
                    str(out_dir),
                    "50",
                    "-5",
                ],
                check=True,
            )

            summary = (out_dir / "orthologs.indexed.summary.tsv").read_text()
            orthologs = (out_dir / "orthologs.indexed.txt").read_text()
            self.assertIn("records\t12\n", summary)
            self.assertIn("ortholog_edges\t2\n", summary)
            self.assertEqual(
                orthologs,
                "a1\tb1\t0.952\n"
                "a2\tb2\t1.048\n",
            )


if __name__ == "__main__":
    unittest.main()
