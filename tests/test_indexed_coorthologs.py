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
class IndexedCoorthologsTest(unittest.TestCase):
    def test_indexed_coortholog_stage_emits_expected_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture = ROOT / "tests" / "fixtures" / "toy_pairs" / "similarSequences.txt"
            compiled_dir = tmp_path / "compiled"
            orth_dir = tmp_path / "orth"
            inpara_dir = tmp_path / "inpara"
            co_dir = tmp_path / "co"
            compile_similarities(fixture, compiled_dir)

            subprocess.run(["make", "build-indexed-orthologs"], cwd=ROOT, check=True)
            subprocess.run(
                [
                    str(ROOT / "build" / "indexed_orthologs"),
                    str(compiled_dir / "similarities.bin"),
                    str(compiled_dir / "proteins.tsv"),
                    str(compiled_dir / "taxa.tsv"),
                    str(orth_dir),
                    "50",
                    "-5",
                ],
                check=True,
            )

            subprocess.run(["make", "build-indexed-inparalogs"], cwd=ROOT, check=True)
            subprocess.run(
                [
                    str(ROOT / "build" / "indexed_inparalogs"),
                    str(compiled_dir / "similarities.bin"),
                    str(compiled_dir / "proteins.tsv"),
                    str(compiled_dir / "taxa.tsv"),
                    str(orth_dir / "orthologs.indexed.txt"),
                    str(inpara_dir),
                    "50",
                    "-5",
                ],
                check=True,
            )

            subprocess.run(["make", "build-indexed-coorthologs"], cwd=ROOT, check=True)
            subprocess.run(
                [
                    str(ROOT / "build" / "indexed_coorthologs"),
                    str(compiled_dir / "similarities.bin"),
                    str(compiled_dir / "proteins.tsv"),
                    str(compiled_dir / "taxa.tsv"),
                    str(orth_dir / "orthologs.indexed.txt"),
                    str(inpara_dir / "inparalogs.indexed.txt"),
                    str(co_dir),
                    "50",
                    "-5",
                ],
                check=True,
            )

            summary = (co_dir / "coorthologs.indexed.summary.tsv").read_text()
            output = (co_dir / "coorthologs.indexed.txt").read_text()
            self.assertIn("records\t12\n", summary)
            self.assertIn("coortholog_edges\t2\n", summary)
            self.assertEqual(
                output,
                "a1\tb2\t1.000\n"
                "a2\tb1\t1.000\n",
            )


if __name__ == "__main__":
    unittest.main()
