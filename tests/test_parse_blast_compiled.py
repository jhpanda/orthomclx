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
from orthomcl.parse_blast_compiled import parse_blast_to_compiled


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


@unittest.skipUnless(shutil.which("cc") and shutil.which("make"), "C toolchain not available")
class ParseBlastCompiledTest(unittest.TestCase):
    def test_c_direct_parser_matches_expected_compiled_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_dir = tmp_path / "fasta"
            blast_path = tmp_path / "blast.tsv"
            expected_similar = tmp_path / "similarSequences.txt"
            expected_compiled = tmp_path / "expected_compiled"
            c_compiled = tmp_path / "c_compiled"

            write_file(
                fasta_dir / "TAXA.fasta",
                ">TAXA|A1\nMPEPTIDE\n>TAXA|A2\nMPEPTIDER\n",
            )
            write_file(
                fasta_dir / "TAXB.fasta",
                ">TAXB|B1\nMPEPTIDE\n>TAXB|B2\nMPEPTIDER\n",
            )
            write_file(
                blast_path,
                "\n".join(
                    [
                        "TAXA|A1\tTAXB|B1\t100.0\t4\t0\t0\t1\t4\t1\t4\t1e-50\t200",
                        "TAXA|A1\tTAXB|B1\t100.0\t4\t0\t0\t5\t8\t5\t8\t1e-50\t200",
                        "TAXB|B1\tTAXA|A1\t100.0\t8\t0\t0\t1\t8\t1\t8\t1e-50\t200",
                        "TAXA|A2\tTAXB|B2\t87.5\t8\t1\t0\t1\t8\t1\t8\t0\t180",
                        "TAXB|B2\tTAXA|A2\t87.5\t8\t1\t0\t1\t8\t1\t8\t0\t180",
                    ]
                )
                + "\n",
            )
            write_file(
                expected_similar,
                "\n".join(
                    [
                        "TAXA|A1\tTAXB|B1\tTAXA\tTAXB\t1\t-50\t100.0\t100.0",
                        "TAXB|B1\tTAXA|A1\tTAXB\tTAXA\t1\t-50\t100.0\t100.0",
                        "TAXA|A2\tTAXB|B2\tTAXA\tTAXB\t0\t0\t87.5\t88.9",
                        "TAXB|B2\tTAXA|A2\tTAXB\tTAXA\t0\t0\t87.5\t88.9",
                    ]
                )
                + "\n",
            )

            expected_summary = compile_similarities(expected_similar, expected_compiled)
            c_summary = parse_blast_to_compiled(blast_path, fasta_dir, c_compiled)

            self.assertEqual(c_summary.record_count, expected_summary.record_count)
            self.assertEqual(c_summary.protein_count, expected_summary.protein_count)
            self.assertEqual(c_summary.taxon_count, expected_summary.taxon_count)
            self.assertEqual(c_summary.proteins_path.read_text(), expected_summary.proteins_path.read_text())
            self.assertEqual(c_summary.taxa_path.read_text(), expected_summary.taxa_path.read_text())
            self.assertEqual(c_summary.binary_path.read_bytes(), expected_summary.binary_path.read_bytes())


if __name__ == "__main__":
    unittest.main()
