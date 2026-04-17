from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.fasta import adjust_fasta, filter_fasta_dir


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


class FastaToolsTest(unittest.TestCase):
    def test_adjust_fasta_rewrites_headers(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_fasta = tmp_path / "input.fa"
            output_fasta = tmp_path / "hsa.fasta"
            write_file(
                input_fasta,
                ">gi|89106888|ref|AP_000668.1|\nMPEP\n>gi|2|ref|AP_000669.1|\nAAAA\n",
            )

            adjust_fasta("hsa", input_fasta, 4, output_fasta)

            self.assertEqual(
                output_fasta.read_text(),
                ">hsa|AP_000668.1\nMPEP\n>hsa|AP_000669.1\nAAAA\n",
            )

    def test_filter_fasta_splits_good_and_poor_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_dir = tmp_path / "compliantFasta"
            good_output = tmp_path / "goodProteins.fasta"
            poor_output = tmp_path / "poorProteins.fasta"

            write_file(
                fasta_dir / "aaa.fasta",
                ">aaa|gene1\nMPEPTIDE\n>aaa|gene2\nAA**AA\n>aaa|gene3\nAAA\n",
            )
            write_file(
                fasta_dir / "bbb.fasta",
                ">bbb|gene4\nMMMMMM\n",
            )

            reject_rates = filter_fasta_dir(
                fasta_dir,
                min_length=5,
                max_percent_stop=20,
                good_output_path=good_output,
                poor_output_path=poor_output,
            )

            self.assertEqual(
                good_output.read_text(),
                ">aaa|gene1\nMPEPTIDE\n>bbb|gene4\nMMMMMM\n",
            )
            self.assertEqual(
                poor_output.read_text(),
                ">aaa|gene2\nAAAA\n>aaa|gene3\nAAA\n",
            )
            self.assertEqual(reject_rates, [("aaa.fasta", 66.66666666666666)])


if __name__ == "__main__":
    unittest.main()
