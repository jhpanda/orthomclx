from __future__ import annotations

import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.blast_parser import parse_blast_m8


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


class BlastParserTest(unittest.TestCase):
    def test_parse_blast_m8_matches_legacy_shape(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_dir = tmp_path / "compliantFasta"
            write_file(
                fasta_dir / "aaa.fasta",
                ">geneA\nAAAAAAAAAA\n>geneB\nMMMMMMMMMMMM\n",
            )
            write_file(
                fasta_dir / "bbb.fasta",
                ">geneC\nCCCCCCCC\n>geneD\nDDDDDDDDDD\n",
            )
            blast_file = tmp_path / "blast.tsv"
            write_file(
                blast_file,
                "\n".join(
                    [
                        "geneA geneC 100 4 0 0 1 4 1 4 1e-20 50",
                        "geneA geneC 50 4 0 0 5 8 3 6 1e-20 40",
                        "geneB geneD 80 5 0 0 2 6 1 5 2e-10 30",
                    ]
                )
                + "\n",
            )

            output_file = tmp_path / "similarSequences.txt"
            with output_file.open("w") as handle:
                parse_blast_m8(blast_file, fasta_dir, handle)

            self.assertEqual(
                output_file.read_text(),
                "geneA\tgeneC\taaa\tbbb\t1\t-20\t75.0\t75.0\n"
                "geneB\tgeneD\taaa\tbbb\t2\t-10\t80.0\t50.0\n",
            )


if __name__ == "__main__":
    unittest.main()
