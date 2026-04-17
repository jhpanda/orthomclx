from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.rbh import read_blast_m8_hits, read_rbh_pairs, reciprocal_best_hits


class RbhTest(unittest.TestCase):
    def test_real_blast_file_matches_uploaded_rbh_pairs(self) -> None:
        blast_path = ROOT / "tests" / "Dmel_Dsub_blast.tsv"
        rbh_path = ROOT / "tests" / "Dsub_Dmel.RBH.longest.tsv"

        calculated = reciprocal_best_hits(read_blast_m8_hits(blast_path))
        expected = read_rbh_pairs(rbh_path)

        self.assertEqual(calculated, expected)
        self.assertEqual(len(calculated), 100)


if __name__ == "__main__":
    unittest.main()
