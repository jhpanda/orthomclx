from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.evalue import parse_evalue_cutoff


class EvalueTest(unittest.TestCase):
    def test_parse_real_evalue_cutoff(self) -> None:
        self.assertEqual(parse_evalue_cutoff(1e-3), (1.0, -3))
        self.assertEqual(parse_evalue_cutoff(0.05), (5.0, -2))

    def test_parse_legacy_exponent_cutoff(self) -> None:
        self.assertEqual(parse_evalue_cutoff(-3), (1.0, -3))


if __name__ == "__main__":
    unittest.main()
