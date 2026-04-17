from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.groups import mcl_to_groups_file


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


class GroupsTest(unittest.TestCase):
    def test_mcl_to_groups_matches_legacy_format(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            mcl_output = tmp_path / "mclOutput"
            groups_output = tmp_path / "groups.txt"
            fixture_dir = ROOT / "tests" / "fixtures" / "toy_pairs"
            write_file(mcl_output, (fixture_dir / "mclOutput").read_text())

            mcl_to_groups_file(mcl_output, groups_output, "OG000", 1000)

            self.assertEqual(
                groups_output.read_text(),
                (fixture_dir / "expected" / "groups.txt").read_text(),
            )


if __name__ == "__main__":
    unittest.main()
