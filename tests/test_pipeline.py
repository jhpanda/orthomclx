from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.pipeline import RunConfig, run_pipeline


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


class PipelineTest(unittest.TestCase):
    def test_run_pipeline_creates_expected_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fixture_dir = ROOT / "tests" / "fixtures" / "toy_pipeline"
            out_dir = tmp_path / "run"

            run_pipeline(
                RunConfig(
                    blast_file=fixture_dir / "blast.tsv",
                    fasta_dir=fixture_dir / "compliantFasta",
                    out_dir=out_dir,
                    percent_match_cutoff=50.0,
                    evalue_exp_cutoff=-5,
                    mcl_output=fixture_dir / "mclOutput",
                    groups_prefix="OG000",
                    start_group_id=1000,
                )
            )

            expected_dir = fixture_dir / "expected"
            self.assertEqual(
                (out_dir / "pairs" / "orthologs.txt").read_text(),
                (expected_dir / "orthologs.txt").read_text(),
            )
            self.assertEqual(
                (out_dir / "pairs" / "inparalogs.txt").read_text(),
                (expected_dir / "inparalogs.txt").read_text(),
            )
            self.assertEqual(
                (out_dir / "pairs" / "coorthologs.txt").read_text(),
                (expected_dir / "coorthologs.txt").read_text(),
            )
            self.assertEqual(
                (out_dir / "mclInput").read_text(),
                (expected_dir / "mclInput").read_text(),
            )
            self.assertEqual(
                (out_dir / "groups.txt").read_text(),
                (expected_dir / "groups.txt").read_text(),
            )


if __name__ == "__main__":
    unittest.main()
