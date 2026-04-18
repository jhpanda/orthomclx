from __future__ import annotations

import shutil
import stat
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orthomcl.mcl import MclRunConfig, run_mcl
from orthomcl.pipeline import IntegratedRunConfig, run_integrated_pipeline


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def make_fake_mcl(path: Path) -> None:
    script = """#!/usr/bin/env python3
import pathlib
import sys

args = sys.argv[1:]
output_path = pathlib.Path(args[args.index("-o") + 1])
output_path.write_text("TAXA|A\\tTAXB|B\\n")
"""
    path.write_text(script)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


class MclTest(unittest.TestCase):
    def test_run_mcl_invokes_external_binary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            mcl_input = tmp_path / "mclInput"
            fake_mcl = tmp_path / "fake_mcl"
            output_path = tmp_path / "mclOutput"
            write_file(mcl_input, "TAXA|A\tTAXB|B\t1.0\n")
            make_fake_mcl(fake_mcl)

            run_mcl(
                MclRunConfig(
                    input_path=mcl_input,
                    output_path=output_path,
                    binary=str(fake_mcl),
                    inflation=2.0,
                    threads=3,
                )
            )

            self.assertEqual(output_path.read_text(), "TAXA|A\tTAXB|B\n")

    @unittest.skipUnless(shutil.which("cc") and shutil.which("make"), "C toolchain not available")
    def test_integrated_pipeline_can_run_mcl_and_write_groups(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_dir = tmp_path / "fasta"
            blast_path = tmp_path / "blast.tsv"
            out_dir = tmp_path / "out"
            fake_mcl = tmp_path / "fake_mcl"

            write_file(fasta_dir / "TAXA.fasta", ">TAXA|A\nMPEPTIDE\n")
            write_file(fasta_dir / "TAXB.fasta", ">TAXB|B\nMPEPTIDE\n")
            write_file(
                blast_path,
                "TAXA|A\tTAXB|B\t99.0\t8\t0\t0\t1\t8\t1\t8\t1e-50\t200\n"
                "TAXB|B\tTAXA|A\t99.0\t8\t0\t0\t1\t8\t1\t8\t1e-50\t200\n",
            )
            make_fake_mcl(fake_mcl)

            summary = run_integrated_pipeline(
                IntegratedRunConfig(
                    blast_file=blast_path,
                    fasta_dir=fasta_dir,
                    out_dir=out_dir,
                    percent_match_cutoff=50.0,
                    evalue_exp_cutoff=-5,
                    threads=1,
                    run_mcl=True,
                    mcl_binary=str(fake_mcl),
                    mcl_inflation=1.5,
                    mcl_threads=2,
                    groups_prefix="OG",
                    start_group_id=1000,
                )
            )

            self.assertEqual(summary.rbh_count, 1)
            self.assertTrue((out_dir / "mclInput").exists())
            self.assertEqual((out_dir / "mclOutput").read_text(), "TAXA|A\tTAXB|B\n")
            self.assertEqual((out_dir / "groups.txt").read_text(), "OG1000: TAXA|A TAXB|B\n")


if __name__ == "__main__":
    unittest.main()
