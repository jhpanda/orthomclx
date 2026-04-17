from __future__ import annotations

import argparse
import sys
from pathlib import Path

from orthomcl.blast_parser import parse_blast_m8
from orthomcl.compile_similarities import compile_similarities
from orthomcl.fasta import adjust_fasta, filter_fasta_dir
from orthomcl.groups import mcl_to_groups_file
from orthomcl.indexed_pairs import run_indexed_pairs
from orthomcl.mcl import MclRunConfig, run_mcl
from orthomcl.pairs import build_pairs
from orthomcl.pipeline import IntegratedRunConfig, RunConfig, run_integrated_pipeline, run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="orthomclx")
    parser.add_argument(
        "--input",
        help="Integrated mode: directory of compliant FASTA files named as taxon.fasta",
    )
    parser.add_argument(
        "--blast",
        help="Integrated mode: BLAST output in m8 format",
    )
    parser.add_argument(
        "--out",
        "--out-dir",
        dest="integrated_out_dir",
        help="Integrated mode: output directory",
    )
    parser.add_argument(
        "--pcut",
        "--percent-match-cutoff",
        dest="integrated_percent_match_cutoff",
        type=float,
        help="Integrated mode: minimum percent match cutoff",
    )
    parser.add_argument(
        "--ecut",
        "--evalue-exp-cutoff",
        dest="integrated_evalue_exp_cutoff",
        type=int,
        help="Integrated mode: maximum allowed e-value exponent",
    )
    parser.add_argument(
        "--jobs",
        dest="integrated_jobs",
        type=int,
        default=1,
        help="Integrated mode: number of parallel shard workers to use",
    )
    parser.add_argument(
        "--run-mcl",
        dest="integrated_run_mcl",
        action="store_true",
        help="Integrated mode: run the external mcl binary on the generated mclInput",
    )
    parser.add_argument(
        "--mcl-output",
        dest="integrated_mcl_output",
        help="Integrated mode: optional MCL label-format output file to convert into groups.txt",
    )
    parser.add_argument(
        "--mcl-binary",
        dest="integrated_mcl_binary",
        default="mcl",
        help="Integrated mode: MCL executable to run when --run-mcl is used",
    )
    parser.add_argument(
        "--inflation",
        dest="integrated_mcl_inflation",
        type=float,
        default=1.5,
        help="Integrated mode: MCL inflation value when --run-mcl is used",
    )
    parser.add_argument(
        "--mcl-threads",
        dest="integrated_mcl_threads",
        type=int,
        help="Integrated mode: MCL thread count when --run-mcl is used",
    )
    parser.add_argument(
        "--groups-prefix",
        dest="integrated_groups_prefix",
        default="OG",
        help="Integrated mode: group ID prefix for groups.txt generation",
    )
    parser.add_argument(
        "--start-group-id",
        dest="integrated_start_group_id",
        type=int,
        default=1000,
        help="Integrated mode: starting numeric suffix for groups.txt generation",
    )

    subparsers = parser.add_subparsers(dest="command", required=False)

    parse_blast = subparsers.add_parser(
        "parse-blast",
        help="Parse BLAST m8 output into OrthoMCL similarity rows",
    )
    parse_blast.add_argument("blast_file", help="BLAST output in m8 format")
    parse_blast.add_argument(
        "fasta_dir",
        help="Directory of compliant FASTA files named as taxon.fasta",
    )
    parse_blast.add_argument(
        "-o",
        "--output",
        help="Write parsed output to this file instead of stdout",
    )

    adjust = subparsers.add_parser(
        "adjust-fasta",
        help="Convert FASTA headers into OrthoMCL-compliant taxon|id headers",
    )
    adjust.add_argument("taxon_code", help="Three or four letter taxon code")
    adjust.add_argument("fasta_file", help="Input FASTA file")
    adjust.add_argument("id_field", type=int, help="1-based field number holding the protein ID")
    adjust.add_argument(
        "-o",
        "--output",
        help="Output FASTA path; defaults to <taxon_code>.fasta in the current directory",
    )

    filter_parser = subparsers.add_parser(
        "filter-fasta",
        help="Filter compliant FASTA files into good and poor protein files",
    )
    filter_parser.add_argument("input_dir", help="Directory containing compliant *.fasta files")
    filter_parser.add_argument("min_length", type=int, help="Minimum allowed protein length")
    filter_parser.add_argument(
        "max_percent_stop",
        type=float,
        help="Maximum allowed percent non-letter residues in a sequence",
    )
    filter_parser.add_argument(
        "--good-output",
        default="goodProteins.fasta",
        help="Path for accepted proteins",
    )
    filter_parser.add_argument(
        "--poor-output",
        default="poorProteins.fasta",
        help="Path for rejected proteins",
    )

    pairs_parser = subparsers.add_parser(
        "pairs",
        help="Build ortholog, inparalog, and coortholog files from parsed similarities",
    )
    pairs_parser.add_argument("similar_sequences", help="Parsed similarity file")
    pairs_parser.add_argument("out_dir", help="Output directory")
    pairs_parser.add_argument(
        "--percent-match-cutoff",
        type=float,
        required=True,
        help="Minimum percent match cutoff",
    )
    pairs_parser.add_argument(
        "--evalue-exp-cutoff",
        type=int,
        required=True,
        help="Maximum allowed e-value exponent",
    )
    pairs_parser.add_argument(
        "--engine",
        choices=["auto", "python", "c"],
        default="python",
        help="Execution engine for the heavy pair-building stage",
    )

    mcl_runner = subparsers.add_parser(
        "mcl",
        help="Run the external MCL binary on an mclInput file",
    )
    mcl_runner.add_argument("mcl_input", help="Path to mclInput in ABC format")
    mcl_runner.add_argument(
        "-o",
        "--output",
        default="mclOutput",
        help="Path to the generated MCL label-format output",
    )
    mcl_runner.add_argument(
        "--binary",
        default="mcl",
        help="MCL executable to run",
    )
    mcl_runner.add_argument(
        "--inflation",
        type=float,
        default=1.5,
        help="MCL inflation value",
    )
    mcl_runner.add_argument(
        "--threads",
        type=int,
        help="Optional MCL thread count",
    )

    mcl_groups = subparsers.add_parser(
        "mcl-to-groups",
        help="Convert MCL label-format output into OrthoMCL groups.txt format",
    )
    mcl_groups.add_argument("mcl_output", help="Path to MCL output in label format")
    mcl_groups.add_argument("prefix", help="Group ID prefix, for example OG000")
    mcl_groups.add_argument("start_id", type=int, help="Starting numeric suffix")
    mcl_groups.add_argument(
        "-o",
        "--output",
        default="groups.txt",
        help="Path to the generated groups file",
    )

    run_parser = subparsers.add_parser(
        "run",
        help="Run parse-blast and pairs, and optionally convert MCL output to groups",
    )
    run_parser.add_argument("--blast", required=True, help="BLAST output in m8 format")
    run_parser.add_argument(
        "--fasta-dir",
        required=True,
        help="Directory of compliant FASTA files named as taxon.fasta",
    )
    run_parser.add_argument("--out-dir", required=True, help="Output directory")
    run_parser.add_argument(
        "--percent-match-cutoff",
        type=float,
        required=True,
        help="Minimum percent match cutoff",
    )
    run_parser.add_argument(
        "--evalue-exp-cutoff",
        type=int,
        required=True,
        help="Maximum allowed e-value exponent",
    )
    run_parser.add_argument(
        "--run-mcl",
        action="store_true",
        help="Run the external mcl binary on the generated mclInput",
    )
    run_parser.add_argument(
        "--mcl-output",
        help="Optional MCL label-format output file to convert into groups.txt",
    )
    run_parser.add_argument(
        "--mcl-binary",
        default="mcl",
        help="MCL executable to run when --run-mcl is used",
    )
    run_parser.add_argument(
        "--inflation",
        type=float,
        default=1.5,
        help="MCL inflation value when --run-mcl is used",
    )
    run_parser.add_argument(
        "--mcl-threads",
        type=int,
        help="Optional MCL thread count when --run-mcl is used",
    )
    run_parser.add_argument(
        "--groups-prefix",
        default="OG",
        help="Group ID prefix for groups.txt generation",
    )
    run_parser.add_argument(
        "--start-group-id",
        type=int,
        default=1000,
        help="Starting numeric suffix for groups.txt generation",
    )
    run_parser.add_argument(
        "--engine",
        choices=["auto", "python", "c"],
        default="python",
        help="Execution engine for the heavy pair-building stage",
    )

    compile_parser = subparsers.add_parser(
        "compile-similarities",
        help="Convert parsed similarities into integer-coded binary and index files",
    )
    compile_parser.add_argument("similar_sequences", help="Parsed similarity file")
    compile_parser.add_argument("out_dir", help="Output directory for compiled artifacts")

    indexed_orthologs = subparsers.add_parser(
        "indexed-orthologs",
        help="Run the first indexed C ortholog pass on compiled similarity binaries",
    )
    indexed_orthologs.add_argument("compiled_dir", help="Directory containing similarities.bin, proteins.tsv, and taxa.tsv")
    indexed_orthologs.add_argument("out_dir", help="Output directory for indexed ortholog artifacts")
    indexed_orthologs.add_argument(
        "--percent-match-cutoff",
        type=float,
        required=True,
        help="Minimum percent match cutoff",
    )
    indexed_orthologs.add_argument(
        "--evalue-exp-cutoff",
        type=int,
        required=True,
        help="Maximum allowed e-value exponent",
    )

    indexed_inparalogs = subparsers.add_parser(
        "indexed-inparalogs",
        help="Run the indexed C inparalog pass on compiled similarity binaries",
    )
    indexed_inparalogs.add_argument("compiled_dir", help="Directory containing similarities.bin, proteins.tsv, and taxa.tsv")
    indexed_inparalogs.add_argument("orthologs_file", help="Ortholog output file from the indexed ortholog stage")
    indexed_inparalogs.add_argument("out_dir", help="Output directory for indexed inparalog artifacts")
    indexed_inparalogs.add_argument(
        "--percent-match-cutoff",
        type=float,
        required=True,
        help="Minimum percent match cutoff",
    )
    indexed_inparalogs.add_argument(
        "--evalue-exp-cutoff",
        type=int,
        required=True,
        help="Maximum allowed e-value exponent",
    )

    indexed_coorthologs = subparsers.add_parser(
        "indexed-coorthologs",
        help="Run the indexed C coortholog pass on compiled similarity binaries",
    )
    indexed_coorthologs.add_argument("compiled_dir", help="Directory containing similarities.bin, proteins.tsv, and taxa.tsv")
    indexed_coorthologs.add_argument("orthologs_file", help="Ortholog output file from the indexed ortholog stage")
    indexed_coorthologs.add_argument("inparalogs_file", help="Inparalog output file from the indexed inparalog stage")
    indexed_coorthologs.add_argument("out_dir", help="Output directory for indexed coortholog artifacts")
    indexed_coorthologs.add_argument(
        "--percent-match-cutoff",
        type=float,
        required=True,
        help="Minimum percent match cutoff",
    )
    indexed_coorthologs.add_argument(
        "--evalue-exp-cutoff",
        type=int,
        required=True,
        help="Maximum allowed e-value exponent",
    )

    indexed_pairs = subparsers.add_parser(
        "indexed-pairs",
        help="Run the indexed multi-stage pairs workflow on compiled similarities",
    )
    indexed_pairs.add_argument("compiled_dir", help="Directory containing similarities.bin, proteins.tsv, and taxa.tsv")
    indexed_pairs.add_argument("out_dir", help="Output directory for final indexed pairs outputs")
    indexed_pairs.add_argument(
        "--percent-match-cutoff",
        type=float,
        required=True,
        help="Minimum percent match cutoff",
    )
    indexed_pairs.add_argument(
        "--evalue-exp-cutoff",
        type=int,
        required=True,
        help="Maximum allowed e-value exponent",
    )
    indexed_pairs.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel shard workers to use per indexed stage",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        if not any(
            value is not None
            for value in (
                args.input,
                args.blast,
                args.integrated_out_dir,
                args.integrated_percent_match_cutoff,
                args.integrated_evalue_exp_cutoff,
            )
        ):
            parser.print_help()
            return 0

        missing: list[str] = []
        if not args.input:
            missing.append("--input")
        if not args.blast:
            missing.append("--blast")
        if args.integrated_percent_match_cutoff is None:
            missing.append("--pcut")
        if args.integrated_evalue_exp_cutoff is None:
            missing.append("--ecut")
        if missing:
            parser.error("Integrated mode is missing required arguments: " + ", ".join(missing))

        out_dir = Path(args.integrated_out_dir) if args.integrated_out_dir else Path("orthomclx_out")
        if args.integrated_run_mcl and args.integrated_mcl_output:
            parser.error("Integrated mode cannot use --run-mcl together with --mcl-output")
        summary = run_integrated_pipeline(
            IntegratedRunConfig(
                blast_file=Path(args.blast),
                fasta_dir=Path(args.input),
                out_dir=out_dir,
                percent_match_cutoff=args.integrated_percent_match_cutoff,
                evalue_exp_cutoff=args.integrated_evalue_exp_cutoff,
                jobs=args.integrated_jobs,
                run_mcl=args.integrated_run_mcl,
                mcl_output=Path(args.integrated_mcl_output) if args.integrated_mcl_output else None,
                mcl_binary=args.integrated_mcl_binary,
                mcl_inflation=args.integrated_mcl_inflation,
                mcl_threads=args.integrated_mcl_threads,
                groups_prefix=args.integrated_groups_prefix,
                start_group_id=args.integrated_start_group_id,
            )
        )
        sys.stdout.write(
            f"integrated run complete: orthologs={summary.ortholog_count}, "
            f"inparalogs={summary.inparalog_count}, "
            f"coorthologs={summary.coortholog_count}, "
            f"rbh_1to1={summary.rbh_count}, "
            f"mcl={summary.mcl_count}\n"
        )
        return 0

    if args.command == "parse-blast":
        if args.output:
            output_path = Path(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w") as handle:
                parse_blast_m8(args.blast_file, args.fasta_dir, handle)
        else:
            parse_blast_m8(args.blast_file, args.fasta_dir, sys.stdout)
        return 0

    if args.command == "adjust-fasta":
        output_path = Path(args.output) if args.output else Path(f"{args.taxon_code}.fasta")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adjust_fasta(args.taxon_code, args.fasta_file, args.id_field, output_path)
        return 0

    if args.command == "filter-fasta":
        good_output = Path(args.good_output)
        poor_output = Path(args.poor_output)
        good_output.parent.mkdir(parents=True, exist_ok=True)
        poor_output.parent.mkdir(parents=True, exist_ok=True)
        reject_rates = filter_fasta_dir(
            args.input_dir,
            args.min_length,
            args.max_percent_stop,
            good_output,
            poor_output,
        )
        if reject_rates:
            sys.stdout.write("\nProteomes with > 10% poor proteins:\n")
            for file_name, reject_percent in reject_rates:
                sys.stdout.write(f"  {file_name}\t{int(reject_percent)}%\n")
        return 0

    if args.command == "pairs":
        build_pairs(
            args.similar_sequences,
            args.out_dir,
            args.percent_match_cutoff,
            args.evalue_exp_cutoff,
            engine=args.engine,
        )
        return 0

    if args.command == "mcl":
        output_path = Path(args.output)
        run_mcl(
            MclRunConfig(
                input_path=Path(args.mcl_input),
                output_path=output_path,
                inflation=args.inflation,
                threads=args.threads,
                binary=args.binary,
            )
        )
        sys.stdout.write(f"mcl output written to {output_path}\n")
        return 0

    if args.command == "mcl-to-groups":
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        mcl_to_groups_file(args.mcl_output, output_path, args.prefix, args.start_id)
        return 0

    if args.command == "run":
        if args.run_mcl and args.mcl_output:
            parser.error("run cannot use --run-mcl together with --mcl-output")
        run_pipeline(
            RunConfig(
                blast_file=Path(args.blast),
                fasta_dir=Path(args.fasta_dir),
                out_dir=Path(args.out_dir),
                percent_match_cutoff=args.percent_match_cutoff,
                evalue_exp_cutoff=args.evalue_exp_cutoff,
                run_mcl=args.run_mcl,
                mcl_output=Path(args.mcl_output) if args.mcl_output else None,
                mcl_binary=args.mcl_binary,
                mcl_inflation=args.inflation,
                mcl_threads=args.mcl_threads,
                groups_prefix=args.groups_prefix,
                start_group_id=args.start_group_id,
                engine=args.engine,
            )
        )
        return 0

    if args.command == "compile-similarities":
        summary = compile_similarities(args.similar_sequences, args.out_dir)
        sys.stdout.write(
            f"compiled {summary.record_count} records, "
            f"{summary.protein_count} proteins, "
            f"{summary.taxon_count} taxa\n"
        )
        sys.stdout.write(f"binary: {summary.binary_path}\n")
        sys.stdout.write(f"proteins: {summary.proteins_path}\n")
        sys.stdout.write(f"taxa: {summary.taxa_path}\n")
        return 0

    if args.command == "indexed-orthologs":
        root = Path(__file__).resolve().parents[2]
        binary = root / "build" / "indexed_orthologs"
        if not binary.exists():
            import subprocess

            subprocess.run(["make", "build-indexed-orthologs"], cwd=root, check=True)
        compiled_dir = Path(args.compiled_dir)
        output_path = Path(args.out_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        import subprocess

        subprocess.run(
            [
                str(binary),
                str(compiled_dir / "similarities.bin"),
                str(compiled_dir / "proteins.tsv"),
                str(compiled_dir / "taxa.tsv"),
                str(output_path),
                str(args.percent_match_cutoff),
                str(args.evalue_exp_cutoff),
            ],
            check=True,
        )
        sys.stdout.write(f"indexed ortholog outputs written to {output_path}\n")
        return 0

    if args.command == "indexed-inparalogs":
        root = Path(__file__).resolve().parents[2]
        binary = root / "build" / "indexed_inparalogs"
        if not binary.exists():
            import subprocess

            subprocess.run(["make", "build-indexed-inparalogs"], cwd=root, check=True)
        compiled_dir = Path(args.compiled_dir)
        output_path = Path(args.out_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        import subprocess

        subprocess.run(
            [
                str(binary),
                str(compiled_dir / "similarities.bin"),
                str(compiled_dir / "proteins.tsv"),
                str(compiled_dir / "taxa.tsv"),
                str(args.orthologs_file),
                str(output_path),
                str(args.percent_match_cutoff),
                str(args.evalue_exp_cutoff),
            ],
            check=True,
        )
        sys.stdout.write(f"indexed inparalog outputs written to {output_path}\n")
        return 0

    if args.command == "indexed-coorthologs":
        root = Path(__file__).resolve().parents[2]
        binary = root / "build" / "indexed_coorthologs"
        if not binary.exists():
            import subprocess

            subprocess.run(["make", "build-indexed-coorthologs"], cwd=root, check=True)
        compiled_dir = Path(args.compiled_dir)
        output_path = Path(args.out_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        import subprocess

        subprocess.run(
            [
                str(binary),
                str(compiled_dir / "similarities.bin"),
                str(compiled_dir / "proteins.tsv"),
                str(compiled_dir / "taxa.tsv"),
                str(args.orthologs_file),
                str(args.inparalogs_file),
                str(output_path),
                str(args.percent_match_cutoff),
                str(args.evalue_exp_cutoff),
            ],
            check=True,
        )
        sys.stdout.write(f"indexed coortholog outputs written to {output_path}\n")
        return 0

    if args.command == "indexed-pairs":
        summary = run_indexed_pairs(
            args.compiled_dir,
            args.out_dir,
            args.percent_match_cutoff,
            args.evalue_exp_cutoff,
            jobs=args.jobs,
        )
        sys.stdout.write(
            f"indexed pairs complete: orthologs={summary.ortholog_count}, "
            f"inparalogs={summary.inparalog_count}, "
            f"coorthologs={summary.coortholog_count}, "
            f"rbh_1to1={summary.rbh_count}, "
            f"mcl={summary.mcl_count}\n"
        )
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
