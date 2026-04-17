from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from orthomcl.blast_parser import parse_blast_m8
from orthomcl.groups import mcl_to_groups_file
from orthomcl.compile_similarities import compile_similarities
from orthomcl.indexed_pairs import IndexedPairsSummary, run_indexed_pairs
from orthomcl.mcl import MclRunConfig, run_mcl
from orthomcl.pairs import build_pairs


@dataclass(slots=True)
class RunConfig:
    blast_file: Path
    fasta_dir: Path
    out_dir: Path
    percent_match_cutoff: float
    evalue_exp_cutoff: int
    engine: str = "python"
    mcl_output: Path | None = None
    run_mcl: bool = False
    mcl_binary: str = "mcl"
    mcl_inflation: float = 1.5
    mcl_threads: int | None = None
    groups_prefix: str = "OG"
    start_group_id: int = 1000


def run_pipeline(config: RunConfig) -> None:
    config.out_dir.mkdir(parents=True, exist_ok=True)

    similar_sequences_path = config.out_dir / "similarSequences.txt"
    with similar_sequences_path.open("w") as handle:
        parse_blast_m8(config.blast_file, config.fasta_dir, handle)

    build_pairs(
        similar_sequences_path,
        config.out_dir,
        config.percent_match_cutoff,
        config.evalue_exp_cutoff,
        engine=config.engine,
    )

    generated_mcl_output = config.out_dir / "mclOutput" if config.run_mcl else None
    if generated_mcl_output is not None:
        run_mcl(
            MclRunConfig(
                input_path=config.out_dir / "mclInput",
                output_path=generated_mcl_output,
                inflation=config.mcl_inflation,
                threads=config.mcl_threads,
                binary=config.mcl_binary,
            )
        )

    groups_input = generated_mcl_output or config.mcl_output
    if groups_input is not None:
        mcl_to_groups_file(
            groups_input,
            config.out_dir / "groups.txt",
            config.groups_prefix,
            config.start_group_id,
        )


@dataclass(slots=True)
class IntegratedRunConfig:
    blast_file: Path
    fasta_dir: Path
    out_dir: Path
    percent_match_cutoff: float
    evalue_exp_cutoff: int
    jobs: int = 1
    mcl_output: Path | None = None
    run_mcl: bool = False
    mcl_binary: str = "mcl"
    mcl_inflation: float = 1.5
    mcl_threads: int | None = None
    groups_prefix: str = "OG"
    start_group_id: int = 1000


def run_integrated_pipeline(config: IntegratedRunConfig) -> IndexedPairsSummary:
    config.out_dir.mkdir(parents=True, exist_ok=True)

    similar_sequences_path = config.out_dir / "similarSequences.txt"
    with similar_sequences_path.open("w") as handle:
        parse_blast_m8(config.blast_file, config.fasta_dir, handle)

    compiled_dir = config.out_dir / "compiled"
    compile_similarities(similar_sequences_path, compiled_dir)

    summary = run_indexed_pairs(
        compiled_dir,
        config.out_dir,
        config.percent_match_cutoff,
        config.evalue_exp_cutoff,
        jobs=config.jobs,
    )

    generated_mcl_output = config.out_dir / "mclOutput" if config.run_mcl else None
    if generated_mcl_output is not None:
        run_mcl(
            MclRunConfig(
                input_path=config.out_dir / "mclInput",
                output_path=generated_mcl_output,
                inflation=config.mcl_inflation,
                threads=config.mcl_threads,
                binary=config.mcl_binary,
            )
        )

    groups_input = generated_mcl_output or config.mcl_output
    if groups_input is not None:
        mcl_to_groups_file(
            groups_input,
            config.out_dir / "groups.txt",
            config.groups_prefix,
            config.start_group_id,
        )
    return summary
