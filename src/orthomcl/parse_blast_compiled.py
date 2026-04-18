from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Union

from orthomcl.compile_similarities import CompiledSimilaritySummary, RECORD_STRUCT
from orthomcl.similarity_indexes import build_similarity_indexes


def ensure_c_blast_compiler_built() -> Path:
    root = Path(__file__).resolve().parents[2]
    binary = root / "build" / "parse_blast_compiled"
    source = root / "src" / "c" / "parse_blast_compiled.c"
    if binary.exists() and source.exists() and binary.stat().st_mtime >= source.stat().st_mtime:
        return binary
    subprocess.run(["make", "build-parse-blast-compiled"], cwd=root, check=True)
    return binary


def parse_blast_to_compiled(
    blast_path: Union[str, Path],
    fasta_dir: Union[str, Path],
    out_dir: Union[str, Path],
) -> CompiledSimilaritySummary:
    binary = ensure_c_blast_compiler_built()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[orthomclx] parsing BLAST to compiled binary: {blast_path}")
    subprocess.run([str(binary), str(blast_path), str(fasta_dir), str(out_dir)], check=True)
    build_similarity_indexes(out_dir, required=True)

    binary_path = out_dir / "similarities.bin"
    proteins_path = out_dir / "proteins.tsv"
    taxa_path = out_dir / "taxa.tsv"
    record_count = binary_path.stat().st_size // RECORD_STRUCT.size
    with proteins_path.open() as handle:
        protein_count = sum(1 for _ in handle)
    with taxa_path.open() as handle:
        taxon_count = sum(1 for _ in handle)
    return CompiledSimilaritySummary(
        record_count=record_count,
        protein_count=protein_count,
        taxon_count=taxon_count,
        binary_path=binary_path,
        proteins_path=proteins_path,
        taxa_path=taxa_path,
    )
