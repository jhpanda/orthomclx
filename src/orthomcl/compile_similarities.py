from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

from orthomcl.io import ensure_directory, read_similarity_records


RECORD_STRUCT = struct.Struct("<IIIIfiff")


@dataclass(slots=True)
class CompiledSimilaritySummary:
    record_count: int
    protein_count: int
    taxon_count: int
    binary_path: Path
    proteins_path: Path
    taxa_path: Path


def compile_similarities(
    similar_sequences_path: str | Path,
    out_dir: str | Path,
) -> CompiledSimilaritySummary:
    records = read_similarity_records(similar_sequences_path)
    min_nonzero_exp = min(
        (record.evalue_exp for record in records if record.evalue_mant != 0),
        default=None,
    )
    adjusted_zero_exp = (min_nonzero_exp - 1) if min_nonzero_exp is not None else 0
    output_dir = ensure_directory(out_dir)
    proteins_path = output_dir / "proteins.tsv"
    taxa_path = output_dir / "taxa.tsv"
    binary_path = output_dir / "similarities.bin"

    protein_index: dict[str, int] = {}
    taxon_index: dict[str, int] = {}

    def get_index(mapping: dict[str, int], value: str) -> int:
        existing = mapping.get(value)
        if existing is not None:
            return existing
        idx = len(mapping)
        mapping[value] = idx
        return idx

    with binary_path.open("wb") as binary_handle:
        for record in records:
            query_idx = get_index(protein_index, record.query_id)
            subject_idx = get_index(protein_index, record.subject_id)
            query_taxon_idx = get_index(taxon_index, record.query_taxon)
            subject_taxon_idx = get_index(taxon_index, record.subject_taxon)
            evalue_exp = int(record.evalue_exp)
            if float(record.evalue_mant) == 0.0 and evalue_exp == 0:
                evalue_exp = adjusted_zero_exp
            binary_handle.write(
                RECORD_STRUCT.pack(
                    query_idx,
                    subject_idx,
                    query_taxon_idx,
                    subject_taxon_idx,
                    float(record.evalue_mant),
                    evalue_exp,
                    float(record.percent_identity),
                    float(record.percent_match),
                )
            )

    write_index_file(proteins_path, protein_index)
    write_index_file(taxa_path, taxon_index)

    return CompiledSimilaritySummary(
        record_count=len(records),
        protein_count=len(protein_index),
        taxon_count=len(taxon_index),
        binary_path=binary_path,
        proteins_path=proteins_path,
        taxa_path=taxa_path,
    )


def write_index_file(path: Path, mapping: dict[str, int]) -> None:
    items = sorted(mapping.items(), key=lambda item: item[1])
    with path.open("w") as handle:
        for value, idx in items:
            handle.write(f"{idx}\t{value}\n")
