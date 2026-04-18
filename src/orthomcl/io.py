from __future__ import annotations

import struct
from pathlib import Path
from typing import Iterable, Union

from orthomcl.models import EdgeRecord, SimilarityRecord


def read_similarity_records(path: Union[str, Path]) -> list[SimilarityRecord]:
    records: list[SimilarityRecord] = []
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 8:
                raise ValueError(
                    f"Similarity line {line_number} in '{path}' does not have 8 columns"
                )
            records.append(
                SimilarityRecord(
                    query_id=fields[0],
                    subject_id=fields[1],
                    query_taxon=fields[2],
                    subject_taxon=fields[3],
                    evalue_mant=float(fields[4]),
                    evalue_exp=int(fields[5]),
                    percent_identity=float(fields[6]),
                    percent_match=float(fields[7]),
                )
            )
    return records


def ensure_directory(path: Union[str, Path]) -> Path:
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def round_score(score: float) -> float:
    # Legacy OrthoMCL stores normalized scores in MySQL FLOAT columns before
    # dumping them, so cast through float32 before the final 3-decimal rounding.
    score32 = struct.unpack("<f", struct.pack("<f", float(score)))[0]
    return int(score32 * 1000 + 0.5) / 1000


def format_score(score: float) -> str:
    rounded = round_score(score)
    return f"{rounded:.3f}".rstrip("0").rstrip(".")


def write_edge_file(path: Union[str, Path], edges: Iterable[EdgeRecord]) -> None:
    with Path(path).open("w") as handle:
        for edge in edges:
            handle.write(f"{edge.seq_a}\t{edge.seq_b}\t{format_score(edge.normalized_score)}\n")


def write_mcl_input(path: Union[str, Path], edges: Iterable[EdgeRecord]) -> None:
    with Path(path).open("w") as handle:
        for edge in edges:
            handle.write(f"{edge.seq_a}\t{edge.seq_b}\t{format_score(edge.normalized_score)}\n")


def read_edge_file(path: Union[str, Path], edge_type: str, same_taxon: bool = False) -> list[EdgeRecord]:
    edges: list[EdgeRecord] = []
    with Path(path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            seq_a, seq_b, score = line.split("\t")
            edges.append(
                EdgeRecord(
                    seq_a=seq_a,
                    seq_b=seq_b,
                    taxon_a="",
                    taxon_b="" if not same_taxon else "",
                    unnormalized_score=float(score),
                    normalized_score=float(score),
                    edge_type=edge_type,
                )
            )
    return edges
