from __future__ import annotations

from pathlib import Path
from typing import Dict, Set, Tuple, Union

from orthomcl.compat import dataclass


@dataclass(frozen=True, slots=True)
class BlastHit:
    query_id: str
    subject_id: str
    evalue: float
    bits: float


def read_blast_m8_hits(path: Union[str, Path]) -> list[BlastHit]:
    hits: list[BlastHit] = []
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 12:
                raise ValueError(f"BLAST line {line_number} in '{path}' does not have 12 columns")
            hits.append(
                BlastHit(
                    query_id=fields[0],
                    subject_id=fields[1],
                    evalue=float(fields[10]),
                    bits=float(fields[11]),
                )
            )
    return hits


def reciprocal_best_hits(hits: list[BlastHit]) -> Set[Tuple[str, str]]:
    best_by_query: Dict[str, BlastHit] = {}
    for hit in hits:
        current = best_by_query.get(hit.query_id)
        if current is None or (hit.evalue, -hit.bits, hit.subject_id) < (
            current.evalue,
            -current.bits,
            current.subject_id,
        ):
            best_by_query[hit.query_id] = hit

    pairs: Set[Tuple[str, str]] = set()
    for query_id, best_hit in best_by_query.items():
        reverse = best_by_query.get(best_hit.subject_id)
        if reverse is None or reverse.subject_id != query_id:
            continue
        pairs.add(tuple(sorted((query_id, best_hit.subject_id))))
    return pairs


def read_rbh_pairs(path: Union[str, Path]) -> Set[Tuple[str, str]]:
    pairs: Set[Tuple[str, str]] = set()
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(f"RBH line {line_number} in '{path}' does not have 2 columns")
            pairs.add(tuple(sorted((fields[0], fields[1]))))
    return pairs
