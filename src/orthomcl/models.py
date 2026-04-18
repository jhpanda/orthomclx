from __future__ import annotations

from orthomcl.compat import dataclass


@dataclass(frozen=True, slots=True)
class SimilarityRecord:
    query_id: str
    subject_id: str
    query_taxon: str
    subject_taxon: str
    evalue_mant: float
    evalue_exp: int
    percent_identity: float
    percent_match: float


@dataclass(frozen=True, slots=True)
class EdgeRecord:
    seq_a: str
    seq_b: str
    taxon_a: str
    taxon_b: str
    unnormalized_score: float
    normalized_score: float
    edge_type: str
