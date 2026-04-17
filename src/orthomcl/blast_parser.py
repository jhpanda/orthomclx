from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(slots=True)
class GeneInfo:
    taxon: str
    length: int


@dataclass(slots=True)
class SubjectAccumulator:
    query_id: str
    subject_id: str
    query_taxon: str
    subject_taxon: str
    query_length: int
    subject_length: int
    query_shorter: bool
    evalue_mant: str
    evalue_exp: int
    hspspans: list[tuple[int, int]] = field(default_factory=list)
    total_identities: float = 0.0
    total_length: int = 0


def get_genes_from_fasta(fasta_dir: str | Path) -> dict[str, GeneInfo]:
    genes: dict[str, GeneInfo] = {}
    fasta_path = Path(fasta_dir)
    for file_path in sorted(fasta_path.iterdir()):
        if file_path.name.startswith("."):
            continue
        if not file_path.name.endswith(".fasta"):
            raise ValueError(f"'{file_path.name}' is not in 'taxon.fasta' format")
        taxon = file_path.stem
        gene_id: str | None = None
        length = 0
        with file_path.open() as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if gene_id is not None:
                        genes[gene_id] = GeneInfo(taxon=taxon, length=length)
                    gene_id = line[1:].split()[0]
                    length = 0
                else:
                    length += len(line)
        if gene_id is not None:
            genes[gene_id] = GeneInfo(taxon=taxon, length=length)
    return genes


def format_evalue(evalue: str) -> tuple[str, int]:
    normalized = f"1{evalue}" if evalue.startswith("e") else evalue
    mantissa, exponent = f"{float(normalized):.3e}".split("e")
    formatted_mantissa = f"{float(mantissa):.2f}".rstrip("0").rstrip(".")
    formatted_exponent = int(exponent)
    return formatted_mantissa, formatted_exponent


def get_start_end(hsp: tuple[int, int]) -> tuple[int, int]:
    start, end = hsp
    if start > end:
        return end, start
    return start, end


def compute_non_overlapping_match_length(hspspans: Iterable[tuple[int, int]]) -> int:
    sorted_hsps = sorted(hspspans, key=lambda hsp: get_start_end(hsp)[0])
    if not sorted_hsps:
        return 0

    start, end = get_start_end(sorted_hsps[0])
    length = 0
    for hsp in sorted_hsps[1:]:
        hsp_start, hsp_end = get_start_end(hsp)
        if hsp_end <= end:
            continue
        if hsp_start <= end:
            end = hsp_end
        else:
            length += end - start + 1
            start, end = hsp_start, hsp_end
    length += end - start + 1
    return length


def build_subject_accumulator(
    query_id: str,
    subject_id: str,
    evalue: str,
    genes: dict[str, GeneInfo],
) -> SubjectAccumulator:
    query_gene = genes.get(query_id)
    if query_gene is None:
        raise KeyError(f"couldn't find taxon for gene '{query_id}'")
    subject_gene = genes.get(subject_id)
    if subject_gene is None:
        raise KeyError(f"couldn't find taxon for gene '{subject_id}'")
    evalue_mant, evalue_exp = format_evalue(evalue)
    return SubjectAccumulator(
        query_id=query_id,
        subject_id=subject_id,
        query_taxon=query_gene.taxon,
        subject_taxon=subject_gene.taxon,
        query_length=query_gene.length,
        subject_length=subject_gene.length,
        query_shorter=query_gene.length < subject_gene.length,
        evalue_mant=evalue_mant,
        evalue_exp=evalue_exp,
    )


def subject_to_row(subject: SubjectAccumulator) -> str:
    non_overlap_match_len = compute_non_overlapping_match_length(subject.hspspans)
    percent_ident = round(subject.total_identities / subject.total_length, 1)
    shorter_length = subject.query_length if subject.query_shorter else subject.subject_length
    percent_match = round(non_overlap_match_len / shorter_length * 100, 1)
    return "\t".join(
        [
            subject.query_id,
            subject.subject_id,
            subject.query_taxon,
            subject.subject_taxon,
            subject.evalue_mant,
            str(subject.evalue_exp),
            f"{percent_ident:.1f}",
            f"{percent_match:.1f}",
        ]
    )


def parse_blast_m8(blast_path: str | Path, fasta_dir: str | Path, out_handle: TextIO) -> None:
    genes = get_genes_from_fasta(fasta_dir)
    current: SubjectAccumulator | None = None

    with Path(blast_path).open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) != 12:
                raise ValueError(
                    f"BLAST line {line_number} in '{blast_path}' does not have 12 columns"
                )
            (
                query_id,
                subject_id,
                percent_identity,
                length,
                _mismatches,
                _ngaps,
                query_start,
                query_end,
                subject_start,
                subject_end,
                evalue,
                _bits,
            ) = fields

            if (
                current is None
                or query_id != current.query_id
                or subject_id != current.subject_id
            ):
                if current is not None:
                    out_handle.write(subject_to_row(current) + "\n")
                current = build_subject_accumulator(query_id, subject_id, evalue, genes)

            span = (
                (int(query_start), int(query_end))
                if current.query_shorter
                else (int(subject_start), int(subject_end))
            )
            current.hspspans.append(span)
            current.total_identities += float(percent_identity) * int(length)
            current.total_length += int(length)

    if current is not None:
        out_handle.write(subject_to_row(current) + "\n")

