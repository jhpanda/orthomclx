from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, Optional, Tuple, Union

from orthomcl.compat import dataclass


@dataclass(slots=True)
class FastaRecord:
    header: str
    sequence: str


@dataclass(slots=True)
class FilterStats:
    total_sequences: int = 0
    rejected_sequences: int = 0


def read_fasta_records(path: Union[str, Path]) -> Iterator[FastaRecord]:
    file_path = Path(path)
    header: Optional[str] = None
    sequence_parts: list[str] = []

    with file_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header=header, sequence="".join(sequence_parts))
                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if header is not None:
        yield FastaRecord(header=header, sequence="".join(sequence_parts))


def write_fasta_record(handle, header: str, sequence: str, wrap: int = 60) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), wrap):
        handle.write(sequence[start : start + wrap] + "\n")


def normalize_header_fields(header: str) -> list[str]:
    normalized = header.strip()
    while "  " in normalized:
        normalized = normalized.replace("  ", " ")
    normalized = normalized.replace(" |", "|").replace("| ", "|")
    fields = [field for field in normalized.replace("|", " ").split(" ") if field]
    return fields


def adjust_fasta(taxon_code: str, input_path: Union[str, Path], id_field: int, output_path: Union[str, Path]) -> None:
    if id_field < 1:
        raise ValueError("id_field must be >= 1")

    seen_ids: set[str] = set()
    with Path(output_path).open("w") as out_handle:
        for record in read_fasta_records(input_path):
            fields = normalize_header_fields(record.header)
            if id_field > len(fields):
                raise ValueError(
                    f"Definition line '{record.header}' does not have field {id_field}"
                )
            protein_id = fields[id_field - 1]
            if protein_id in seen_ids:
                raise ValueError(
                    f"Fasta file '{input_path}' contains a duplicate id: {protein_id}"
                )
            seen_ids.add(protein_id)
            write_fasta_record(out_handle, f"{taxon_code}|{protein_id}", record.sequence)


def iter_compliant_fastas(input_dir: Union[str, Path]) -> Iterable[Path]:
    directory = Path(input_dir)
    for path in sorted(directory.iterdir()):
        if path.name.startswith("."):
            continue
        if path.suffix == ".fasta":
            yield path


def validate_compliant_header(header: str, expected_taxon: str) -> None:
    prefix = header.split("|", 1)[0]
    if prefix != expected_taxon:
        raise ValueError(
            f"The ID on def line '>{header}' is missing the prefix '{expected_taxon}|'"
        )


def sanitize_sequence_for_filter(sequence: str) -> Tuple[str, int, int]:
    raw_length = len(sequence)
    cleaned = "".join(char for char in sequence if char.isalpha())
    removed_count = raw_length - len(cleaned)
    return cleaned, raw_length, removed_count


def filter_fasta_dir(
    input_dir: Union[str, Path],
    min_length: int,
    max_percent_stop: float,
    good_output_path: Union[str, Path],
    poor_output_path: Union[str, Path],
) -> list[tuple[str, float]]:
    reject_rates: list[tuple[str, float]] = []

    with Path(good_output_path).open("w") as good_handle, Path(poor_output_path).open(
        "w"
    ) as poor_handle:
        fasta_files = list(iter_compliant_fastas(input_dir))
        if not fasta_files:
            raise ValueError(f"Input directory {input_dir} does not contain any .fasta files")

        for fasta_file in fasta_files:
            abbrev = fasta_file.stem
            stats = FilterStats()
            for record in read_fasta_records(fasta_file):
                validate_compliant_header(record.header, abbrev)
                cleaned_sequence, raw_length, removed_count = sanitize_sequence_for_filter(
                    record.sequence
                )
                if raw_length == 0:
                    raise ValueError(
                        f"Error: zero length sequence in file {fasta_file.name} for header '{record.header}'"
                    )

                stats.total_sequences += 1
                stop_percent = removed_count / raw_length * 100
                target_handle = (
                    poor_handle
                    if raw_length < min_length or stop_percent > max_percent_stop
                    else good_handle
                )
                if target_handle is poor_handle:
                    stats.rejected_sequences += 1
                write_fasta_record(target_handle, record.header, cleaned_sequence)

            if stats.total_sequences and stats.rejected_sequences:
                reject_percent = stats.rejected_sequences / stats.total_sequences * 100
                if reject_percent > 10:
                    reject_rates.append((fasta_file.name, reject_percent))

    return sorted(reject_rates, key=lambda item: item[1], reverse=True)
