from __future__ import annotations

from pathlib import Path
from typing import TextIO, Union


def mcl_to_groups_stream(input_handle: TextIO, output_handle: TextIO, prefix: str, start_id: int) -> None:
    current_id = start_id
    for raw_line in input_handle:
        line = raw_line.rstrip("\n")
        normalized = line.replace("\t", " ")
        output_handle.write(f"{prefix}{current_id}: {normalized}\n")
        current_id += 1


def mcl_to_groups_file(input_path: Union[str, Path], output_path: Union[str, Path], prefix: str, start_id: int) -> None:
    with Path(input_path).open() as input_handle, Path(output_path).open("w") as output_handle:
        mcl_to_groups_stream(input_handle, output_handle, prefix, start_id)
