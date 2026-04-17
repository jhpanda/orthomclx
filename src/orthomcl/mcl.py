from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class MclRunConfig:
    input_path: Path
    output_path: Path
    inflation: float = 1.5
    threads: int | None = None
    binary: str = "mcl"


def resolve_mcl_binary(binary: str) -> str:
    if binary != "mcl":
        resolved = shutil.which(binary)
        if resolved is None:
            raise FileNotFoundError(
                f"Could not find MCL binary '{binary}' on PATH or as an executable path."
            )
        return resolved

    root = Path(__file__).resolve().parents[2]
    bundled = root / "build" / "mcl"
    if bundled.exists() and bundled.is_file():
        return str(bundled)

    resolved = shutil.which(binary)
    if resolved is None:
        raise FileNotFoundError(
            "Could not find MCL binary 'mcl' on PATH and no bundled build/mcl was found. "
            "Install MCL, vendor it under src/mcl and build it, or pass --mcl-binary."
        )
    return resolved


def run_mcl(config: MclRunConfig) -> Path:
    binary = resolve_mcl_binary(config.binary)

    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        binary,
        str(config.input_path),
        "--abc",
        "-I",
        str(config.inflation),
        "-o",
        str(config.output_path),
    ]
    if config.threads is not None and config.threads > 0:
        command.extend(["-te", str(config.threads)])

    subprocess.run(command, check=True)
    return config.output_path
