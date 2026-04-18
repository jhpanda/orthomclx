from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Union


def ensure_similarity_index_builder_built() -> Path:
    root = Path(__file__).resolve().parents[2]
    binary = root / "build" / "build_similarity_indexes"
    source = root / "src" / "c" / "build_similarity_indexes.c"
    if binary.exists() and source.exists() and binary.stat().st_mtime >= source.stat().st_mtime:
        return binary
    subprocess.run(["make", "build-similarity-indexes"], cwd=root, check=True)
    return binary


def build_similarity_indexes(compiled_dir: Union[str, Path], *, required: bool = False) -> None:
    compiled_dir = Path(compiled_dir)
    if compiled_dir.joinpath("refs.query_taxon_best.bin").exists() and compiled_dir.joinpath("refs.query_subject.bin").exists() and compiled_dir.joinpath("refs.query_evalue.bin").exists():
        return
    if shutil.which("make") is None:
        if required:
            raise RuntimeError("make is required to build reusable similarity indexes")
        return
    try:
        binary = ensure_similarity_index_builder_built()
        print(f"[orthomclx] building reusable similarity indexes in {compiled_dir}")
        subprocess.run([str(binary), str(compiled_dir / "similarities.bin"), str(compiled_dir)], check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        if required:
            raise
