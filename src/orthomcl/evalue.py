from __future__ import annotations

from typing import Tuple


def parse_evalue_cutoff(value: float) -> Tuple[float, int]:
    if value < 0:
        return 1.0, int(value)
    if value == 0:
        raise ValueError("e-value cutoff must be non-zero")
    mantissa_text, exponent_text = f"{float(value):.3e}".split("e")
    mantissa = float(mantissa_text)
    mantissa = float(f"{mantissa:.2f}".rstrip("0").rstrip("."))
    exponent = int(exponent_text)
    return mantissa, exponent


def evalue_passes(record_mant: float, record_exp: int, cutoff_mant: float, cutoff_exp: int) -> bool:
    if record_mant == 0.0:
        return True
    if record_exp != cutoff_exp:
        return record_exp < cutoff_exp
    return record_mant <= cutoff_mant
