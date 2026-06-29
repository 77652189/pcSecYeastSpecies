from __future__ import annotations

from typing import Literal


CompatibilityMode = Literal["matlab_compat", "corrected"]
GlycosylationMode = Literal["native", "humanized"]
ResultStatus = Literal["draft", "matlab_aligned", "corrected_condition"]

DEFAULT_COMPATIBILITY_MODE: CompatibilityMode = "corrected"
DEFAULT_GLYCOSYLATION_MODE: GlycosylationMode = "native"
DEFAULT_RESULT_STATUS: ResultStatus = "draft"


def validate_compatibility_mode(value: str) -> CompatibilityMode:
    if value not in ("matlab_compat", "corrected"):
        raise ValueError(f"Unsupported compatibility mode: {value}")
    return value  # type: ignore[return-value]


def validate_glycosylation_mode(value: str) -> GlycosylationMode:
    if value not in ("native", "humanized"):
        raise ValueError(f"Unsupported glycosylation mode: {value}")
    return value  # type: ignore[return-value]
