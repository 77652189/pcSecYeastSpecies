from __future__ import annotations

from pathlib import Path

from pcsec_pichia.adapters.soplex_output import parse_soplex_text
from app.core.models import SoplexSummary


def parse_soplex_output(text: str) -> SoplexSummary:
    parsed = parse_soplex_text(text)
    return SoplexSummary(
        optimal=parsed.is_optimal,
        objective_value=parsed.objective_text,
        status_line=parsed.status,
        solution_type=parsed.solution_type,
        diagnostic=parsed.diagnostic,
        condition_number=parsed.condition_number,
        max_bound_violation=parsed.max_bound_violation,
        max_row_violation=parsed.max_row_violation,
        termination_despite_violations=parsed.termination_despite_violations,
    )


def parse_soplex_file(path: Path) -> SoplexSummary:
    return parse_soplex_output(path.read_text(encoding="utf-8", errors="replace"))
