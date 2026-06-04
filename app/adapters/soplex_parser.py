from __future__ import annotations

import re
from pathlib import Path

from app.core.models import SoplexSummary


STATUS_PATTERN = re.compile(r"SoPlex status\s*:\s*(.+)")
OBJECTIVE_PATTERN = re.compile(r"Objective value\s*:?\s*(.+)")


def parse_soplex_output(text: str) -> SoplexSummary:
    status_matches = STATUS_PATTERN.findall(text)
    objective_matches = OBJECTIVE_PATTERN.findall(text)
    status_line = status_matches[-1].strip() if status_matches else None
    objective_value = objective_matches[-1].strip() if objective_matches else None
    return SoplexSummary(
        optimal=bool(status_line and "problem is solved [optimal]" in status_line),
        objective_value=objective_value,
        status_line=status_line,
    )


def parse_soplex_file(path: Path) -> SoplexSummary:
    return parse_soplex_output(path.read_text(encoding="utf-8", errors="replace"))
