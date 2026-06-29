from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


LpSection = Literal["header", "objective", "constraints", "bounds", "end"]

CONSTRAINT_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*):\s*(.*)$")
VARIABLE_RE = re.compile(r"\bX(\d+)\b")
OBJECTIVE_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*):\s*(.*)$")


@dataclass(frozen=True)
class LpConstraintSummary:
    name: str
    prefix: str
    sense: str
    variable_count: int
    max_variable_index: int | None = None


@dataclass(frozen=True)
class LpFileSummary:
    path: Path
    optimization_sense: str
    objective_name: str
    objective_variable: str | None
    constraint_count: int
    bounds_count: int
    distinct_variable_count: int
    max_variable_index: int | None
    constraint_prefix_counts: dict[str, int]
    constraint_sense_counts: dict[str, int]
    constraints: list[LpConstraintSummary] = field(default_factory=list)

    def to_dict(self) -> dict[str, object]:
        return {
            "path": str(self.path),
            "optimization_sense": self.optimization_sense,
            "objective_name": self.objective_name,
            "objective_variable": self.objective_variable,
            "constraint_count": self.constraint_count,
            "bounds_count": self.bounds_count,
            "distinct_variable_count": self.distinct_variable_count,
            "max_variable_index": self.max_variable_index,
            "constraint_prefix_counts": self.constraint_prefix_counts,
            "constraint_sense_counts": self.constraint_sense_counts,
        }


def parse_lp_file(path: Path) -> LpFileSummary:
    optimization_sense = ""
    objective_name = ""
    objective_variable: str | None = None
    bounds_count = 0
    constraints: list[LpConstraintSummary] = []
    all_variables: set[int] = set()
    current_constraint_name: str | None = None
    current_constraint_lines: list[str] = []
    section: LpSection = "header"

    def flush_constraint() -> None:
        nonlocal current_constraint_name, current_constraint_lines
        if current_constraint_name is None:
            return
        text = " ".join(current_constraint_lines)
        variables = [int(value) for value in VARIABLE_RE.findall(text)]
        all_variables.update(variables)
        constraints.append(
            LpConstraintSummary(
                name=current_constraint_name,
                prefix=_constraint_prefix(current_constraint_name),
                sense=_constraint_sense(text),
                variable_count=len(set(variables)),
                max_variable_index=max(variables) if variables else None,
            )
        )
        current_constraint_name = None
        current_constraint_lines = []

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            if stripped in {"Maximize", "Minimize"}:
                flush_constraint()
                optimization_sense = stripped
                section = "objective"
                continue
            if stripped == "Subject To":
                flush_constraint()
                section = "constraints"
                continue
            if stripped == "Bounds":
                flush_constraint()
                section = "bounds"
                continue
            if stripped == "End":
                flush_constraint()
                section = "end"
                continue

            if section == "objective":
                match = OBJECTIVE_RE.match(line)
                if match:
                    objective_name = match.group(1)
                    objective_text = match.group(2)
                    objective_match = VARIABLE_RE.search(objective_text)
                    objective_variable = f"X{objective_match.group(1)}" if objective_match else None
                    if objective_match:
                        all_variables.add(int(objective_match.group(1)))
                continue

            if section == "constraints":
                match = CONSTRAINT_RE.match(line)
                if match:
                    flush_constraint()
                    current_constraint_name = match.group(1)
                    current_constraint_lines = [match.group(2)]
                elif current_constraint_name is not None:
                    current_constraint_lines.append(line)
                continue

            if section == "bounds":
                variables = [int(value) for value in VARIABLE_RE.findall(line)]
                if variables:
                    bounds_count += 1
                    all_variables.update(variables)

    flush_constraint()
    prefix_counts = Counter(item.prefix for item in constraints)
    sense_counts = Counter(item.sense for item in constraints)
    return LpFileSummary(
        path=path,
        optimization_sense=optimization_sense,
        objective_name=objective_name,
        objective_variable=objective_variable,
        constraint_count=len(constraints),
        bounds_count=bounds_count,
        distinct_variable_count=len(all_variables),
        max_variable_index=max(all_variables) if all_variables else None,
        constraint_prefix_counts=dict(prefix_counts),
        constraint_sense_counts=dict(sense_counts),
        constraints=constraints,
    )


def write_lp_summary_json(summary: LpFileSummary, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(summary.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path


def _constraint_prefix(name: str) -> str:
    if name == "Cmito":
        return "Cmito"
    match = re.match(r"([A-Za-z_]+)", name)
    return match.group(1) if match else name


def _constraint_sense(text: str) -> str:
    if "<=" in text:
        return "<="
    if ">=" in text:
        return ">="
    if "=" in text:
        return "="
    return "unknown"
