from __future__ import annotations

import json
import re
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


STATUS_RE = re.compile(r"SoPlex status\s*:\s*(.+)")
OBJECTIVE_RE = re.compile(r"Objective value\s*:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")
SOLVING_TIME_RE = re.compile(r"Solving time \(sec\)\s*:\s*([-+]?\d+(?:\.\d+)?)")
ITERATIONS_RE = re.compile(r"Iterations\s*:\s*(\d+)")
SOLUTION_TYPE_RE = re.compile(r"Solution \(([^)]+)\)\s*:")
CONDITION_NUMBER_RE = re.compile(r"Condition Number\s*:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")
MAX_BOUND_VIOLATION_RE = re.compile(r"Max\. bound violation\s*=\s*([^\s]+)")
MAX_ROW_VIOLATION_RE = re.compile(r"Max\. row violation\s*=\s*([^\s]+)")
TERMINATION_DESPITE_VIOLATIONS_RE = re.compile(r"termination despite violations", re.IGNORECASE)


@dataclass(frozen=True)
class SoplexOutputSummary:
    path: Path
    status: str | None
    is_optimal: bool
    objective_value: float | None
    objective_text: str | None = None
    solving_time_seconds: float | None = None
    iterations: int | None = None
    solution_type: str | None = None
    condition_number: float | None = None
    max_bound_violation: float | None = None
    max_row_violation: float | None = None
    termination_despite_violations: bool = False
    diagnostic: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "path": str(self.path),
            "status": self.status,
            "is_optimal": self.is_optimal,
            "objective_value": self.objective_value,
            "objective_text": self.objective_text,
            "solving_time_seconds": self.solving_time_seconds,
            "iterations": self.iterations,
            "solution_type": self.solution_type,
            "condition_number": self.condition_number,
            "max_bound_violation": self.max_bound_violation,
            "max_row_violation": self.max_row_violation,
            "termination_despite_violations": self.termination_despite_violations,
            "diagnostic": self.diagnostic,
        }


@dataclass(frozen=True)
class SoplexObjectiveComparison:
    reference_path: Path
    candidate_path: Path
    reference_objective: float | None
    candidate_objective: float | None
    absolute_difference: float | None
    tolerance: float
    both_optimal: bool
    within_tolerance: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "reference_path": str(self.reference_path),
            "candidate_path": str(self.candidate_path),
            "reference_objective": self.reference_objective,
            "candidate_objective": self.candidate_objective,
            "absolute_difference": self.absolute_difference,
            "tolerance": self.tolerance,
            "both_optimal": self.both_optimal,
            "within_tolerance": self.within_tolerance,
        }


def parse_soplex_output(path: Path) -> SoplexOutputSummary:
    text = path.read_text(encoding="utf-8", errors="replace")
    return parse_soplex_text(text, path=path)


def parse_soplex_text(text: str, path: Path | None = None) -> SoplexOutputSummary:
    statuses = STATUS_RE.findall(text)
    objectives = OBJECTIVE_RE.findall(text)
    solving_times = SOLVING_TIME_RE.findall(text)
    iterations = ITERATIONS_RE.findall(text)
    solution_types = SOLUTION_TYPE_RE.findall(text)
    condition_numbers = CONDITION_NUMBER_RE.findall(text)
    max_bound_violations = MAX_BOUND_VIOLATION_RE.findall(text)
    max_row_violations = MAX_ROW_VIOLATION_RE.findall(text)
    status = statuses[-1].strip() if statuses else None
    objective_text = objectives[-1] if objectives else None
    objective_value = float(objectives[-1]) if objectives else None
    solution_type = solution_types[-1].strip().lower() if solution_types else None
    termination_despite_violations = bool(TERMINATION_DESPITE_VIOLATIONS_RE.search(text))
    is_optimal = bool(status and "problem is solved [optimal]" in status)
    return SoplexOutputSummary(
        path=path or Path("<memory>"),
        status=status,
        is_optimal=is_optimal,
        objective_value=objective_value,
        objective_text=objective_text,
        solving_time_seconds=float(solving_times[-1]) if solving_times else None,
        iterations=int(iterations[-1]) if iterations else None,
        solution_type=solution_type,
        condition_number=float(condition_numbers[-1]) if condition_numbers else None,
        max_bound_violation=_parse_violation_value(max_bound_violations[-1]) if max_bound_violations else None,
        max_row_violation=_parse_violation_value(max_row_violations[-1]) if max_row_violations else None,
        termination_despite_violations=termination_despite_violations,
        diagnostic=_diagnostic(status, is_optimal, solution_type, termination_despite_violations),
    )


def compare_soplex_objectives(
    reference_path: Path,
    candidate_path: Path,
    tolerance: float = 1e-8,
) -> SoplexObjectiveComparison:
    reference = parse_soplex_output(reference_path)
    candidate = parse_soplex_output(candidate_path)
    difference = None
    within_tolerance = False
    if reference.objective_value is not None and candidate.objective_value is not None:
        difference = abs(reference.objective_value - candidate.objective_value)
        within_tolerance = difference <= tolerance
    return SoplexObjectiveComparison(
        reference_path=reference.path,
        candidate_path=candidate.path,
        reference_objective=reference.objective_value,
        candidate_objective=candidate.objective_value,
        absolute_difference=difference,
        tolerance=tolerance,
        both_optimal=reference.is_optimal and candidate.is_optimal,
        within_tolerance=within_tolerance,
    )


def write_soplex_summary_json(summary: SoplexOutputSummary | SoplexObjectiveComparison, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(summary.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path


def _parse_violation_value(text: str) -> float | None:
    if text == "-":
        return None
    if "/" in text:
        return float(Fraction(text))
    return float(text)


def _diagnostic(
    status: str | None,
    is_optimal: bool,
    solution_type: str | None,
    termination_despite_violations: bool,
) -> str:
    if is_optimal:
        return "optimal"
    if status and "optimal with unscaled violations" in status:
        return "optimal_with_unscaled_violations"
    if status and status.startswith("error"):
        if solution_type == "rational" or termination_despite_violations:
            return "rational_numerical_difficulty"
        return "solver_error"
    if status:
        return "non_optimal"
    return "missing_status"
