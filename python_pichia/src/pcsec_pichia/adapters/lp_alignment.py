from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pcsec_pichia.adapters.lp_parser import LpFileSummary, parse_lp_file


DEFAULT_LP_ALIGNMENT_FIELDS: tuple[str, ...] = (
    "optimization_sense",
    "objective_variable",
    "constraint_count",
    "bounds_count",
    "distinct_variable_count",
    "max_variable_index",
    "constraint_prefix_counts",
    "constraint_sense_counts",
)


@dataclass(frozen=True)
class LpFieldDifference:
    reference: Any
    candidate: Any


@dataclass(frozen=True)
class LpAlignmentComparison:
    reference_path: Path
    candidate_path: Path
    compared_fields: tuple[str, ...]
    matching_fields: tuple[str, ...]
    differing_fields: dict[str, LpFieldDifference] = field(default_factory=dict)

    @property
    def is_match(self) -> bool:
        return not self.differing_fields

    def to_dict(self) -> dict[str, object]:
        return {
            "reference_path": str(self.reference_path),
            "candidate_path": str(self.candidate_path),
            "is_match": self.is_match,
            "compared_fields": list(self.compared_fields),
            "matching_fields": list(self.matching_fields),
            "differing_fields": {
                field_name: {
                    "reference": difference.reference,
                    "candidate": difference.candidate,
                }
                for field_name, difference in self.differing_fields.items()
            },
        }


@dataclass(frozen=True)
class LpConstraintMathDifference:
    name: str
    kind: str
    reference_sense: str | None = None
    candidate_sense: str | None = None
    reference_rhs: float | None = None
    candidate_rhs: float | None = None
    variable: str | None = None
    reference_coefficient: float | None = None
    candidate_coefficient: float | None = None
    absolute_difference: float | None = None
    reference_variable_count: int | None = None
    candidate_variable_count: int | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "kind": self.kind,
            "reference_sense": self.reference_sense,
            "candidate_sense": self.candidate_sense,
            "reference_rhs": self.reference_rhs,
            "candidate_rhs": self.candidate_rhs,
            "variable": self.variable,
            "reference_coefficient": self.reference_coefficient,
            "candidate_coefficient": self.candidate_coefficient,
            "absolute_difference": self.absolute_difference,
            "reference_variable_count": self.reference_variable_count,
            "candidate_variable_count": self.candidate_variable_count,
        }


@dataclass(frozen=True)
class LpBoundDifference:
    variable: str
    reference_lower: float | None
    reference_upper: float | None
    candidate_lower: float | None
    candidate_upper: float | None

    def to_dict(self) -> dict[str, object]:
        return {
            "variable": self.variable,
            "reference_lower": self.reference_lower,
            "reference_upper": self.reference_upper,
            "candidate_lower": self.candidate_lower,
            "candidate_upper": self.candidate_upper,
        }


@dataclass(frozen=True)
class LpMathAlignmentComparison:
    reference_path: Path
    candidate_path: Path
    coefficient_tolerance: float
    rhs_tolerance: float
    bound_tolerance: float
    constraint_count: int
    bound_count: int
    constraint_differences: tuple[LpConstraintMathDifference, ...] = ()
    bound_differences: tuple[LpBoundDifference, ...] = ()
    max_abs_coefficient_difference: float = 0.0
    max_abs_coefficient_difference_item: dict[str, object] | None = None

    @property
    def is_match(self) -> bool:
        return not self.constraint_differences and not self.bound_differences

    @property
    def constraint_difference_count(self) -> int:
        return len(self.constraint_differences)

    @property
    def bound_difference_count(self) -> int:
        return len(self.bound_differences)

    def to_dict(self) -> dict[str, object]:
        return {
            "reference_path": str(self.reference_path),
            "candidate_path": str(self.candidate_path),
            "is_match": self.is_match,
            "coefficient_tolerance": self.coefficient_tolerance,
            "rhs_tolerance": self.rhs_tolerance,
            "bound_tolerance": self.bound_tolerance,
            "constraint_count": self.constraint_count,
            "bound_count": self.bound_count,
            "constraint_difference_count": self.constraint_difference_count,
            "bound_difference_count": self.bound_difference_count,
            "max_abs_coefficient_difference": self.max_abs_coefficient_difference,
            "max_abs_coefficient_difference_item": self.max_abs_coefficient_difference_item,
            "constraint_differences": [difference.to_dict() for difference in self.constraint_differences],
            "bound_differences": [difference.to_dict() for difference in self.bound_differences],
        }


def compare_lp_files(
    reference_path: Path,
    candidate_path: Path,
    fields: tuple[str, ...] = DEFAULT_LP_ALIGNMENT_FIELDS,
) -> LpAlignmentComparison:
    return compare_lp_summaries(parse_lp_file(reference_path), parse_lp_file(candidate_path), fields=fields)


def compare_lp_summaries(
    reference: LpFileSummary,
    candidate: LpFileSummary,
    fields: tuple[str, ...] = DEFAULT_LP_ALIGNMENT_FIELDS,
) -> LpAlignmentComparison:
    matching: list[str] = []
    differing: dict[str, LpFieldDifference] = {}
    for field_name in fields:
        reference_value = getattr(reference, field_name)
        candidate_value = getattr(candidate, field_name)
        if reference_value == candidate_value:
            matching.append(field_name)
        else:
            differing[field_name] = LpFieldDifference(reference=reference_value, candidate=candidate_value)
    return LpAlignmentComparison(
        reference_path=reference.path,
        candidate_path=candidate.path,
        compared_fields=fields,
        matching_fields=tuple(matching),
        differing_fields=differing,
    )


def compare_lp_math(
    reference_path: Path,
    candidate_path: Path,
    coefficient_tolerance: float = 1e-12,
    rhs_tolerance: float = 1e-12,
    bound_tolerance: float = 1e-12,
) -> LpMathAlignmentComparison:
    reference = _parse_math_lp(reference_path)
    candidate = _parse_math_lp(candidate_path)

    constraint_differences: list[LpConstraintMathDifference] = []
    for name in _natural_sorted_names(set(reference.constraints) | set(candidate.constraints)):
        reference_constraint = reference.constraints.get(name)
        candidate_constraint = candidate.constraints.get(name)
        if reference_constraint is None:
            constraint_differences.append(
                LpConstraintMathDifference(
                    name=name,
                    kind="extra_constraint",
                    candidate_sense=candidate_constraint.sense,
                    candidate_rhs=candidate_constraint.rhs,
                    candidate_variable_count=len(candidate_constraint.coefficients),
                )
            )
            continue
        if candidate_constraint is None:
            constraint_differences.append(
                LpConstraintMathDifference(
                    name=name,
                    kind="missing_constraint",
                    reference_sense=reference_constraint.sense,
                    reference_rhs=reference_constraint.rhs,
                    reference_variable_count=len(reference_constraint.coefficients),
                )
            )
            continue
        if (
            reference_constraint.sense != candidate_constraint.sense
            or set(reference_constraint.coefficients) != set(candidate_constraint.coefficients)
            or not math.isclose(reference_constraint.rhs, candidate_constraint.rhs, rel_tol=0.0, abs_tol=rhs_tolerance)
        ):
            constraint_differences.append(
                LpConstraintMathDifference(
                    name=name,
                    kind="structure_or_rhs",
                    reference_sense=reference_constraint.sense,
                    candidate_sense=candidate_constraint.sense,
                    reference_rhs=reference_constraint.rhs,
                    candidate_rhs=candidate_constraint.rhs,
                    reference_variable_count=len(reference_constraint.coefficients),
                    candidate_variable_count=len(candidate_constraint.coefficients),
                )
            )
            continue
        largest_variable, largest_difference = _largest_coefficient_difference(
            reference_constraint.coefficients,
            candidate_constraint.coefficients,
        )
        if largest_variable is not None and largest_difference > coefficient_tolerance:
            constraint_differences.append(
                LpConstraintMathDifference(
                    name=name,
                    kind="coefficient",
                    variable=f"X{largest_variable}",
                    reference_coefficient=reference_constraint.coefficients[largest_variable],
                    candidate_coefficient=candidate_constraint.coefficients[largest_variable],
                    absolute_difference=largest_difference,
                    reference_variable_count=len(reference_constraint.coefficients),
                    candidate_variable_count=len(candidate_constraint.coefficients),
                )
            )

    bound_differences: list[LpBoundDifference] = []
    for variable in _sorted_variables(set(reference.bounds) | set(candidate.bounds)):
        reference_bound = reference.bounds.get(variable)
        candidate_bound = candidate.bounds.get(variable)
        if reference_bound is None:
            bound_differences.append(
                LpBoundDifference(
                    variable=f"X{variable}",
                    reference_lower=None,
                    reference_upper=None,
                    candidate_lower=candidate_bound[0],
                    candidate_upper=candidate_bound[1],
                )
            )
            continue
        if candidate_bound is None:
            bound_differences.append(
                LpBoundDifference(
                    variable=f"X{variable}",
                    reference_lower=reference_bound[0],
                    reference_upper=reference_bound[1],
                    candidate_lower=None,
                    candidate_upper=None,
                )
            )
            continue
        if not _bounds_close(reference_bound, candidate_bound, bound_tolerance):
            bound_differences.append(
                LpBoundDifference(
                    variable=f"X{variable}",
                    reference_lower=reference_bound[0],
                    reference_upper=reference_bound[1],
                    candidate_lower=candidate_bound[0],
                    candidate_upper=candidate_bound[1],
                )
            )

    max_difference = 0.0
    max_item: dict[str, object] | None = None
    for name in _natural_sorted_names(set(reference.constraints) & set(candidate.constraints)):
        reference_constraint = reference.constraints[name]
        candidate_constraint = candidate.constraints[name]
        if set(reference_constraint.coefficients) != set(candidate_constraint.coefficients):
            continue
        variable, difference = _largest_coefficient_difference(
            reference_constraint.coefficients,
            candidate_constraint.coefficients,
        )
        if variable is not None and difference > max_difference:
            max_difference = difference
            max_item = {
                "constraint": name,
                "variable": f"X{variable}",
                "reference_coefficient": reference_constraint.coefficients[variable],
                "candidate_coefficient": candidate_constraint.coefficients[variable],
            }

    return LpMathAlignmentComparison(
        reference_path=reference_path,
        candidate_path=candidate_path,
        coefficient_tolerance=coefficient_tolerance,
        rhs_tolerance=rhs_tolerance,
        bound_tolerance=bound_tolerance,
        constraint_count=len(reference.constraints),
        bound_count=len(reference.bounds),
        constraint_differences=tuple(constraint_differences),
        bound_differences=tuple(bound_differences),
        max_abs_coefficient_difference=max_difference,
        max_abs_coefficient_difference_item=max_item,
    )


def write_lp_math_alignment_json(comparison: LpMathAlignmentComparison, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(comparison.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path


@dataclass(frozen=True)
class _MathConstraint:
    sense: str
    coefficients: dict[int, float]
    rhs: float


@dataclass(frozen=True)
class _MathLpData:
    constraints: dict[str, _MathConstraint]
    bounds: dict[int, tuple[float, float]]


_TERM_RE = re.compile(r"([+-]?)\s*(?:((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*)?X(\d+)")
_BOUND_RE = re.compile(r"^\s*(\S+)\s*<=\s*X(\d+)\s*<=\s*(\S+)\s*$")


def _parse_math_lp(path: Path) -> _MathLpData:
    section = "header"
    constraints: dict[str, str] = {}
    bounds: dict[int, tuple[float, float]] = {}
    current_name: str | None = None
    current_lines: list[str] = []

    def flush_constraint() -> None:
        nonlocal current_name, current_lines
        if current_name is not None:
            constraints[current_name] = " ".join(line.strip() for line in current_lines)
            current_name = None
            current_lines = []

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped in {"Maximize", "Minimize"}:
                flush_constraint()
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
            if section == "constraints":
                match = re.match(r"^\s*([A-Za-z][A-Za-z0-9_]*):\s*(.*)$", line)
                if match:
                    flush_constraint()
                    current_name = match.group(1)
                    current_lines = [match.group(2)]
                elif current_name is not None:
                    current_lines.append(line)
            elif section == "bounds":
                match = _BOUND_RE.match(stripped)
                if match:
                    bounds[int(match.group(2))] = (_bound_value(match.group(1)), _bound_value(match.group(3)))

    flush_constraint()
    return _MathLpData(
        constraints={name: _parse_constraint_expression(expression) for name, expression in constraints.items()},
        bounds=bounds,
    )


def _parse_constraint_expression(expression: str) -> _MathConstraint:
    sense = ""
    lhs = ""
    rhs_text = ""
    for candidate_sense in ("<=", ">=", "="):
        if candidate_sense in expression:
            sense = candidate_sense
            lhs, rhs_text = expression.rsplit(candidate_sense, 1)
            break
    if not sense:
        raise ValueError(f"Could not parse LP constraint sense: {expression[:120]}")
    coefficients: dict[int, float] = {}
    for sign, coefficient_text, variable_text in _TERM_RE.findall(lhs):
        coefficient = float(coefficient_text) if coefficient_text else 1.0
        if sign == "-":
            coefficient = -coefficient
        variable = int(variable_text)
        coefficients[variable] = coefficients.get(variable, 0.0) + coefficient
    return _MathConstraint(sense=sense, coefficients=coefficients, rhs=float(rhs_text.strip()))


def _bound_value(text: str) -> float:
    normalized = text.strip().lower()
    if normalized in {"+infinity", "infinity", "+inf", "inf"}:
        return math.inf
    if normalized in {"-infinity", "-inf"}:
        return -math.inf
    return float(normalized)


def _largest_coefficient_difference(
    reference_coefficients: dict[int, float],
    candidate_coefficients: dict[int, float],
) -> tuple[int | None, float]:
    largest_variable: int | None = None
    largest_difference = 0.0
    for variable in reference_coefficients:
        difference = abs(reference_coefficients[variable] - candidate_coefficients[variable])
        if difference > largest_difference:
            largest_variable = variable
            largest_difference = difference
    return largest_variable, largest_difference


def _bounds_close(reference: tuple[float, float], candidate: tuple[float, float], tolerance: float) -> bool:
    return math.isclose(reference[0], candidate[0], rel_tol=0.0, abs_tol=tolerance) and math.isclose(
        reference[1],
        candidate[1],
        rel_tol=0.0,
        abs_tol=tolerance,
    )


def _natural_sorted_names(names: set[str]) -> list[str]:
    def key(name: str) -> tuple[str, int, str]:
        match = re.match(r"([A-Za-z_]+)(\d*)", name)
        if not match:
            return name, 0, name
        return match.group(1), int(match.group(2) or 0), name

    return sorted(names, key=key)


def _sorted_variables(variables: set[int]) -> list[int]:
    return sorted(variables)
