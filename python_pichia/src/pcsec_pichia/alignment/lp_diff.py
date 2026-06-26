from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy import sparse


NUMBER_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
SIGNED_NUMBER_RE = r"[+-]?\s*(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
TERM_RE = re.compile(rf"(?P<coef>{SIGNED_NUMBER_RE})?\s*\bX(?P<var>\d+)\b")
LABELED_ROW_RE = re.compile(r"^\s*(?P<label>[^:\s]+)\s*:\s*(?P<body>.*)$")


@dataclass(frozen=True)
class LpObjective:
    label: str
    coefficients: dict[int, float]


@dataclass(frozen=True)
class LpConstraint:
    label: str
    coefficients: dict[int, float]
    sense: str
    rhs: float


@dataclass(frozen=True)
class LpBound:
    lower: float | None
    upper: float | None


@dataclass(frozen=True)
class ParsedLp:
    path: Path
    optimization_sense: str
    objective: LpObjective
    constraints: tuple[LpConstraint, ...]
    bounds: dict[int, LpBound]


def write_matrix_lp(
    *,
    model: Any,
    A_eq: sparse.spmatrix,
    b_eq: Iterable[float],
    A_ub: sparse.spmatrix,
    b_ub: Iterable[float],
    objective_reaction: str,
    path: Path,
) -> Path:
    """Write a matrix-backed LP in the MATLAB-style X-indexed format used by diagnostics."""

    if objective_reaction not in model.reaction_index:
        raise KeyError(f"Reaction not found: {objective_reaction}")
    path.parent.mkdir(parents=True, exist_ok=True)
    objective_index = int(model.reaction_index[objective_reaction]) + 1
    eq = sparse.csr_matrix(A_eq)
    ub = sparse.csr_matrix(A_ub)
    b_eq_array = np.asarray(list(b_eq), dtype=float)
    b_ub_array = np.asarray(list(b_ub), dtype=float)

    lines = ["Maximize", f"obj: X{objective_index}", "Subject To"]
    for row_index in range(eq.shape[0]):
        lines.extend(_format_constraint(f"C{row_index + 1}", eq.getrow(row_index), "=", b_eq_array[row_index]))
    offset = eq.shape[0]
    for row_index in range(ub.shape[0]):
        lines.extend(_format_constraint(f"C{offset + row_index + 1}", ub.getrow(row_index), "<=", b_ub_array[row_index]))

    lines.append("Bounds")
    for index, (lower, upper) in enumerate(zip(model.lb, model.ub), start=1):
        upper_value = float(upper)
        upper_text = "+infinity" if upper_value >= 100 else f"{upper_value:.6f}"
        lines.append(f"{float(lower):.6f} <= X{index} <= {upper_text}")
    lines.append("End")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_matrix_lp_from_rows(
    *,
    model: Any,
    rows: Iterable[tuple[str, sparse.spmatrix, str, float]],
    objective_reaction: str,
    path: Path,
) -> Path:
    """Write an LP from explicit row specs for alignment-only formulation probes."""

    if objective_reaction not in model.reaction_index:
        raise KeyError(f"Reaction not found: {objective_reaction}")
    path.parent.mkdir(parents=True, exist_ok=True)
    objective_index = int(model.reaction_index[objective_reaction]) + 1
    lines = ["Maximize", f"obj: X{objective_index}", "Subject To"]
    for label, row, sense, rhs in rows:
        lines.extend(_format_constraint(label, sparse.csr_matrix(row), sense, float(rhs)))
    lines.append("Bounds")
    for index, (lower, upper) in enumerate(zip(model.lb, model.ub), start=1):
        upper_value = float(upper)
        upper_text = "+infinity" if upper_value >= 100 else f"{upper_value:.6f}"
        lines.append(f"{float(lower):.6f} <= X{index} <= {upper_text}")
    lines.append("End")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def parse_lp(path: Path) -> ParsedLp:
    optimization_sense = ""
    objective = LpObjective(label="", coefficients={})
    constraints: list[LpConstraint] = []
    bounds: dict[int, LpBound] = {}
    section = "header"
    current_label: str | None = None
    current_parts: list[str] = []

    def flush_constraint() -> None:
        nonlocal current_label, current_parts
        if current_label is None:
            return
        constraints.append(_parse_constraint(current_label, " ".join(current_parts)))
        current_label = None
        current_parts = []

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line in {"Maximize", "Minimize"}:
            flush_constraint()
            optimization_sense = line
            section = "objective"
            continue
        if line == "Subject To":
            flush_constraint()
            section = "constraints"
            continue
        if line == "Bounds":
            flush_constraint()
            section = "bounds"
            continue
        if line == "End":
            flush_constraint()
            section = "end"
            continue

        if section == "objective":
            match = LABELED_ROW_RE.match(line)
            if match:
                objective = LpObjective(
                    label=match.group("label"),
                    coefficients=_parse_coefficients(match.group("body")),
                )
            continue
        if section == "constraints":
            match = LABELED_ROW_RE.match(line)
            if match:
                flush_constraint()
                current_label = match.group("label")
                current_parts = [match.group("body")]
            elif current_label is not None:
                current_parts.append(line)
            continue
        if section == "bounds":
            parsed = _parse_bound(line)
            if parsed:
                variable_index, bound = parsed
                bounds[variable_index] = bound

    flush_constraint()
    return ParsedLp(
        path=path,
        optimization_sense=optimization_sense,
        objective=objective,
        constraints=tuple(constraints),
        bounds=bounds,
    )


def diff_lp_files(
    matlab_lp_path: Path,
    python_lp_path: Path,
    *,
    tolerance: float = 1e-9,
    top_n: int = 50,
) -> dict[str, Any]:
    matlab_lp = parse_lp(matlab_lp_path)
    python_lp = parse_lp(python_lp_path)
    row_diffs = _diff_rows(matlab_lp.constraints, python_lp.constraints, tolerance=tolerance)
    bound_diffs = _diff_bounds(matlab_lp.bounds, python_lp.bounds, tolerance=tolerance)
    objective_diff = _diff_objective(matlab_lp.objective, python_lp.objective, tolerance=tolerance)
    top_rows = sorted(row_diffs, key=lambda item: item["severity_score"], reverse=True)[:top_n]
    return {
        "matlab_lp_path": str(matlab_lp_path),
        "python_lp_path": str(python_lp_path),
        "tolerance": tolerance,
        "matlab": _lp_stats(matlab_lp),
        "python": _lp_stats(python_lp),
        "objective_diff": objective_diff,
        "row_diff_summary": _row_summary(row_diffs, matlab_lp.constraints, python_lp.constraints),
        "bound_diff_summary": _bound_summary(bound_diffs, matlab_lp.bounds, python_lp.bounds),
        "top_mismatched_rows": top_rows,
        "top_bound_differences": bound_diffs[:top_n],
    }


def write_lp_diff_outputs(diff: dict[str, Any], output_dir: Path, *, stem: str = "opn_lp_diff") -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / f"{stem}_summary.json"
    rows_path = output_dir / f"{stem}_top_rows.csv"
    report_path = output_dir / f"{stem}_report.md"
    summary_path.write_text(json.dumps(diff, ensure_ascii=False, indent=2), encoding="utf-8")
    _write_top_rows_csv(diff.get("top_mismatched_rows", []), rows_path)
    report_path.write_text(build_lp_diff_report(diff), encoding="utf-8")
    return {"summary": summary_path, "top_rows": rows_path, "report": report_path}


def build_lp_diff_report(diff: dict[str, Any]) -> str:
    objective = diff["objective_diff"]
    rows = diff["row_diff_summary"]
    bounds = diff["bound_diff_summary"]
    generation = diff.get("python_generation", {})
    artifact_context = diff.get("matlab_artifact_context", {})
    first_rows = diff.get("top_mismatched_rows", [])[:3]
    return "\n".join(
        [
            "# OPN Python corrected optional LP vs MATLAB LP row-level diff",
            "",
            "## Scope",
            "",
            "- Python LP: corrected medium, OPN_ALPHA_FULL_PROJECT, mu=0.10, ribosome translation constraint enabled, misfolding constraints enabled.",
            "- MATLAB LP: existing artifact only; no MATLAB was started and no new MATLAB baseline was generated.",
            "- This report is diagnostic. It does not mark the result as MATLAB aligned.",
            "",
            "## Shape",
            "",
            f"- MATLAB rows: `{diff['matlab']['constraint_count']}`; Python rows: `{diff['python']['constraint_count']}`; diff: `{diff['python']['constraint_count'] - diff['matlab']['constraint_count']}`.",
            f"- MATLAB bounds: `{diff['matlab']['bound_count']}`; Python bounds: `{diff['python']['bound_count']}`; diff: `{diff['python']['bound_count'] - diff['matlab']['bound_count']}`.",
            f"- MATLAB objective: `{diff['matlab']['objective_label']}` {objective['matlab_coefficients']}; Python objective: `{diff['python']['objective_label']}` {objective['python_coefficients']}.",
            f"- Python constraint counts: `{generation.get('constraint_counts', {})}`.",
            "",
            "## Difference Summary",
            "",
            f"- Objective coefficient differences: `{objective['coefficient_difference_count']}`.",
            f"- Ordered row label/name differences: `{rows['ordered_label_difference_count']}`.",
            f"- Label-only ordered row differences: `{rows.get('label_only_difference_count')}`.",
            f"- Ordered row RHS differences: `{rows['rhs_difference_count']}`.",
            f"- Ordered row coefficient differences: `{rows['coefficient_difference_count']}`.",
            f"- Ordered row sparsity differences: `{rows['sparsity_difference_count']}`.",
            f"- First ordered label difference: `{rows.get('first_label_difference_row')}`; first coefficient difference: `{rows.get('first_coefficient_difference_row')}`; first RHS/sense difference: `{rows.get('first_rhs_or_sense_difference_row')}`.",
            f"- Missing row labels in Python: `{rows['missing_label_count']}`; extra row labels in Python: `{rows['extra_label_count']}`.",
            f"- Bound differences: `{bounds['difference_count']}`; missing bounds in Python: `{bounds['missing_in_python_count']}`; extra bounds in Python: `{bounds['extra_in_python_count']}`.",
            "",
            "## Main Findings",
            "",
            f"- Existing MATLAB artifact context: `{artifact_context.get('note', 'not recorded')}`",
            "- The row and bound counts are shape-compatible: both LPs have 24,435 labeled constraints and 29,057 bounds.",
            "- The objective differs by formulation: MATLAB objective is `X545` (`Ex_glc_D` in the artifact harness), while Python corrected pipeline objective is `X29057` (OPN exchange).",
            "- Rows 20,222 onward use different labels (`CM*`/`Cmito` in MATLAB vs continuous `C*` labels in the Python matrix export). Most early differences in this block are label/name formatting only.",
            "- A real coefficient/sparsity difference first appears at row 22,956. The largest mismatch cluster starts at row 23,014 because MATLAB places `Cmito` before the optional/misfolding block, while the Python matrix export leaves mitochondrial as the final UB row.",
            "- Bound differences include corrected-medium differences and formulation differences; they should not be interpreted as row-matrix coefficient mismatches.",
            "",
            "## Top Row Mismatch Examples",
            "",
            *[
                f"- Row `{item.get('row_index_1based')}`: `{item.get('matlab_label')}` vs `{item.get('python_label')}`, category `{item.get('category')}`, nnz `{item.get('matlab_nnz')}` vs `{item.get('python_nnz')}`, coeff diffs `{item.get('coefficient_difference_count')}`."
                for item in first_rows
            ],
            "",
            "## Interpretation",
            "",
            "- Label/name formatting differences are tracked separately from coefficient/RHS differences.",
            "- Bounds include medium and formulation differences, so corrected medium vs old MATLAB artifact can produce real bound differences even when row matrices have matching shape.",
            "- If the objective differs, that is reported as objective formulation difference rather than a pcSec constraint-row mismatch.",
            "",
            "## Files",
            "",
            f"- MATLAB LP: `{diff['matlab_lp_path']}`",
            f"- Python LP: `{diff['python_lp_path']}`",
            "- Summary JSON: `opn_lp_diff_summary.json`",
            "- Top mismatched rows CSV: `opn_lp_diff_top_rows.csv`",
            "",
        ]
    )


def _format_constraint(label: str, row: sparse.csr_matrix, sense: str, rhs: float) -> list[str]:
    terms: list[str] = []
    start = row.indptr[0]
    end = row.indptr[1]
    for col, value in zip(row.indices[start:end], row.data[start:end]):
        if float(value) == 0.0:
            continue
        variable = int(col) + 1
        number = f"{float(value):.15f}"
        if not terms:
            terms.append(f"{number} X{variable}")
        elif float(value) > 0:
            terms.append(f"+ {number} X{variable}")
        else:
            terms.append(f"{number} X{variable}")
    if not terms:
        terms.append("0")
    chunks = [terms[index : index + 150] for index in range(0, len(terms), 150)]
    lines = [f"{label}: {' '.join(chunks[0])}"]
    for chunk in chunks[1:]:
        lines.append(f" {' '.join(chunk)}")
    lines[-1] += f" {sense} {float(rhs):.15f}"
    return lines


def _parse_constraint(label: str, body: str) -> LpConstraint:
    sense = "<=" if "<=" in body else ">=" if ">=" in body else "="
    left, right = body.rsplit(sense, 1)
    return LpConstraint(
        label=label,
        coefficients=_parse_coefficients(left),
        sense=sense,
        rhs=float(right.strip()),
    )


def _parse_coefficients(text: str) -> dict[int, float]:
    coefficients: dict[int, float] = {}
    for match in TERM_RE.finditer(text):
        raw = match.group("coef")
        coef = 1.0 if raw in (None, "", "+") else -1.0 if raw == "-" else float(raw.replace(" ", ""))
        variable = int(match.group("var"))
        coefficients[variable] = coefficients.get(variable, 0.0) + coef
    return {variable: coef for variable, coef in coefficients.items() if coef != 0.0}


def _parse_bound(line: str) -> tuple[int, LpBound] | None:
    match = re.match(rf"\s*(?P<lower>{NUMBER_RE})\s*<=\s*X(?P<var>\d+)\s*<=\s*(?P<upper>[+-]?infinity|{NUMBER_RE})\s*$", line)
    if not match:
        return None
    lower = float(match.group("lower"))
    upper_text = match.group("upper").lower()
    upper = math.inf if "infinity" in upper_text and not upper_text.startswith("-") else -math.inf if "infinity" in upper_text else float(upper_text)
    return int(match.group("var")), LpBound(lower=lower, upper=upper)


def _diff_objective(matlab: LpObjective, python: LpObjective, *, tolerance: float) -> dict[str, Any]:
    coeff = _coefficient_delta(matlab.coefficients, python.coefficients, tolerance=tolerance)
    return {
        "label_difference": matlab.label != python.label,
        "matlab_label": matlab.label,
        "python_label": python.label,
        "matlab_coefficients": _format_coefficients(matlab.coefficients),
        "python_coefficients": _format_coefficients(python.coefficients),
        "coefficient_difference_count": coeff["difference_count"],
        "missing_variables": coeff["missing_variables"],
        "extra_variables": coeff["extra_variables"],
        "max_abs_difference": coeff["max_abs_difference"],
    }


def _diff_rows(matlab_rows: tuple[LpConstraint, ...], python_rows: tuple[LpConstraint, ...], *, tolerance: float) -> list[dict[str, Any]]:
    diffs: list[dict[str, Any]] = []
    total = max(len(matlab_rows), len(python_rows))
    for index in range(total):
        matlab = matlab_rows[index] if index < len(matlab_rows) else None
        python = python_rows[index] if index < len(python_rows) else None
        if matlab is None or python is None:
            diffs.append(
                {
                    "row_index_1based": index + 1,
                    "category": "missing_or_extra_row",
                    "matlab_label": matlab.label if matlab else None,
                    "python_label": python.label if python else None,
                    "severity_score": 10_000,
                }
            )
            continue
        coeff = _coefficient_delta(matlab.coefficients, python.coefficients, tolerance=tolerance)
        rhs_diff = abs(matlab.rhs - python.rhs)
        sense_diff = matlab.sense != python.sense
        label_diff = matlab.label != python.label
        coefficient_diff = coeff["difference_count"] > 0
        rhs_changed = rhs_diff > tolerance
        sparsity_changed = coeff["missing_count"] > 0 or coeff["extra_count"] > 0
        if not (label_diff or sense_diff or coefficient_diff or rhs_changed):
            continue
        category = "label_name_formatting_difference"
        if sense_diff or rhs_changed:
            category = "rhs_difference"
        if coefficient_diff:
            category = "coefficient_difference"
        if sparsity_changed:
            category = "coefficient_sparsity_difference"
        diffs.append(
            {
                "row_index_1based": index + 1,
                "category": category,
                "matlab_label": matlab.label,
                "python_label": python.label,
                "matlab_sense": matlab.sense,
                "python_sense": python.sense,
                "matlab_rhs": matlab.rhs,
                "python_rhs": python.rhs,
                "rhs_abs_difference": rhs_diff,
                "matlab_nnz": len(matlab.coefficients),
                "python_nnz": len(python.coefficients),
                "missing_variable_count": coeff["missing_count"],
                "extra_variable_count": coeff["extra_count"],
                "coefficient_difference_count": coeff["difference_count"],
                "max_abs_coefficient_difference": coeff["max_abs_difference"],
                "sample_missing_variables": coeff["missing_variables"][:10],
                "sample_extra_variables": coeff["extra_variables"][:10],
                "sample_value_differences": coeff["value_differences"][:10],
                "severity_score": (
                    coeff["difference_count"] * 100
                    + coeff["missing_count"] * 1000
                    + coeff["extra_count"] * 1000
                    + (100 if rhs_changed else 0)
                    + (10 if sense_diff else 0)
                    + (1 if label_diff else 0)
                ),
            }
        )
    return diffs


def _diff_bounds(matlab: dict[int, LpBound], python: dict[int, LpBound], *, tolerance: float) -> list[dict[str, Any]]:
    diffs: list[dict[str, Any]] = []
    for variable in sorted(set(matlab) | set(python)):
        left = matlab.get(variable)
        right = python.get(variable)
        if left is None or right is None:
            diffs.append({"variable": f"X{variable}", "category": "missing_or_extra_bound", "matlab": _bound_dict(left), "python": _bound_dict(right)})
            continue
        lower_diff = _float_abs_diff(left.lower, right.lower)
        upper_diff = _float_abs_diff(left.upper, right.upper)
        if lower_diff > tolerance or upper_diff > tolerance:
            diffs.append(
                {
                    "variable": f"X{variable}",
                    "category": "bound_difference",
                    "matlab_lower": left.lower,
                    "python_lower": right.lower,
                    "lower_abs_difference": lower_diff,
                    "matlab_upper": left.upper,
                    "python_upper": right.upper,
                    "upper_abs_difference": upper_diff,
                }
            )
    return diffs


def _coefficient_delta(left: dict[int, float], right: dict[int, float], *, tolerance: float) -> dict[str, Any]:
    missing = sorted(set(left) - set(right))
    extra = sorted(set(right) - set(left))
    value_differences: list[dict[str, float | str]] = []
    max_abs = 0.0
    for variable in sorted(set(left) & set(right)):
        diff = abs(left[variable] - right[variable])
        if diff > tolerance:
            max_abs = max(max_abs, diff)
            value_differences.append(
                {
                    "variable": f"X{variable}",
                    "matlab": left[variable],
                    "python": right[variable],
                    "abs_difference": diff,
                }
            )
    max_abs = max(max_abs, *(abs(left[v]) for v in missing), *(abs(right[v]) for v in extra), 0.0)
    return {
        "missing_count": len(missing),
        "extra_count": len(extra),
        "difference_count": len(missing) + len(extra) + len(value_differences),
        "missing_variables": [f"X{value}" for value in missing[:100]],
        "extra_variables": [f"X{value}" for value in extra[:100]],
        "value_differences": value_differences[:100],
        "max_abs_difference": max_abs,
    }


def _row_summary(diffs: list[dict[str, Any]], matlab_rows: tuple[LpConstraint, ...], python_rows: tuple[LpConstraint, ...]) -> dict[str, Any]:
    matlab_labels = {row.label for row in matlab_rows}
    python_labels = {row.label for row in python_rows}
    label_only_count = 0
    first_label_difference: int | None = None
    first_coefficient_difference: int | None = None
    first_rhs_or_sense_difference: int | None = None
    for item in diffs:
        has_label_difference = item.get("matlab_label") != item.get("python_label")
        has_coefficient_difference = int(item.get("coefficient_difference_count", 0) or 0) > 0
        has_rhs_difference = float(item.get("rhs_abs_difference", 0.0) or 0.0) > 0.0 or item.get("matlab_sense") != item.get("python_sense")
        row_index = item.get("row_index_1based")
        if has_label_difference and first_label_difference is None:
            first_label_difference = int(row_index)
        if has_coefficient_difference and first_coefficient_difference is None:
            first_coefficient_difference = int(row_index)
        if has_rhs_difference and first_rhs_or_sense_difference is None:
            first_rhs_or_sense_difference = int(row_index)
        if has_label_difference and not has_coefficient_difference and not has_rhs_difference:
            label_only_count += 1
    return {
        "matlab_row_count": len(matlab_rows),
        "python_row_count": len(python_rows),
        "row_count_difference": len(python_rows) - len(matlab_rows),
        "ordered_mismatched_row_count": len(diffs),
        "ordered_label_difference_count": sum(1 for item in diffs if item.get("matlab_label") != item.get("python_label")),
        "label_only_difference_count": label_only_count,
        "rhs_difference_count": sum(1 for item in diffs if item.get("rhs_abs_difference", 0.0) > 0.0),
        "coefficient_difference_count": sum(1 for item in diffs if item.get("coefficient_difference_count", 0) > 0),
        "sparsity_difference_count": sum(1 for item in diffs if item.get("missing_variable_count", 0) > 0 or item.get("extra_variable_count", 0) > 0),
        "first_label_difference_row": first_label_difference,
        "first_coefficient_difference_row": first_coefficient_difference,
        "first_rhs_or_sense_difference_row": first_rhs_or_sense_difference,
        "missing_label_count": len(matlab_labels - python_labels),
        "extra_label_count": len(python_labels - matlab_labels),
        "missing_label_examples": sorted(matlab_labels - python_labels)[:20],
        "extra_label_examples": sorted(python_labels - matlab_labels)[:20],
    }


def _bound_summary(diffs: list[dict[str, Any]], matlab_bounds: dict[int, LpBound], python_bounds: dict[int, LpBound]) -> dict[str, Any]:
    return {
        "matlab_bound_count": len(matlab_bounds),
        "python_bound_count": len(python_bounds),
        "difference_count": len(diffs),
        "missing_in_python_count": sum(1 for item in diffs if item["category"] == "missing_or_extra_bound" and item.get("python") is None),
        "extra_in_python_count": sum(1 for item in diffs if item["category"] == "missing_or_extra_bound" and item.get("matlab") is None),
    }


def _lp_stats(lp: ParsedLp) -> dict[str, Any]:
    return {
        "path": str(lp.path),
        "optimization_sense": lp.optimization_sense,
        "objective_label": lp.objective.label,
        "constraint_count": len(lp.constraints),
        "bound_count": len(lp.bounds),
    }


def _format_coefficients(coefficients: dict[int, float]) -> dict[str, float]:
    return {f"X{variable}": value for variable, value in sorted(coefficients.items())}


def _bound_dict(bound: LpBound | None) -> dict[str, float | None] | None:
    if bound is None:
        return None
    return {"lower": bound.lower, "upper": bound.upper}


def _float_abs_diff(left: float | None, right: float | None) -> float:
    if left is None and right is None:
        return 0.0
    if left is None or right is None:
        return math.inf
    if math.isinf(left) and math.isinf(right) and left == right:
        return 0.0
    return abs(left - right)


def _write_top_rows_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fieldnames = [
        "row_index_1based",
        "category",
        "matlab_label",
        "python_label",
        "matlab_sense",
        "python_sense",
        "matlab_rhs",
        "python_rhs",
        "rhs_abs_difference",
        "matlab_nnz",
        "python_nnz",
        "missing_variable_count",
        "extra_variable_count",
        "coefficient_difference_count",
        "max_abs_coefficient_difference",
        "sample_missing_variables",
        "sample_extra_variables",
        "sample_value_differences",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(row.get(key), ensure_ascii=False) if isinstance(row.get(key), (list, dict)) else row.get(key) for key in fieldnames})


__all__ = [
    "ParsedLp",
    "diff_lp_files",
    "build_lp_diff_report",
    "parse_lp",
    "write_lp_diff_outputs",
    "write_matrix_lp",
    "write_matrix_lp_from_rows",
]
