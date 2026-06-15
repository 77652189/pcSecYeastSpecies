"""Check pcSec LP feasibility with SciPy/HiGHS.

This is a diagnostic helper for generated local LP files. It is not a
replacement for SoPlex validation; it gives a fast pass/fail signal and, for
infeasible files, bisects the first constraint prefix that becomes infeasible.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


LABEL_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*\d*):\s*(.*)$")
TERM_RE = re.compile(
    r"([+-]?\s*(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+X(\d+)"
)
REL_RE = re.compile(
    r"\s*(.*?)\s*(<=|>=|=)\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$"
)
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[+-]?infinity"
BOUND_RE = re.compile(rf"({NUMBER})\s*<=\s*X(\d+)\s*<=\s*({NUMBER})", re.I)


def _bound_value(text: str) -> float | None:
    lowered = text.lower()
    if "inf" in lowered:
        return None if not lowered.startswith("-") else -np.inf
    return float(text)


def parse_lp(path: Path):
    rows: list[tuple[str, str]] = []
    bounds_lines: list[str] = []
    section: str | None = None
    current_label: str | None = None
    current_expr: list[str] = []

    def flush_current() -> None:
        nonlocal current_label, current_expr
        if current_label is not None:
            rows.append((current_label, " ".join(current_expr)))
        current_label = None
        current_expr = []

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        lowered = stripped.lower()
        if lowered == "subject to":
            section = "rows"
            continue
        if lowered == "bounds":
            flush_current()
            section = "bounds"
            continue
        if lowered == "end":
            flush_current()
            section = None
            continue

        if section == "rows":
            match = LABEL_RE.match(line)
            if match:
                flush_current()
                current_label = match.group(1)
                current_expr = [match.group(2).strip()]
            else:
                current_expr.append(stripped)
        elif section == "bounds":
            bounds_lines.append(stripped)

    parsed_rows = []
    max_var = 0
    for label, expr in rows:
        relation_match = REL_RE.match(expr)
        if not relation_match:
            raise ValueError(f"Could not parse row {label}: {expr[-200:]}")
        lhs, relation, rhs = relation_match.groups()
        terms = []
        for coef_text, var_text in TERM_RE.findall(lhs):
            var = int(var_text)
            max_var = max(max_var, var)
            terms.append((var - 1, float(coef_text.replace(" ", ""))))
        parsed_rows.append((label, relation, float(rhs), terms, expr))

    bounds: list[tuple[float | None, float | None]] = [(0.0, None) for _ in range(max_var)]
    for line in bounds_lines:
        match = BOUND_RE.fullmatch(line)
        if not match:
            raise ValueError(f"Could not parse bound: {line}")
        lower, var_text, upper = match.groups()
        bounds[int(var_text) - 1] = (_bound_value(lower), _bound_value(upper))

    return parsed_rows, bounds


def solve_prefix(parsed_rows, bounds, prefix: int | None = None):
    if prefix is None:
        prefix = len(parsed_rows)

    num_vars = len(bounds)
    eq_i: list[int] = []
    eq_j: list[int] = []
    eq_v: list[float] = []
    b_eq: list[float] = []
    ub_i: list[int] = []
    ub_j: list[int] = []
    ub_v: list[float] = []
    b_ub: list[float] = []
    eq_row = 0
    ub_row = 0

    for _label, relation, rhs, terms, _expr in parsed_rows[:prefix]:
        if relation == "=":
            for var, coef in terms:
                eq_i.append(eq_row)
                eq_j.append(var)
                eq_v.append(coef)
            b_eq.append(rhs)
            eq_row += 1
        elif relation == "<=":
            for var, coef in terms:
                ub_i.append(ub_row)
                ub_j.append(var)
                ub_v.append(coef)
            b_ub.append(rhs)
            ub_row += 1
        else:
            for var, coef in terms:
                ub_i.append(ub_row)
                ub_j.append(var)
                ub_v.append(-coef)
            b_ub.append(-rhs)
            ub_row += 1

    a_eq = coo_matrix((eq_v, (eq_i, eq_j)), shape=(eq_row, num_vars)).tocsr() if eq_row else None
    a_ub = coo_matrix((ub_v, (ub_i, ub_j)), shape=(ub_row, num_vars)).tocsr() if ub_row else None
    return linprog(
        np.zeros(num_vars),
        A_ub=a_ub,
        b_ub=np.array(b_ub) if ub_row else None,
        A_eq=a_eq,
        b_eq=np.array(b_eq) if eq_row else None,
        bounds=bounds,
        method="highs",
        options={"presolve": True, "time_limit": 60},
    )


def first_bad_constraint(parsed_rows, bounds):
    low = 0
    high = len(parsed_rows)
    while low + 1 < high:
        mid = (low + high) // 2
        result = solve_prefix(parsed_rows, bounds, mid)
        if result.status == 0:
            low = mid
        else:
            high = mid
    label, _relation, _rhs, _terms, expr = parsed_rows[high - 1]
    return high, label, expr


def check_path(path: Path, show_first_bad: bool) -> int:
    parsed_rows, bounds = parse_lp(path)
    result = solve_prefix(parsed_rows, bounds)
    print(f"{path}\trows={len(parsed_rows)}\tvars={len(bounds)}\tstatus={result.status}\t{result.message}")
    if result.status != 0 and show_first_bad:
        index, label, expr = first_bad_constraint(parsed_rows, bounds)
        print(f"  first_bad={index}:{label}:{expr[:240]}")
    return result.status


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("lp", nargs="+", type=Path)
    parser.add_argument("--first-bad", action="store_true")
    args = parser.parse_args()

    statuses = [check_path(path, args.first_bad) for path in args.lp]
    return 0 if all(status == 0 for status in statuses) else 1


if __name__ == "__main__":
    raise SystemExit(main())
