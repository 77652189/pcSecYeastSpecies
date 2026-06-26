"""Parse and solve the MATLAB baseline LP with scipy HiGHS.

This tells us whether the 81% objective gap is in the LP itself or in how
the probe constructs the LP. If scipy gets -1.07282263 from the MATLAB LP,
the probe's LP construction is wrong. If scipy gets -0.195, the issue is
in the solver or LP parsing.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

from pcsec_pichia.probe import repo_root

MATLAB_LP = (
    "local_runs/OPN_PPA_glc_smoke/"
    "Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_"
    "misfolddefault_noMisfoldEq_noRiboEq_PP.lp"
)
MATLAB_OBJECTIVE = -1.07282263e00


def parse_lp_file(path: Path) -> dict[str, object]:
    """Parse an LP file into objective, constraints, and bounds."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    # Find sections
    sense = None
    obj_var = None
    subject_to_start = None
    bounds_start = None
    end_start = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "Maximize" or stripped == "Minimize":
            sense = "Maximize" if stripped == "Maximize" else "Minimize"
        elif stripped.startswith("obj:"):
            m = re.match(r"obj:\s*X(\d+)", stripped)
            if m:
                obj_var = int(m.group(1))
        elif stripped == "Subject To":
            subject_to_start = i
        elif stripped == "Bounds":
            bounds_start = i
        elif stripped == "End":
            end_start = i
            break

    # Parse constraints (may span multiple lines)
    # A constraint line starts with a label like "C1:" or "CM1:" or "Cmito:"
    constraint_lines = lines[subject_to_start + 1 : bounds_start]

    # Join multi-line constraints
    constraints: list[tuple[str, str]] = []  # (label, full_text)
    current_label = None
    current_text = ""

    for line in constraint_lines:
        stripped = line.strip()
        if not stripped:
            continue
        m = re.match(r"^(C\d+|CM\d+|Cmito)\s*:\s*(.*)", stripped)
        if m:
            # New constraint
            if current_label is not None:
                constraints.append((current_label, current_text))
            current_label = m.group(1)
            current_text = m.group(2)
        else:
            # Continuation of previous constraint
            current_text += " " + stripped

    if current_label is not None:
        constraints.append((current_label, current_text))

    # Parse each constraint into (coefficients, sense, rhs)
    # Format: "coef X%d + coef X%d ... (=|<=|>=) rhs"
    eq_rows: list[list[tuple[int, float]]] = []
    eq_rhs: list[float] = []
    ub_rows: list[list[tuple[int, float]]] = []
    ub_rhs: list[float] = []

    for label, text in constraints:
        # Split into LHS and RHS
        for op in (" = ", " <= ", " >= "):
            if op in text:
                lhs, rhs_str = text.rsplit(op, 1)
                rhs = float(rhs_str.strip())
                sense_char = op.strip()
                break
        else:
            continue

        # Parse LHS coefficients
        terms: list[tuple[int, float]] = []
        # Tokenize: split by spaces, handle "+ " and "- " prefixes
        # Replace "+ " with "\n+" and "- " with "\n-" for easier parsing
        normalized = lhs.replace("+ ", "\n+ ").replace("- ", "\n- ")
        tokens = [t.strip() for t in normalized.split("\n") if t.strip()]
        for token in tokens:
            # Try explicit coefficient: "coef X%d"
            m = re.match(r"([+-]?[\d.eE+-]+)\s+X(\d+)", token)
            if m:
                coef = float(m.group(1))
                var = int(m.group(2))
                terms.append((var - 1, coef))  # 0-based
                continue
            # Try implied coefficient 1.0: "X%d"
            m = re.match(r"X(\d+)", token)
            if m:
                var = int(m.group(1))
                terms.append((var - 1, 1.0))  # 0-based, implied 1.0

        if sense_char == "=":
            eq_rows.append(terms)
            eq_rhs.append(rhs)
        elif sense_char == "<=":
            ub_rows.append(terms)
            ub_rhs.append(rhs)
        elif sense_char == ">=":
            # Convert >= to <= by negating
            ub_rows.append([(v, -c) for v, c in terms])
            ub_rhs.append(-rhs)

    # Parse bounds
    bound_lines = lines[bounds_start + 1 : end_start]
    # Determine n_vars from the highest variable index in constraints AND bounds
    max_var_in_constraints = 0
    for terms in eq_rows + ub_rows:
        for v, _ in terms:
            if v > max_var_in_constraints:
                max_var_in_constraints = v
    # Also check bounds for max var
    max_var_in_bounds = 0
    for line in bound_lines:
        stripped = line.strip()
        m = re.match(r"([-\d.eE+]+)\s*<=\s*X(\d+)\s*<=\s*(\+?infinity|[-\d.eE+]+)", stripped)
        if m:
            var = int(m.group(2)) - 1
            if var > max_var_in_bounds:
                max_var_in_bounds = var
    n_vars = max(max_var_in_constraints, max_var_in_bounds) + 1
    lb = np.full(n_vars, -np.inf)
    ub = np.full(n_vars, np.inf)

    for line in bound_lines:
        stripped = line.strip()
        if not stripped:
            continue
        # Format: "lb <= X%d <= ub" or "lb <= X%d <= +infinity"
        m = re.match(r"([-\d.eE+]+)\s*<=\s*X(\d+)\s*<=\s*(\+?infinity|[-\d.eE+]+)", stripped)
        if m:
            lower = float(m.group(1))
            var = int(m.group(2)) - 1  # 0-based
            upper_str = m.group(3)
            if "infinity" in upper_str:
                upper = np.inf
            else:
                upper = float(upper_str)
            lb[var] = lower
            ub[var] = upper

    # Build sparse matrices
    def build_matrix(rows: list[list[tuple[int, float]]], n_cols: int) -> sparse.csr_matrix:
        data: list[float] = []
        row_ind: list[int] = []
        col_ind: list[int] = []
        for i, terms in enumerate(rows):
            for col, val in terms:
                data.append(val)
                row_ind.append(i)
                col_ind.append(col)
        return sparse.csr_matrix((data, (row_ind, col_ind)), shape=(len(rows), n_cols))

    A_eq = build_matrix(eq_rows, n_vars)
    b_eq = np.array(eq_rhs, dtype=float)
    A_ub = build_matrix(ub_rows, n_vars)
    b_ub = np.array(ub_rhs, dtype=float)

    return {
        "sense": sense,
        "obj_var": obj_var - 1,  # 0-based
        "n_vars": n_vars,
        "A_eq": A_eq,
        "b_eq": b_eq,
        "A_ub": A_ub,
        "b_ub": b_ub,
        "lb": lb,
        "ub": ub,
        "n_eq": len(eq_rows),
        "n_ub": len(ub_rows),
    }


def main() -> int:
    root = repo_root()
    lp_path = root / MATLAB_LP

    print(f"Parsing MATLAB LP: {lp_path}")
    lp = parse_lp_file(lp_path)

    print(f"  sense: {lp['sense']}")
    print(f"  obj_var: X{lp['obj_var']+1} (0-based: {lp['obj_var']})")
    print(f"  n_vars: {lp['n_vars']}")
    print(f"  n_eq: {lp['n_eq']}")
    print(f"  n_ub: {lp['n_ub']}")

    # Build objective vector (maximize → minimize negative)
    c = np.zeros(lp["n_vars"], dtype=float)
    if lp["sense"] == "Maximize":
        c[lp["obj_var"]] = -1.0
    else:
        c[lp["obj_var"]] = 1.0

    print(f"\nSolving with scipy HiGHS...")
    result = linprog(
        c=c,
        A_eq=lp["A_eq"],
        b_eq=lp["b_eq"],
        A_ub=lp["A_ub"],
        b_ub=lp["b_ub"],
        bounds=list(zip(lp["lb"], lp["ub"])),
        method="highs",
        options={"presolve": True},
    )

    print(f"  status: {result.status}")
    print(f"  success: {result.success}")
    print(f"  message: {result.message}")
    if result.success:
        obj_val = -result.fun if lp["sense"] == "Maximize" else result.fun
        print(f"  objective_value: {obj_val}")
        print(f"  MATLAB SoPlex objective: {MATLAB_OBJECTIVE}")
        print(f"  diff: {abs(obj_val - MATLAB_OBJECTIVE):.6e}")
        print(f"  relative diff: {abs(obj_val - MATLAB_OBJECTIVE) / abs(MATLAB_OBJECTIVE):.4%}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
