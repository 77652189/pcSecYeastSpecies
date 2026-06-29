"""Diff the probe LP against the MATLAB LP constraint by constraint.

Both LPs have the same dimensions (23016 rows, 29057 cols) but different
objective values. This script parses both into normalized (variable, coef)
dicts per constraint and reports which constraints differ.
"""

from __future__ import annotations

import re
from pathlib import Path

from pcsec_pichia.probe import repo_root

MATLAB_LP = (
    "local_runs/OPN_PPA_glc_smoke/"
    "Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_"
    "misfolddefault_noMisfoldEq_noRiboEq_PP.lp"
)
PROBE_LP = "local_runs/pichia_python/alignment/probe_opn.lp"


def parse_lp_constraints(path: Path) -> dict[str, dict[int, float] | tuple[str, float]]:
    """Parse an LP file into {label: (terms_dict, sense, rhs)}."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    subject_to_start = None
    bounds_start = None
    for i, line in enumerate(lines):
        if line.strip() == "Subject To":
            subject_to_start = i
        elif line.strip() == "Bounds":
            bounds_start = i
            break

    constraint_lines = lines[subject_to_start + 1 : bounds_start]

    constraints: dict[str, tuple[dict[int, float], str, float]] = {}
    current_label = None
    current_text = ""

    for line in constraint_lines:
        stripped = line.strip()
        if not stripped:
            continue
        m = re.match(r"^(C\d+|CM\d+|Cmito)\s*:\s*(.*)", stripped)
        if m:
            if current_label is not None:
                constraints[current_label] = _parse_constraint_text(current_text)
            current_label = m.group(1)
            current_text = m.group(2)
        else:
            current_text += " " + stripped

    if current_label is not None:
        constraints[current_label] = _parse_constraint_text(current_text)

    return constraints


def _parse_constraint_text(text: str) -> tuple[dict[int, float], str, float]:
    """Parse constraint text into (terms, sense, rhs)."""
    sense = None
    rhs = 0.0
    lhs = text

    for op in (" = ", " <= ", " >= "):
        if op in text:
            lhs, rhs_str = text.rsplit(op, 1)
            rhs = float(rhs_str.strip())
            sense = op.strip()
            break

    terms: dict[int, float] = {}
    normalized = lhs.replace("+ ", "\n+ ").replace("- ", "\n- ")
    tokens = [t.strip() for t in normalized.split("\n") if t.strip()]
    for token in tokens:
        m = re.match(r"([+-]?[\d.eE+-]+)\s+X(\d+)", token)
        if m:
            coef = float(m.group(1))
            var = int(m.group(2)) - 1  # 0-based
            terms[var] = terms.get(var, 0.0) + coef
            continue
        m = re.match(r"X(\d+)", token)
        if m:
            var = int(m.group(1)) - 1
            terms[var] = terms.get(var, 0.0) + 1.0

    return (terms, sense or "=", rhs)


def main() -> int:
    root = repo_root()
    ml = parse_lp_constraints(root / MATLAB_LP)
    py = parse_lp_constraints(root / PROBE_LP)

    # Convert to ordered lists by position (not label)
    # MATLAB order: C1..C20221 (stoich), CM1..CM2732 (met), CM2733..CM2790 (sec),
    #               CM2791 (protein_mass), CM2792 (unmodeled_ER), Cmito,
    #               CM2793 (proteasome), CM2796 (ribo_assembly)
    # Probe order: C1..C23016 (all C-prefixed, sequential)

    def sort_key(label: str) -> int:
        m = re.match(r"[A-Za-z]+(\d+)", label)
        return int(m.group(1)) if m else 0

    # For MATLAB: CM constraints come after C constraints in writeLPGlc order
    # But in the LP file, they appear in write order: C1..C20221, CM1..CM2790, CM2791, CM2792, Cmito, CM2793, CM2796
    # For probe: all are C1..C23016 in order

    # Build ordered lists
    ml_ordered = sorted(ml.items(), key=lambda x: _matlab_order(x[0]))
    py_ordered = sorted(py.items(), key=lambda x: sort_key(x[0]))

    print(f"MATLAB constraints: {len(ml_ordered)}")
    print(f"Probe constraints:  {len(py_ordered)}")

    if len(ml_ordered) != len(py_ordered):
        print("ERROR: constraint count mismatch!")
        return 1

    diffs = []
    for i, ((ml_label, (ml_terms, ml_sense, ml_rhs)), (py_label, (py_terms, py_sense, py_rhs))) in enumerate(zip(ml_ordered, py_ordered)):
        # Compare RHS
        rhs_diff = abs(ml_rhs - py_rhs)
        if rhs_diff > 1e-10:
            diffs.append((i, ml_label, py_label, "rhs", ml_rhs, py_rhs, rhs_diff))
            continue

        # Compare sense
        if ml_sense != py_sense:
            diffs.append((i, ml_label, py_label, "sense", ml_sense, py_sense, None))
            continue

        # Compare terms
        all_vars = set(ml_terms.keys()) | set(py_terms.keys())
        term_diffs = []
        for var in sorted(all_vars):
            ml_coef = ml_terms.get(var, 0.0)
            py_coef = py_terms.get(var, 0.0)
            diff = abs(ml_coef - py_coef)
            if diff > 1e-10:
                term_diffs.append((var, ml_coef, py_coef, diff))

        if term_diffs:
            diffs.append((i, ml_label, py_label, "terms", len(term_diffs), term_diffs[:5], None))

    print(f"\nTotal differing constraints: {len(diffs)}")
    print(f"Total matching constraints:  {len(ml_ordered) - len(diffs)}")

    if diffs:
        print("\n=== Differences (first 30) ===")
        for row_idx, ml_label, py_label, diff_type, *rest in diffs[:30]:
            if diff_type == "rhs":
                print(f"  row {row_idx} ({ml_label} vs {py_label}): RHS diff  MATLAB={rest[0]:.6g}  probe={rest[1]:.6g}  |diff|={rest[2]:.6e}")
            elif diff_type == "sense":
                print(f"  row {row_idx} ({ml_label} vs {py_label}): SENSE diff  MATLAB={rest[0]}  probe={rest[1]}")
            elif diff_type == "terms":
                n_diffs = rest[0]
                examples = rest[1]
                print(f"  row {row_idx} ({ml_label} vs {py_label}): {n_diffs} term coefficient differences")
                for var, ml_c, py_c, d in examples:
                    print(f"    X{var+1}: MATLAB={ml_c:.6g}  probe={py_c:.6g}  |diff|={d:.6e}")

    # Section summary
    print("\n=== Section summary ===")
    n_stoich_diffs = sum(1 for d in diffs if d[0] < 20221)
    n_met_diffs = sum(1 for d in diffs if 20221 <= d[0] < 20221 + 2732)
    n_sec_diffs = sum(1 for d in diffs if 20221 + 2732 <= d[0] < 20221 + 2732 + 58)
    n_protein_diffs = sum(1 for d in diffs if d[0] == 20221 + 2732 + 58)
    n_unmodeler_diffs = sum(1 for d in diffs if d[0] == 20221 + 2732 + 58 + 1)
    n_mito_diffs = sum(1 for d in diffs if d[0] == 20221 + 2732 + 58 + 2)
    n_prot_diffs = sum(1 for d in diffs if d[0] == 20221 + 2732 + 58 + 3)
    n_ribo_diffs = sum(1 for d in diffs if d[0] == 20221 + 2732 + 58 + 4)
    print(f"  stoichiometric:    {n_stoich_diffs} diffs (out of 20221)")
    print(f"  metabolic:         {n_met_diffs} diffs (out of 2732)")
    print(f"  secretory:         {n_sec_diffs} diffs (out of 58)")
    print(f"  protein_mass:      {n_protein_diffs} diffs (out of 1)")
    print(f"  unmodeled_ER:      {n_unmodeler_diffs} diffs (out of 1)")
    print(f"  mitochondrial:     {n_mito_diffs} diffs (out of 1)")
    print(f"  proteasome:        {n_prot_diffs} diffs (out of 1)")
    print(f"  ribosome_assembly: {n_ribo_diffs} diffs (out of 1)")

    return 0


def _matlab_order(label: str) -> int:
    """Return the write order index for a MATLAB LP label."""
    if label == "Cmito":
        return 20221 + 2732 + 58 + 2  # after protein_mass + unmodeled_ER
    m = re.match(r"([A-Za-z]+)(\d+)", label)
    if not m:
        return 99999
    prefix, num = m.group(1), int(m.group(2))
    if prefix == "C":
        return num  # C1..C20221
    if prefix == "CM":
        # CM1..CM2732 = metabolic (20222..22953)
        # CM2733..CM2790 = secretory (22954..23011)
        # CM2791 = protein_mass (23012)
        # CM2792 = unmodeled_ER (23013)
        # CM2793 = proteasome (23014)
        # CM2796 = ribosome_assembly (23015)
        if num <= 2732:
            return 20221 + num
        elif num <= 2790:
            return 20221 + 2732 + (num - 2732)
        elif num == 2791:
            return 20221 + 2732 + 58 + 1
        elif num == 2792:
            return 20221 + 2732 + 58 + 2
        elif num == 2793:
            return 20221 + 2732 + 58 + 3
        elif num == 2796:
            return 20221 + 2732 + 58 + 4
    return 99999


if __name__ == "__main__":
    raise SystemExit(main())
