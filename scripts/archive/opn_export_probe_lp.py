"""Export the probe's LP to a file in the same format as MATLAB writeLPGlc.

This lets us do a line-by-line diff between the probe LP and the MATLAB LP
to pinpoint exact coefficient differences.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import sparse

from pcsec_pichia.probe import (
    build_pcsec_constraint_matrices,
    build_supported_target_model,
    build_target_enzymedata,
    load_aa_stoichiometry,
    load_combined_enzymedata,
    load_metabolic_enzymedata,
    load_pcsec_pichia_model,
    load_secretory_enzymedata,
    load_targets,
    prepare_glucose_model,
    repo_root,
)


def write_lp_file(model, A_eq, b_eq, A_ub, b_ub, objective_rxn, path: Path) -> None:
    """Write an LP file in the same format as MATLAB writeLPGlc."""
    reaction_index = model.reaction_index
    n = len(model.rxns)
    obj_idx = reaction_index[objective_rxn] + 1  # 1-based

    lines: list[str] = []
    lines.append("Maximize")
    lines.append(f"obj: X{obj_idx}")
    lines.append("Subject To")

    def format_constraint(label: str, terms: list[str], sense: str, rhs: float) -> list[str]:
        """Format a constraint, breaking long lines every 150 terms (like MATLAB)."""
        result_lines = []
        # Split terms into chunks of 150
        chunk_size = 150
        chunks = [terms[k:k+chunk_size] for k in range(0, len(terms), chunk_size)]
        if not chunks:
            chunks = [[]]
        # First chunk gets the label
        result_lines.append(f"{label}: {' '.join(chunks[0])}")
        for chunk in chunks[1:]:
            result_lines.append(f" {' '.join(chunk)}")
        # Append sense + rhs to the last line
        result_lines[-1] += f" {sense} {rhs:.15f}"
        return result_lines

    # Equality constraints
    A_eq_csr = sparse.csr_matrix(A_eq)
    for i in range(A_eq.shape[0]):
        row = A_eq_csr.getrow(i)
        terms = []
        for j in range(row.indptr[0], row.indptr[1]):
            col = row.indices[j] + 1  # 1-based
            val = row.data[j]
            if val == 0:
                continue
            if not terms:
                terms.append(f"{val:.15f} X{col}")
            else:
                if val > 0:
                    terms.append(f"+ {val:.15f} X{col}")
                else:
                    terms.append(f"{val:.15f} X{col}")
        rhs = b_eq[i]
        lines.extend(format_constraint(f"C{i+1}", terms, "=", rhs))

    # Inequality constraints (<=)
    A_ub_csr = sparse.csr_matrix(A_ub)
    offset = A_eq.shape[0]
    for i in range(A_ub.shape[0]):
        row = A_ub_csr.getrow(i)
        terms = []
        for j in range(row.indptr[0], row.indptr[1]):
            col = row.indices[j] + 1
            val = row.data[j]
            if val == 0:
                continue
            if not terms:
                terms.append(f"{val:.15f} X{col}")
            else:
                if val > 0:
                    terms.append(f"+ {val:.15f} X{col}")
                else:
                    terms.append(f"{val:.15f} X{col}")
        rhs = b_ub[i]
        lines.extend(format_constraint(f"C{offset+i+1}", terms, "<=", rhs))

    # Bounds
    lines.append("Bounds")
    for i in range(n):
        lb = float(model.lb[i])
        ub = float(model.ub[i])
        if ub >= 100:
            lines.append(f"{lb:.6f} <= X{i+1} <= +infinity")
        else:
            lines.append(f"{lb:.6f} <= X{i+1} <= {ub:.6f}")

    lines.append("End")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    root = repo_root()
    model = prepare_glucose_model(load_pcsec_pichia_model(root), media_type=4)
    amino_acids = load_aa_stoichiometry(root)
    metabolic = load_metabolic_enzymedata(root)
    secretory = load_secretory_enzymedata(root)
    combined = load_combined_enzymedata(root)
    targets = load_targets(root, targets_json=None)
    opn = next(t for t in targets if "OPN" in t.protein_id)

    build = build_supported_target_model(model, opn, amino_acids)
    target_enzyme = build_target_enzymedata(opn, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzyme.reaction_coefficients)
    target_combined = combined.with_target(target_enzyme)

    dilution_blocks = {rxn: (None, 0.0) for rxn in build.model.rxns if "dilution_misfolding" in rxn}
    fixed_model = build.model.with_bounds({
        "BIOMASS": (0.10, 0.10),
        build.exchange_reaction_id: (1e-8, 1e-8),
        **dilution_blocks,
    })

    A_eq, b_eq, A_ub, b_ub, counts = build_pcsec_constraint_matrices(
        fixed_model,
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu=0.10,
        total_protein_content=0.37,
        unmodeled_er_protein_fraction=0.040,
        mitochondrial_protein_fraction=0.05,
        write_ribosome_translation_constraint=False,
        write_misfolding_constraints=False,
    )

    out_path = root / "local_runs" / "pichia_python" / "alignment" / "probe_opn.lp"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_lp_file(fixed_model, A_eq, b_eq, A_ub, b_ub, "Ex_glc_D", out_path)

    print(f"Probe LP written to: {out_path}")
    print(f"  rows_eq={A_eq.shape[0]}  rows_ub={A_ub.shape[0]}  cols={A_eq.shape[1]}")
    print(f"  MATLAB LP: local_runs/OPN_PPA_glc_smoke/Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_misfolddefault_noMisfoldEq_noRiboEq_PP.lp")

    # Quick diff: compare first 5 metabolic constraints
    ml_path = root / "local_runs" / "OPN_PPA_glc_smoke" / "Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_misfolddefault_noMisfoldEq_noRiboEq_PP.lp"
    ml_lines = ml_path.read_text(encoding="utf-8", errors="replace").splitlines()
    py_lines = out_path.read_text(encoding="utf-8", errors="replace").splitlines()

    print(f"\nMATLAB LP: {len(ml_lines)} lines")
    print(f"Probe LP:  {len(py_lines)} lines")

    # Find first difference
    print("\n=== First differences (skipping header) ===")
    diffs_found = 0
    for i, (ml_line, py_line) in enumerate(zip(ml_lines[3:], py_lines[3:])):
        if ml_line.strip() != py_line.strip():
            diffs_found += 1
            print(f"  line {i+4}:")
            print(f"    MATLAB: {ml_line[:150]}")
            print(f"    probe:  {py_line[:150]}")
            if diffs_found >= 5:
                break

    if diffs_found == 0:
        print("  (no differences found in overlapping lines)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
