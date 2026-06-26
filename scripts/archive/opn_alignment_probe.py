"""OPN MATLAB-vs-Python alignment probe.

This script runs the pcSecPichia probe prototype under the same setup as the
MATLAB ``local_opn_pichia_glc`` baseline and prints a side-by-side comparison
so we can quantify where the two LPs diverge.

MATLAB baseline (already on disk):
    local_runs/OPN_PPA_glc_smoke/
        Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_misfolddefault_noMisfoldEq_noRiboEq_PP.lp
        Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_misfolddefault_noMisfoldEq_noRiboEq_PP.lp.float.out

Setup:
    candidate = OPN_ALPHA_FULL_PROJECT
    mu = 0.10
    mediaType = 4
    productionRatio = 1e-8
    blockMisfoldDilution = true (default)
    objective = Maximize Ex_glc_D
"""

from __future__ import annotations

import json
from pathlib import Path

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
    solve_pcsec_maximize,
)

MATLAB_BASELINE_LP = (
    "local_runs/OPN_PPA_glc_smoke/"
    "Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_"
    "misfolddefault_noMisfoldEq_noRiboEq_PP.lp"
)
MATLAB_BASELINE_OUT = MATLAB_BASELINE_LP + ".float.out"
MATLAB_OBJECTIVE = -1.07282263e00
MATLAB_ROWS = 23016
MATLAB_COLS = 29057


def _extract_matlab_dimensions() -> dict[str, object]:
    """Pull the LP row/column counts from the SoPlex .out file."""
    out_path = repo_root() / MATLAB_BASELINE_OUT
    text = out_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("LP has") and "rows" in stripped and "columns" in stripped:
            parts = stripped.replace("LP has", "").replace("rows", "").replace("columns", "")
            tokens = [t for t in parts.split() if t.isdigit()]
            if len(tokens) >= 2:
                return {"rows": int(tokens[0]), "cols": int(tokens[1]), "nonzeros": int(tokens[2]) if len(tokens) > 2 else None}
    return {"rows": MATLAB_ROWS, "cols": MATLAB_COLS, "nonzeros": None}


def _extract_matlab_objective() -> float | None:
    out_path = repo_root() / MATLAB_BASELINE_OUT
    text = out_path.read_text(encoding="utf-8", errors="replace")
    for marker in ("Objective value", "  Objective value"):
        for line in text.splitlines():
            if marker in line and ":" in line:
                try:
                    return float(line.split(":", 1)[1].strip())
                except ValueError:
                    continue
    return None


def run_python_alignment() -> dict[str, object]:
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

    result, _ = solve_pcsec_maximize(
        fixed_model,
        "Ex_glc_D",
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu=0.10,
        key_reactions=("BIOMASS", "Ex_glc_D", "Ex_o2", build.exchange_reaction_id),
    )

    return {
        "python": {
            "rows_eq": int(A_eq.shape[0]),
            "rows_ub": int(A_ub.shape[0]),
            "rows_total": int(A_eq.shape[0] + A_ub.shape[0]),
            "cols": int(A_eq.shape[1]),
            "nonzeros_eq": int(A_eq.nnz),
            "nonzeros_ub": int(A_ub.nnz),
            "objective_value": result.objective_value,
            "status": result.status,
            "success": result.success,
            "fluxes": result.fluxes,
            "constraint_counts": counts,
            "dilution_misfolding_blocked": len(dilution_blocks),
        },
        "matlab": {
            **_extract_matlab_dimensions(),
            "objective_value": _extract_matlab_objective() or MATLAB_OBJECTIVE,
            "baseline_lp": MATLAB_BASELINE_LP,
            "baseline_out": MATLAB_BASELINE_OUT,
        },
    }


def main() -> int:
    report = run_python_alignment()
    py = report["python"]
    ml = report["matlab"]

    print("=" * 70)
    print("OPN pcSecPichia alignment probe")
    print("=" * 70)
    print(f"{'metric':<32} {'python':>20} {'matlab':>20}")
    print("-" * 70)
    print(f"{'rows (eq+ub)':<32} {py['rows_total']:>20} {ml['rows']:>20}")
    print(f"{'columns':<32} {py['cols']:>20} {ml['cols']:>20}")
    print(f"{'objective (maximize Ex_glc_D)':<32} {py['objective_value']:>20.8g} {ml['objective_value']:>20.8g}")
    print(f"{'dilution_misfolding blocked':<32} {py['dilution_misfolding_blocked']:>20} {'(default true)':>20}")
    print("-" * 70)

    if py["objective_value"] is not None and ml["objective_value"] is not None:
        diff = abs(py["objective_value"] - ml["objective_value"])
        rel = diff / abs(ml["objective_value"]) if ml["objective_value"] else None
        print(f"absolute objective diff: {diff:.6e}")
        if rel is not None:
            print(f"relative objective diff: {rel:.4%}")

    rows_diff = py["rows_total"] - ml["rows"]
    cols_diff = py["cols"] - ml["cols"]
    print(f"row count diff: {rows_diff:+d}  (python - matlab)")
    print(f"column count diff: {cols_diff:+d}  (python - matlab)")
    print()
    print("constraint counts (python):")
    for key, value in py["constraint_counts"].items():
        print(f"  {key:<28} {value:>8}")

    out_path = repo_root() / "local_runs" / "pichia_python" / "alignment" / "opn_alignment_probe.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nreport written to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
