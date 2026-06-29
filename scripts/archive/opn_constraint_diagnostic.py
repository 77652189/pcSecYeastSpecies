"""Diagnostic: dump probe constraint details and compare with MATLAB LP.

Reads the MATLAB baseline LP file and the probe's constraint matrices, then
reports per-section coefficient statistics so we can pinpoint where the 81%
objective gap comes from.
"""

from __future__ import annotations

import re
from collections import Counter
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
    metabolic_coupling_rows,
    prepare_glucose_model,
    proteasome_rows,
    protein_mass_rows,
    ribosome_assembly_rows,
    repo_root,
    secretory_coupling_rows,
)

MATLAB_LP = (
    "local_runs/OPN_PPA_glc_smoke/"
    "Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_"
    "misfolddefault_noMisfoldEq_noRiboEq_PP.lp"
)


def parse_matlab_lp(path: Path) -> dict[str, object]:
    """Parse the MATLAB LP file and extract per-section statistics."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    # Find section boundaries
    subject_to_idx = None
    bounds_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Subject To":
            subject_to_idx = i
        elif line.strip() == "Bounds":
            bounds_idx = i
            break

    constraint_lines = lines[subject_to_idx + 1 : bounds_idx] if subject_to_idx and bounds_idx else []

    # Categorize constraints by prefix
    sections: dict[str, list[str]] = {
        "stoichiometric_C": [],
        "metabolic_CM": [],
        "secretory_CM": [],
        "protein_mass_CM": [],
        "mito_Cmito": [],
        "proteasome_CM": [],
        "ribo_assembly_CM": [],
        "other": [],
    }

    for line in constraint_lines:
        line = line.strip()
        if not line:
            continue
        # Constraints can span multiple lines; we only look at lines starting with a label
        m = re.match(r"^(C\d+|CM\d+|Cmito)\s*:", line)
        if not m:
            continue
        label = m.group(1)
        if label.startswith("Cmito"):
            sections["mito_Cmito"].append(line)
        elif label.startswith("CM"):
            # CM constraints are numbered sequentially. We need to figure out which section.
            # MATLAB writeLPGlc order: stoich (C1..Cn), met (CM1..CM_nmet), sec (CM_nmet+1..),
            # protein_mass (CM_k+1, CM_k+2), mito (Cmito), proteasome (CM_k+3), ribo_assembly (CM_k+6)
            sections["metabolic_CM"].append(line)  # placeholder, will refine below
        elif label.startswith("C"):
            sections["stoichiometric_C"].append(line)
        else:
            sections["other"].append(line)

    # Refine CM categorization by parsing the constraint structure
    # Metabolic: "X%d - %.15f X%d = 0" (2 terms, equality)
    # Secretory: multiple terms + "- coef X%d = 0" (equality)
    # Protein mass: "sum = %.15f" (equality, no X_formation subtraction pattern)
    # Proteasome: "sum - coef X%d = 0"
    # Ribosome assembly: "X%d - coef X%d = 0" (2 terms, but specific indices)

    # Count stoichiometric rows
    n_stoich = len(sections["stoichiometric_C"])

    # Parse all CM constraints
    cm_lines = sections["metabolic_CM"]
    n_cm = len(cm_lines)

    # The last few CM constraints are special:
    # protein_mass (1 eq), unmodeled_ER (1 eq), proteasome (1 eq), ribo_assembly (1 eq)
    # = 4 special CM constraints at the end
    # But mito is "Cmito" not "CM"

    # Let's just report counts
    return {
        "total_constraints": len(constraint_lines),
        "stoichiometric_C": n_stoich,
        "CM_total": n_cm,
        "Cmito": len(sections["mito_Cmito"]),
        "sample_stoich": sections["stoichiometric_C"][:2],
        "sample_cm_first": cm_lines[:2] if cm_lines else [],
        "sample_cm_last": cm_lines[-6:] if len(cm_lines) >= 6 else cm_lines,
        "sample_mito": sections["mito_Cmito"][:1],
    }


def dump_probe_constraints() -> dict[str, object]:
    """Run the probe and dump per-section constraint statistics."""
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

    # Per-section row counts and coefficient stats
    met_rows, met_rhs = metabolic_coupling_rows(fixed_model, metabolic, mu=0.10)
    sec_rows, sec_rhs = secretory_coupling_rows(fixed_model, target_secretory, mu=0.10)
    protein_rows, protein_rhs = protein_mass_rows(
        fixed_model, target_combined, mu=0.10,
        total_protein_content=0.37, unmodeled_er_protein_fraction=0.040,
    )
    mito_rows, mito_rhs = [], []
    from pcsec_pichia.probe import mitochondrial_rows
    mito_rows, mito_rhs = mitochondrial_rows(fixed_model, target_combined, mu=0.10, mitochondrial_protein_fraction=0.05)
    prot_rows, prot_rhs = proteasome_rows(fixed_model, target_combined, mu=0.10)
    ribo_rows, ribo_rhs = ribosome_assembly_rows(fixed_model, target_combined, mu=0.10)

    def row_stats(rows: list, rhs: list, name: str) -> dict[str, object]:
        if not rows:
            return {"name": name, "count": 0}
        stacked = sparse.vstack(rows, format="csr")
        return {
            "name": name,
            "count": len(rows),
            "nnz": int(stacked.nnz),
            "coef_abs_min": float(np.min(np.abs(stacked.data))),
            "coef_abs_max": float(np.max(np.abs(stacked.data))),
            "coef_abs_mean": float(np.mean(np.abs(stacked.data))),
            "rhs_abs_max": float(np.max(np.abs(rhs))) if rhs else 0.0,
        }

    return {
        "metabolic": row_stats(met_rows, met_rhs, "metabolic"),
        "secretory": row_stats(sec_rows, sec_rhs, "secretory"),
        "protein_mass": row_stats(protein_rows, protein_rhs, "protein_mass"),
        "mitochondrial": row_stats(mito_rows, mito_rhs, "mitochondrial"),
        "proteasome": row_stats(prot_rows, prot_rhs, "proteasome"),
        "ribosome_assembly": row_stats(ribo_rows, ribo_rhs, "ribosome_assembly"),
        "metabolic_enzyme_count": len(metabolic.enzymes),
        "metabolic_rows_generated": len(met_rows),
        "metabolic_rows_skipped": len(metabolic.enzymes) - len(met_rows),
        "secretory_complex_count": len(secretory.complexes),
        "secretory_unique_complexes": len(secretory.unique_complex_entries()),
        "secretory_rows_generated": len(sec_rows),
    }


def main() -> int:
    root = repo_root()
    ml = parse_matlab_lp(root / MATLAB_LP)
    py = dump_probe_constraints()

    print("=" * 72)
    print("Constraint diagnostic: probe vs MATLAB LP")
    print("=" * 72)
    print()
    print("MATLAB LP sections:")
    print(f"  stoichiometric (C*):  {ml['stoichiometric_C']}")
    print(f"  CM total:             {ml['CM_total']}")
    print(f"  Cmito:                {ml['Cmito']}")
    print(f"  total:                {ml['total_constraints']}")
    print()
    print("Probe constraint sections:")
    for key in ("metabolic", "secretory", "protein_mass", "mitochondrial", "proteasome", "ribosome_assembly"):
        s = py[key]
        print(f"  {s['name']:<20} count={s['count']:>5}  nnz={s.get('nnz',0):>8}  "
              f"coef_abs=[{s.get('coef_abs_min',0):.4g}, {s.get('coef_abs_max',0):.4g}]  "
              f"rhs_abs_max={s.get('rhs_abs_max',0):.4g}")
    print()
    print(f"metabolic enzymes loaded:   {py['metabolic_enzyme_count']}")
    print(f"metabolic rows generated:   {py['metabolic_rows_generated']}")
    print(f"metabolic rows SKIPPED:     {py['metabolic_rows_skipped']}")
    print(f"secretory complexes total:  {py['secretory_complex_count']}")
    print(f"secretory unique complexes: {py['secretory_unique_complexes']}")
    print(f"secretory rows generated:   {py['secretory_rows_generated']}")
    print()
    print("MATLAB sample stoichiometric:")
    for s in ml["sample_stoich"][:1]:
        print(f"  {s[:120]}...")
    print("MATLAB sample CM (first):")
    for s in ml["sample_cm_first"][:1]:
        print(f"  {s[:120]}...")
    print("MATLAB sample CM (last 6 — protein_mass, mito, proteasome, ribo):")
    for s in ml["sample_cm_last"]:
        print(f"  {s[:120]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
