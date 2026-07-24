"""Stacked strain modifications for a "modified-strain re-solve" (ADR-004 #1 迭代候选).

The base secretion solve (`solve_secretion_capacity`) and the OE dose-response sweep both solve
the **wildtype** target strain: KO/OE only ever enter as single-candidate perturbations tested
*against* that wildtype baseline (the genome/reaction screens), never stacked into one solved
strain. To answer "given a strain I've *already* modified, what's the next bottleneck / next OE
lever" the modifications have to be applied to the model **simultaneously and then solved once** —
so the binding bottleneck actually shifts as the strain is iterated (松开 A 之后 B 顶上来).

This module is the one shared, pure place that turns a `StrainModifications` spec into a perturbed
`(model, secretory, combined)`, reusing the *exact* per-reaction OE mechanics of
`run_pcsec_oe_screen` (sec_ complex → secretory kcat ×factor; Mach_ complex → combined enzyme kcat
×factor; other reaction → upper bound ×factor) plus KO as a hard `(0, 0)` reaction bound. Both
`solve_secretion_capacity` and `run_oe_dose_response_sweep` take an **opt-in** `strain_modifications`
argument that defaults to no-op; an empty/None spec leaves the solve byte-identical (glucose
corrected_reference regression stays locked).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pcsec_pichia.probe import CobraModel, CombinedEnzymeData, SecretoryEnzymeData

_OE_FORMATION_SUFFIX = "_formation"


@dataclass(frozen=True)
class StrainModifications:
    """Reaction/complex-level modifications defining one already-modified strain.

    Complex-level by design (reaction_id is the OE target = the complex-formation reaction whose
    capacity ceiling is relaxed), matching the "下一步 OE 候选" readout. `oe_factor` is a single
    capacity multiplier shared by every OE'd complex (a relative capacity dose, not a measured
    expression level).
    """

    ko_reaction_ids: tuple[str, ...] = ()
    oe_reaction_ids: tuple[str, ...] = ()
    oe_factor: float = 2.0

    def is_empty(self) -> bool:
        return not self.ko_reaction_ids and not self.oe_reaction_ids


def apply_strain_modifications(
    model: CobraModel,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    modifications: StrainModifications,
) -> tuple[CobraModel, SecretoryEnzymeData, CombinedEnzymeData, tuple[dict[str, Any], ...], tuple[str, ...]]:
    """Apply stacked KO/OE modifications to a prepared model + enzyme data (pure; returns copies).

    Returns `(model, secretory, combined, applied, warnings)`. `applied` records each resolved
    change (reaction_id / kind / capacity_basis) so the caller can surface *what actually happened*;
    `warnings` flags unresolved reaction ids and KO/OE conflicts. Modifications that don't resolve
    are skipped (never silently faked). An empty spec returns the inputs unchanged.
    """

    applied: list[dict[str, Any]] = []
    warnings: list[str] = []
    if modifications.is_empty():
        return model, secretory, combined, (), ()

    factor = float(modifications.oe_factor)
    reaction_index = model.reaction_index
    bound_changes: dict[str, tuple[float | None, float | None]] = {}
    secretory_eff = secretory
    combined_eff = combined

    ko_ids = tuple(dict.fromkeys(str(r) for r in modifications.ko_reaction_ids if str(r)))
    oe_ids = tuple(dict.fromkeys(str(r) for r in modifications.oe_reaction_ids if str(r)))
    ko_set = set(ko_ids)

    for reaction_id in ko_ids:
        if reaction_id not in reaction_index:
            warnings.append(f"KO reaction not found in model, skipped: {reaction_id}")
            continue
        bound_changes[reaction_id] = (0.0, 0.0)
        applied.append({"reaction_id": reaction_id, "kind": "KO", "capacity_basis": "reaction_bounds_zero"})

    if oe_ids and (factor <= 0.0 or abs(factor - 1.0) < 1e-9):
        warnings.append(
            f"OE factor {factor} has no capacity effect (must be > 0 and != 1); OE modifications skipped."
        )
        oe_ids = ()

    for reaction_id in oe_ids:
        if reaction_id in ko_set:
            warnings.append(f"Reaction requested for both KO and OE; kept as KO, dropped OE: {reaction_id}")
            continue
        if reaction_id not in reaction_index:
            warnings.append(f"OE reaction not found in model, skipped: {reaction_id}")
            continue
        complex_id = reaction_id.replace(_OE_FORMATION_SUFFIX, "")
        if complex_id.startswith("sec_"):
            secretory_eff = secretory_eff.with_complex_kcat_multiplier(complex_id, factor)
            applied.append(
                {"reaction_id": reaction_id, "kind": "OE", "capacity_basis": "secretory_complex_kcat_multiplier", "factor": factor}
            )
        elif complex_id.startswith("Mach_"):
            combined_eff = combined_eff.with_enzyme_kcat_multiplier(complex_id, factor)
            applied.append(
                {"reaction_id": reaction_id, "kind": "OE", "capacity_basis": "machine_complex_kcat_multiplier", "factor": factor}
            )
        else:
            old_upper = float(model.ub[reaction_index[reaction_id]])
            new_upper = old_upper * factor if old_upper > 0 else factor
            bound_changes[reaction_id] = (None, new_upper)
            applied.append(
                {"reaction_id": reaction_id, "kind": "OE", "capacity_basis": "reaction_upper_bound", "factor": factor}
            )

    model_eff = model.with_bounds(bound_changes) if bound_changes else model
    return model_eff, secretory_eff, combined_eff, tuple(applied), tuple(warnings)


__all__ = ["StrainModifications", "apply_strain_modifications"]
