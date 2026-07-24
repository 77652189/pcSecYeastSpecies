"""Engine orchestration for the "modified-strain → 下一步 OE 候选" two-pass (ADR-004 #1 迭代候选).

This is the engine half of #1: given an *already-modified* strain (stacked KO/OE), re-solve it and
find its **new** binding bottleneck, then quantify the bounded OE dose-response of those bottleneck
complexes — so the readout shifts as the strain is iterated instead of always returning the same
wildtype #1 lever.

Two passes, both on the modified strain (see `strain_modifications`):
1. `solve_secretion_capacity(..., strain_modifications=mods)` → LP attribution → `oe_actionable_bottlenecks`
   (binding upper-bound complexes = OE-actionable leads for *this* strain).
2. `run_oe_dose_response_sweep(reactions=top-N bottlenecks, strain_modifications=mods)` → per-complex
   shape (linear vs saturating, how much / where it saturates).

It returns the two raw pieces (`oe_actionable_bottlenecks` + `dose_response`) that the app-layer
assembler `app/services/per_strain_oe_candidates.build_next_oe_candidates_readout` (C1) ranks; the
ranking/assembly stays app-side and pure. Relative signal only, complex-level, no absolute titer.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from pcsec_pichia.analysis import (
    analyze_target_protein_lp_attribution,
    classify_oe_dose_response_sweep,
    summarize_oe_dose_response_shape,
    summarize_protein_lp_attribution,
)
from pcsec_pichia.constraints import build_pcsec_constraints
from pcsec_pichia.engines.base import PichiaSimulationRequest
from pcsec_pichia.loading import load_pcsec_pichia_inputs
from pcsec_pichia.pipeline import _resolve_target
from pcsec_pichia.screens import DEFAULT_OE_DOSE_RESPONSE_FACTORS, run_oe_dose_response_sweep
from pcsec_pichia.secretion_plan import build_secretion_plan
from pcsec_pichia.simulation import solve_secretion_capacity
from pcsec_pichia.strain_modifications import StrainModifications

DEFAULT_NEXT_OE_TOP_N = 6


def analyze_next_oe_candidates(
    *,
    target_id: str,
    ko_reaction_ids: Sequence[str] = (),
    oe_reaction_ids: Sequence[str] = (),
    oe_factor: float = 2.0,
    mu: float = 0.10,
    media_type: int = 4,
    carbon_source_id: str = "glucose",
    compatibility_mode: str = "corrected",
    enable_ribosome_translation_constraint: bool = False,
    enable_misfolding_constraint: bool = False,
    target_input: Any | None = None,
    leader_candidate: Any | None = None,
    top_n: int = DEFAULT_NEXT_OE_TOP_N,
    dose_response_factors: Sequence[float] = (),
    root: Path | None = None,
) -> dict[str, Any]:
    """Two-pass modified-strain bottleneck + bounded dose-response for the next-OE-candidate readout.

    `ko_reaction_ids` / `oe_reaction_ids` define the already-applied strain modifications (complex/
    reaction level; `oe_factor` is the shared OE capacity multiplier). Returns a dict with the raw
    `oe_actionable_bottlenecks` and `dose_response` (`shapes_by_reaction`) for C1 to rank, plus
    solve status, applied modifications, and honest warnings. Never fabricates a candidate: if the
    modified solve is infeasible, bottlenecks/dose-response come back empty with the reason.
    """

    root = root or Path(__file__).resolve().parents[3]
    ko_ids = tuple(str(r).strip() for r in ko_reaction_ids if str(r).strip())
    oe_ids = tuple(str(r).strip() for r in oe_reaction_ids if str(r).strip())

    request = PichiaSimulationRequest(
        target_id=target_id,
        candidate_id=target_id,
        target_input=target_input,
        leader_candidate=leader_candidate,
        mu=mu,
        media_type=media_type,
        carbon_source_id=carbon_source_id,
        compatibility_mode=compatibility_mode,
    )
    inputs = load_pcsec_pichia_inputs(
        root,
        media_type=media_type,
        compatibility_mode=compatibility_mode,
        carbon_source_id=carbon_source_id,
    )
    target = _resolve_target(request, root)
    plan = build_secretion_plan(target)
    constraint_result = build_pcsec_constraints(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        mu=mu,
        write_ribosome_translation_constraint=enable_ribosome_translation_constraint,
        write_misfolding_constraints=enable_misfolding_constraint,
    )

    # Honest validity check against the base model index: warn (don't fake) for modifications that
    # cannot resolve to a model reaction. Target augmentation only *adds* reactions, so a base-model
    # hit is also present in the solved (augmented) model.
    reaction_index = inputs.prepared_model.reaction_index
    modification_warnings = [
        f"KO reaction not found in model, skipped: {reaction_id}"
        for reaction_id in ko_ids
        if reaction_id not in reaction_index
    ] + [
        f"OE reaction not found in model, skipped: {reaction_id}"
        for reaction_id in oe_ids
        if reaction_id not in reaction_index
    ]

    modifications = StrainModifications(ko_reaction_ids=ko_ids, oe_reaction_ids=oe_ids, oe_factor=float(oe_factor))
    strain_modifications = None if modifications.is_empty() else modifications

    simulation = solve_secretion_capacity(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_rate=mu,
        write_ribosome_translation_constraint=enable_ribosome_translation_constraint,
        write_misfolding_constraints=enable_misfolding_constraint,
        strain_modifications=strain_modifications,
    )
    reaction_ids = tuple(inputs.prepared_model.rxns)
    lp_attribution = analyze_target_protein_lp_attribution(
        target,
        plan,
        constraint_result.constraint_counts,
        simulation,
        reaction_ids=reaction_ids,
    )
    summary = summarize_protein_lp_attribution(lp_attribution)
    bottlenecks = list(summary.get("oe_actionable_bottlenecks") or [])

    top_reactions = [
        str(entry["reaction_id"])
        for entry in bottlenecks[: max(0, int(top_n))]
        if isinstance(entry, dict) and entry.get("reaction_id")
    ]

    dose_response: dict[str, Any] | None = None
    if simulation.success and top_reactions:
        sweep = run_oe_dose_response_sweep(
            inputs.prepared_model,
            target,
            inputs.amino_acids,
            inputs.metabolic,
            inputs.secretory,
            inputs.combined,
            reactions=top_reactions,
            factors=tuple(dose_response_factors) or DEFAULT_OE_DOSE_RESPONSE_FACTORS,
            growth_rate=mu,
            write_ribosome_translation_constraint=enable_ribosome_translation_constraint,
            write_misfolding_constraints=enable_misfolding_constraint,
            strain_modifications=strain_modifications,
        )
        if sweep.success:
            shape_dicts = [
                summarize_oe_dose_response_shape(shape)
                for shape in classify_oe_dose_response_sweep(sweep.reaction_points, sweep.baseline_objective)
            ]
            dose_response = {
                "baseline_objective": sweep.baseline_objective,
                "tested_factors": list(sweep.tested_factors),
                "shapes_by_reaction": {
                    str(shape.get("reaction_id")): shape
                    for shape in shape_dicts
                    if isinstance(shape, dict) and shape.get("reaction_id")
                },
                "warnings": list(sweep.warnings),
            }

    return {
        "target_id": target.target_id,
        "modified_solve_success": bool(simulation.success),
        "modified_objective_value": simulation.objective_value,
        "lp_attribution_status": summary.get("result_status"),
        "oe_actionable_bottlenecks": bottlenecks,
        "floor_constraints_not_oe_addressable": list(summary.get("floor_constraints_not_oe_addressable") or []),
        "dose_response": dose_response,
        "applied_modifications": {
            "ko_reaction_ids": list(ko_ids),
            "oe_reaction_ids": list(oe_ids),
            "oe_factor": float(oe_factor),
        },
        "modification_warnings": modification_warnings,
        "carbon_source_id": inputs.carbon_source_id,
        "medium_condition_id": inputs.medium_condition_id,
        "mu": float(mu),
        "top_n": int(top_n),
    }


__all__ = ["DEFAULT_NEXT_OE_TOP_N", "analyze_next_oe_candidates"]
