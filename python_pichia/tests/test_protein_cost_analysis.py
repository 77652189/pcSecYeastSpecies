from __future__ import annotations

from pathlib import Path

from pcsec_pichia.analysis import (
    ProteinLpAttributionResult,
    _lp_reaction_process,
    analyze_target_protein_lp_attribution,
    classify_oe_dose_response_shape,
    compare_ranking_robustness,
    compare_solver_robustness,
    prioritize_value_of_information,
    summarize_oe_dose_response_shape,
    summarize_protein_cost_slope_compatibility,
    summarize_protein_lp_attribution,
    summarize_ranking_robustness,
    summarize_solver_robustness,
    summarize_value_of_information,
)
from pcsec_pichia.simulation import ProteinCostSlopeCompatibilityResult
from pcsec_pichia.secretion_plan import build_secretion_plan
from pcsec_pichia.simulation import SecretionSimulationResult
from pcsec_pichia.targets import load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


def _builtin(target_id: str):
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}[target_id]


def test_lp_attribution_summarizes_top_n_without_full_marginal_arrays() -> None:
    target = _builtin("OPN_ALPHA_FULL_PROJECT")
    plan = build_secretion_plan(target)
    simulation = SecretionSimulationResult(
        success=True,
        target_id=target.target_id,
        objective_value=0.006,
        growth_rate=0.10,
        secretion_flux=0.006,
        status="0",
        message="ok",
        constraint_counts={
            "stoichiometric": 2,
            "secretory_coupling": 2,
            "protein_mass": 1,
            "ribosome_translation": 1,
            "misfolding": 1,
            "mitochondrial": 1,
            "eq_total": 7,
            "ub_total": 1,
        },
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
        exchange_reaction_id="r_OPN_ALPHA_FULL_PROJECT_exchange",
        build_status="supported",
        lp_sensitivity={
            "eq_marginals": (0.0, 0.2, -0.5, 0.1, 0.7, -0.9, 0.3),
            "ub_marginals": (0.4,),
            "lower_marginals": (0.0, -0.8, 0.0, 0.0),
            "upper_marginals": (0.0, 0.0, 0.6, 0.0),
        },
        key_fluxes={
            "BIOMASS": 0.10,
            "r_OPN_ALPHA_FULL_PROJECT_exchange": 0.006,
        },
    )

    result = analyze_target_protein_lp_attribution(
        target,
        plan,
        simulation.constraint_counts,
        simulation,
        reaction_ids=("BIOMASS", "r_OPN_ALPHA_FULL_PROJECT_exchange", "Ex_glc_D", "OTHER"),
        top_n=3,
    )
    summary = summarize_protein_lp_attribution(result)

    assert summary["result_status"] == "draft_lp_sensitivity"
    assert len(summary["top_constraint_marginals"]) == 3
    assert summary["top_constraint_marginals"][0]["block"] == "ribosome_translation"
    assert len(summary["top_bound_marginals"]) == 2
    assert summary["top_bound_marginals"][0]["reaction_id"] == "r_OPN_ALPHA_FULL_PROJECT_exchange"
    assert summary["active_bound_counts"]["total_bound_marginal_nonzero"] == 2
    assert "eq_marginals" not in summary
    assert "lower_marginals" not in summary


def test_oe_actionable_bottlenecks_exclude_lower_bound_floor() -> None:
    # R1 (ADR-004): the biggest bound marginal here is a LOWER bound (floor, abs 0.8) that
    # OE cannot relax; the OE-addressable one is a smaller UPPER bound (abs 0.6). The derived
    # oe_actionable list must contain only the upper bound and never the larger lower one -
    # this is the PDI1-alone / ribosome false-positive guard made programmatic.
    target = _builtin("OPN_ALPHA_FULL_PROJECT")
    plan = build_secretion_plan(target)
    simulation = SecretionSimulationResult(
        success=True,
        target_id=target.target_id,
        objective_value=0.006,
        growth_rate=0.10,
        secretion_flux=0.006,
        status="0",
        message="ok",
        constraint_counts={"stoichiometric": 7, "eq_total": 7, "ub_total": 1},
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
        exchange_reaction_id="r_OPN_ALPHA_FULL_PROJECT_exchange",
        build_status="supported",
        lp_sensitivity={
            "eq_marginals": (0.0,) * 7,
            "ub_marginals": (0.0,),
            "lower_marginals": (0.0, -0.8, 0.0, 0.0),
            "upper_marginals": (0.0, 0.0, 0.6, 0.0),
        },
        key_fluxes={"BIOMASS": 0.10},
    )
    reaction_ids = ("BIOMASS", "r_OPN_ALPHA_FULL_PROJECT_exchange", "Ex_glc_D", "OTHER")

    result = analyze_target_protein_lp_attribution(
        target, plan, simulation.constraint_counts, simulation, reaction_ids=reaction_ids, top_n=5
    )

    oe_ids = [row["reaction_id"] for row in result.oe_actionable_bottlenecks]
    floor_ids = [row["reaction_id"] for row in result.floor_constraints_not_oe_addressable]
    assert oe_ids == ["Ex_glc_D"]
    assert all(row["bound_type"] == "upper" and row["oe_actionable"] is True for row in result.oe_actionable_bottlenecks)
    # the larger lower-bound floor is segregated, never promoted into the OE-actionable list
    assert "r_OPN_ALPHA_FULL_PROJECT_exchange" in floor_ids
    assert "r_OPN_ALPHA_FULL_PROJECT_exchange" not in oe_ids
    assert all(row["bound_type"] == "lower" for row in result.floor_constraints_not_oe_addressable)
    summary = summarize_protein_lp_attribution(result)
    assert summary["oe_actionable_bottlenecks"] == list(result.oe_actionable_bottlenecks)
    assert summary["floor_constraints_not_oe_addressable"] == list(result.floor_constraints_not_oe_addressable)


def _attr_with_top_bottleneck(reaction_id: str, block: str) -> ProteinLpAttributionResult:
    return ProteinLpAttributionResult(
        target_id="hLF",
        result_status="draft_lp_sensitivity",
        objective_evidence={},
        dominant_constraint_blocks=({"block": block, "sum_abs_marginal": 1.0},),
        top_constraint_marginals=(),
        top_bound_marginals=(),
        target_related_fluxes=(),
        active_bound_counts={},
        warnings=(),
        oe_actionable_bottlenecks=({"reaction_id": reaction_id, "bound_type": "upper", "abs_marginal": 1.0},),
    )


def test_compare_solver_robustness_flags_solver_dependent_bottleneck() -> None:
    stable = compare_solver_robustness(
        "hLF",
        {
            "highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
            "highs-ds": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
        },
    )
    assert stable.classification == "ranking-insensitive-to-solver"
    assert stable.top_bottleneck_stable is True

    flipped = compare_solver_robustness(
        "hLF",
        {
            "highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
            "highs-ipm": _attr_with_top_bottleneck("sec_Y_complex_formation", "secretory_coupling"),
        },
    )
    assert flipped.classification == "ranking-sensitive-to-solver"
    assert flipped.top_bottleneck_stable is False
    assert summarize_solver_robustness(flipped)["classification"] == "ranking-sensitive-to-solver"

    inconclusive = compare_solver_robustness(
        "hLF", {"highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling")}
    )
    assert inconclusive.classification == "inconclusive"

    # a re-solve that errored is fed in as an unavailable attribution: even with a good
    # method present, the verdict degrades to inconclusive rather than pretending stability.
    errored = ProteinLpAttributionResult(
        target_id="hLF",
        result_status="draft_lp_sensitivity_unavailable",
        objective_evidence={},
        dominant_constraint_blocks=(),
        top_constraint_marginals=(),
        top_bound_marginals=(),
        target_related_fluxes=(),
        active_bound_counts={},
        warnings=(),
    )
    with_error = compare_solver_robustness(
        "hLF",
        {
            "highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
            "highs-ipm": errored,
        },
    )
    assert with_error.classification == "inconclusive"
    assert any(row["result_status"] == "draft_lp_sensitivity_unavailable" for row in with_error.per_method)


def test_classify_oe_dose_response_shape_distinguishes_curve_shapes() -> None:
    # Concave: most of the gain arrives early, then plateaus -> saturating (modest OE is enough).
    saturating = classify_oe_dose_response_shape(
        "sec_X_complex_formation", [(1.0, 1.0), (2.0, 1.08), (4.0, 1.11), (8.0, 1.12)]
    )
    assert saturating.shape == "saturating"
    assert saturating.normalized_auc is not None and saturating.normalized_auc > 0.6
    assert saturating.monotonic_non_decreasing is True
    # relative gain is expressed against the baseline, never as an absolute capacity
    assert abs(saturating.max_relative_gain - 0.12) < 1e-9
    assert saturating.half_gain_factor is not None and saturating.half_gain_factor < 2.0

    # Straight line: gain grows proportionally, no plateau reached -> linear (push expression higher).
    linear = classify_oe_dose_response_shape(
        "sec_X_complex_formation", [(1.0, 1.0), (2.0, 1.02), (4.0, 1.06), (8.0, 1.14)]
    )
    assert linear.shape == "linear"

    # Convex: little until a higher dose -> threshold (a minimum dose is required).
    threshold = classify_oe_dose_response_shape(
        "sec_X_complex_formation", [(1.0, 1.0), (2.0, 1.005), (4.0, 1.02), (8.0, 1.14)]
    )
    assert threshold.shape == "threshold"
    assert threshold.half_gain_factor is not None and threshold.half_gain_factor > 4.0


def test_classify_oe_dose_response_shape_flat_and_artifact_guards() -> None:
    # All gains under the noise floor -> no detectable response, no AUC computed.
    flat = classify_oe_dose_response_shape(
        "sec_flat", [(1.0, 1.0), (2.0, 1.0000005), (4.0, 1.0000008), (8.0, 1.0000009)]
    )
    assert flat.shape == "flat_no_response"
    assert flat.normalized_auc is None

    # OE relaxes a ceiling (enlarges the feasible region), so a beyond-noise decrease cannot be
    # a real dose-response: it must be reported as a numerical artifact, not a shape to trust.
    artifact = classify_oe_dose_response_shape(
        "sec_degenerate", [(1.0, 1.0), (2.0, 1.1), (4.0, 1.05), (8.0, 1.2)]
    )
    assert artifact.shape == "non_monotonic_numerical_artifact"
    assert artifact.monotonic_non_decreasing is False
    assert artifact.normalized_auc is None

    # Fewer than two usable points cannot yield a shape.
    insufficient = classify_oe_dose_response_shape("sec_X", [(2.0, 1.1)])
    assert insufficient.shape == "insufficient_points"
    assert insufficient.result_status == "draft_oe_dose_response_insufficient"

    # A failed solve (None objective) is dropped, not treated as zero.
    with_failed = classify_oe_dose_response_shape(
        "sec_X", [(1.0, 1.0), (2.0, None), (4.0, 1.06), (8.0, 1.14)]
    )
    assert with_failed.shape == "linear"
    assert [p["factor"] for p in with_failed.point_deltas] == [1.0, 4.0, 8.0]

    # The summary carries the relative-signal caveats and never fabricates an absolute capacity.
    summary = summarize_oe_dose_response_shape(with_failed)
    assert summary["shape"] == "linear"
    assert any("RELATIVE signal" in w for w in summary["warnings"])
    assert "capacity" not in {k for k in summary}  # no absolute-capacity field is emitted


def test_compare_ranking_robustness_covers_bandwidth_and_solver_classes() -> None:
    # R3 (ADR-004): a candidate ranking that survives capacity-bandwidth and solver perturbations is
    # a trustworthy relative signal; one that flips is an artifact. Absolute status stays unavailable
    # in every class, and no 'capacity-robust'-style label is used.
    baseline = ["A", "B", "C", "D"]

    # class 1: bandwidth-stable -> ranking-insensitive-to-capacity
    stable = compare_ranking_robustness(
        "hLF",
        baseline,
        capacity_bandwidth_rankings={"x0.7": ["A", "B", "C", "Z"], "x1.3": ["A", "B", "C", "Q"]},
    )
    assert stable.capacity_classification == "ranking-insensitive-to-capacity"
    assert stable.absolute_status == "unavailable"

    # class 2: bandwidth-flip -> ranking-sensitive-to-capacity
    flip = compare_ranking_robustness(
        "hLF",
        baseline,
        capacity_bandwidth_rankings={"x0.7": ["A", "B", "C"], "x1.3": ["B", "A", "C"]},
    )
    assert flip.capacity_classification == "ranking-sensitive-to-capacity"
    assert flip.absolute_status == "unavailable"

    # class 3: solver-flip -> ranking-sensitive-to-solver
    solver_flip = compare_ranking_robustness("hLF", baseline, solver_rankings={"highs-ipm": ["B", "A", "C"]})
    assert solver_flip.solver_classification == "ranking-sensitive-to-solver"
    assert solver_flip.absolute_status == "unavailable"

    # invariant: no 'capacity-robust'-style naming in any label, absolute stays unavailable everywhere
    labels = {stable.capacity_classification, flip.capacity_classification, solver_flip.solver_classification}
    assert not any("capacity-robust" in label for label in labels)
    for result in (stable, flip, solver_flip):
        assert summarize_ranking_robustness(result)["absolute_status"] == "unavailable"

    # no perturbations supplied -> inconclusive, not a false verdict
    empty = compare_ranking_robustness("hLF", baseline)
    assert empty.capacity_classification == "inconclusive"
    assert empty.solver_classification == "inconclusive"


def test_prioritize_value_of_information_ranks_measurements_without_predicting_yield() -> None:
    # R4 (ADR-004): flag where the model cannot confidently order candidates and prioritize the
    # minimal measurement that resolves it; never predict an outcome or promote a candidate.
    near_tie = prioritize_value_of_information("hLF", [("A", 1.0), ("B", 0.98), ("C", 0.5)])
    assert near_tie.has_actionable_ambiguity is True
    top_item = near_tie.information_items[0]
    assert top_item["ambiguity_kind"] == "near_tie"
    assert top_item["candidates"] == ["A", "B"]
    assert top_item["resolves_top_of_ranking"] is True
    assert "target-specific" in top_item["suggested_measurement"]
    # the ranking (with scores) is carried so the UI can chart near-ties as near-equal bars
    assert [row["candidate_id"] for row in near_tie.ranked_candidates] == ["A", "B", "C"]

    # clearly separated scores -> no actionable ambiguity
    separated = prioritize_value_of_information("hLF", [("A", 1.0), ("B", 0.5), ("C", 0.2)])
    assert separated.has_actionable_ambiguity is False

    # an R3 ranking flip becomes the highest-priority information item (fed straight from R3)
    robustness = compare_ranking_robustness(
        "hLF", ["A", "B", "C"], capacity_bandwidth_rankings={"x1.3": ["B", "A", "C"]}
    )
    with_flip = prioritize_value_of_information(
        "hLF", [("A", 1.0), ("B", 0.5), ("C", 0.2)], ranking_robustness=robustness
    )
    assert with_flip.information_items[0]["ambiguity_kind"] == "capacity_flip"
    assert with_flip.information_items[0]["priority_rank"] == 1

    # invariants: no absolute yield prediction / promotion; absolute status unavailable
    summary = summarize_value_of_information(with_flip)
    assert summary["absolute_status"] == "unavailable"
    assert any("does not predict" in w for w in summary["warnings"])
    assert any("never promotes" in w for w in summary["warnings"])


def test_lp_attribution_handles_missing_sensitivity_without_crashing() -> None:
    target = _builtin("hLF")
    plan = build_secretion_plan(target)
    simulation = SecretionSimulationResult(
        success=False,
        target_id=target.target_id,
        objective_value=None,
        growth_rate=0.10,
        secretion_flux=None,
        status="2",
        message="infeasible",
        constraint_counts={},
        result_status="draft",
        target_parameter_status="draft_matlab_alignment_pending",
        matlab_alignment_status="pending",
        exchange_reaction_id=None,
        build_status="failed",
    )

    result = analyze_target_protein_lp_attribution(target, plan, {}, simulation)

    assert result.result_status == "draft_lp_sensitivity_unavailable"
    assert result.top_constraint_marginals == ()
    assert any("unavailable" in warning for warning in result.warnings)


def test_cost_slope_compatibility_summary_is_disabled_by_default_payload() -> None:
    summary = summarize_protein_cost_slope_compatibility(None)

    assert summary["enabled"] is False
    assert summary["result_status"] == "disabled"


def test_cost_slope_compatibility_summary_keeps_matlab_style_definition() -> None:
    result = ProteinCostSlopeCompatibilityResult(
        target_id="OPN_ALPHA_FULL_PROJECT",
        enabled=True,
        success=True,
        growth_rates=(0.05,),
        secretion_ratios=(5e-7, 1e-6),
        rows=(
            {
                "mu": 0.05,
                "target_exchange_ratio": 5e-7,
                "success": True,
                "glucose_cost": 1.0,
            },
        ),
        glucose_cost_slopes=(
            {
                "mu": 0.05,
                "cost_key": "glucose_cost",
                "success": True,
                "slope": 10.0,
                "point_count": 2,
                "status": "slope_estimated",
            },
        ),
        ribosome_cost_slopes=(),
        result_status="draft_matlab_compatible_cost_slope",
        warnings=("draft",),
    )

    summary = summarize_protein_cost_slope_compatibility(result)

    assert summary["enabled"] is True
    assert summary["result_status"] == "draft_matlab_compatible_cost_slope"
    assert summary["medium_compatibility_mode"] == "corrected"
    assert summary["medium_bound_overrides"] == []
    assert summary["comparison_scope"]["not_default_pipeline"] is True
    assert summary["comparison_scope"]["medium_compatibility"] == "corrected"
    assert summary["comparison_scope"]["current_default_definition"] == (
        "fixed growth rate, corrected medium, maximize target secretion flux"
    )
    assert summary["comparison_scope"]["ratio_policy"] == "explicit_absolute_ratios"
    assert summary["secretion_ratio_policy"] == "explicit_absolute_ratios"
    assert "fixed target exchange" in summary["comparison_scope"]["definition"]


def test_lp_reaction_process_classifies_host_sec_complexes_by_gene_token() -> None:
    # Root fix: the host sec_* secretory-machine complexes used to fall through to "unknown"; the
    # engine now delegates their classification to the shared gene-token classifier so the payload
    # carries a real process label. These are host-shared reactions (not this target's own handles
    # and not embedding its protein id), so they exercise the classifier, not the target-aware prefix.
    target = _builtin("hLF")
    plan = build_secretion_plan(target)
    assert "sec_Pdi1p_complex_formation" not in plan.reaction_ids  # host-shared, not a target handle

    def process(reaction_id: str) -> str:
        return _lp_reaction_process(reaction_id, target, plan)

    # High-confidence, one representative per process layer (verified against the yeast secretory
    # pathway): the dominant hLF bottleneck sec_Pdi1p -> disulfide_folding is the headline case.
    assert process("sec_Pdi1p_complex_formation") == "disulfide_folding"
    assert process("sec_Kar2p_complex_formation") == "chaperone_folding"
    assert process("sec_OSTC_complex_formation") == "n_glycan_processing"
    assert process("sec_Pmt_complex_formation") == "o_glycan_processing"
    assert process("sec_SEC61SEC63C_complex_formation") == "er_translocation"
    assert process("Mach_Ribosome_complex_formation") == "ribosome"
    # Conservative bucket: sec_* machinery whose gene token is genuinely ambiguous (GPI-anchor
    # remodeling here) stays in the generic secretory_capacity layer rather than a wrong guess.
    assert process("sec_Bst1p_complex_formation") == "secretory_capacity"


def test_lp_reaction_process_preserves_target_aware_cases() -> None:
    # The target-aware cases are resolved before delegation and must not regress: they are specific
    # to this attribution and have no general-classifier equivalent.
    target = _builtin("hLF")
    plan = build_secretion_plan(target)
    own_reaction = next(iter(plan.reaction_ids))

    assert _lp_reaction_process("", target, plan) == "unknown"
    assert _lp_reaction_process(own_reaction, target, plan) == "target_secretory_reaction"
    assert _lp_reaction_process(f"r_{target.protein_id}_exchange", target, plan) == "target_exchange"
    assert _lp_reaction_process(f"{target.protein_id}_exchange", target, plan) == "target_exchange"
    # a reaction whose id embeds the target's protein id (but is neither a plan handle nor the
    # exchange) is target_related, never leaked into a shared secretory-process bucket
    assert _lp_reaction_process(f"probe_{target.protein_id}_custom_reaction", target, plan) == "target_related"


def test_lp_attribution_payload_tags_sec_complexes_not_unknown() -> None:
    # End-to-end (the "payload is correct for all consumers" goal): with sec_* complexes as binding
    # bounds, the LP-attribution payload's secretory_process must carry real process labels, and the
    # relative-signal invariant (lower-bound floors segregated from OE-actionable upper bounds) holds.
    target = _builtin("hLF")
    plan = build_secretion_plan(target)
    reaction_ids = (
        "Mach_Ribosome_complex_formation",
        "sec_Pdi1p_complex_formation",
        "sec_OSTC_complex_formation",
        "sec_Bst1p_complex_formation",
    )
    simulation = SecretionSimulationResult(
        success=True,
        target_id=target.target_id,
        objective_value=0.01,
        growth_rate=0.10,
        secretion_flux=0.01,
        status="0",
        message="ok",
        constraint_counts={"stoichiometric": 7, "eq_total": 7, "ub_total": 1},
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
        exchange_reaction_id="r_hLF_exchange",
        build_status="supported",
        lp_sensitivity={
            "eq_marginals": (0.0,) * 7,
            "ub_marginals": (0.0,),
            # index 0 (ribosome) and index 1 (PDI) are binding lower-bound floors; index 2 (OSTC) and
            # index 3 (Bst1) are binding upper-bound ceilings (OE-actionable). PDI carries the ~5074
            # shadow price that motivated this fix.
            "lower_marginals": (180.8, 5073.9, 0.0, 0.0),
            "upper_marginals": (0.0, 0.0, 12.0, 3.0),
        },
        key_fluxes={"BIOMASS": 0.10},
    )

    result = analyze_target_protein_lp_attribution(
        target, plan, simulation.constraint_counts, simulation, reaction_ids=reaction_ids, top_n=5
    )

    floors = {row["reaction_id"]: row["secretory_process"] for row in result.floor_constraints_not_oe_addressable}
    oe = {row["reaction_id"]: row["secretory_process"] for row in result.oe_actionable_bottlenecks}
    assert floors["sec_Pdi1p_complex_formation"] == "disulfide_folding"
    assert floors["Mach_Ribosome_complex_formation"] == "ribosome"
    assert oe["sec_OSTC_complex_formation"] == "n_glycan_processing"
    assert oe["sec_Bst1p_complex_formation"] == "secretory_capacity"
    # relative-signal invariants: floors are lower bounds only, OE-actionable are upper bounds only
    assert all(row["bound_type"] == "lower" for row in result.floor_constraints_not_oe_addressable)
    assert all(row["bound_type"] == "upper" for row in result.oe_actionable_bottlenecks)
    assert "sec_Pdi1p_complex_formation" not in oe
