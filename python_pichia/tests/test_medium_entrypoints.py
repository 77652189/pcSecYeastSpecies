from __future__ import annotations

import pytest

from pcsec_pichia.loading import load_pcsec_pichia_inputs, load_pcsec_pichia_model, prepare_carbon_source_model, repo_root
from pcsec_pichia.media import (
    MATLAB_LEGACY_COST_MEDIUM_EXCHANGES,
    compose_medium_condition_bounds,
    diff_medium_conditions,
    list_carbon_source_specs,
    list_carbon_source_formulations,
    list_medium_condition_specs,
    load_base_medium_spec,
    load_carbon_source_formulation,
    load_compatibility_overlay_spec,
    load_medium_condition_spec,
    medium_condition_warnings,
    summarize_medium_condition,
)


def _actual_bound(model, reaction_id: str):
    index = model.reaction_index[reaction_id]
    return float(model.lb[index]), float(model.ub[index])


def _bound_map(condition_id: str):
    return {bound.reaction_id: bound for bound in compose_medium_condition_bounds(condition_id)}


def test_current_default_medium_condition_is_glucose_core_aa_corrected() -> None:
    condition = load_medium_condition_spec("glucose_ynb_core_aa_corrected")
    bounds = _bound_map(condition.condition_id)

    assert condition.status == "active"
    assert condition.carbon_source_id == "glucose"
    assert condition.base_medium_id == "ynb_core_aa"
    assert condition.compatibility_overlay_id == "none"
    assert condition.legacy_media_type == 4

    assert bounds["Ex_glc_D"].lower_bound == -1000.0
    assert bounds["Ex_o2"].lower_bound == -1000.0
    assert bounds["Ex_nh4"].lower_bound == -1000.0
    assert bounds["Ex_pi"].lower_bound == -1000.0
    assert bounds["Ex_glyc"].lower_bound == 0.0
    assert bounds["Ex_meoh"].lower_bound == 0.0
    assert bounds["BIOMASS_glyc"].upper_bound == 0.0
    assert bounds["BIOMASS_meoh"].upper_bound == 0.0
    assert bounds["Ex_btn"].lower_bound == -2.0
    assert bounds["Ex_arg_L"].lower_bound == -0.08
    assert bounds["Ex_ala_L"].lower_bound == 0.0


def test_matlab_legacy_cost_overlay_is_explicit_and_limited_to_cost_scope() -> None:
    overlay = load_compatibility_overlay_spec("matlab_legacy_cost")
    assert overlay.scope == "protein_cost_slope_only"
    assert tuple(bound.reaction_id for bound in overlay.bounds) == MATLAB_LEGACY_COST_MEDIUM_EXCHANGES
    assert all(bound.lower_bound == 0.0 for bound in overlay.bounds)
    assert all(bound.confidence == "legacy" for bound in overlay.bounds)

    diff = diff_medium_conditions(
        "glucose_ynb_core_aa_corrected",
        "glucose_ynb_core_aa_matlab_legacy_cost",
    )
    assert diff["difference_count"] == len(MATLAB_LEGACY_COST_MEDIUM_EXCHANGES)
    assert tuple(row["reaction_id"] for row in diff["differences"]) == tuple(sorted(MATLAB_LEGACY_COST_MEDIUM_EXCHANGES))
    for row in diff["differences"]:
        assert row["left_lower_bound"] == -1000.0
        assert row["right_lower_bound"] == 0.0


def test_base_medium_and_carbon_source_are_separate_axes() -> None:
    base = load_base_medium_spec("ynb_core_aa")
    carbon_sources = {item.carbon_source_id: item for item in list_carbon_source_specs()}
    formulations = {item.carbon_source_id: item for item in list_carbon_source_formulations()}

    assert base.legacy_media_type == 4
    assert "glucose" in carbon_sources
    assert carbon_sources["glucose"].status == "active"
    assert carbon_sources["glycerol"].status == "active"
    assert carbon_sources["methanol"].status == "active"
    assert carbon_sources["glucose_glycerol"].status == "active"
    assert carbon_sources["glycerol_methanol"].status == "active"
    assert set(formulations) == set(carbon_sources)


@pytest.mark.parametrize(
    ("carbon_source_id", "selected_growth", "objective_weights"),
    (
        ("glucose", "BIOMASS", {"Ex_glc_D": 1.0}),
        ("glycerol", "BIOMASS_glyc", {"Ex_glyc": 1.0}),
        ("methanol", "BIOMASS_meoh", {"Ex_meoh": 1.0}),
        ("glucose_glycerol", "BIOMASS", {"Ex_glc_D": 1.0, "Ex_glyc": 1.0}),
        ("glycerol_methanol", "BIOMASS_glyc", {"Ex_glyc": 1.0, "Ex_meoh": 1.0}),
    ),
)
def test_carbon_source_formulation_centralizes_growth_and_objective_selection(
    carbon_source_id: str,
    selected_growth: str,
    objective_weights: dict[str, float],
) -> None:
    formulation = load_carbon_source_formulation(carbon_source_id)

    assert formulation.selected_growth_reaction_id == selected_growth
    assert formulation.carbon_objective_weights == objective_weights
    assert set(formulation.active_uptake_reaction_ids) == set(objective_weights)
    assert selected_growth in formulation.candidate_growth_reaction_ids


def test_medium_summary_exposes_carbon_source_formulation_for_reports_and_harnesses() -> None:
    summary = summarize_medium_condition("glycerol_methanol_ynb_core_aa_corrected")
    formulation = summary["carbon_source_formulation"]

    assert formulation["carbon_source_id"] == "glycerol_methanol"
    assert formulation["selected_growth_reaction_id"] == "BIOMASS_glyc"
    assert formulation["carbon_objective_weights"] == {"Ex_glyc": 1.0, "Ex_meoh": 1.0}
    assert "Ex_glc_D" in formulation["blocked_uptake_reaction_ids"]


def test_medium_condition_summary_is_auditable_without_running_solver() -> None:
    summary = summarize_medium_condition("glucose_ynb_core_aa_matlab_legacy_cost")

    assert summary["condition_id"] == "glucose_ynb_core_aa_matlab_legacy_cost"
    assert summary["carbon_source_id"] == "glucose"
    assert summary["base_medium_id"] == "ynb_core_aa"
    assert summary["compatibility_overlay_id"] == "matlab_legacy_cost"
    assert "Ex_glc_D" in summary["active_uptake_reactions"]
    assert "Ex_nh4" in summary["closed_uptake_reactions"]
    assert summary["scientific_status"] == "matlab_legacy_artifact_compatibility"
    assert any("历史 artifact 对齐" in item for item in summary["warnings"])


def test_loaded_inputs_expose_runtime_carbon_source_formulation() -> None:
    inputs = load_pcsec_pichia_inputs(repo_root(), carbon_source_id="methanol")
    summary = summarize_medium_condition(inputs.medium_condition_id)

    assert inputs.carbon_source_formulation.carbon_source_id == "methanol"
    assert inputs.carbon_source_formulation.selected_growth_reaction_id == "BIOMASS_meoh"
    assert summary["carbon_source_formulation"]["carbon_objective_weights"] == {"Ex_meoh": 1.0}


def test_non_glucose_medium_conditions_expose_scientific_warnings() -> None:
    glucose = summarize_medium_condition("glucose_ynb_core_aa_corrected")
    mixed = summarize_medium_condition("glucose_glycerol_ynb_core_aa_corrected")

    assert glucose["scientific_status"] == "default_python_corrected_reference"
    assert glucose["warnings"] == []
    assert mixed["scientific_status"] == "draft_co_carbon_boundary_requires_promoter_context"
    assert any("共碳源摄取" in item for item in mixed["warnings"])
    assert any("旧 MATLAB baseline fully aligned" in item for item in medium_condition_warnings(mixed["condition_id"]))


def test_core_medium_conditions_are_active_for_supported_carbon_sources() -> None:
    conditions = {item.condition_id: item for item in list_medium_condition_specs()}

    assert conditions["glucose_ynb_core_aa_corrected"].status == "active"
    assert conditions["glucose_ynb_core_aa_matlab_legacy_cost"].status == "active"
    assert conditions["glycerol_ynb_core_aa_corrected"].status == "active"
    assert conditions["methanol_ynb_core_aa_corrected"].status == "active"
    assert conditions["glucose_glycerol_ynb_core_aa_corrected"].status == "active"
    assert conditions["glycerol_methanol_ynb_core_aa_corrected"].status == "active"


def test_minimal_and_all_aa_conditions_exist_for_supported_carbon_sources() -> None:
    conditions = {item.condition_id: item for item in list_medium_condition_specs()}

    for carbon_source_id in ("glucose", "glycerol", "methanol", "glucose_glycerol", "glycerol_methanol"):
        minimal = conditions[f"{carbon_source_id}_ynb_minimal_corrected"]
        all_aa = conditions[f"{carbon_source_id}_ynb_all_aa_corrected"]
        assert minimal.status == "active"
        assert minimal.legacy_media_type == 2
        assert all_aa.status == "active"
        assert all_aa.legacy_media_type == 5


@pytest.mark.parametrize(
    ("carbon_source_id", "expected"),
    (
        (
            "glucose",
            {
                "BIOMASS": (0.0, 1000.0),
                "BIOMASS_glyc": (0.0, 0.0),
                "BIOMASS_meoh": (0.0, 0.0),
                "Ex_glc_D": (-1000.0, 1000.0),
                "Ex_glyc": (0.0, 1000.0),
                "Ex_meoh": (0.0, 1000.0),
            },
        ),
        (
            "glycerol",
            {
                "BIOMASS": (0.0, 0.0),
                "BIOMASS_glyc": (0.0, 1000.0),
                "BIOMASS_meoh": (0.0, 0.0),
                "Ex_glc_D": (0.0, 1000.0),
                "Ex_glyc": (-1000.0, 1000.0),
                "Ex_meoh": (0.0, 1000.0),
            },
        ),
        (
            "methanol",
            {
                "BIOMASS": (0.0, 0.0),
                "BIOMASS_glyc": (0.0, 0.0),
                "BIOMASS_meoh": (0.0, 1000.0),
                "Ex_glc_D": (0.0, 1000.0),
                "Ex_glyc": (0.0, 1000.0),
                "Ex_meoh": (-1000.0, 1000.0),
            },
        ),
        (
            "glucose_glycerol",
            {
                "BIOMASS": (0.0, 1000.0),
                "BIOMASS_glyc": (0.0, 1000.0),
                "BIOMASS_meoh": (0.0, 0.0),
                "Ex_glc_D": (-1000.0, 1000.0),
                "Ex_glyc": (-1000.0, 1000.0),
                "Ex_meoh": (0.0, 1000.0),
            },
        ),
        (
            "glycerol_methanol",
            {
                "BIOMASS": (0.0, 0.0),
                "BIOMASS_glyc": (0.0, 1000.0),
                "BIOMASS_meoh": (0.0, 1000.0),
                "Ex_glc_D": (0.0, 1000.0),
                "Ex_glyc": (-1000.0, 1000.0),
                "Ex_meoh": (-1000.0, 1000.0),
            },
        ),
    ),
)
def test_prepared_model_uses_medium_condition_bounds_for_carbon_sources(
    carbon_source_id: str,
    expected: dict[str, tuple[float, float]],
) -> None:
    model = load_pcsec_pichia_model(repo_root())
    prepared = prepare_carbon_source_model(model, carbon_source_id=carbon_source_id)

    for reaction_id, bounds in expected.items():
        assert _actual_bound(prepared, reaction_id) == bounds
