from __future__ import annotations

from dataclasses import replace

import pytest
from types import SimpleNamespace

from pcsec_pichia.oe_capacity import (
    ConfidenceLevel,
    EvidenceSourceType,
    GeneCapacityCatalog,
    GeneCapacityParameterSet,
    GeneEnzymeReactionMapping,
    GPRRole,
    OEDoseMode,
    OECapacityParameterConflictError,
    OECapacityValidationError,
    OEExecutionMode,
    OEExecutionStatus,
    ParameterEstimate,
    ParameterPolicy,
    ParameterScenario,
    build_gene_capacity_specs,
    build_current_model_parameter_policy,
    build_oe_dose_spec,
    fingerprint_oe_capacity_model,
    plan_gene_level_overexpression,
)
from pcsec_pichia.oe_capacity.parameters import load_capacity_anchor_catalog
from pcsec_pichia.oe_capacity.schema import CapacityAnchor, CapacityAnchorCatalog


def test_explicit_and_categorical_dose_inputs_never_infer_unreviewed_multiplier() -> None:
    explicit = build_oe_dose_spec(
        {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        }
    )
    assert explicit.dose_mode is OEDoseMode.EXPLICIT_MULTIPLIER
    assert explicit.expression_multiplier == 2.0

    categorical = build_oe_dose_spec(
        {
            "dose_id": "pGAP-unknown-copy",
            "promoter": "pGAP",
        }
    )
    assert categorical.dose_mode is OEDoseMode.CATEGORICAL_ONLY
    assert categorical.expression_multiplier is None

    with pytest.raises(OECapacityValidationError, match="positive"):
        build_oe_dose_spec(
            {
                "dose_id": "invalid",
                "dose_mode": "explicit_multiplier",
                "expression_multiplier": -1.0,
            }
        )

    with pytest.raises(OECapacityValidationError, match="categorical_only"):
        build_oe_dose_spec(
            {
                "dose_id": "contradictory",
                "dose_mode": "categorical_only",
                "promoter": "pGAP",
                "expression_multiplier": 2.0,
            }
        )


def test_reviewed_promoter_copy_mapping_is_traceable() -> None:
    mapped = build_oe_dose_spec(
        {
            "dose_id": "pGAP-copy2",
            "promoter": "pGAP",
            "copy_number": 2,
        },
        dose_mapping={
            "pGAP|2": {
                "expression_multiplier": 1.8,
                "mapping_source": "reviewed-dose-table-v1",
            }
        },
    )

    assert mapped.dose_mode is OEDoseMode.PROMOTER_COPY_MAPPING
    assert mapped.expression_multiplier == 1.8
    assert mapped.mapping_source == "reviewed-dose-table-v1"


def test_parameter_policy_prefers_current_sources_and_builds_three_scenarios() -> None:
    catalog = GeneCapacityCatalog(
        model_fingerprint="model-v1",
        mappings=(_mapping(),),
    )
    dose = build_oe_dose_spec(
        {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        }
    )
    smoke = _parameter_set(
        "smoke",
        EvidenceSourceType.SMOKE_FIXTURE,
        kcat=90.0,
    )
    local = _parameter_set(
        "local",
        EvidenceSourceType.LOCAL_ENZYME_DATA,
        kcat=120.0,
    )
    policy = ParameterPolicy(
        parameter_sets=(smoke, local),
        test_only_allow_smoke_fixture=True,
    )

    specs = build_gene_capacity_specs("G1", catalog, dose, policy)

    assert tuple(spec.parameter_scenario for spec in specs) == (
        ParameterScenario.LOW,
        ParameterScenario.NOMINAL,
        ParameterScenario.HIGH,
    )
    assert all(spec.kcat is not None and spec.kcat.nominal_value == 120.0 for spec in specs)


def test_current_model_policy_requires_reviewed_absolute_capacity_anchor() -> None:
    combined = SimpleNamespace(
        exact_enzyme_kcat=lambda enzyme_id: 7200.0,
        exact_enzyme_mw=lambda enzyme_id: 60000.0,
    )
    catalog = GeneCapacityCatalog("model-v1", (_mapping(),))
    assert not build_current_model_parameter_policy(catalog, combined).parameter_sets

    policy = build_current_model_parameter_policy(
        catalog,
        combined,
        capacity_anchors=CapacityAnchorCatalog(
            model_fingerprint="model-v1",
            anchors=(_anchor(),),
        ),
        target_id="hLF",
        context_id="ctx-hlf",
        relative_uncertainty=0.1,
    )

    parameter_set = policy.parameter_sets[0]
    assert parameter_set.kcat is not None
    assert parameter_set.kcat.unit == "1/h"
    assert parameter_set.kcat.lower_bound == pytest.approx(6480.0)
    assert parameter_set.molecular_weight is not None
    assert parameter_set.molecular_weight.unit == "g/mol"
    assert parameter_set.baseline_enzyme_amount is not None
    assert parameter_set.baseline_enzyme_amount.nominal_value == 2.5
    assert parameter_set.baseline_enzyme_amount.unit == "model_flux"
    assert parameter_set.baseline_enzyme_amount.source_ref == "reviewed-capacity/v1"


def test_capacity_anchor_asset_loader_validates_review_and_model_identity() -> None:
    loaded = load_capacity_anchor_catalog(
        {
            "schema_version": 1,
            "model_fingerprint": "model-v1",
            "anchors": [
                {
                    "anchor_id": "anchor-G1",
                    "target_id": "hLF",
                    "context_id": "ctx-hlf",
                    "gene_id": "G1",
                    "enzyme_id": "R1_complex",
                    "formation_or_dilution_reaction_id": "R1_complex_formation",
                    "model_fingerprint": "model-v1",
                    "baseline_capacity": 2.5,
                    "unit": "model_flux",
                    "source_ref": "reviewed-capacity/v1",
                    "source_version": "v1",
                    "reviewed_by": "capacity-review-board",
                    "reviewed_at": "2026-07-14",
                }
            ],
        }
    )
    assert loaded.anchors == (_anchor(),)

    second_model = CapacityAnchor(
        **{
            **_anchor().__dict__,
            "anchor_id": "anchor-G1-model-v2",
            "model_fingerprint": "model-v2",
        }
    )
    CapacityAnchorCatalog(
        model_fingerprint="registry-v1",
        anchors=(_anchor(), second_model),
    ).validate()


def test_smoke_fixture_parameters_require_explicit_test_mode() -> None:
    smoke = _parameter_set("smoke", EvidenceSourceType.SMOKE_FIXTURE, kcat=90.0)
    with pytest.raises(OECapacityValidationError, match="test-only"):
        ParameterPolicy(parameter_sets=(smoke,)).validate()

    ParameterPolicy(
        parameter_sets=(smoke,),
        test_only_allow_smoke_fixture=True,
    ).validate()


def test_same_priority_parameter_conflict_is_not_silently_resolved() -> None:
    catalog = GeneCapacityCatalog("model-v1", (_mapping(),))
    dose = build_oe_dose_spec(
        {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        }
    )
    policy = ParameterPolicy(
        parameter_sets=(
            _parameter_set("local-a", EvidenceSourceType.LOCAL_ENZYME_DATA, kcat=100.0),
            _parameter_set("local-b", EvidenceSourceType.LOCAL_ENZYME_DATA, kcat=140.0),
        )
    )

    with pytest.raises(OECapacityParameterConflictError, match="conflicting"):
        build_gene_capacity_specs("G1", catalog, dose, policy)


def test_parameter_set_cannot_cross_fill_a_different_mapping_handle() -> None:
    first = _mapping()
    second = replace(
        first,
        mapping_id="map-G1-R2",
        reaction_id="R2",
        enzyme_variable_id="R2_complex_formation",
        formation_or_dilution_reaction_id="R2_complex_formation",
    )
    catalog = GeneCapacityCatalog("model-v1", (first, second))
    dose = build_oe_dose_spec(
        {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        }
    )

    specs = build_gene_capacity_specs(
        "G1",
        catalog,
        dose,
        ParameterPolicy(
            parameter_sets=(
                _parameter_set(
                    "only-first",
                    EvidenceSourceType.LOCAL_ENZYME_DATA,
                    kcat=100.0,
                ),
            )
        ),
    )

    assert {spec.mapping.mapping_id for spec in specs} == {"map-G1-R1"}


def test_complete_lower_priority_parameters_can_fill_missing_current_values() -> None:
    catalog = GeneCapacityCatalog("model-v1", (_mapping(),))
    dose = build_oe_dose_spec(
        {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        }
    )
    incomplete_current = GeneCapacityParameterSet(
        parameter_set_id="current-incomplete",
        mapping_id="map-G1-R1",
        gene_id="G1",
        enzyme_id="R1_complex",
        kcat=None,
        molecular_weight=None,
        baseline_enzyme_amount=ParameterEstimate(
            "baseline_enzyme_amount",
            1.0,
            1.0,
            1.0,
            "relative_capacity",
            EvidenceSourceType.CURRENT_MODEL,
            "R1_complex_formation",
            "model-v1",
            ConfidenceLevel.HIGH,
        ),
    )
    smoke = _parameter_set("smoke", EvidenceSourceType.SMOKE_FIXTURE, kcat=90.0)

    specs = build_gene_capacity_specs(
        "G1",
        catalog,
        dose,
        ParameterPolicy(
            parameter_sets=(incomplete_current, smoke),
            test_only_allow_smoke_fixture=True,
        ),
    )

    assert specs
    assert all(spec.kcat is not None and spec.kcat.nominal_value == 90.0 for spec in specs)


def test_plan_exposes_comparison_proxy_fallback_and_categorical_boundary() -> None:
    model = SimpleNamespace(
        rxns=["R1", "R1_complex_formation"],
        rules=["x(1)", ""],
        gr_rules=["G1", ""],
        genes=["G1"],
        gene_index={"G1": 0},
    )
    model_fingerprint = fingerprint_oe_capacity_model(model)
    catalog = GeneCapacityCatalog(
        model_fingerprint,
        (_mapping(model_fingerprint=model_fingerprint),),
    )
    numeric_dose = build_oe_dose_spec(
        {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        }
    )
    complete = ParameterPolicy(
        parameter_sets=(
            _parameter_set("local", EvidenceSourceType.LOCAL_ENZYME_DATA, kcat=120.0),
        )
    )

    comparison = plan_gene_level_overexpression(
        model,
        "G1",
        "hLF",
        "ctx-hlf",
        numeric_dose,
        catalog,
        complete,
    )
    assert comparison.execution_mode is OEExecutionMode.COMPARISON
    assert comparison.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE
    assert comparison.proxy_reaction_ids == ("R1",)
    assert len(comparison.executable_capacity_specs) == 3

    external_trace = GeneEnzymeReactionMapping(
        mapping_id="external-trace",
        model_fingerprint="",
        gene_id="G1",
        enzyme_id="external-enzyme",
        reaction_id="external-reaction",
        gpr_rule="G1",
        gpr_role=GPRRole.SINGLE_GENE,
        mapping_source=EvidenceSourceType.EXTERNAL_PICHIA_MODEL,
        mapping_confidence=ConfidenceLevel.LOW,
        execution_status=OEExecutionStatus.EXTERNAL_EVIDENCE_ONLY,
        missing_information=("current_model_gene_enzyme_reaction_mapping",),
    )
    partial_catalog = GeneCapacityCatalog(
        model_fingerprint,
        (*catalog.mappings, external_trace),
    )
    partial = plan_gene_level_overexpression(
        model,
        "G1",
        "hLF",
        "ctx-hlf",
        numeric_dose,
        partial_catalog,
        complete,
    )
    assert partial.execution_status is OEExecutionStatus.PARTIAL_MAPPING

    proxy_only = plan_gene_level_overexpression(
        model,
        "G1",
        "hLF",
        "ctx-hlf",
        numeric_dose,
        catalog,
        ParameterPolicy(parameter_sets=()),
    )
    assert proxy_only.execution_mode is OEExecutionMode.REACTION_PROXY
    assert proxy_only.execution_status is OEExecutionStatus.PROXY_ONLY
    assert "capacity_parameters" in proxy_only.missing_information

    categorical = plan_gene_level_overexpression(
        model,
        "G1",
        "hLF",
        "ctx-hlf",
        build_oe_dose_spec({"dose_id": "pGAP", "promoter": "pGAP"}),
        catalog,
        complete,
    )
    assert categorical.execution_mode is OEExecutionMode.NOT_EXECUTABLE
    assert categorical.execution_status is OEExecutionStatus.CATEGORICAL_DOSE_ONLY


def test_any_incomplete_expected_mapping_downgrades_entire_plan_to_partial() -> None:
    model = SimpleNamespace(
        rxns=["R1", "R2", "R1_complex_formation", "R2_complex_formation"],
        rules=["x(1)", "x(1)", "", ""],
        gr_rules=["G1", "G1", "", ""],
        genes=["G1"],
        gene_index={"G1": 0},
    )
    fingerprint = fingerprint_oe_capacity_model(model)
    first = _mapping(model_fingerprint=fingerprint)
    second = GeneEnzymeReactionMapping(
        **{
            **first.__dict__,
            "mapping_id": "map-G1-R2",
            "enzyme_id": "R2_complex",
            "reaction_id": "R2",
            "enzyme_variable_id": "R2_complex_formation",
            "formation_or_dilution_reaction_id": "R2_complex_formation",
        }
    )
    plan = plan_gene_level_overexpression(
        model,
        "G1",
        "hLF",
        "ctx-hlf",
        build_oe_dose_spec(
            {
                "dose_id": "2x",
                "dose_mode": "explicit_multiplier",
                "expression_multiplier": 2.0,
            }
        ),
        GeneCapacityCatalog(fingerprint, (first, second)),
        ParameterPolicy(
            parameter_sets=(
                _parameter_set(
                    "local",
                    EvidenceSourceType.LOCAL_ENZYME_DATA,
                    kcat=120.0,
                ),
            )
        ),
    )
    assert plan.execution_status is OEExecutionStatus.PARTIAL_MAPPING
    assert "capacity_parameters" in plan.missing_information


def _mapping(*, model_fingerprint: str = "model-v1") -> GeneEnzymeReactionMapping:
    return GeneEnzymeReactionMapping(
        mapping_id="map-G1-R1",
        model_fingerprint=model_fingerprint,
        gene_id="G1",
        enzyme_id="R1_complex",
        reaction_id="R1",
        gpr_rule="G1",
        gpr_role=GPRRole.SINGLE_GENE,
        enzyme_variable_id="R1_complex_formation",
        formation_or_dilution_reaction_id="R1_complex_formation",
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=ConfidenceLevel.HIGH,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
    )


def _anchor() -> CapacityAnchor:
    return CapacityAnchor(
        anchor_id="anchor-G1",
        target_id="hLF",
        context_id="ctx-hlf",
        gene_id="G1",
        enzyme_id="R1_complex",
        formation_or_dilution_reaction_id="R1_complex_formation",
        model_fingerprint="model-v1",
        baseline_capacity=2.5,
        unit="model_flux",
        source_ref="reviewed-capacity/v1",
        source_version="v1",
        reviewed_by="capacity-review-board",
        reviewed_at="2026-07-14",
    )


def _parameter_set(
    parameter_set_id: str,
    source_type: EvidenceSourceType,
    *,
    kcat: float,
) -> GeneCapacityParameterSet:
    def estimate(name: str, value: float, unit: str) -> ParameterEstimate:
        return ParameterEstimate(
            name,
            value,
            value * 0.8,
            value * 1.2,
            unit,
            source_type,
            parameter_set_id,
            "v1",
            ConfidenceLevel.HIGH,
        )

    return GeneCapacityParameterSet(
        parameter_set_id=parameter_set_id,
        mapping_id="map-G1-R1",
        gene_id="G1",
        enzyme_id="R1_complex",
        kcat=estimate("kcat", kcat, "1/s"),
        molecular_weight=estimate("molecular_weight", 60.0, "kDa"),
        baseline_enzyme_amount=estimate(
            "baseline_enzyme_amount",
            1.0,
            "model_flux",
        ),
    )
