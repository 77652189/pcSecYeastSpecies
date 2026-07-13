from __future__ import annotations

import json
from pathlib import Path

import pytest

from pcsec_pichia.experimental_feedback import (
    ConditionContext,
    ExperimentBundle,
    ExperimentImportManifest,
    ExperimentRecord,
    HostContext,
    InterventionRecord,
    InterventionType,
    MeasurementRecord,
    MeasurementStatus,
    PredictionLinkRecord,
    PredictionLinkStatus,
    SchemaValidationError,
    QualityStatus,
    UnitValidationError,
    build_calibration_summary,
    build_prediction_index,
    link_experiments_to_predictions,
    load_experiment_bundle,
    validate_experiment_bundle,
    write_experiment_feedback_cache,
)


FIXTURE_PATH = Path(__file__).parent / "fixtures" / "experimental_feedback" / "round0_cases.json"


def test_bundle_preserves_combination_intervention_components() -> None:
    experiment = ExperimentRecord(
        experiment_id="HLF-B01-R1",
        target_id="hLF",
        host=HostContext(
            species="Komagataella phaffii",
            strain="X33",
            parent_strain="X33",
        ),
        batch_id="B01",
        condition=ConditionContext(
            medium="BMMY",
            carbon_source="methanol",
            culture_mode="shake_flask",
            temperature_c=30.0,
            ph=6.0,
            oxygen_or_agitation="250 rpm",
            sampling_time_h=72.0,
        ),
    )
    interventions = (
        InterventionRecord(
            experiment_id=experiment.experiment_id,
            intervention_id="KO-1",
            component_index=1,
            intervention_type=InterventionType.KO,
            gene_id="PAS_chr1-1_0001",
            construction_method="CRISPR-Cas9",
        ),
        InterventionRecord(
            experiment_id=experiment.experiment_id,
            intervention_id="OE-1",
            component_index=2,
            intervention_type=InterventionType.OE,
            gene_id="PAS_chr2-1_0140",
            construct_id="pGAP-KAR2",
            promoter="pGAP",
            induction_mode="constitutive",
            warnings=("copy_number_unknown",),
        ),
    )

    bundle = ExperimentBundle(experiments=(experiment,), interventions=interventions)

    bundle.validate()
    assert [item.component_index for item in bundle.interventions] == [1, 2]
    assert [item.intervention_type.value for item in bundle.interventions] == ["KO", "OE"]


def test_assay_failure_remains_missing_not_zero() -> None:
    measurement = MeasurementRecord(
        experiment_id="OPN-B02-R1",
        measurement_id="TITER-T2",
        assay_type="titer",
        assay_method="ELISA",
        compartment="extracellular",
        raw_value=None,
        raw_unit="mg/L",
        canonical_value=None,
        canonical_unit="mg/L",
        status=MeasurementStatus.ASSAY_FAILED,
        technical_replicate_id="T2",
        status_reason="plate control failed",
    )

    measurement.validate()
    assert measurement.raw_value is None
    assert measurement.canonical_value is None
    assert measurement.status.value == "assay_failed"


def test_common_name_only_cannot_form_a_matched_prediction_link() -> None:
    link = PredictionLinkRecord(
        experiment_id="HLF-B01-R1",
        intervention_id="OE-1",
        prediction_run_id="screen-hlf-001",
        evidence_id="evidence-17",
        target_id="hLF",
        gene_id="",
        common_name="KAR2",
        intervention_type=InterventionType.OE,
        status=PredictionLinkStatus.MATCHED,
    )

    with pytest.raises(SchemaValidationError, match="gene_id"):
        link.validate()


def test_unregistered_titer_unit_is_rejected() -> None:
    measurement = MeasurementRecord(
        experiment_id="HLF-B01-R1",
        measurement_id="TITER-T1",
        assay_type="titer",
        assay_method="ELISA",
        compartment="extracellular",
        raw_value=0.12,
        raw_unit="g/L",
        canonical_value=0.12,
        canonical_unit="g/L",
        status=MeasurementStatus.VALID,
    )

    with pytest.raises(UnitValidationError, match="mg/L"):
        measurement.validate()


def test_round0_fixture_covers_required_sanitized_cases() -> None:
    payload = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))

    assert set(payload["case_labels"]) == {
        "hLF",
        "OPN",
        "control",
        "KO",
        "OE",
        "combination_components",
        "technical_replicates",
        "bad_unit",
        "missing_control",
        "ambiguous_link",
    }
    assert all(str(row["operator_id"]).startswith("operator-") for row in payload["experiments"])
    assert len({row["measurement_id"] for row in payload["technical_replicates"]}) == 2


def test_public_api_exposes_structured_validation_result() -> None:
    invalid = ExperimentBundle(
        experiments=(),
        interventions=(
            InterventionRecord(
                experiment_id="missing",
                intervention_id="KO-1",
                component_index=1,
                intervention_type=InterventionType.KO,
                gene_id="PAS_chr1-1_0001",
                construction_method="CRISPR-Cas9",
            ),
        ),
    )

    result = validate_experiment_bundle(invalid)

    assert result.is_valid is False
    assert result.errors[0].code == "schema_validation_error"
    assert all(
        callable(api)
        for api in (
            load_experiment_bundle,
            validate_experiment_bundle,
            write_experiment_feedback_cache,
            build_prediction_index,
            link_experiments_to_predictions,
            build_calibration_summary,
        )
    )


def test_import_manifest_and_quality_status_are_explicit() -> None:
    manifest = ExperimentImportManifest(
        source_file="sanitized_round0.csv",
        source_sha256="a" * 64,
        imported_at="2026-07-13T00:00:00+00:00",
        record_count=2,
        warnings=("copy_number_unknown",),
    )

    manifest.validate()
    assert manifest.schema_version == 1
    assert QualityStatus.INVALID.value == "invalid"


def test_custom_target_name_and_unknown_oe_copy_number_are_not_silently_defaulted() -> None:
    custom = ExperimentRecord(
        experiment_id="CUSTOM-B01-R1",
        target_id="custom-target-001",
        target_name="sanitized custom target",
        host=HostContext("Komagataella phaffii", "X33", "X33"),
        batch_id="B01",
        condition=ConditionContext("BMMY", "methanol", "shake_flask", 30.0, 6.0, "250 rpm", 72.0),
    )
    oe_without_warning = InterventionRecord(
        experiment_id=custom.experiment_id,
        intervention_id="OE-1",
        component_index=1,
        intervention_type=InterventionType.OE,
        gene_id="PAS_chr2-1_0140",
        construct_id="pGAP-KAR2",
        promoter="pGAP",
        induction_mode="constitutive",
    )

    custom.validate()
    with pytest.raises(SchemaValidationError, match="copy_number_unknown"):
        oe_without_warning.validate()


def test_canonical_schema_rejects_raw_enum_strings_and_empty_bundles() -> None:
    raw_string_intervention = InterventionRecord(
        experiment_id="HLF-B01-R1",
        intervention_id="OE-RAW",
        component_index=1,
        intervention_type="OE",  # type: ignore[arg-type]
    )

    with pytest.raises(SchemaValidationError, match="InterventionType"):
        raw_string_intervention.validate()
    with pytest.raises(SchemaValidationError, match="at least one experiment"):
        ExperimentBundle(experiments=(), interventions=()).validate()
