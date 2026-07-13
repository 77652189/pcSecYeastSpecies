from __future__ import annotations

import json
from dataclasses import replace

import pytest

from pcsec_pichia.experimental_feedback import (
    CalibrationConfig,
    ConditionContext,
    ExperimentBundle,
    ExperimentRecord,
    HostContext,
    InterventionRecord,
    InterventionType,
    MeasurementRecord,
    MeasurementStatus,
    QualityStatus,
    SchemaValidationError,
    build_calibration_summary,
    build_prediction_index,
    link_experiments_to_predictions,
    validate_experiment_bundle,
    write_calibration_outputs,
)


def test_control_matching_builds_auditable_eligible_effect_record() -> None:
    bundle = _single_candidate_bundle(candidate_value=15.0, control_value=10.0)
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "calibration-run",
                "evidence_items": [
                    {
                        "evidence_id": "hLF-KO-1",
                        "target_id": "hLF",
                        "gene_id": "G1",
                        "intervention_type": "KO",
                        "context_id": "ctx-hlf",
                        "rank": 1,
                        "predicted_direction": "increase",
                        "evidence_tier": "high",
                    }
                ],
            },
        )
    )
    linkage = link_experiments_to_predictions(bundle, index)
    config = CalibrationConfig(
        primary_assay_type="titer",
        increase_threshold_ratio=1.10,
        decrease_threshold_ratio=0.90,
        top_k=(1, 2),
        baseline_hit_rate=0.25,
        minimum_rank_pairs=2,
    )

    summary = build_calibration_summary(validate_experiment_bundle(bundle), linkage, config)

    record = summary.records[0]
    assert record.eligibility_status == "eligible"
    assert record.control_experiment_ids == ("HLF-CONTROL",)
    assert record.observed_ratio == pytest.approx(1.5)
    assert record.observed_direction == "increase"
    assert record.direction_consistent is True
    assert record.hit is True


def test_calibration_metrics_preserve_negative_failed_and_missing_control_records(tmp_path) -> None:
    bundle, index = _metric_fixture()
    linkage = link_experiments_to_predictions(bundle, index)
    config = CalibrationConfig(
        primary_assay_type="titer",
        increase_threshold_ratio=1.10,
        decrease_threshold_ratio=0.90,
        top_k=(1, 2),
        baseline_hit_rate=0.25,
        minimum_rank_pairs=2,
    )

    summary = build_calibration_summary(validate_experiment_bundle(bundle), linkage, config)

    by_target = {item.target_id: item for item in summary.targets}
    hlf = by_target["hLF"]
    assert (hlf.eligible_count, hlf.ineligible_count) == (2, 2)
    assert hlf.direction_consistency_rate == pytest.approx(0.5)
    assert hlf.rank_correlation == pytest.approx(1.0)
    assert [(item.k, item.tested_count, item.hit_rate) for item in hlf.top_k_metrics] == [
        (1, 1, 1.0),
        (2, 2, 0.5),
    ]
    assert [item.relative_baseline_enrichment for item in hlf.top_k_metrics] == [4.0, 2.0]
    tier_rates = {item.evidence_tier: item.hit_rate for item in hlf.evidence_tier_metrics}
    assert tier_rates == {"high": 1.0, "medium": 0.0}
    assert by_target["OPN"].eligible_count == 1
    assert by_target["OPN"].rank_correlation is None
    by_experiment = {record.experiment_id: record for record in summary.records}
    assert by_experiment["HLF-NEGATIVE"].hit is False
    assert by_experiment["HLF-FAILED"].measurement_statuses == ("assay_failed",)
    assert by_experiment["HLF-FAILED"].eligibility_status == "ineligible"
    assert by_experiment["HLF-NO-CONTROL"].ineligibility_reasons == ("control_match_missing",)
    outputs = write_calibration_outputs(summary, tmp_path)
    manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    assert manifest["config"]["top_k"] == [1, 2]
    assert manifest["config"]["increase_threshold_ratio"] == 1.10
    assert manifest["config"]["decrease_threshold_ratio"] == 0.90
    assert manifest["config"]["baseline_hit_rate"] == 0.25
    records_text = outputs.records_path.read_text(encoding="utf-8")
    assert "HLF-FAILED" in records_text
    assert "HLF-NO-CONTROL" in records_text


def test_control_biological_replicates_are_not_weighted_by_technical_replicate_count() -> None:
    base = _single_candidate_bundle(candidate_value=15.0, control_value=10.0)
    host = base.experiments[0].host
    condition = base.experiments[0].condition
    bundle = ExperimentBundle(
        experiments=(
            *base.experiments,
            ExperimentRecord("HLF-CONTROL-2", "hLF", host, "B01", condition, context_id="ctx-hlf"),
        ),
        interventions=(
            *base.interventions,
            InterventionRecord("HLF-CONTROL-2", "CONTROL-1", 1, InterventionType.CONTROL),
        ),
        measurements=(
            *base.measurements,
            _measurement("HLF-CONTROL", "CONTROL-T2", 10.0),
            _measurement("HLF-CONTROL-2", "CONTROL-T1", 20.0),
        ),
    )
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "calibration-run",
                "evidence_items": [
                    _prediction_row("hLF-KO-1", "hLF", "G1", 1, "high", "ctx-hlf")
                ],
            },
        )
    )
    linkage = link_experiments_to_predictions(bundle, index)

    summary = build_calibration_summary(
        validate_experiment_bundle(bundle),
        linkage,
        CalibrationConfig(increase_threshold_ratio=1.10, decrease_threshold_ratio=0.90),
    )

    record = summary.records[0]
    assert record.control_experiment_ids == ("HLF-CONTROL", "HLF-CONTROL-2")
    assert record.control_value == pytest.approx(15.0)
    assert record.observed_ratio == pytest.approx(1.0)
    assert record.hit is False


def test_calibration_config_rejects_non_finite_thresholds_and_bool_top_k() -> None:
    with pytest.raises(SchemaValidationError, match="top_k"):
        CalibrationConfig(top_k=(True,)).validate()
    with pytest.raises(SchemaValidationError, match="finite"):
        CalibrationConfig(baseline_hit_rate=float("nan")).validate()


def test_incomplete_condition_and_invalid_quality_cannot_enter_calibration() -> None:
    base = _single_candidate_bundle(candidate_value=15.0, control_value=10.0)
    incomplete = replace(base.experiments[0].condition, sampling_time_h=None)
    bundle = replace(
        base,
        experiments=(
            replace(base.experiments[0], condition=incomplete),
            replace(base.experiments[1], condition=incomplete, quality_status=QualityStatus.INVALID),
        ),
    )
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "calibration-run",
                "evidence_items": [
                    _prediction_row("hLF-KO-1", "hLF", "G1", 1, "high", "ctx-hlf")
                ],
            },
        )
    )
    linkage = link_experiments_to_predictions(bundle, index)

    summary = build_calibration_summary(
        validate_experiment_bundle(bundle),
        linkage,
        CalibrationConfig(),
    )

    assert summary.records[0].eligibility_status == "ineligible"
    assert "experiment_quality_status:invalid" in summary.records[0].ineligibility_reasons


def test_ranking_assessment_counts_only_records_with_rank_and_observed_ratio() -> None:
    base = _single_candidate_bundle(candidate_value=11.0, control_value=10.0)
    host = base.experiments[0].host
    condition = base.experiments[0].condition
    bundle = replace(
        base,
        experiments=(
            *base.experiments,
            ExperimentRecord("HLF-CANDIDATE-2", "hLF", host, "B01", condition, context_id="ctx-hlf"),
            ExperimentRecord("HLF-CANDIDATE-3", "hLF", host, "B01", condition, context_id="ctx-hlf"),
        ),
        interventions=(
            *base.interventions,
            InterventionRecord(
                "HLF-CANDIDATE-2",
                "KO-1",
                1,
                InterventionType.KO,
                gene_id="G2",
                construction_method="CRISPR-Cas9",
                prediction_run_id="calibration-run",
                evidence_id="hLF-KO-2",
            ),
            InterventionRecord(
                "HLF-CANDIDATE-3",
                "KO-1",
                1,
                InterventionType.KO,
                gene_id="G3",
                construction_method="CRISPR-Cas9",
                prediction_run_id="calibration-run",
                evidence_id="hLF-KO-3",
            ),
        ),
        measurements=(
            *base.measurements,
            _measurement("HLF-CANDIDATE-2", "CANDIDATE-2-T1", 12.0),
            _measurement("HLF-CANDIDATE-3", "CANDIDATE-3-T1", 13.0),
        ),
    )
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "calibration-run",
                "evidence_items": [
                    _prediction_row("hLF-KO-1", "hLF", "G1", 1, "high", "ctx-hlf"),
                    _prediction_row("hLF-KO-2", "hLF", "G2", None, "medium", "ctx-hlf"),
                    _prediction_row("hLF-KO-3", "hLF", "G3", None, "low", "ctx-hlf"),
                ],
            },
        )
    )
    linkage = link_experiments_to_predictions(bundle, index)

    summary = build_calibration_summary(
        validate_experiment_bundle(bundle),
        linkage,
        CalibrationConfig(minimum_rank_pairs=2),
    )

    target = summary.targets[0]
    assert target.eligible_count == 3
    assert target.comparable_rank_pair_count == 1
    assert target.ranking_assessment == "insufficient_evidence"
    assert target.rank_correlation is None


def _single_candidate_bundle(*, candidate_value: float, control_value: float) -> ExperimentBundle:
    host = HostContext("Komagataella phaffii", "X33", "X33")
    condition = ConditionContext("BMMY", "methanol", "shake_flask", 30.0, 6.0, "250 rpm", 72.0)
    experiments = (
        ExperimentRecord("HLF-CONTROL", "hLF", host, "B01", condition, context_id="ctx-hlf"),
        ExperimentRecord("HLF-CANDIDATE", "hLF", host, "B01", condition, context_id="ctx-hlf"),
    )
    interventions = (
        InterventionRecord("HLF-CONTROL", "CONTROL-1", 1, InterventionType.CONTROL),
        InterventionRecord(
            "HLF-CANDIDATE",
            "KO-1",
            1,
            InterventionType.KO,
            gene_id="G1",
            construction_method="CRISPR-Cas9",
            prediction_run_id="calibration-run",
            evidence_id="hLF-KO-1",
        ),
    )
    measurements = (
        _measurement("HLF-CONTROL", "CONTROL-T1", control_value),
        _measurement("HLF-CANDIDATE", "CANDIDATE-T1", candidate_value),
    )
    return ExperimentBundle(experiments=experiments, interventions=interventions, measurements=measurements)


def _measurement(experiment_id: str, measurement_id: str, value: float) -> MeasurementRecord:
    return MeasurementRecord(
        experiment_id=experiment_id,
        measurement_id=measurement_id,
        assay_type="titer",
        assay_method="ELISA",
        compartment="extracellular",
        raw_value=value,
        raw_unit="mg/L",
        canonical_value=value,
        canonical_unit="mg/L",
        status=MeasurementStatus.VALID,
    )


def _metric_fixture() -> tuple[ExperimentBundle, object]:
    host = HostContext("Komagataella phaffii", "X33", "X33")
    condition = ConditionContext("BMMY", "methanol", "shake_flask", 30.0, 6.0, "250 rpm", 72.0)
    specs = (
        ("HLF-CONTROL-M", "hLF", "B01", "ctx-hlf"),
        ("HLF-HIT", "hLF", "B01", "ctx-hlf"),
        ("HLF-NEGATIVE", "hLF", "B01", "ctx-hlf"),
        ("HLF-FAILED", "hLF", "B01", "ctx-hlf"),
        ("HLF-NO-CONTROL", "hLF", "B99", "ctx-hlf"),
        ("OPN-CONTROL-M", "OPN", "B02", "ctx-opn"),
        ("OPN-HIT", "OPN", "B02", "ctx-opn"),
    )
    experiments = tuple(
        ExperimentRecord(experiment_id, target_id, host, batch, condition, context_id=context_id)
        for experiment_id, target_id, batch, context_id in specs
    )
    interventions = (
        InterventionRecord("HLF-CONTROL-M", "CONTROL", 1, InterventionType.CONTROL),
        _linked_ko("HLF-HIT", "G1", "HLF-E1"),
        _linked_ko("HLF-NEGATIVE", "G2", "HLF-E2"),
        _linked_ko("HLF-FAILED", "G3", "HLF-E3"),
        _linked_ko("HLF-NO-CONTROL", "G4", "HLF-E4"),
        InterventionRecord("OPN-CONTROL-M", "CONTROL", 1, InterventionType.CONTROL),
        _linked_ko("OPN-HIT", "G5", "OPN-E1"),
    )
    measurements = (
        _measurement("HLF-CONTROL-M", "H-C", 10.0),
        _measurement("HLF-HIT", "H-1", 15.0),
        _measurement("HLF-NEGATIVE", "H-2", 9.0),
        MeasurementRecord(
            "HLF-FAILED",
            "H-3",
            "titer",
            "ELISA",
            "extracellular",
            None,
            "mg/L",
            None,
            "mg/L",
            MeasurementStatus.ASSAY_FAILED,
            status_reason="sanitized assay failure",
        ),
        _measurement("HLF-NO-CONTROL", "H-4", 12.0),
        _measurement("OPN-CONTROL-M", "O-C", 20.0),
        _measurement("OPN-HIT", "O-1", 22.0),
    )
    prediction_rows = [
        _prediction_row("HLF-E1", "hLF", "G1", 1, "high", "ctx-hlf"),
        _prediction_row("HLF-E2", "hLF", "G2", 2, "medium", "ctx-hlf"),
        _prediction_row("HLF-E3", "hLF", "G3", 3, "high", "ctx-hlf"),
        _prediction_row("HLF-E4", "hLF", "G4", 4, "low", "ctx-hlf"),
        _prediction_row("OPN-E1", "OPN", "G5", 1, "high", "ctx-opn"),
    ]
    bundle = ExperimentBundle(experiments=experiments, interventions=interventions, measurements=measurements)
    index = build_prediction_index(
        ({"prediction_run_id": "metric-run", "evidence_items": prediction_rows},)
    )
    return bundle, index


def _linked_ko(experiment_id: str, gene_id: str, evidence_id: str) -> InterventionRecord:
    return InterventionRecord(
        experiment_id,
        "KO-1",
        1,
        InterventionType.KO,
        gene_id=gene_id,
        construction_method="CRISPR-Cas9",
        prediction_run_id="metric-run",
        evidence_id=evidence_id,
    )


def _prediction_row(
    evidence_id: str,
    target_id: str,
    gene_id: str,
    rank: int | None,
    evidence_tier: str,
    context_id: str,
) -> dict[str, object]:
    return {
        "evidence_id": evidence_id,
        "target_id": target_id,
        "gene_id": gene_id,
        "intervention_type": "KO",
        "context_id": context_id,
        "rank": rank,
        "predicted_direction": "increase",
        "evidence_tier": evidence_tier,
    }
