from __future__ import annotations

import argparse
from pathlib import Path

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
    build_calibration_summary,
    build_prediction_index,
    link_experiments_to_predictions,
    validate_experiment_bundle,
    write_calibration_outputs,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("local_runs/experiment_feedback/round3_calibration"),
    )
    args = parser.parse_args()
    bundle = _bundle()
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "sanitized-calibration-run",
                "evidence_items": [
                    _prediction("HLF-E1", "hLF", "G1", 1, "high", "ctx-hlf"),
                    _prediction("HLF-E2", "hLF", "G2", 2, "medium", "ctx-hlf"),
                    _prediction("OPN-E1", "OPN", "G3", 1, "high", "ctx-opn"),
                ],
            },
        )
    )
    linkage = link_experiments_to_predictions(bundle, index)
    config = CalibrationConfig(
        increase_threshold_ratio=1.10,
        decrease_threshold_ratio=0.90,
        top_k=(1, 2),
        baseline_hit_rate=0.25,
        minimum_rank_pairs=2,
    )
    summary = build_calibration_summary(validate_experiment_bundle(bundle), linkage, config)
    outputs = write_calibration_outputs(summary, args.output_dir)
    by_target = {target.target_id: target for target in summary.targets}
    if set(by_target) != {"hLF", "OPN"}:
        return 1
    if not any(record.measurement_statuses == ("assay_failed",) for record in summary.records):
        return 1
    print(outputs.manifest_path)
    return 0


def _bundle() -> ExperimentBundle:
    host = HostContext("Komagataella phaffii", "sanitized-strain", "sanitized-parent")
    condition = ConditionContext("sanitized-medium", "methanol", "shake_flask", 30.0, 6.0, "sanitized", 72.0)
    experiments = tuple(
        ExperimentRecord(experiment_id, target_id, host, batch, condition, context_id=context_id)
        for experiment_id, target_id, batch, context_id in (
            ("CAL-H-C", "hLF", "B1", "ctx-hlf"),
            ("CAL-H-HIT", "hLF", "B1", "ctx-hlf"),
            ("CAL-H-FAIL", "hLF", "B1", "ctx-hlf"),
            ("CAL-O-C", "OPN", "B2", "ctx-opn"),
            ("CAL-O-NEG", "OPN", "B2", "ctx-opn"),
        )
    )
    interventions = (
        InterventionRecord("CAL-H-C", "CONTROL", 1, InterventionType.CONTROL),
        _ko("CAL-H-HIT", "G1", "HLF-E1"),
        _ko("CAL-H-FAIL", "G2", "HLF-E2"),
        InterventionRecord("CAL-O-C", "CONTROL", 1, InterventionType.CONTROL),
        _ko("CAL-O-NEG", "G3", "OPN-E1"),
    )
    measurements = (
        _measurement("CAL-H-C", "H-C", 10.0),
        _measurement("CAL-H-HIT", "H-HIT", 15.0),
        MeasurementRecord(
            "CAL-H-FAIL",
            "H-FAIL",
            "titer",
            "sanitized-assay",
            "extracellular",
            None,
            "mg/L",
            None,
            "mg/L",
            MeasurementStatus.ASSAY_FAILED,
            status_reason="sanitized assay failure",
        ),
        _measurement("CAL-O-C", "O-C", 20.0),
        _measurement("CAL-O-NEG", "O-NEG", 18.0),
    )
    return ExperimentBundle(experiments=experiments, interventions=interventions, measurements=measurements)


def _ko(experiment_id: str, gene_id: str, evidence_id: str) -> InterventionRecord:
    return InterventionRecord(
        experiment_id,
        "KO-1",
        1,
        InterventionType.KO,
        gene_id=gene_id,
        construction_method="CRISPR-Cas9",
        prediction_run_id="sanitized-calibration-run",
        evidence_id=evidence_id,
    )


def _measurement(experiment_id: str, measurement_id: str, value: float) -> MeasurementRecord:
    return MeasurementRecord(
        experiment_id,
        measurement_id,
        "titer",
        "sanitized-assay",
        "extracellular",
        value,
        "mg/L",
        value,
        "mg/L",
        MeasurementStatus.VALID,
    )


def _prediction(
    evidence_id: str,
    target_id: str,
    gene_id: str,
    rank: int,
    tier: str,
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
        "evidence_tier": tier,
    }


if __name__ == "__main__":
    raise SystemExit(main())
