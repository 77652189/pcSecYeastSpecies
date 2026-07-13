from __future__ import annotations

import csv
import json
from dataclasses import asdict
from enum import Enum
from pathlib import Path

from app.services.pichia_experiment_feedback_service import (
    export_experiment_feedback_issues,
    export_experiment_feedback_report,
    load_experiment_feedback_run,
    submit_experiment_feedback_import,
)
from pcsec_pichia.experimental_feedback import (
    ConditionContext,
    ExperimentRecord,
    HostContext,
    InterventionRecord,
    InterventionType,
    MeasurementRecord,
    MeasurementStatus,
)


def test_service_imports_links_calibrates_and_exposes_issue_exports(tmp_path) -> None:
    csv_bytes = _experiment_csv_bytes(tmp_path)
    fact_pack = {
        "prediction_run_id": "ui-screen-run",
        "evidence_items": [
            {
                "evidence_id": "hLF-KO-UI-1",
                "target_id": "hLF",
                "gene_id": "G1",
                "intervention_type": "KO",
                "context_id": "ctx-hlf-ui",
                "rank": 1,
                "predicted_direction": "increase",
                "evidence_tier": "high",
            }
        ],
    }

    result = submit_experiment_feedback_import(
        experiment_filename="sanitized_ui.csv",
        experiment_bytes=csv_bytes,
        prediction_filename="fact_pack.json",
        prediction_bytes=json.dumps(fact_pack).encode("utf-8"),
        run_name="ui-contract",
        output_root=tmp_path / "runs",
    )

    assert result["validation"]["is_valid"] is True
    assert result["linkage"]["matched_count"] == 1
    assert result["calibration"]["targets"][0]["target_id"] == "hLF"
    assert result["calibration"]["targets"][0]["eligible_count"] == 1
    run_dir = tmp_path / "runs" / "ui-contract"
    assert (run_dir / "inbox" / "sanitized_ui.csv").exists()
    assert (run_dir / "validated" / "manifest.json").exists()
    assert (run_dir / "linkage" / "linkage_summary.json").exists()
    assert (run_dir / "calibration" / "calibration_manifest.json").exists()
    assert (run_dir / "report" / "prediction_experiment_manifest.json").exists()
    assert (run_dir / "report" / "prediction_experiment_summary.json").exists()
    assert (run_dir / "report" / "prediction_experiment_report.md").exists()
    report_manifest = json.loads(
        (run_dir / "report" / "prediction_experiment_manifest.json").read_text(encoding="utf-8")
    )
    assert report_manifest["source_classification"] == "local_unreviewed_input"
    assert report_manifest["uses_real_experiment_data"] is None
    assert report_manifest["data_approval_status"] == "unreviewed"
    assert Path(result["paths"]["report_path"]).name == "prediction_experiment_report.md"
    assert export_experiment_feedback_report(run_dir).startswith(
        b"# Prediction vs Experiment"
    )
    assert load_experiment_feedback_run(run_dir)["run_name"] == "ui-contract"
    assert export_experiment_feedback_issues(run_dir, issue_kind="conflicts") == b""


def _experiment_csv_bytes(tmp_path) -> bytes:
    host = HostContext("Komagataella phaffii", "X33", "X33")
    condition = ConditionContext("BMMY", "methanol", "shake_flask", 30.0, 6.0, "250 rpm", 72.0)
    rows = (
        ("experiment", ExperimentRecord("UI-CONTROL", "hLF", host, "B01", condition, context_id="ctx-hlf-ui")),
        ("experiment", ExperimentRecord("UI-CANDIDATE", "hLF", host, "B01", condition, context_id="ctx-hlf-ui")),
        ("intervention", InterventionRecord("UI-CONTROL", "CONTROL", 1, InterventionType.CONTROL)),
        (
            "intervention",
            InterventionRecord(
                "UI-CANDIDATE",
                "KO-1",
                1,
                InterventionType.KO,
                gene_id="G1",
                construction_method="CRISPR-Cas9",
                prediction_run_id="ui-screen-run",
                evidence_id="hLF-KO-UI-1",
            ),
        ),
        ("measurement", _measurement("UI-CONTROL", "CONTROL-T1", 10.0)),
        ("measurement", _measurement("UI-CANDIDATE", "CANDIDATE-T1", 15.0)),
    )
    path = tmp_path / "source.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["record_type", "payload_json"])
        writer.writeheader()
        for record_type, record in rows:
            writer.writerow(
                {
                    "record_type": record_type,
                    "payload_json": json.dumps(asdict(record), default=_enum_value),
                }
            )
    return path.read_bytes()


def _measurement(experiment_id: str, measurement_id: str, value: float) -> MeasurementRecord:
    return MeasurementRecord(
        experiment_id,
        measurement_id,
        "titer",
        "ELISA",
        "extracellular",
        value,
        "mg/L",
        value,
        "mg/L",
        MeasurementStatus.VALID,
    )


def _enum_value(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    raise TypeError(type(value).__name__)
