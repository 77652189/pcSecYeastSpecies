from __future__ import annotations

import csv
import json
from dataclasses import asdict
from enum import Enum
from pathlib import Path

from openpyxl import Workbook

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
    assert result["calibration"]["targets"][0]["comparable_rank_pair_count"] == 1
    assert result["calibration"]["targets"][0]["ranking_assessment"] == "insufficient_evidence"
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


def test_service_accepts_xlsx_experiment_uploads(tmp_path) -> None:
    csv_bytes = _experiment_csv_bytes(tmp_path)
    reader = csv.DictReader(csv_bytes.decode("utf-8").splitlines())
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "records"
    worksheet.append(("record_type", "payload_json"))
    for row in reader:
        worksheet.append((row["record_type"], row["payload_json"]))
    source_path = tmp_path / "sanitized_ui.xlsx"
    workbook.save(source_path)

    result = submit_experiment_feedback_import(
        experiment_filename=source_path.name,
        experiment_bytes=source_path.read_bytes(),
        run_name="ui-xlsx-contract",
        output_root=tmp_path / "runs",
    )

    assert result["validation"]["is_valid"] is True
    assert result["linkage"]["missing_prediction_count"] == 1
    run_dir = tmp_path / "runs" / "ui-xlsx-contract"
    assert (run_dir / "inbox" / source_path.name).exists()
    assert (run_dir / "report" / "prediction_experiment_report.md").exists()


def test_service_passes_form_metadata_to_wide_template_core_adapter(tmp_path) -> None:
    fixture = (
        Path(__file__).resolve().parents[1]
        / "python_pichia"
        / "tests"
        / "fixtures"
        / "experimental_feedback"
        / "fermentation_template_sanitized.csv"
    )

    result = submit_experiment_feedback_import(
        experiment_filename=fixture.name,
        experiment_bytes=fixture.read_bytes(),
        run_name="wide-template-contract",
        output_root=tmp_path / "runs",
        experiment_metadata={"target_id": "hLF", "batch_id": "B01"},
    )

    assert result["validation"]["is_valid"] is True
    assert result["calibration"]["available"] is True
    assert len(result["calibration"]["records"]) == 6
    run_dir = tmp_path / "runs" / "wide-template-contract"
    manifest = json.loads(
        (run_dir / "validated" / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["adapter_id"] == "pcsec_pichia.fermentation_template.v1"
    assert manifest["import_metadata"] == {"batch_id": "B01", "target_id": "hLF"}
    assert (run_dir / "inbox" / fixture.name).read_bytes() == fixture.read_bytes()


def test_service_remains_a_facade_without_template_science_rules() -> None:
    service_source = (
        Path(__file__).resolve().parents[1]
        / "app"
        / "services"
        / "pichia_experiment_feedback_service.py"
    ).read_text(encoding="utf-8")

    assert "experiment_metadata" in service_source
    for core_rule in ("改造方案", "亲本对照组编号", "培养失败", "contamination"):
        assert core_rule not in service_source


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
