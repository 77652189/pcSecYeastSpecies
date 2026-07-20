from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

from openpyxl import Workbook
import pytest

from pcsec_pichia.experimental_feedback import (
    CalibrationConfig,
    FERMENTATION_TEMPLATE_ADAPTER_ID,
    FermentationDataStatus,
    ExperimentValidationResult,
    MeasurementStatus,
    ModificationConfirmationStatus,
    PredictionLinkStatus,
    SchemaValidationError,
    QualityStatus,
    build_calibration_summary,
    build_prediction_index,
    link_experiments_to_predictions,
    load_experiment_bundle,
    run_experiment_feedback_replay,
    validate_experiment_bundle,
    write_experiment_feedback_cache,
)


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "experimental_feedback"
WIDE_CSV = FIXTURE_ROOT / "fermentation_template_sanitized.csv"
IMPORT_METADATA = {"target_id": "hLF", "batch_id": "B01"}


def _extracellular_titer(bundle, experiment_id: str):
    return next(
        item
        for item in bundle.measurements
        if item.experiment_id == experiment_id
        and item.assay_type == "titer"
        and item.compartment == "extracellular"
    )


def test_wide_csv_import_preserves_direction1_fields_statuses_and_raw_trace() -> None:
    bundle = load_experiment_bundle(WIDE_CSV, metadata=IMPORT_METADATA)
    validation = validate_experiment_bundle(bundle)

    assert validation.is_valid is True
    assert bundle.import_manifest is not None
    assert bundle.import_manifest.adapter_id == FERMENTATION_TEMPLATE_ADAPTER_ID
    assert json.loads(bundle.import_manifest.metadata_json) == IMPORT_METADATA
    assert len(bundle.experiments) == len(bundle.interventions) == 7
    # each of the 14 real fields maps into exactly 4 measurement slots per row:
    # OD600, UPR, extracellular titer, intracellular titer.
    assert len(bundle.measurements) == 28

    clone_a = [item for item in bundle.experiments if item.clone_id == "CLONE-A"]
    assert {item.biological_replicate_id for item in clone_a} == {"R1", "R2"}
    assert len({item.experiment_id for item in clone_a}) == 2
    assert {item.fermentation_data_status for item in bundle.experiments} == {
        FermentationDataStatus.NORMAL,
        FermentationDataStatus.CONTAMINATION,
        FermentationDataStatus.CULTURE_FAILED,
        FermentationDataStatus.ASSAY_FAILED,
        FermentationDataStatus.OTHER_EXCLUDED,
    }

    by_status = {
        experiment.fermentation_data_status: _extracellular_titer(bundle, experiment.experiment_id)
        for experiment in bundle.experiments
        if experiment.fermentation_data_status is not FermentationDataStatus.NORMAL
    }
    assert by_status[FermentationDataStatus.CONTAMINATION].raw_value == 14.0
    assert by_status[FermentationDataStatus.CONTAMINATION].status is MeasurementStatus.EXCLUDED
    assert by_status[FermentationDataStatus.CULTURE_FAILED].raw_value == 12.0
    assert by_status[FermentationDataStatus.ASSAY_FAILED].status is MeasurementStatus.ASSAY_FAILED
    assert by_status[FermentationDataStatus.OTHER_EXCLUDED].raw_value == 9.0
    assert all(item.source_row_number is not None for item in bundle.measurements)
    raw_fields = json.loads(by_status[FermentationDataStatus.CONTAMINATION].raw_fields_json)
    assert raw_fields["研发备注"] == "同克隆独立培养"
    assert "unmapped_template_column:研发备注" in bundle.warnings

    contaminated_experiment = next(
        item
        for item in bundle.experiments
        if item.fermentation_data_status is FermentationDataStatus.CONTAMINATION
    )
    assert contaminated_experiment.notes == "疑似轻微异味"
    confirmed_intervention = next(
        item for item in bundle.interventions if item.experiment_id == contaminated_experiment.experiment_id
    )
    assert confirmed_intervention.confirmation_status is ModificationConfirmationStatus.CONFIRMED_SUCCESS
    assert confirmed_intervention.confirmation_method == "测序"


def test_wide_xlsx_uses_file_metadata_and_rejects_form_conflict(tmp_path) -> None:
    xlsx_path = _write_wide_xlsx(tmp_path / "fermentation.xlsx")

    bundle = load_experiment_bundle(xlsx_path)

    assert validate_experiment_bundle(bundle).is_valid is True
    assert bundle.import_manifest is not None
    assert json.loads(bundle.import_manifest.metadata_json) == IMPORT_METADATA
    assert {item.target_id for item in bundle.experiments} == {"hLF"}
    assert {item.batch_id for item in bundle.experiments} == {"B01"}
    assert {item.source_sheet for item in bundle.measurements} == {"records"}

    with pytest.raises(SchemaValidationError, match="metadata conflict: target_id"):
        load_experiment_bundle(xlsx_path, metadata={"target_id": "OPN"})

    duplicate_path = _write_wide_xlsx(
        tmp_path / "duplicate-metadata.xlsx",
        metadata_rows=[("key", "value"), ("target_id", "hLF"), ("target_id", "OPN")],
    )
    with pytest.raises(SchemaValidationError, match="metadata conflict in XLSX sheet"):
        load_experiment_bundle(duplicate_path)


def test_wide_jsonl_cache_roundtrip_preserves_source_row_and_unknown_columns(tmp_path) -> None:
    rows = _fixture_rows()
    jsonl_path = tmp_path / "fermentation.jsonl"
    jsonl_path.write_text(
        "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows),
        encoding="utf-8",
    )
    bundle = load_experiment_bundle(jsonl_path, metadata=IMPORT_METADATA)
    outputs = write_experiment_feedback_cache(bundle, tmp_path / "cache")

    reloaded = load_experiment_bundle(outputs.validated_records_path)

    assert validate_experiment_bundle(reloaded).is_valid is True
    contaminated = next(
        item
        for item in reloaded.experiments
        if item.fermentation_data_status is FermentationDataStatus.CONTAMINATION
    )
    measurement = next(
        item for item in reloaded.measurements if item.experiment_id == contaminated.experiment_id
    )
    assert measurement.source_row_number == 3
    assert json.loads(measurement.raw_fields_json)["研发备注"] == "同克隆独立培养"
    manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    assert manifest["adapter_id"] == FERMENTATION_TEMPLATE_ADAPTER_ID
    assert manifest["import_metadata"] == IMPORT_METADATA


def test_bad_status_is_rejected_but_missing_gene_is_retained_for_audit(tmp_path) -> None:
    bad_status = _minimal_row(
        clone_id="BAD-STATUS",
        plan="KO:G1",
        status="未知状态",
        group="C1",
        repeat="R1",
        value="1",
    )
    bad_status_path = _write_wide_csv(tmp_path / "bad-status.csv", [bad_status])
    with pytest.raises(SchemaValidationError, match="invalid data_status"):
        load_experiment_bundle(bad_status_path, metadata=IMPORT_METADATA)

    missing_gene = _minimal_row(
        clone_id="MISSING-GENE",
        plan="KO",
        status="正常",
        group="C1",
        repeat="R1",
        value="1",
    )
    missing_gene_path = _write_wide_csv(tmp_path / "missing-gene.csv", [missing_gene])
    bundle = load_experiment_bundle(missing_gene_path, metadata=IMPORT_METADATA)
    validation = validate_experiment_bundle(bundle)
    linkage = link_experiments_to_predictions(bundle, build_prediction_index(()))

    assert validation.is_valid is False
    assert any("gene_id" in issue.message for issue in validation.errors)
    assert bundle.interventions[0].gene_id == ""
    assert json.loads(bundle.measurements[0].raw_fields_json)["改造方案"] == "KO"
    assert linkage.links[0].status is PredictionLinkStatus.AMBIGUOUS
    assert linkage.links[0].reason == "gene_id_missing"


def test_duplicate_generated_identity_becomes_conflict_instead_of_overwrite(tmp_path) -> None:
    first = _minimal_row(
        clone_id="DUPLICATE",
        plan="KO:G1",
        status="正常",
        group="C1",
        repeat="R1",
        value="10",
    )
    second = dict(first, **{"72h-胞外产量mg/L": "20"})
    path = _write_wide_csv(tmp_path / "duplicate.csv", [first, second])

    bundle = load_experiment_bundle(path, metadata=IMPORT_METADATA)
    validation = validate_experiment_bundle(bundle)

    assert validation.is_valid is False
    # raw_fields_json captures the whole source row, so the changed extracellular-titer
    # value makes all 4 of the row's measurement slots differ from their row-1 counterparts.
    assert len(bundle.import_conflicts) == 4
    assert all(item.record_type == "measurement" for item in bundle.import_conflicts)
    extracellular_conflict = next(
        item
        for item in bundle.import_conflicts
        if item.record_id.endswith("titer-extracellular")
    )
    assert "10.0" in extracellular_conflict.first_payload_json
    assert "20.0" in extracellular_conflict.conflicting_payload_json
    # the first row's full set of 4 measurements survives; the conflicting second-row
    # measurements are recorded as conflicts, not appended.
    assert len(bundle.measurements) == 4


def test_identical_duplicate_source_row_is_warning_not_content_conflict(tmp_path) -> None:
    row = _minimal_row("DUPLICATE", "KO:G1", "正常", "C1", "R1", "10")
    path = _write_wide_csv(tmp_path / "identical-duplicate.csv", [row, dict(row)])

    bundle = load_experiment_bundle(path, metadata=IMPORT_METADATA)

    assert bundle.import_conflicts == ()
    assert len(bundle.experiments) == len(bundle.interventions) == 1
    assert len(bundle.measurements) == 4
    assert any("duplicate_record_ignored:measurement" in item for item in bundle.warnings)


def test_generated_identity_includes_control_group_and_condition_scope(tmp_path) -> None:
    first = _minimal_row("LOCAL-REPEAT", "KO:G1", "正常", "C1", "R1", "10")
    second = dict(first, **{"亲本对照组编号": "C2"})
    path = _write_wide_csv(tmp_path / "identity-scope.csv", [first, second])

    bundle = load_experiment_bundle(path, metadata=IMPORT_METADATA)

    assert len(bundle.experiments) == 2
    assert len({item.experiment_id for item in bundle.experiments}) == 2
    assert bundle.import_conflicts == ()


def test_non_normal_canonical_record_cannot_bypass_calibration_gate(tmp_path) -> None:
    rows = [
        _minimal_row("CONTROL", "亲本对照", "正常", "C1", "R1", "10"),
        _minimal_row(
            "CONTAMINATED",
            "KO:G1",
            "污染",
            "C1",
            "R1",
            "15",
            prediction_run_id="template-run",
            evidence_id="HLF-G1",
        ),
    ]
    path = _write_wide_csv(tmp_path / "canonical-bypass.csv", rows)
    bundle = load_experiment_bundle(path, metadata=IMPORT_METADATA)
    contaminated = next(
        item for item in bundle.experiments if item.clone_id == "CONTAMINATED"
    )
    forged_bundle = replace(
        bundle,
        experiments=tuple(
            replace(item, quality_status=QualityStatus.VALID)
            if item.experiment_id == contaminated.experiment_id
            else item
            for item in bundle.experiments
        ),
        measurements=tuple(
            replace(
                item,
                status=MeasurementStatus.VALID,
                excluded=False,
                exclusion_reason="",
                status_reason="",
                canonical_value=item.raw_value,
            )
            if item.experiment_id == contaminated.experiment_id
            else item
            for item in bundle.measurements
        ),
    )
    validation = validate_experiment_bundle(forged_bundle)
    linkage = link_experiments_to_predictions(
        forged_bundle,
        build_prediction_index((_prediction_run(("G1",)),)),
    )
    forced_validation = ExperimentValidationResult(
        bundle=forged_bundle,
        is_valid=True,
    )
    calibration = build_calibration_summary(
        forced_validation,
        linkage,
        CalibrationConfig(minimum_rank_pairs=2),
    )

    assert validation.is_valid is False
    assert any("non-normal fermentation_data_status" in item.message for item in validation.errors)
    assert calibration.records[0].eligibility_status == "ineligible"
    assert calibration.records[0].ineligibility_reasons == (
        "fermentation_data_status:contamination",
    )


def test_parent_control_group_prevents_pooling_and_wide_xlsx_replays_end_to_end(tmp_path) -> None:
    rows = [
        _minimal_row("CONTROL-C1", "亲本对照", "正常", "C1", "R1", "10"),
        _minimal_row("CONTROL-C2", "亲本对照", "正常", "C2", "R1", "100"),
        _minimal_row(
            "CANDIDATE-C1",
            "KO:G1",
            "正常",
            "C1",
            "R1",
            "15",
            prediction_run_id="template-run",
            evidence_id="HLF-G1",
        ),
    ]
    path = _write_wide_xlsx(tmp_path / "control-groups.xlsx", rows=rows)
    bundle = load_experiment_bundle(path)
    prediction_run = _prediction_run(("G1",))
    linkage = link_experiments_to_predictions(bundle, build_prediction_index((prediction_run,)))
    linked = replace(bundle, prediction_links=linkage.links)
    validation = validate_experiment_bundle(linked)
    calibration = build_calibration_summary(validation, linkage, CalibrationConfig(minimum_rank_pairs=2))

    assert validation.is_valid is True
    assert len(calibration.records) == 1
    assert calibration.records[0].control_value == 10.0
    assert calibration.records[0].observed_ratio == 1.5
    assert len(calibration.records[0].control_experiment_ids) == 1

    replay_path = _write_wide_xlsx(tmp_path / "full-replay.xlsx")
    replay_prediction = _prediction_run(("G1", "G2", "G3", "G4", "G5"))
    result = run_experiment_feedback_replay(
        experiment_path=replay_path,
        prediction_runs=(replay_prediction,),
        output_dir=tmp_path / "replay",
        source_classification="sanitized_fixture",
        config=CalibrationConfig(minimum_rank_pairs=2),
    )

    assert result.validation.is_valid is True
    assert result.linkage.matched_count == 6
    assert sum(record.eligibility_status == "eligible" for record in result.calibration.records) == 2
    summary = json.loads(result.outputs.summary_path.read_text(encoding="utf-8"))
    for status in ("contamination", "culture_failed", "assay_failed", "other_excluded"):
        assert summary["preserved_status_counts"][f"fermentation_data_status:{status}"] == 1
    negative = next(
        item
        for item in summary["targets"]["hLF"]["direction_discordant_candidates"]
        if item["gene_id"] == "G5"
    )
    assert negative["observed_ratio"] == 0.8
    report = result.outputs.report_path.read_text(encoding="utf-8")
    assert "contamination" in report
    assert "culture_failed" in report
    assert "assay_failed" in report
    assert "other_excluded" in report
    assert "未使用真实实验数据" in report


def test_failed_control_is_preserved_in_summary_and_markdown(tmp_path) -> None:
    rows = [
        _minimal_row("CONTROL-FAILED", "亲本对照", "污染", "C1", "R1", "10"),
        _minimal_row(
            "CANDIDATE",
            "KO:G1",
            "正常",
            "C1",
            "R1",
            "15",
            prediction_run_id="template-run",
            evidence_id="HLF-G1",
        ),
    ]
    path = _write_wide_xlsx(tmp_path / "failed-control.xlsx", rows=rows)
    result = run_experiment_feedback_replay(
        experiment_path=path,
        prediction_runs=(_prediction_run(("G1",)),),
        output_dir=tmp_path / "failed-control-replay",
        source_classification="sanitized_fixture",
        config=CalibrationConfig(minimum_rank_pairs=2),
    )

    summary = json.loads(result.outputs.summary_path.read_text(encoding="utf-8"))
    preserved = summary["preserved_experiment_evidence"]
    failed_control = next(item for item in preserved if item["role"] == "control")
    assert failed_control["fermentation_data_status"] == "contamination"
    extracellular_titer = next(
        item for item in failed_control["measurements"] if item["measurement_id"] == "titer-extracellular"
    )
    assert extracellular_titer["raw_value"] == 10.0
    assert result.calibration.records[0].ineligibility_reasons == ("control_match_missing",)
    report = result.outputs.report_path.read_text(encoding="utf-8")
    assert "## 保留的失败与排除实验" in report
    assert "CONTROL-FAILED" in report
    assert "control" in report
    assert "contamination" in report


def test_excluded_zero_value_is_preserved_only_as_raw_evidence(tmp_path) -> None:
    row = _minimal_row("ZERO-EXCLUDED", "KO:G1", "污染", "C1", "R1", "0")
    path = _write_wide_csv(tmp_path / "excluded-zero.csv", [row])

    bundle = load_experiment_bundle(path, metadata=IMPORT_METADATA)
    validation = validate_experiment_bundle(bundle)
    measurement = _extracellular_titer(bundle, bundle.experiments[0].experiment_id)

    assert validation.is_valid is True
    assert measurement.raw_value == 0.0
    assert measurement.canonical_value is None
    assert measurement.status is MeasurementStatus.EXCLUDED


def test_above_range_flag_excludes_canonical_value_and_preserves_audit_trail(tmp_path) -> None:
    row = _minimal_row("ABOVE-RANGE", "KO:G1", "正常", "C1", "R1", "9999")
    row["胞外是否超标曲"] = "是"
    row["胞外ELISA稀释倍数"] = "200"
    path = _write_wide_csv(tmp_path / "above-range.csv", [row])

    bundle = load_experiment_bundle(path, metadata=IMPORT_METADATA)
    validation = validate_experiment_bundle(bundle)
    measurement = _extracellular_titer(bundle, bundle.experiments[0].experiment_id)

    assert validation.is_valid is True
    assert measurement.status is MeasurementStatus.ABOVE_RANGE
    assert measurement.raw_value == 9999.0
    assert measurement.canonical_value is None
    assert measurement.dilution_factor == 200.0
    assert measurement.status_reason


def _fixture_rows() -> list[dict[str, str]]:
    with WIDE_CSV.open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_wide_xlsx(
    path: Path,
    *,
    rows: list[dict[str, str]] | None = None,
    metadata_rows: list[tuple[str, str]] | None = None,
) -> Path:
    source_rows = rows or _fixture_rows()
    workbook = Workbook()
    records = workbook.active
    records.title = "records"
    headers = list(source_rows[0])
    records.append(headers)
    for row in source_rows:
        records.append([row.get(header) for header in headers])
    metadata = workbook.create_sheet("metadata")
    for metadata_row in metadata_rows or [
        ("key", "value"),
        ("target_id", IMPORT_METADATA["target_id"]),
        ("batch_id", IMPORT_METADATA["batch_id"]),
    ]:
        metadata.append(metadata_row)
    workbook.save(path)
    return path


def _write_wide_csv(path: Path, rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _minimal_row(
    clone_id: str,
    plan: str,
    status: str,
    group: str,
    repeat: str,
    value: str,
    *,
    prediction_run_id: str = "",
    evidence_id: str = "",
) -> dict[str, str]:
    return {
        "克隆编号": clone_id,
        "宿主物种": "Komagataella phaffii",
        "发酵菌株": "X33",
        "本底菌株": "X33",
        "发酵条件": "BMMY+methanol+shake_flask+250rpm",
        "取样时间_h": "72",
        "72h-OD600": "50",
        "72h-UPR": "1.0",
        "72h-胞外产量mg/L": value,
        "胞外ELISA稀释倍数": "5",
        "胞外是否超标曲": "否",
        "72h-胞内产量mg/L": "5",
        "胞内ELISA稀释倍数": "5",
        "胞内是否超标曲": "否",
        "备注": "",
        "改造方案": plan,
        "改造确认": "",
        "确认方式": "",
        "数据状态": status,
        "亲本对照组编号": group,
        "重复编号": repeat,
        "预测条件ID": "ctx-hlf-template",
        "预测Run编号": prediction_run_id,
        "证据编号": evidence_id,
        "状态原因": "",
        "研发备注": "脱敏测试",
    }


def _prediction_run(genes: tuple[str, ...]) -> dict[str, object]:
    return {
        "prediction_run_id": "template-run",
        "evidence_items": [
            {
                "evidence_id": f"HLF-{gene}",
                "target_id": "hLF",
                "gene_id": gene,
                "intervention_type": "KO",
                "context_id": "ctx-hlf-template",
                "rank": index,
                "predicted_direction": "increase",
                "evidence_tier": "sanitized",
            }
            for index, gene in enumerate(genes, start=1)
        ],
    }
