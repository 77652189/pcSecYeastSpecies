from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Sequence

from pcsec_pichia.experimental_feedback.calibration import (
    CalibrationConfig,
    CalibrationOutputs,
    CalibrationRecord,
    CalibrationSummary,
    build_calibration_summary,
    write_calibration_outputs,
)
from pcsec_pichia.experimental_feedback.io import (
    ExperimentFeedbackOutputs,
    load_experiment_bundle,
    write_experiment_feedback_cache,
)
from pcsec_pichia.experimental_feedback.linkage import (
    PredictionLinkageResult,
    build_prediction_index,
    link_experiments_to_predictions,
)
from pcsec_pichia.experimental_feedback.quality import (
    ExperimentValidationResult,
    validate_experiment_bundle,
)
from pcsec_pichia.experimental_feedback.schema import SCHEMA_VERSION, SchemaValidationError


_SOURCE_CLASSIFICATIONS = {
    "sanitized_fixture",
    "approved_local_data",
    "local_unreviewed_input",
}


@dataclass(frozen=True)
class PredictionExperimentReportOutputs:
    manifest_path: Path
    summary_path: Path
    report_path: Path


@dataclass(frozen=True)
class ExperimentReplayOutputs:
    validated: ExperimentFeedbackOutputs
    calibration: CalibrationOutputs
    linkage_path: Path
    manifest_path: Path
    summary_path: Path
    report_path: Path


@dataclass(frozen=True)
class ExperimentReplayResult:
    validation: ExperimentValidationResult
    linkage: PredictionLinkageResult
    calibration: CalibrationSummary
    outputs: ExperimentReplayOutputs


def run_experiment_feedback_replay(
    *,
    experiment_path: str | Path,
    prediction_runs: Sequence[Mapping[str, Any]],
    output_dir: str | Path,
    source_classification: str,
    config: CalibrationConfig | None = None,
    prediction_sources: Sequence[str | Path] = (),
) -> ExperimentReplayResult:
    """Run an auditable import-to-report replay without promoting curated data."""

    if source_classification not in _SOURCE_CLASSIFICATIONS:
        raise SchemaValidationError(
            f"source_classification must be one of {sorted(_SOURCE_CLASSIFICATIONS)}."
        )
    if prediction_sources and len(prediction_sources) != len(prediction_runs):
        raise SchemaValidationError("prediction_sources must align with prediction_runs.")

    resolved_output = Path(output_dir)
    resolved_output.mkdir(parents=True, exist_ok=True)
    bundle = load_experiment_bundle(experiment_path)
    prediction_index = build_prediction_index(tuple(prediction_runs))
    linkage = link_experiments_to_predictions(bundle, prediction_index)
    linked_bundle = replace(bundle, prediction_links=linkage.links)
    validation = validate_experiment_bundle(linked_bundle)
    if not validation.is_valid:
        write_experiment_feedback_cache(linked_bundle, resolved_output / "validated")
        raise SchemaValidationError("replay requires a schema-valid imported experiment bundle.")

    validated_outputs = write_experiment_feedback_cache(
        linked_bundle,
        resolved_output / "validated",
    )
    linkage_dir = resolved_output / "linkage"
    linkage_dir.mkdir(parents=True, exist_ok=True)
    linkage_path = linkage_dir / "linkage_summary.json"
    linkage_path.write_text(
        json.dumps(_linkage_payload(linkage), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )

    calibration = build_calibration_summary(
        validation,
        linkage,
        config or CalibrationConfig(),
    )
    calibration_outputs = write_calibration_outputs(
        calibration,
        resolved_output / "calibration",
    )
    report_outputs = write_prediction_experiment_report_outputs(
        calibration=calibration,
        validation=validation,
        linkage=linkage,
        output_dir=resolved_output / "report",
        experiment_path=Path(experiment_path),
        prediction_runs=prediction_runs,
        prediction_sources=prediction_sources,
        source_classification=source_classification,
        supporting_files={
            "validated_records": validated_outputs.validated_records_path.name,
            "import_manifest": validated_outputs.manifest_path.name,
            "linkage": linkage_path.name,
            "calibration_records": calibration_outputs.records_path.name,
            "calibration_summary": calibration_outputs.summary_path.name,
            "calibration_manifest": calibration_outputs.manifest_path.name,
        },
    )
    return ExperimentReplayResult(
        validation=validation,
        linkage=linkage,
        calibration=calibration,
        outputs=ExperimentReplayOutputs(
            validated=validated_outputs,
            calibration=calibration_outputs,
            linkage_path=linkage_path,
            manifest_path=report_outputs.manifest_path,
            summary_path=report_outputs.summary_path,
            report_path=report_outputs.report_path,
        ),
    )


def write_prediction_experiment_report_outputs(
    *,
    calibration: CalibrationSummary,
    validation: ExperimentValidationResult,
    linkage: PredictionLinkageResult,
    output_dir: str | Path,
    experiment_path: str | Path,
    prediction_runs: Sequence[Mapping[str, Any]],
    source_classification: str,
    prediction_sources: Sequence[str | Path] = (),
    supporting_files: Mapping[str, str] | None = None,
) -> PredictionExperimentReportOutputs:
    if source_classification not in _SOURCE_CLASSIFICATIONS:
        raise SchemaValidationError(
            f"source_classification must be one of {sorted(_SOURCE_CLASSIFICATIONS)}."
        )
    if prediction_sources and len(prediction_sources) != len(prediction_runs):
        raise SchemaValidationError("prediction_sources must align with prediction_runs.")
    report_dir = Path(output_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    summary_path = report_dir / "prediction_experiment_summary.json"
    report_path = report_dir / "prediction_experiment_report.md"
    manifest_path = report_dir / "prediction_experiment_manifest.json"
    summary_payload = _summary_payload(calibration)
    summary_path.write_text(
        json.dumps(summary_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    report_path.write_text(
        _render_report(summary_payload, source_classification=source_classification),
        encoding="utf-8",
    )
    manifest_path.write_text(
        json.dumps(
            _manifest_payload(
                experiment_path=Path(experiment_path),
                prediction_runs=prediction_runs,
                prediction_sources=prediction_sources,
                source_classification=source_classification,
                validation=validation,
                linkage=linkage,
                calibration=calibration,
                supporting_files=supporting_files or {},
                summary_path=summary_path,
                report_path=report_path,
            ),
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return PredictionExperimentReportOutputs(
        manifest_path=manifest_path,
        summary_path=summary_path,
        report_path=report_path,
    )


def _summary_payload(calibration: CalibrationSummary) -> dict[str, Any]:
    status_counts: Counter[str] = Counter()
    targets: dict[str, dict[str, Any]] = {}
    for record in calibration.records:
        for status in record.measurement_statuses:
            if status != "valid":
                status_counts[status] += 1
        if record.eligibility_status == "eligible" and record.hit is False:
            status_counts["negative_observation"] += 1
        for reason in record.ineligibility_reasons:
            status_counts[reason] += 1

    for target in calibration.targets:
        records = [record for record in calibration.records if record.target_id == target.target_id]
        consistent = sum(record.direction_consistent is True for record in records)
        discordant = sum(record.direction_consistent is False for record in records)
        targets[target.target_id] = {
            "eligible_count": target.eligible_count,
            "ineligible_count": target.ineligible_count,
            "comparable_rank_pair_count": target.comparable_rank_pair_count,
            "minimum_rank_pairs": calibration.config.minimum_rank_pairs,
            "direction_consistent_count": consistent,
            "direction_discordant_count": discordant,
            "direction_consistency_rate": target.direction_consistency_rate,
            "rank_correlation": target.rank_correlation,
            "ranking_assessment": target.ranking_assessment,
            "top_k_metrics": _json_ready([asdict(item) for item in target.top_k_metrics]),
            "evidence_tier_metrics": _json_ready(
                [asdict(item) for item in target.evidence_tier_metrics]
            ),
            "confirmed_direction_candidates": [
                _candidate_payload(record)
                for record in records
                if record.direction_consistent is True
            ],
            "direction_discordant_candidates": [
                _candidate_payload(record)
                for record in records
                if record.direction_consistent is False
            ],
            "unresolved_candidates": [
                _candidate_payload(record)
                for record in records
                if record.eligibility_status != "eligible"
            ],
        }
    return {
        "schema_version": SCHEMA_VERSION,
        "targets": targets,
        "preserved_status_counts": dict(sorted(status_counts.items())),
        "record_count": len(calibration.records),
        "eligible_count": sum(
            record.eligibility_status == "eligible" for record in calibration.records
        ),
        "ineligible_count": sum(
            record.eligibility_status != "eligible" for record in calibration.records
        ),
        "next_round_guidance": (
            "Retain direction-consistent candidates for human review; inspect discordant and "
            "ineligible records before any ranking change. No tier is changed automatically."
        ),
    }


def _candidate_payload(record: CalibrationRecord) -> dict[str, Any]:
    return {
        "experiment_id": record.experiment_id,
        "intervention_id": record.intervention_id,
        "gene_id": record.gene_id,
        "intervention_type": record.intervention_type,
        "evidence_id": record.evidence_id,
        "prediction_run_id": record.prediction_run_id,
        "prediction_rank": record.prediction_rank,
        "recommendation_tier": record.recommendation_tier,
        "predicted_direction": record.predicted_direction,
        "observed_ratio": record.observed_ratio,
        "observed_direction": record.observed_direction,
        "eligibility_status": record.eligibility_status,
        "ineligibility_reasons": list(record.ineligibility_reasons),
        "measurement_statuses": list(record.measurement_statuses),
    }


def _manifest_payload(
    *,
    experiment_path: Path,
    prediction_runs: Sequence[Mapping[str, Any]],
    prediction_sources: Sequence[str | Path],
    source_classification: str,
    validation: ExperimentValidationResult,
    linkage: PredictionLinkageResult,
    calibration: CalibrationSummary,
    supporting_files: Mapping[str, str],
    summary_path: Path,
    report_path: Path,
) -> dict[str, Any]:
    sources: list[dict[str, Any]] = []
    for index, prediction_run in enumerate(prediction_runs):
        path = Path(prediction_sources[index]) if prediction_sources else None
        sources.append(
            {
                "prediction_run_id": str(prediction_run.get("prediction_run_id") or ""),
                "source_file": str(path) if path else "",
                "source_sha256": _sha256(path) if path else "",
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_classification": source_classification,
        "uses_real_experiment_data": (
            True
            if source_classification == "approved_local_data"
            else False
            if source_classification == "sanitized_fixture"
            else None
        ),
        "data_approval_status": {
            "approved_local_data": "approved",
            "sanitized_fixture": "not_applicable",
            "local_unreviewed_input": "unreviewed",
        }[source_classification],
        "experiment_source": {
            "source_file": str(experiment_path),
            "source_sha256": _sha256(experiment_path),
        },
        "prediction_sources": sources,
        "validation": {
            "is_valid": validation.is_valid,
            "error_count": len(validation.errors),
            "warning_count": len(validation.warnings),
        },
        "linkage": {
            "matched_count": linkage.matched_count,
            "ambiguous_count": linkage.ambiguous_count,
            "missing_prediction_count": linkage.missing_prediction_count,
            "context_mismatch_count": linkage.context_mismatch_count,
        },
        "calibration": {
            "record_count": len(calibration.records),
            "eligible_count": sum(
                record.eligibility_status == "eligible" for record in calibration.records
            ),
            "ineligible_count": sum(
                record.eligibility_status != "eligible" for record in calibration.records
            ),
            "config": _json_ready(asdict(calibration.config)),
        },
        "files": {
            **dict(supporting_files),
            "summary": summary_path.name,
            "report": report_path.name,
        },
        "promotes_curated_data": False,
        "mutates_recommendation_tier": False,
        "mutates_model_constraints": False,
        "sends_raw_experiment_records_to_llm": False,
    }


def _render_report(payload: Mapping[str, Any], *, source_classification: str) -> str:
    lines = [
        "# Prediction vs Experiment 回放报告",
        "",
        "本报告由程序从结构化实验反馈生成，用于研发复核，不是绝对产量预测或实验成功率承诺。",
        "",
    ]
    if source_classification == "sanitized_fixture":
        lines.extend(("本次回放未使用真实实验数据，仅使用脱敏 fixture。", ""))
    elif source_classification == "local_unreviewed_input":
        lines.extend(("本次输入仅在本地处理，数据审批与脱敏状态尚未复核。", ""))
    for target_id in ("hLF", "OPN"):
        target = (payload.get("targets") or {}).get(target_id)
        lines.extend((f"## {target_id}", ""))
        if not target:
            lines.extend(("无可报告记录。", ""))
            continue
        lines.extend(
            (
                f"- Eligible: {target['eligible_count']}",
                f"- 不可校准: {target['ineligible_count']}",
                f"- 方向一致: {target['direction_consistent_count']}",
                f"- 方向不一致: {target['direction_discordant_count']}",
                f"- 排序证据状态: {target['ranking_assessment']}",
                f"- 可比较排名对: {target['comparable_rank_pair_count']}/{target['minimum_rank_pairs']}",
                "",
                "### 模型推荐与实验观察",
                "",
                "| experiment | gene | intervention | evidence | rank | predicted | observed | ratio | eligibility | reasons/status |",
                "| --- | --- | --- | --- | ---: | --- | --- | ---: | --- | --- |",
            )
        )
        candidates = (
            target["confirmed_direction_candidates"]
            + target["direction_discordant_candidates"]
            + target["unresolved_candidates"]
        )
        for item in candidates:
            reasons = list(item["ineligibility_reasons"])
            reasons.extend(status for status in item["measurement_statuses"] if status != "valid")
            ratio = "N/A" if item["observed_ratio"] is None else f"{item['observed_ratio']:.3f}"
            lines.append(
                "| {experiment_id} | {gene_id} | {intervention_type} | {evidence_id} | {rank} | {predicted} | {observed} | "
                "{ratio} | {eligibility} | {reasons} |".format(
                    experiment_id=item["experiment_id"],
                    gene_id=item["gene_id"],
                    intervention_type=item["intervention_type"],
                    evidence_id=item["evidence_id"] or "N/A",
                    rank=item["prediction_rank"] if item["prediction_rank"] is not None else "N/A",
                    predicted=item["predicted_direction"] or "N/A",
                    observed=item["observed_direction"] or "N/A",
                    ratio=ratio,
                    eligibility=item["eligibility_status"],
                    reasons=", ".join(reasons) or "-",
                )
            )
        lines.extend(("", "### 下一轮建议", ""))
        if target["confirmed_direction_candidates"]:
            lines.append("- 保留方向一致候选进入人工复核，不自动提升 recommendation tier。")
        if target["direction_discordant_candidates"]:
            lines.append("- 复核方向不一致候选及对照条件，再决定是否人工调整后续排序。")
        if target["unresolved_candidates"]:
            lines.append("- 先解决不可校准记录的 linkage、assay 或 control 问题，再纳入统计。")
        lines.append("")
    lines.extend(
        (
            "## 明确边界",
            "",
            "- 回放结果不会自动修改 recommendation tier 或模型约束。",
            "- 原始实验记录未发送给 LLM，也未提升为 curated scientific asset。",
            "- 样本量不足时只报告描述性结果，不声称排序已被充分验证。",
            "",
        )
    )
    return "\n".join(lines)


def _linkage_payload(linkage: PredictionLinkageResult) -> dict[str, Any]:
    return {
        "matched_count": linkage.matched_count,
        "ambiguous_count": linkage.ambiguous_count,
        "missing_prediction_count": linkage.missing_prediction_count,
        "context_mismatch_count": linkage.context_mismatch_count,
        "control_count": linkage.control_count,
        "links": [_json_ready(asdict(link)) for link in linkage.links],
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_ready(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "ExperimentReplayOutputs",
    "ExperimentReplayResult",
    "PredictionExperimentReportOutputs",
    "run_experiment_feedback_replay",
    "write_prediction_experiment_report_outputs",
]
