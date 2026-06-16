from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pandas as pd

from app.adapters.uspnet import USPNetAdapter
from app.core.signal_peptides import UniProtCandidateLibraryResult
from app.core.paths import ProjectPaths
from app.services.opn import OpnCandidateCatalog
from app.services.opn_signal_peptides import OpnSignalPeptideCandidateSource
from app.services.signal_peptide_exports import (
    write_candidate_fasta,
    write_csv,
    write_json,
    write_signal_peptide_fasta,
)
from app.services.signal_peptide_library import SignalPeptideLibraryService
from app.services.signal_peptide_rules import score_signal_peptide


UNIPROT_CANDIDATES_CSV = "uniprot_candidates.csv"
UNIPROT_DUPLICATES_CSV = "uniprot_duplicate_candidates.csv"
UNIPROT_DISCOVERY_SUMMARY_JSON = "uniprot_candidate_discovery_summary.json"
METHOD_INPUT_FASTA = "method_comparison_input.fasta"
METHOD_COMPARISON_CSV = "signal_peptide_method_comparison.csv"
RECOMMENDED_FASTA = "method_recommended_candidates.fasta"
METHOD_SUMMARY_JSON = "signal_peptide_method_comparison_summary.json"


@dataclass(frozen=True)
class SignalPeptideScreeningResult:
    available: bool
    success: bool
    message: str
    summary: dict[str, object]
    rows: list[dict[str, object]]
    output_dir: Path
    uniprot_csv: Path | None = None
    duplicate_csv: Path | None = None
    input_fasta: Path | None = None
    comparison_csv: Path | None = None
    recommended_fasta: Path | None = None
    uspnet_raw_dir: Path | None = None
    summary_json: Path | None = None
    errors: list[str] | None = None


class SignalPeptideScreeningService:
    def __init__(
        self,
        paths: ProjectPaths,
        *,
        library_service: SignalPeptideLibraryService | None = None,
        uspnet_adapter: USPNetAdapter | None = None,
    ) -> None:
        self.paths = paths
        self.library_service = library_service or SignalPeptideLibraryService(
            OpnSignalPeptideCandidateSource(OpnCandidateCatalog(paths)).list_candidates()
        )
        self.uspnet_adapter = uspnet_adapter or USPNetAdapter(repo_dir=paths.uspnet_repo)

    def discover_and_persist_uniprot_candidates(
        self,
        *,
        taxon_id: int = 4922,
        max_records: int = 300,
        reviewed_only: bool = False,
        exclude_existing: bool = True,
    ) -> UniProtCandidateLibraryResult:
        discovery = self.library_service.discover_uniprot_candidate_library(
            taxon_id=taxon_id,
            max_records=max_records,
            reviewed_only=reviewed_only,
            exclude_existing=exclude_existing,
        )
        paths = self._output_paths()
        paths["output_dir"].mkdir(parents=True, exist_ok=True)
        self._persist_uniprot_discovery(
            discovery,
            taxon_id=taxon_id,
            max_records=max_records,
            reviewed_only=reviewed_only,
            exclude_existing=exclude_existing,
        )
        return discovery

    def load_persisted_screening_result(self) -> SignalPeptideScreeningResult | None:
        paths = self._output_paths()
        summary_json = paths["summary_json"]
        comparison_csv = paths["comparison_csv"]
        if not summary_json.exists() or not comparison_csv.exists():
            return None
        try:
            payload = json.loads(summary_json.read_text(encoding="utf-8"))
            rows = [_ensure_screening_row_defaults(row) for row in _read_csv_rows(comparison_csv)]
        except (OSError, ValueError, pd.errors.ParserError):
            return None
        summary = {key: value for key, value in payload.items() if key not in {"message", "errors"}}
        summary.setdefault("uniprot_candidate_source", "上次保存的方法比较结果")
        summary.setdefault("uniprot_reused_from_disk", True)
        for key, value in _rules_score_distribution(rows).items():
            summary.setdefault(key, value)
        return SignalPeptideScreeningResult(
            available=bool(payload.get("uspnet_available", False)),
            success=bool(payload.get("success", True)),
            message=str(payload.get("message", "已加载上次保存的筛选结果。")),
            summary=summary,
            rows=rows,
            output_dir=paths["output_dir"],
            uniprot_csv=paths["uniprot_csv"],
            duplicate_csv=paths["duplicate_csv"],
            input_fasta=paths["input_fasta"],
            comparison_csv=comparison_csv,
            recommended_fasta=paths["recommended_fasta"],
            uspnet_raw_dir=Path(str(payload["uspnet_raw_dir"])) if payload.get("uspnet_raw_dir") else None,
            summary_json=summary_json,
            errors=list(payload.get("errors", [])),
        )

    def screen_uniprot_candidates(
        self,
        *,
        taxon_id: int = 4922,
        max_records: int = 300,
        reviewed_only: bool = False,
        timeout_seconds: int = 3600,
    ) -> SignalPeptideScreeningResult:
        paths = self._output_paths()
        output_dir = paths["output_dir"]
        output_dir.mkdir(parents=True, exist_ok=True)
        uniprot_csv = paths["uniprot_csv"]
        duplicate_csv = paths["duplicate_csv"]
        input_fasta = paths["input_fasta"]
        comparison_csv = paths["comparison_csv"]
        recommended_fasta = paths["recommended_fasta"]
        summary_json = paths["summary_json"]

        persisted_discovery = self._load_persisted_uniprot_discovery(
            taxon_id=taxon_id,
            max_records=max_records,
            reviewed_only=reviewed_only,
            exclude_existing=True,
        )
        reused_uniprot = persisted_discovery is not None
        discovery = persisted_discovery or self.discover_and_persist_uniprot_candidates(
            taxon_id=taxon_id,
            max_records=max_records,
            reviewed_only=reviewed_only,
            exclude_existing=True,
        )
        candidate_rows = discovery.rows
        write_candidate_fasta(input_fasta, candidate_rows)

        errors = list(discovery.errors)
        summary: dict[str, object] = {
            "taxon_id": taxon_id,
            "reviewed_only": reviewed_only,
            "max_records": max_records,
            "uniprot_initial_hits": discovery.initial_hit_count,
            "uniprot_fetched_records": discovery.fetched_record_count,
            "uniprot_extracted_signal_count": discovery.extracted_signal_count,
            "uniprot_duplicate_count": discovery.duplicate_count,
            "deduplicated_candidates": discovery.deduplicated_count,
            "uniprot_candidate_source": "已复用本地保存的 UniProt 候选" if reused_uniprot else "UniProt API 实时查询",
            "uniprot_reused_from_disk": reused_uniprot,
            "uniprot_source_url": discovery.source_url,
            "rules_passed": 0,
            "rules_high_priority": 0,
            "rules_score_95_plus": 0,
            "rules_score_80_to_94": 0,
            "rules_score_65_to_79": 0,
            "rules_score_below_65": 0,
            "uspnet_available": False,
            "uspnet_success": False,
            "uspnet_completed": 0,
            "uspnet_passed": 0,
            "consensus_passed": 0,
            "needs_external_review": 0,
        }

        if not candidate_rows:
            message = "UniProt 没有返回可用于比较的候选信号肽。"
            write_json(summary_json, {**summary, "message": message, "errors": errors})
            return SignalPeptideScreeningResult(
                available=False,
                success=False,
                message=message,
                summary=summary,
                rows=[],
                output_dir=output_dir,
                uniprot_csv=uniprot_csv,
                duplicate_csv=duplicate_csv,
                input_fasta=input_fasta,
                summary_json=summary_json,
                errors=errors,
            )

        screened_rows = [_add_rule_screening(row) for row in candidate_rows]
        summary["rules_passed"] = sum(1 for row in screened_rows if row["rules_pass"])
        summary["rules_high_priority"] = sum(1 for row in screened_rows if row["rules_high_priority"])
        summary.update(_rules_score_distribution(screened_rows))

        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        uspnet_raw_dir = output_dir / "uspnet_raw" / run_id
        uspnet_result = self.uspnet_adapter.run(input_fasta, uspnet_raw_dir, timeout_seconds=timeout_seconds)
        summary["uspnet_available"] = uspnet_result.available
        summary["uspnet_success"] = uspnet_result.success
        if not uspnet_result.available:
            errors.append(uspnet_result.message)
        else:
            prediction_by_id = {prediction.candidate_id: prediction for prediction in uspnet_result.predictions}
            screened_rows = [_merge_uspnet_screening(row, prediction_by_id) for row in screened_rows]
            summary["uspnet_completed"] = sum(1 for row in screened_rows if row["uspnet_completed"])
            summary["uspnet_passed"] = sum(1 for row in screened_rows if row["uspnet_pass"])
            if not uspnet_result.success:
                errors.append(uspnet_result.message)

        uspnet_results_usable = bool(summary["uspnet_success"])
        screened_rows = [_finalize_recommendation(row, uspnet_results_usable) for row in screened_rows]
        summary["consensus_passed"] = sum(1 for row in screened_rows if row["consensus_pass"])
        summary["needs_external_review"] = sum(1 for row in screened_rows if row["screening_status"] == "规则高优先级，待 USPNet 复核")

        recommended_rows = [
            row
            for row in screened_rows
            if row["consensus_pass"] or row["screening_status"] == "规则高优先级，待 USPNet 复核"
        ]
        write_csv(comparison_csv, screened_rows)
        write_signal_peptide_fasta(recommended_fasta, recommended_rows)

        if summary["uspnet_success"]:
            message = (
                f"多方法比较完成：UniProt 初始命中 {summary['uniprot_initial_hits']} 条，"
                f"去重后 {summary['deduplicated_candidates']} 条，重复 {summary['uniprot_duplicate_count']} 条，"
                f"规则高优先级 {summary['rules_high_priority']} 条，"
                f"USPNet 通过 {summary['uspnet_passed']} 条，一致通过 {summary['consensus_passed']} 条。"
                f"候选来源：{summary['uniprot_candidate_source']}。"
            )
        elif summary["uspnet_available"]:
            message = (
                f"UniProt + 自研规则预筛完成：初始命中 {summary['uniprot_initial_hits']} 条，"
                f"去重后 {summary['deduplicated_candidates']} 条，重复 {summary['uniprot_duplicate_count']} 条，"
                f"规则高优先级 {summary['rules_high_priority']} 条。"
                "USPNet 已检测到，但本次运行未完成，因此没有多方法一致通过结论。"
                f"候选来源：{summary['uniprot_candidate_source']}。"
            )
        else:
            message = (
                f"UniProt + 自研规则预筛完成：初始命中 {summary['uniprot_initial_hits']} 条，"
                f"去重后 {summary['deduplicated_candidates']} 条，重复 {summary['uniprot_duplicate_count']} 条，"
                f"规则高优先级 {summary['rules_high_priority']} 条。"
                "USPNet 尚未安装，因此没有多方法一致通过结论。"
                f"候选来源：{summary['uniprot_candidate_source']}。"
            )

        write_json(
            summary_json,
            {
                **summary,
                "success": True,
                "message": message,
                "errors": errors,
                "uspnet_raw_dir": str(uspnet_raw_dir),
            },
        )
        return SignalPeptideScreeningResult(
            available=bool(summary["uspnet_available"]),
            success=True,
            message=message,
            summary=summary,
            rows=screened_rows,
            output_dir=output_dir,
            uniprot_csv=uniprot_csv,
            duplicate_csv=duplicate_csv,
            input_fasta=input_fasta,
            comparison_csv=comparison_csv,
            recommended_fasta=recommended_fasta,
            uspnet_raw_dir=uspnet_raw_dir,
            summary_json=summary_json,
            errors=errors,
        )

    def _output_paths(self) -> dict[str, Path]:
        output_dir = self.paths.opn_signal_peptides_dir
        return {
            "output_dir": output_dir,
            "uniprot_csv": output_dir / UNIPROT_CANDIDATES_CSV,
            "duplicate_csv": output_dir / UNIPROT_DUPLICATES_CSV,
            "discovery_summary_json": output_dir / UNIPROT_DISCOVERY_SUMMARY_JSON,
            "input_fasta": output_dir / METHOD_INPUT_FASTA,
            "comparison_csv": output_dir / METHOD_COMPARISON_CSV,
            "recommended_fasta": output_dir / RECOMMENDED_FASTA,
            "summary_json": output_dir / METHOD_SUMMARY_JSON,
        }

    def _persist_uniprot_discovery(
        self,
        discovery: UniProtCandidateLibraryResult,
        *,
        taxon_id: int,
        max_records: int,
        reviewed_only: bool,
        exclude_existing: bool,
    ) -> None:
        paths = self._output_paths()
        write_csv(paths["uniprot_csv"], discovery.rows)
        write_csv(paths["duplicate_csv"], discovery.duplicate_rows)
        write_candidate_fasta(paths["input_fasta"], discovery.rows)
        write_json(
            paths["discovery_summary_json"],
            {
                "taxon_id": taxon_id,
                "max_records": max_records,
                "reviewed_only": reviewed_only,
                "exclude_existing": exclude_existing,
                "source_url": discovery.source_url,
                "initial_hit_count": discovery.initial_hit_count,
                "fetched_record_count": discovery.fetched_record_count,
                "extracted_signal_count": discovery.extracted_signal_count,
                "deduplicated_count": discovery.deduplicated_count,
                "duplicate_count": discovery.duplicate_count,
                "errors": discovery.errors,
                "saved_at": datetime.now().isoformat(timespec="seconds"),
            },
        )

    def _load_persisted_uniprot_discovery(
        self,
        *,
        taxon_id: int,
        max_records: int,
        reviewed_only: bool,
        exclude_existing: bool,
    ) -> UniProtCandidateLibraryResult | None:
        paths = self._output_paths()
        if not paths["uniprot_csv"].exists():
            return None
        rows = _read_csv_rows(paths["uniprot_csv"])
        if not rows:
            return None

        duplicate_rows = _read_csv_rows(paths["duplicate_csv"])
        summary = _read_json_dict(paths["discovery_summary_json"])
        if not summary and not (
            taxon_id == 4922 and max_records <= 300 and reviewed_only is False and exclude_existing is True
        ):
            return None
        if summary and not _discovery_summary_matches(
            summary,
            taxon_id=taxon_id,
            max_records=max_records,
            reviewed_only=reviewed_only,
            exclude_existing=exclude_existing,
        ):
            return None

        inferred_initial_count = len(rows) + len(duplicate_rows)
        return UniProtCandidateLibraryResult(
            rows=rows,
            source_url=str(summary.get("source_url", "local persisted uniprot_candidates.csv")),
            errors=list(summary.get("errors", [])),
            initial_hit_count=int(summary.get("initial_hit_count", inferred_initial_count)),
            fetched_record_count=int(summary.get("fetched_record_count", inferred_initial_count)),
            extracted_signal_count=int(summary.get("extracted_signal_count", inferred_initial_count)),
            deduplicated_count=int(summary.get("deduplicated_count", len(rows))),
            duplicate_count=int(summary.get("duplicate_count", len(duplicate_rows))),
            duplicate_rows=duplicate_rows,
        )


def _add_rule_screening(row: dict[str, object]) -> dict[str, object]:
    sequence = str(row.get("signal_peptide_sequence", ""))
    result = score_signal_peptide(sequence)
    high_priority = result.passed and result.score >= 90 and not result.risks
    return {
        **row,
        "uniprot_signal_annotated": bool(sequence),
        "rules_score": result.score,
        "rules_pass": result.passed,
        "rules_high_priority": high_priority,
        "rules_priority": "高" if high_priority else ("中" if result.passed else "低"),
        "rules_score_note": _rules_score_note(result.score, result.passed),
        "rules_tier": result.tier,
        "rules_reasons": "；".join(result.reasons),
        "rules_risks": "；".join(result.risks),
        "rules_n_region_positive_count": result.n_region_positive_count,
        "rules_h_region_max_hydrophobicity": result.h_region_max_hydrophobicity,
        "rules_c_region_small_neutral": result.c_region_small_neutral_rule,
        "uspnet_completed": False,
        "uspnet_prediction": "",
        "uspnet_prediction_label": "未运行",
        "uspnet_interpretation": "尚未得到 USPNet 预测结果。",
        "uspnet_cleavage_sequence": "",
        "uspnet_pass": False,
    }


def _merge_uspnet_screening(row: dict[str, object], prediction_by_id: dict[str, object]) -> dict[str, object]:
    candidate_id = str(row["candidate_id"])
    prediction = prediction_by_id.get(candidate_id)
    if prediction is None:
        return row
    return {
        **row,
        "uspnet_completed": True,
        "uspnet_prediction": prediction.predicted_type,
        "uspnet_prediction_label": _uspnet_prediction_label(prediction.predicted_type),
        "uspnet_interpretation": _uspnet_interpretation(prediction.predicted_type, prediction.passed),
        "uspnet_cleavage_sequence": prediction.predicted_cleavage,
        "uspnet_pass": prediction.passed,
    }


def _finalize_recommendation(row: dict[str, object], uspnet_results_usable: bool) -> dict[str, object]:
    rules_pass = bool(row.get("rules_pass"))
    rules_high_priority = bool(row.get("rules_high_priority"))
    uspnet_completed = bool(row.get("uspnet_completed"))
    uspnet_pass = bool(row.get("uspnet_pass"))
    has_uspnet_judgement = uspnet_results_usable and uspnet_completed
    consensus_pass = rules_high_priority and has_uspnet_judgement and uspnet_pass
    if consensus_pass:
        status = "多方法一致通过"
    elif rules_high_priority and not has_uspnet_judgement:
        status = "规则高优先级，待 USPNet 复核"
    elif rules_pass and not has_uspnet_judgement:
        status = "规则基础通过，待人工复核"
    elif rules_high_priority and has_uspnet_judgement and not uspnet_pass:
        status = "规则与 USPNet 冲突"
    else:
        status = "暂不推荐"
    return {
        **row,
        "consensus_pass": consensus_pass,
        "screening_status": status,
        "recommended_for_draft_library": consensus_pass or status == "规则高优先级，待 USPNet 复核",
    }


def _read_json_dict(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _read_csv_rows(path: Path) -> list[dict[str, object]]:
    if not path.exists() or path.stat().st_size <= 1:
        return []
    try:
        frame = pd.read_csv(path, keep_default_na=False)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, OSError):
        return []
    return [_coerce_row_types(row) for row in frame.to_dict("records")]


def _coerce_row_types(row: dict[str, object]) -> dict[str, object]:
    bool_columns = {
        "already_in_formal_library",
        "uniprot_signal_annotated",
        "rules_pass",
        "rules_high_priority",
        "rules_c_region_small_neutral",
        "uspnet_completed",
        "uspnet_pass",
        "consensus_pass",
        "recommended_for_draft_library",
    }
    coerced = dict(row)
    for column in bool_columns:
        if column in coerced:
            coerced[column] = _coerce_bool(coerced[column])
    return coerced


def _coerce_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _ensure_screening_row_defaults(row: dict[str, object]) -> dict[str, object]:
    updated = dict(row)
    score = int(float(str(updated.get("rules_score", 0) or 0)))
    rules_pass = bool(updated.get("rules_pass", False))
    predicted_type = str(updated.get("uspnet_prediction", "") or "")
    uspnet_pass = bool(updated.get("uspnet_pass", False))
    updated.setdefault("rules_score_note", _rules_score_note(score, rules_pass))
    updated.setdefault(
        "uspnet_prediction_label",
        _uspnet_prediction_label(predicted_type) if predicted_type else "未运行",
    )
    updated.setdefault(
        "uspnet_interpretation",
        _uspnet_interpretation(predicted_type, uspnet_pass) if predicted_type else "尚未得到 USPNet 预测结果。",
    )
    return updated


def _discovery_summary_matches(
    summary: dict[str, object],
    *,
    taxon_id: int,
    max_records: int,
    reviewed_only: bool,
    exclude_existing: bool,
) -> bool:
    return (
        int(summary.get("taxon_id", taxon_id)) == taxon_id
        and int(summary.get("max_records", max_records)) >= max_records
        and bool(summary.get("reviewed_only", reviewed_only)) == reviewed_only
        and bool(summary.get("exclude_existing", exclude_existing)) == exclude_existing
    )


def _rules_score_distribution(rows: list[dict[str, object]]) -> dict[str, int]:
    scores = [int(row.get("rules_score", 0)) for row in rows]
    return {
        "rules_score_95_plus": sum(1 for score in scores if score >= 95),
        "rules_score_80_to_94": sum(1 for score in scores if 80 <= score < 95),
        "rules_score_65_to_79": sum(1 for score in scores if 65 <= score < 80),
        "rules_score_below_65": sum(1 for score in scores if score < 65),
    }


def _rules_score_note(score: int, passed: bool) -> str:
    if score >= 95:
        return "典型信号肽特征完整；只说明像 signal peptide，不代表产量更高。"
    if passed:
        return "具备基础 signal peptide 特征，但仍需要模型和实验复核。"
    return "规则特征不足，暂不建议直接进入草案库。"


def _uspnet_prediction_label(predicted_type: str) -> str:
    if predicted_type == "SP":
        return "SP：USPNet 判断为信号肽"
    if predicted_type == "NO_SP":
        return "NO_SP：USPNet 未判断为信号肽"
    return f"{predicted_type}：USPNet 原始类别"


def _uspnet_interpretation(predicted_type: str, passed: bool) -> str:
    if passed:
        return "机器学习复核支持该序列具有信号肽特征。"
    if predicted_type == "NO_SP":
        return "机器学习复核不支持该序列作为信号肽，建议降级或人工复核。"
    return "USPNet 给出非 SP 类别，需结合规则和来源证据人工判断。"
