from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from app.adapters.uspnet import USPNetAdapter
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

    def screen_uniprot_candidates(
        self,
        *,
        taxon_id: int = 4922,
        max_records: int = 300,
        reviewed_only: bool = False,
        timeout_seconds: int = 3600,
    ) -> SignalPeptideScreeningResult:
        output_dir = self.paths.opn_signal_peptides_dir
        output_dir.mkdir(parents=True, exist_ok=True)
        uniprot_csv = output_dir / UNIPROT_CANDIDATES_CSV
        input_fasta = output_dir / METHOD_INPUT_FASTA
        comparison_csv = output_dir / METHOD_COMPARISON_CSV
        recommended_fasta = output_dir / RECOMMENDED_FASTA
        summary_json = output_dir / METHOD_SUMMARY_JSON

        discovery = self.library_service.discover_uniprot_candidate_library(
            taxon_id=taxon_id,
            max_records=max_records,
            reviewed_only=reviewed_only,
            exclude_existing=False,
        )
        candidate_rows = discovery.rows
        write_csv(uniprot_csv, candidate_rows)
        write_candidate_fasta(input_fasta, candidate_rows)

        errors = list(discovery.errors)
        summary: dict[str, object] = {
            "taxon_id": taxon_id,
            "reviewed_only": reviewed_only,
            "max_records": max_records,
            "uniprot_initial_hits": discovery.initial_hit_count,
            "uniprot_fetched_records": discovery.fetched_record_count,
            "uniprot_extracted_signal_count": discovery.extracted_signal_count,
            "deduplicated_candidates": discovery.deduplicated_count,
            "rules_passed": 0,
            "rules_high_priority": 0,
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
                input_fasta=input_fasta,
                summary_json=summary_json,
                errors=errors,
            )

        screened_rows = [_add_rule_screening(row) for row in candidate_rows]
        summary["rules_passed"] = sum(1 for row in screened_rows if row["rules_pass"])
        summary["rules_high_priority"] = sum(1 for row in screened_rows if row["rules_high_priority"])

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
                f"去重后 {summary['deduplicated_candidates']} 条，规则高优先级 {summary['rules_high_priority']} 条，"
                f"USPNet 通过 {summary['uspnet_passed']} 条，一致通过 {summary['consensus_passed']} 条。"
            )
        elif summary["uspnet_available"]:
            message = (
                f"UniProt + 自研规则预筛完成：初始命中 {summary['uniprot_initial_hits']} 条，"
                f"去重后 {summary['deduplicated_candidates']} 条，规则高优先级 {summary['rules_high_priority']} 条。"
                "USPNet 已检测到，但本次运行未完成，因此没有多方法一致通过结论。"
            )
        else:
            message = (
                f"UniProt + 自研规则预筛完成：初始命中 {summary['uniprot_initial_hits']} 条，"
                f"去重后 {summary['deduplicated_candidates']} 条，规则高优先级 {summary['rules_high_priority']} 条。"
                "USPNet 尚未安装，因此没有多方法一致通过结论。"
            )

        write_json(summary_json, {**summary, "message": message, "errors": errors})
        return SignalPeptideScreeningResult(
            available=bool(summary["uspnet_available"]),
            success=True,
            message=message,
            summary=summary,
            rows=screened_rows,
            output_dir=output_dir,
            uniprot_csv=uniprot_csv,
            input_fasta=input_fasta,
            comparison_csv=comparison_csv,
            recommended_fasta=recommended_fasta,
            uspnet_raw_dir=uspnet_raw_dir,
            summary_json=summary_json,
            errors=errors,
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
        "rules_tier": result.tier,
        "rules_reasons": "；".join(result.reasons),
        "rules_risks": "；".join(result.risks),
        "rules_n_region_positive_count": result.n_region_positive_count,
        "rules_h_region_max_hydrophobicity": result.h_region_max_hydrophobicity,
        "rules_c_region_small_neutral": result.c_region_small_neutral_rule,
        "uspnet_completed": False,
        "uspnet_prediction": "",
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
