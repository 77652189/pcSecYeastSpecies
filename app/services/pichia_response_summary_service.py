from __future__ import annotations

from typing import Any

from app.services.pichia_secretion_schema import SecretionRunResponse


def response_to_summary(response: SecretionRunResponse) -> dict[str, Any]:
    return {
        "success": response.success,
        "target_id": response.target_id,
        "result_status": response.result_status,
        "matlab_alignment_status": response.matlab_alignment_status,
        "objective_value": response.objective_value,
        "secretion_flux": response.secretion_flux,
        "constraint_counts": response.constraint_counts,
        "warnings": response.warnings,
        "errors": response.errors,
        "output_dir": str(response.output_dir) if response.output_dir else None,
        "summary_path": str(response.summary_path) if response.summary_path else None,
        "report_path": str(response.report_path) if response.report_path else None,
        "candidate_table_path": str(response.candidate_table_path) if response.candidate_table_path else None,
        "tradeoff_path": str(response.tradeoff_path) if response.tradeoff_path else None,
        "alignment_summary": response.alignment_summary,
        "target_metadata": response.target_metadata,
        "target_warnings": response.target_warnings,
    }


__all__ = ["response_to_summary"]
