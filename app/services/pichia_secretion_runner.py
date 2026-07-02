from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any

from app.services.pichia_secretion_schema import SecretionRunRequest, SecretionRunResponse


def run_pichia_pipeline_draft(
    request: SecretionRunRequest,
    *,
    output_dir: Path,
    warnings: list[str],
    engine_target_id: str,
    target_input: Any | None,
    sequence_contract: dict[str, Any],
) -> SecretionRunResponse:
    """Run the formal python_pichia pipeline for the app service facade."""
    _ensure_pcsec_pichia_analysis_api()
    from pcsec_pichia.engines.base import PichiaSimulationRequest
    from pcsec_pichia.pipeline import run_pichia_secretion_simulation

    try:
        engine_request = PichiaSimulationRequest(
            target_id=engine_target_id,
            candidate_id=engine_target_id,
            target_input=target_input,
            compatibility_mode="corrected",
            glycosylation_mode="native",
            mu=request.mu,
            media_type=request.media_type,
            carbon_source_id=request.carbon_source_id,
            enable_ribosome_translation_constraint=request.enable_ribosome_translation_constraint,
            enable_misfolding_constraint=request.enable_misfolding_constraint,
            growth_points=tuple(request.growth_points),
            ko_gene_ids=tuple(request.ko_gene_ids or request.ko_candidates),
            ko_reaction_ids=tuple(request.ko_reaction_ids),
            oe_gene_ids=tuple(request.oe_gene_ids),
            oe_reaction_ids=tuple(request.oe_reaction_ids or request.oe_candidates),
            screen_candidate_limit=int(request.screen_candidate_limit),
            enable_gene_rule_overlay=bool(request.enable_gene_rule_overlay),
            enable_cost_slope_compatibility=bool(request.enable_cost_slope_compatibility),
            cost_slope_growth_rates=tuple(request.cost_slope_growth_rates),
            cost_slope_secretion_ratios=tuple(request.cost_slope_secretion_ratios),
            cost_slope_capacity_fractions=tuple(request.cost_slope_capacity_fractions),
            cost_slope_medium_compatibility_mode=str(request.cost_slope_medium_compatibility_mode),
            **sequence_contract,
        )
        result = run_pichia_secretion_simulation(engine_request, output_dir=output_dir)
    except Exception as exc:  # keep UI/API facade stable; engine details remain in warnings/errors.
        return SecretionRunResponse(
            success=False,
            target_id=request.target_id,
            result_status="error",
            matlab_alignment_status="pending",
            warnings=warnings,
            errors=[f"{type(exc).__name__}: {exc}"],
            output_dir=output_dir,
        )

    summary_payload = _result_summary_payload(result.summary_path)
    summary_warnings = [str(item) for item in summary_payload.get("screen_warnings") or []]
    medium_warnings = [str(item) for item in summary_payload.get("medium_condition_warnings") or []]
    target_warnings = [str(item) for item in summary_payload.get("target_warnings") or []]
    target_metadata = summary_payload.get("target_metadata") if isinstance(summary_payload.get("target_metadata"), dict) else {}
    protein_cost = (
        summary_payload.get("protein_cost_analysis")
        if isinstance(summary_payload.get("protein_cost_analysis"), dict)
        else {}
    )
    target_growth = (
        summary_payload.get("target_growth_analysis")
        if isinstance(summary_payload.get("target_growth_analysis"), dict)
        else {}
    )
    yield_recommendations = (
        summary_payload.get("yield_improvement_recommendations")
        if isinstance(summary_payload.get("yield_improvement_recommendations"), dict)
        else {}
    )
    medium_condition = (
        summary_payload.get("medium_condition")
        if isinstance(summary_payload.get("medium_condition"), dict)
        else {}
    )

    return SecretionRunResponse(
        success=result.success,
        target_id=result.target_id,
        result_status=result.result_status,
        matlab_alignment_status=result.matlab_alignment_status,
        objective_value=result.objective_value,
        secretion_flux=result.objective_value,
        constraint_counts=dict(result.constraint_counts),
        warnings=[*warnings, *medium_warnings, *summary_warnings],
        output_dir=output_dir,
        summary_path=result.summary_path,
        report_path=result.report_path,
        candidate_table_path=result.candidate_table_path,
        tradeoff_path=result.tradeoff_path,
        alignment_summary=dict(result.alignment_summary),
        target_metadata=dict(target_metadata),
        target_warnings=target_warnings,
        protein_cost_analysis=dict(protein_cost),
        target_growth_analysis=dict(target_growth),
        yield_improvement_recommendations=dict(yield_recommendations),
        medium_condition=dict(medium_condition),
    )


def _result_summary_payload(summary_path: Path | None) -> dict[str, Any]:
    if summary_path is None or not summary_path.exists():
        return {}
    try:
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _ensure_pcsec_pichia_analysis_api() -> None:
    """Refresh long-lived Streamlit imports after engine modules change on disk."""
    runtime = importlib.import_module("sys")
    module = runtime.modules.get("pcsec_pichia.analysis")
    required = {
        "analyze_target_growth_impact",
        "analyze_yield_improvement_candidates",
        "summarize_protein_cost_slope_compatibility",
        "summarize_yield_improvement_recommendations",
    }
    if module is not None and not required.issubset(set(dir(module))):
        module = importlib.reload(module)
    if module is None:
        module = importlib.import_module("pcsec_pichia.analysis")
    missing = sorted(name for name in required if not hasattr(module, name))
    if missing:
        source = getattr(module, "__file__", "<unknown>")
        raise ImportError(
            f"pcsec_pichia.analysis is missing {', '.join(missing)} "
            f"after reload from {source}"
        )


__all__ = ["run_pichia_pipeline_draft"]
