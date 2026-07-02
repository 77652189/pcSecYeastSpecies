from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from pcsec_pichia.secretion_plan import SecretionPlanResult, build_secretion_plan
from pcsec_pichia.simulation import (
    GrowthTradeoffResult,
    ProteinCostSlopeCompatibilityResult,
    SecretionSimulationResult,
)
from pcsec_pichia.targets import TargetSpec, target_features
from pcsec_pichia.services.gene_evidence import recommendation_tier_for_candidate


@dataclass(frozen=True)
class ProteinCostItem:
    category: str
    label: str
    basis: str
    raw_value: float
    relative_score: float
    interpretation: str


@dataclass(frozen=True)
class ProteinCostAnalysisResult:
    target_id: str
    protein_id: str
    route_kind: str
    sequence_features: dict[str, Any]
    ptm_counts: dict[str, int]
    stage_counts: dict[str, int]
    cost_items: tuple[ProteinCostItem, ...]
    total_relative_score: float
    dominant_cost_categories: tuple[str, ...]
    warnings: tuple[str, ...]
    result_status: str = "draft_explanatory"


@dataclass(frozen=True)
class GrowthTradeoffPoint:
    mu: float
    success: bool
    secretion_flux: float | None
    secretion_per_biomass: float | None
    status: str
    interpretation: str


@dataclass(frozen=True)
class TargetGrowthAnalysisResult:
    target_id: str
    growth_points: tuple[float, ...]
    valid_point_count: int
    best_secretion_point: dict[str, object] | None
    best_secretion_per_biomass_point: dict[str, object] | None
    growth_sensitivity_label: str
    growth_sensitivity_reason: str
    tradeoff_points: tuple[GrowthTradeoffPoint, ...]
    warnings: tuple[str, ...]
    result_status: str = "draft_explanatory"


@dataclass(frozen=True)
class ProteinLpAttributionResult:
    target_id: str
    result_status: str
    objective_evidence: dict[str, object]
    dominant_constraint_blocks: tuple[dict[str, object], ...]
    top_constraint_marginals: tuple[dict[str, object], ...]
    top_bound_marginals: tuple[dict[str, object], ...]
    target_related_fluxes: tuple[dict[str, object], ...]
    active_bound_counts: dict[str, int]
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class YieldImprovementCandidateRecommendation:
    candidate_id: str
    display_name: str
    model_gene_id: str
    locus_tag: str
    intervention_type: str
    execution_mode: str
    delta_objective: float | None
    effect_label: str
    growth_risk_label: str
    secretory_process: str
    evidence_tier: str
    wet_lab_readiness: str
    database_annotation_sources: tuple[str, ...]
    database_annotation_confidence: str
    model_gpr_executable: bool
    oe_reaction_proxy: bool
    phenotype_evidence: dict[str, object]
    recommendation_tier: str
    recommendation_tier_reason: str
    recommendation_label: str
    recommendation_score: float
    rationale: str
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class YieldImprovementRecommendationResult:
    target_id: str
    medium_condition: dict[str, object]
    baseline_objective: float | None
    baseline_secretion_flux: float | None
    candidate_count: int
    recommended_candidates: tuple[YieldImprovementCandidateRecommendation, ...]
    not_recommended_candidates: tuple[YieldImprovementCandidateRecommendation, ...]
    unresolved_candidates: tuple[YieldImprovementCandidateRecommendation, ...]
    warnings: tuple[str, ...]
    result_status: str = "draft_model_recommendation"


def analyze_target_protein_cost(
    target: TargetSpec,
    plan: SecretionPlanResult | None = None,
) -> ProteinCostAnalysisResult:
    resolved_plan = plan or build_secretion_plan(target)
    features = _target_feature_payload(target)
    raw_items = _cost_items_from_target_and_plan(target, resolved_plan, features)
    cost_items = _normalise_cost_scores(raw_items)
    return ProteinCostAnalysisResult(
        target_id=target.target_id,
        protein_id=target.protein_id,
        route_kind=resolved_plan.route_kind,
        sequence_features=features,
        ptm_counts=dict(resolved_plan.ptm_counts),
        stage_counts=dict(resolved_plan.stage_counts),
        cost_items=cost_items,
        total_relative_score=round(sum(item.relative_score for item in cost_items), 3),
        dominant_cost_categories=_dominant_categories(cost_items),
        warnings=tuple(_cost_warnings(target)),
    )


def analyze_target_protein_lp_attribution(
    target: TargetSpec,
    plan: SecretionPlanResult,
    constraint_counts: dict[str, int],
    simulation_result: SecretionSimulationResult,
    reaction_ids: tuple[str, ...] = (),
    top_n: int = 10,
) -> ProteinLpAttributionResult:
    sensitivity = simulation_result.lp_sensitivity or {}
    if not simulation_result.success or not sensitivity:
        return ProteinLpAttributionResult(
            target_id=target.target_id,
            result_status="draft_lp_sensitivity_unavailable",
            objective_evidence=_objective_evidence(simulation_result),
            dominant_constraint_blocks=(),
            top_constraint_marginals=(),
            top_bound_marginals=(),
            target_related_fluxes=_target_related_fluxes(target, plan, simulation_result, reaction_ids),
            active_bound_counts={},
            warnings=(
                "LP sensitivity is unavailable because the Python solve failed or did not return HiGHS marginals.",
                "This does not change the corrected pipeline objective or constraints.",
            ),
        )

    eq_marginals = tuple(float(value) for value in sensitivity.get("eq_marginals", ()))
    ub_marginals = tuple(float(value) for value in sensitivity.get("ub_marginals", ()))
    lower_marginals = tuple(float(value) for value in sensitivity.get("lower_marginals", ()))
    upper_marginals = tuple(float(value) for value in sensitivity.get("upper_marginals", ()))
    blocks = _constraint_blocks(constraint_counts, len(eq_marginals), len(ub_marginals))
    return ProteinLpAttributionResult(
        target_id=target.target_id,
        result_status="draft_lp_sensitivity",
        objective_evidence=_objective_evidence(simulation_result),
        dominant_constraint_blocks=_dominant_constraint_blocks(blocks, eq_marginals, ub_marginals, top_n=top_n),
        top_constraint_marginals=_top_constraint_marginals(blocks, eq_marginals, ub_marginals, top_n=top_n),
        top_bound_marginals=_top_bound_marginals(
            lower_marginals,
            upper_marginals,
            reaction_ids=reaction_ids,
            target=target,
            plan=plan,
            top_n=top_n,
        ),
        target_related_fluxes=_target_related_fluxes(target, plan, simulation_result, reaction_ids),
        active_bound_counts={
            "lower_marginal_nonzero": _nonzero_count(lower_marginals),
            "upper_marginal_nonzero": _nonzero_count(upper_marginals),
            "total_bound_marginal_nonzero": _nonzero_count(lower_marginals) + _nonzero_count(upper_marginals),
        },
        warnings=(
            "LP sensitivity is a Python draft based on SciPy HiGHS marginals; it is not MATLAB/SoPlex fully aligned shadow pricing.",
            "The maximization problem is solved through SciPy minimization, so signs should be interpreted as draft sensitivity evidence.",
            "Only compressed top-N attribution rows are written to reports and summaries.",
        ),
    )


def analyze_yield_improvement_candidates(
    candidate_rows: tuple[dict[str, Any], ...] | list[dict[str, Any]],
    *,
    target_id: str = "",
    medium_condition: dict[str, object] | None = None,
    baseline_objective: float | None = None,
    baseline_secretion_flux: float | None = None,
    curated_gene_evidence: tuple[dict[str, object], ...] | list[dict[str, object]] = (),
    top_n: int = 10,
) -> YieldImprovementRecommendationResult:
    rows = tuple(dict(row) for row in candidate_rows)
    evidence_index = _curated_evidence_index(curated_gene_evidence)
    resolved_target_id = target_id or (str(rows[0].get("target_id") or "") if rows else "")
    recommendations = tuple(
        _yield_candidate_recommendation(row, evidence_index, resolved_target_id)
        for row in rows
    )
    recommended = tuple(
        item
        for item in sorted(recommendations, key=_yield_recommendation_sort_key)
        if item.recommendation_label in {"strong_model_candidate", "promising_but_proxy_only"}
    )[: max(0, int(top_n))]
    unresolved = tuple(
        item
        for item in recommendations
        if item.recommendation_label == "unresolved_not_actionable"
    )
    not_recommended = tuple(
        item
        for item in sorted(recommendations, key=_yield_recommendation_sort_key)
        if item.recommendation_label
        not in {"strong_model_candidate", "promising_but_proxy_only", "unresolved_not_actionable"}
    )
    warnings = [
        "当前结果是 Python corrected draft 模型内推荐，不代表真实发酵产量或实验成功率。",
        "KO 表示 gene-level GPR knockout；OE 仍是 reaction-level proxy，不能直接等同于基因过表达湿实验。",
        "第一版只解释当前小批量候选，不做全模型 KO/OE 筛选。",
    ]
    if not rows:
        warnings.append("没有可用候选行；推荐结果为空。")
    return YieldImprovementRecommendationResult(
        target_id=resolved_target_id,
        medium_condition=dict(medium_condition or {}),
        baseline_objective=baseline_objective,
        baseline_secretion_flux=baseline_secretion_flux,
        candidate_count=len(rows),
        recommended_candidates=recommended,
        not_recommended_candidates=not_recommended,
        unresolved_candidates=unresolved,
        warnings=tuple(warnings),
    )


def analyze_target_growth_impact(
    tradeoff: GrowthTradeoffResult,
    baseline_growth_rate: float | None = None,
) -> TargetGrowthAnalysisResult:
    points = tuple(_growth_tradeoff_points(tradeoff))
    valid_points = tuple(point for point in points if point.success and point.secretion_flux is not None)
    best_flux = _best_point(valid_points, key="secretion_flux")
    best_per_biomass = _best_point(valid_points, key="secretion_per_biomass")
    interpreted = tuple(
        _with_growth_point_interpretation(point, best_flux, best_per_biomass, baseline_growth_rate)
        for point in points
    )
    interpreted_valid_points = tuple(point for point in interpreted if point.success and point.secretion_flux is not None)
    interpreted_best_flux = _best_point(interpreted_valid_points, key="secretion_flux")
    interpreted_best_per_biomass = _best_point(interpreted_valid_points, key="secretion_per_biomass")
    growth_sensitivity_label, growth_sensitivity_reason = _growth_sensitivity(points)
    return TargetGrowthAnalysisResult(
        target_id=tradeoff.target_id,
        growth_points=tuple(float(value) for value in tradeoff.growth_points),
        valid_point_count=len(valid_points),
        best_secretion_point=asdict(interpreted_best_flux) if interpreted_best_flux else None,
        best_secretion_per_biomass_point=asdict(interpreted_best_per_biomass) if interpreted_best_per_biomass else None,
        growth_sensitivity_label=growth_sensitivity_label,
        growth_sensitivity_reason=growth_sensitivity_reason,
        tradeoff_points=interpreted,
        warnings=(
            "当前生长分析是 Python draft explanatory tradeoff，不代表真实发酵生长预测。",
            "该分析只解释已有 small-grid tradeoff rows，不新增求解算法或长网格任务。",
        ),
    )


def summarize_protein_cost_analysis(result: ProteinCostAnalysisResult) -> dict[str, object]:
    return {
        "target_id": result.target_id,
        "protein_id": result.protein_id,
        "route_kind": result.route_kind,
        "sequence_features": result.sequence_features,
        "ptm_counts": result.ptm_counts,
        "stage_counts": result.stage_counts,
        "cost_items": build_cost_item_table(result),
        "total_relative_score": result.total_relative_score,
        "dominant_cost_categories": list(result.dominant_cost_categories),
        "warnings": list(result.warnings),
        "result_status": result.result_status,
    }


def summarize_protein_lp_attribution(result: ProteinLpAttributionResult) -> dict[str, object]:
    return {
        "target_id": result.target_id,
        "result_status": result.result_status,
        "objective_evidence": result.objective_evidence,
        "dominant_constraint_blocks": list(result.dominant_constraint_blocks),
        "top_constraint_marginals": list(result.top_constraint_marginals),
        "top_bound_marginals": list(result.top_bound_marginals),
        "target_related_fluxes": list(result.target_related_fluxes),
        "active_bound_counts": result.active_bound_counts,
        "warnings": list(result.warnings),
    }


def summarize_target_growth_analysis(result: TargetGrowthAnalysisResult) -> dict[str, object]:
    return {
        "target_id": result.target_id,
        "growth_points": list(result.growth_points),
        "valid_point_count": result.valid_point_count,
        "best_secretion_point": result.best_secretion_point,
        "best_secretion_per_biomass_point": result.best_secretion_per_biomass_point,
        "growth_sensitivity_label": result.growth_sensitivity_label,
        "growth_sensitivity_reason": result.growth_sensitivity_reason,
        "tradeoff_points": build_growth_tradeoff_item_table(result),
        "warnings": list(result.warnings),
        "result_status": result.result_status,
    }


def summarize_yield_improvement_recommendations(
    result: YieldImprovementRecommendationResult,
) -> dict[str, object]:
    return {
        "target_id": result.target_id,
        "medium_condition": result.medium_condition,
        "baseline_objective": result.baseline_objective,
        "baseline_secretion_flux": result.baseline_secretion_flux,
        "candidate_count": result.candidate_count,
        "recommended_candidates": build_yield_recommendation_table(result, bucket="recommended"),
        "not_recommended_candidates": build_yield_recommendation_table(result, bucket="not_recommended"),
        "unresolved_candidates": build_yield_recommendation_table(result, bucket="unresolved"),
        "summary_counts": {
            "recommended": len(result.recommended_candidates),
            "not_recommended": len(result.not_recommended_candidates),
            "unresolved": len(result.unresolved_candidates),
        },
        "warnings": list(result.warnings),
        "result_status": result.result_status,
    }


def summarize_protein_cost_slope_compatibility(
    result: ProteinCostSlopeCompatibilityResult | None,
) -> dict[str, object]:
    if result is None:
        return {
            "enabled": False,
            "result_status": "disabled",
            "warnings": [
                "MATLAB-compatible protein cost slope mode is disabled by default.",
            ],
        }
    return {
        "target_id": result.target_id,
        "enabled": result.enabled,
        "success": result.success,
        "growth_rates": list(result.growth_rates),
        "secretion_ratios": list(result.secretion_ratios),
        "secretion_ratio_policy": result.secretion_ratio_policy,
        "capacity_reference": result.capacity_reference,
        "capacity_fractions": list(result.capacity_fractions),
        "rows": list(result.rows),
        "glucose_cost_slopes": list(result.glucose_cost_slopes),
        "ribosome_cost_slopes": list(result.ribosome_cost_slopes),
        "medium_compatibility_mode": result.medium_compatibility_mode,
        "medium_bound_overrides": list(result.medium_bound_overrides),
        "result_status": result.result_status,
        "warnings": list(result.warnings),
        "comparison_scope": {
            "matlab_reference": "Code/pcSecPichia/Simulation/SimulateProteinCost.m and Results/Protein_cost_TP/all_proteincost_gluPP.mat",
            "definition": "fixed target exchange ratios, fixed growth rates, optimize Ex_glc_D, then estimate cost slopes",
            "current_default_definition": "fixed growth rate, corrected medium, maximize target secretion flux",
            "ratio_policy": result.secretion_ratio_policy,
            "not_default_pipeline": True,
            "medium_compatibility": result.medium_compatibility_mode,
            "medium_note": "corrected keeps Python corrected medium; matlab_legacy_cost applies 9 historical exchange lower bounds only for MATLAB artifact comparison.",
        },
    }


def build_cost_item_table(result: ProteinCostAnalysisResult) -> tuple[dict[str, object], ...]:
    return tuple(asdict(item) for item in result.cost_items)


def build_growth_tradeoff_item_table(result: TargetGrowthAnalysisResult) -> tuple[dict[str, object], ...]:
    return tuple(asdict(item) for item in result.tradeoff_points)


def build_yield_recommendation_table(
    result: YieldImprovementRecommendationResult,
    bucket: str = "recommended",
) -> tuple[dict[str, object], ...]:
    if bucket == "not_recommended":
        candidates = result.not_recommended_candidates
    elif bucket == "unresolved":
        candidates = result.unresolved_candidates
    else:
        candidates = result.recommended_candidates
    return tuple(asdict(item) for item in candidates)


def _curated_evidence_index(rows: tuple[dict[str, object], ...] | list[dict[str, object]]) -> dict[str, dict[str, object]]:
    index: dict[str, dict[str, object]] = {}
    for row in rows:
        record = dict(row)
        keys = {
            row.get("common_name"),
            row.get("declared_model_gene_id"),
            row.get("mapped_model_gene_id"),
            row.get("ko_reaction_id"),
            row.get("oe_reaction_id"),
        }
        for item in row.get("reaction_evidence") or ():
            if isinstance(item, dict):
                keys.add(item.get("reaction_id"))
        for key in keys:
            text = str(key or "").strip()
            if text:
                index[text.lower()] = record
    return index


def _yield_candidate_recommendation(
    row: dict[str, Any],
    evidence_index: dict[str, dict[str, object]],
    target_protein_context: str,
) -> YieldImprovementCandidateRecommendation:
    evidence = _evidence_for_candidate(row, evidence_index)
    intervention_type = str(row.get("intervention_type") or row.get("screen_type") or "")
    execution_mode = _yield_execution_mode(row, evidence)
    status = str(row.get("status") or "")
    unresolved = status.startswith("unresolved") or str(row.get("mapping_level") or "") == "unresolved"
    model_gpr_executable = execution_mode == "gene_level_ko" and not unresolved
    oe_reaction_proxy = (
        execution_mode == "reaction_level_oe_proxy"
        and not unresolved
        and status
        not in {"not_run_complex_subunit_limited", "not_run_gene_oe_proxy", "not_run_no_gpr_effect"}
    )
    phenotype = row.get("phenotype_evidence") if isinstance(row.get("phenotype_evidence"), dict) else {}
    recommendation_tier, tier_reason, matched_phenotype = _candidate_recommendation_tier(
        row,
        evidence,
        intervention_type,
        target_protein_context,
        model_gpr_executable=model_gpr_executable,
        oe_reaction_proxy=oe_reaction_proxy,
    )
    phenotype_payload = phenotype or (matched_phenotype.to_dict() if matched_phenotype else {})
    delta = _optional_float(row.get("delta_objective"))
    success = _truthy(row.get("success"))
    effect_label = str(row.get("effect_label") or "")
    label = _yield_recommendation_label(
        success=success,
        unresolved=unresolved,
        delta_objective=delta,
        execution_mode=execution_mode,
        evidence=evidence,
        recommendation_tier=recommendation_tier,
    )
    score = _yield_recommendation_score(label, delta, execution_mode, evidence)
    display_name = _yield_display_name(row, evidence)
    evidence_tier = _yield_evidence_tier(row, evidence, execution_mode)
    warnings = tuple(_yield_candidate_warnings(label, execution_mode, evidence))
    return YieldImprovementCandidateRecommendation(
        candidate_id=str(row.get("candidate_id") or row.get("input_gene_id") or row.get("reaction_id") or display_name),
        display_name=display_name,
        model_gene_id=str(row.get("canonical_gene_id") or row.get("gene_id") or evidence.get("mapped_model_gene_id") or evidence.get("declared_model_gene_id") or ""),
        locus_tag=str(row.get("canonical_gene_id") or row.get("gene_id") or evidence.get("mapped_model_gene_id") or evidence.get("declared_model_gene_id") or ""),
        intervention_type=intervention_type,
        execution_mode=execution_mode,
        delta_objective=delta,
        effect_label=effect_label,
        growth_risk_label=_growth_risk_label(row, recommendation_tier),
        secretory_process=str(row.get("secretory_process") or ""),
        evidence_tier=evidence_tier,
        wet_lab_readiness=str(row.get("wet_lab_readiness") or _wet_lab_readiness_from_evidence(evidence)),
        database_annotation_sources=tuple(str(item) for item in row.get("database_annotation_sources") or row.get("evidence_sources") or ()),
        database_annotation_confidence=str(row.get("database_annotation_confidence") or row.get("evidence_confidence") or ""),
        model_gpr_executable=model_gpr_executable,
        oe_reaction_proxy=oe_reaction_proxy,
        phenotype_evidence=phenotype_payload,
        recommendation_tier=recommendation_tier,
        recommendation_tier_reason=tier_reason,
        recommendation_label=label,
        recommendation_score=score,
        rationale=_yield_rationale(label, delta, execution_mode, evidence_tier, row),
        warnings=warnings,
    )


def _evidence_for_candidate(
    row: dict[str, Any],
    evidence_index: dict[str, dict[str, object]],
) -> dict[str, object]:
    keys = (
        row.get("input_gene_id"),
        row.get("canonical_gene_id"),
        row.get("gene_id"),
        row.get("candidate_id"),
        row.get("resolved_reaction_id"),
        row.get("reaction_id"),
    )
    for key in keys:
        text = str(key or "").strip().lower()
        if text and text in evidence_index:
            return evidence_index[text]
    return {}


def _candidate_recommendation_tier(
    row: dict[str, Any],
    evidence: dict[str, object],
    intervention_type: str,
    target_protein_context: str,
    *,
    model_gpr_executable: bool,
    oe_reaction_proxy: bool,
) -> tuple[str, str, Any | None]:
    row_tier = str(row.get("recommendation_tier") or "").strip()
    if row_tier:
        return row_tier, str(row.get("recommendation_tier_reason") or ""), None
    gene_id = str(
        row.get("canonical_gene_id")
        or row.get("gene_id")
        or evidence.get("mapped_model_gene_id")
        or evidence.get("declared_model_gene_id")
        or row.get("input_gene_id")
        or evidence.get("common_name")
        or ""
    )
    aliases = tuple(
        str(value)
        for value in (
            row.get("input_gene_id"),
            row.get("display_name"),
            row.get("standard_gene_symbol"),
            evidence.get("common_name"),
        )
        if str(value or "").strip()
    )
    return recommendation_tier_for_candidate(
        gene_id=gene_id,
        intervention_type=intervention_type,
        target_protein_context=target_protein_context,
        model_gpr_executable=model_gpr_executable,
        oe_reaction_proxy=oe_reaction_proxy,
        resolved=not str(row.get("status") or "").startswith("unresolved"),
        database_annotation_available=bool(row.get("evidence_sources") or evidence),
        aliases=aliases,
    )


def _yield_execution_mode(row: dict[str, Any], evidence: dict[str, object]) -> str:
    intervention_type = str(row.get("intervention_type") or "")
    ko_status = str(row.get("ko_support_status") or "")
    if intervention_type == "KO" and ko_status == "ko_runnable_gpr_gene_deletion":
        return "gene_level_ko"
    if intervention_type == "KO" and (row.get("gene_id") or row.get("canonical_gene_id")) and str(row.get("mapping_level") or "") != "unresolved":
        return "gene_level_ko"
    if intervention_type in {"OE_gene_proxy", "OE_reaction"}:
        return "reaction_level_oe_proxy"
    if intervention_type == "KO_reaction":
        return "reaction_level_ko_proxy"
    if evidence.get("reaction_proxy_ready") or evidence.get("ko_reaction_id") or evidence.get("oe_reaction_id"):
        return "reaction_level_proxy"
    if evidence:
        return "evidence_only_manual_review"
    return "unresolved_or_unknown"


def _yield_recommendation_label(
    *,
    success: bool,
    unresolved: bool,
    delta_objective: float | None,
    execution_mode: str,
    evidence: dict[str, object],
    recommendation_tier: str = "",
) -> str:
    if recommendation_tier == "not_recommended_growth_risk":
        return "not_recommended_growth_risk"
    if unresolved or execution_mode == "unresolved_or_unknown":
        return "unresolved_not_actionable"
    if not success:
        return "not_recommended_solver_failed"
    if delta_objective is None or delta_objective <= 0:
        return "not_recommended_no_model_gain"
    if execution_mode == "gene_level_ko":
        return "strong_model_candidate"
    if execution_mode in {"reaction_level_oe_proxy", "reaction_level_ko_proxy", "reaction_level_proxy"}:
        return "promising_but_proxy_only"
    if evidence:
        return "biology_interesting_manual_review"
    return "unresolved_not_actionable"


def _yield_recommendation_score(
    label: str,
    delta_objective: float | None,
    execution_mode: str,
    evidence: dict[str, object],
) -> float:
    if label in {
        "not_recommended_solver_failed",
        "not_recommended_no_model_gain",
        "not_recommended_growth_risk",
        "unresolved_not_actionable",
    }:
        return 0.0
    delta_score = min(40.0, max(0.0, float(delta_objective or 0.0) * 1_000_000.0))
    mode_score = {
        "gene_level_ko": 35.0,
        "reaction_level_oe_proxy": 22.0,
        "reaction_level_ko_proxy": 18.0,
        "reaction_level_proxy": 15.0,
        "evidence_only_manual_review": 8.0,
    }.get(execution_mode, 0.0)
    evidence_score = 15.0 if evidence.get("mapping_status") == "model_gpr_gene_available" else 10.0 if evidence else 0.0
    return round(delta_score + mode_score + evidence_score, 3)


def _yield_recommendation_sort_key(item: YieldImprovementCandidateRecommendation) -> tuple[int, float, float, str]:
    priority = {
        "strong_model_candidate": 0,
        "promising_but_proxy_only": 1,
        "biology_interesting_manual_review": 2,
        "not_recommended_no_model_gain": 3,
        "not_recommended_solver_failed": 4,
        "not_recommended_growth_risk": 4,
        "unresolved_not_actionable": 5,
    }.get(item.recommendation_label, 9)
    return (priority, -float(item.recommendation_score), -float(item.delta_objective or 0.0), item.display_name)


def _yield_display_name(row: dict[str, Any], evidence: dict[str, object]) -> str:
    for value in (
        evidence.get("common_name"),
        row.get("display_name"),
        row.get("standard_gene_symbol"),
        row.get("input_gene_id"),
        row.get("gene_id"),
        row.get("candidate_id"),
        row.get("reaction_id"),
    ):
        text = str(value or "").strip()
        if text:
            return text
    return "unknown_candidate"


def _yield_evidence_tier(row: dict[str, Any], evidence: dict[str, object], execution_mode: str) -> str:
    if execution_mode == "gene_level_ko":
        return "模型可执行 GPR + curated 证据" if evidence else "模型可执行 GPR"
    if execution_mode in {"reaction_level_oe_proxy", "reaction_level_ko_proxy", "reaction_level_proxy"}:
        return "模型反应代理 + curated 证据" if evidence else "模型反应代理"
    confidence = str(row.get("evidence_confidence") or "")
    if confidence:
        return confidence
    if evidence:
        return "curated 证据；需人工确认"
    return "无可用证据"


def _wet_lab_readiness_from_evidence(evidence: dict[str, object]) -> str:
    if evidence.get("mapping_status") == "model_gpr_gene_available":
        return "database_supported_experiment_candidate"
    if evidence:
        return "manual_review_required"
    return "model_only_not_experiment_ready"


def _yield_candidate_warnings(label: str, execution_mode: str, evidence: dict[str, object]) -> list[str]:
    warnings: list[str] = []
    if execution_mode in {"reaction_level_oe_proxy", "reaction_level_ko_proxy", "reaction_level_proxy"}:
        warnings.append("该候选是 reaction-level proxy，不是严格 gene-level 湿实验操作。")
    if label == "not_recommended_solver_failed":
        warnings.append("当前固定生长率和约束组合下求解失败或不可行。")
    if label == "not_recommended_growth_risk":
        warnings.append("该 KO 候选有 essentiality 表型证据，优先标记为生长风险。")
    if label == "unresolved_not_actionable":
        warnings.append("该候选未解析到可执行模型基因或反应。")
    if evidence and evidence.get("mapping_status") != "model_gpr_gene_available":
        warnings.append("湿实验前需要确认 K. phaffii locus/GPR 证据。")
    return warnings


def _yield_rationale(
    label: str,
    delta_objective: float | None,
    execution_mode: str,
    evidence_tier: str,
    row: dict[str, Any],
) -> str:
    delta_text = "NA" if delta_objective is None else f"{delta_objective:.6g}"
    process = str(row.get("secretory_process") or "未标注环节")
    if label == "strong_model_candidate":
        return f"模型显示分泌目标提升 Δobjective={delta_text}，且操作为 gene-level KO；关联环节：{process}。"
    if label == "promising_but_proxy_only":
        return f"模型显示分泌目标提升 Δobjective={delta_text}，但操作为 {execution_mode}；关联环节：{process}。"
    if label == "biology_interesting_manual_review":
        return f"候选有生物学/证据意义，但当前只适合作人工复核；证据等级：{evidence_tier}。"
    if label == "not_recommended_no_model_gain":
        return f"当前模型未显示提升（Δobjective={delta_text}），暂不作为优先候选。"
    if label == "not_recommended_solver_failed":
        return "当前候选求解失败或约束不可行，不能作为模型内提升推荐。"
    if label == "not_recommended_growth_risk":
        return "该 KO 候选有 essentiality 表型证据，优先标记为生长风险，不作为推荐扰动。"
    return "候选未解析到可执行模型基因或反应，暂不可操作。"


def _growth_risk_label(row: dict[str, Any], recommendation_tier: str = "") -> str:
    if str(recommendation_tier or row.get("recommendation_tier") or "") == "not_recommended_growth_risk":
        return "essential_gene_ko"
    status = str(row.get("status") or "")
    if status.startswith("unresolved"):
        return "unresolved"
    if not _truthy(row.get("success")):
        return "solver_failed_or_infeasible"
    delta = _optional_float(row.get("delta_objective"))
    if delta is not None and delta < 0:
        return "secretion_decrease"
    return "not_assessed"


def _truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"1", "true", "yes", "y", "success", "optimal"}


def _optional_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _target_feature_payload(target: TargetSpec) -> dict[str, Any]:
    features = target_features(target)
    return {
        "full_sequence_length": int(features.get("full_length") or 0),
        "mature_sequence_length": int(features.get("mature_length") or 0),
        "leader_sequence_length": int(features.get("leader_length") or 0),
        "signal_peptide_length": int(features.get("signal_peptide_length") or 0),
        "protein_mw": float(features.get("protein_mw") or 0.0),
        "cysteine_count": int(features.get("cysteine_count") or 0),
        "inferred_disulfide_pairs_from_cysteines": int(
            features.get("inferred_disulfide_pairs_from_cysteines") or 0
        ),
        "n_glycosylation_motif_positions": list(features.get("n_glycosylation_motif_positions") or []),
        "ser_thr_count": int(features.get("ser_thr_count") or 0),
        "valid_mature_sequence": bool(features.get("valid_mature_sequence")),
    }


def _cost_items_from_target_and_plan(
    target: TargetSpec,
    plan: SecretionPlanResult,
    features: dict[str, Any],
) -> tuple[ProteinCostItem, ...]:
    full_len = float(features["full_sequence_length"])
    mature_len = float(features["mature_sequence_length"])
    leader_len = float(features["leader_sequence_length"])
    signal_len = float(features["signal_peptide_length"])
    stage_counts = plan.stage_counts
    dsb = int(target.disulfide_sites)
    ng = int(target.n_glycosylation_sites)
    og = int(target.o_glycosylation_sites)
    motif_count = len(features["n_glycosylation_motif_positions"])
    ser_thr_count = int(features["ser_thr_count"])
    ptm_burden = dsb + ng + og + int(target.gpi_sites) + int(target.transmembrane)

    return (
        ProteinCostItem(
            "translation",
            "翻译负担",
            "full sequence length",
            full_len,
            0.0,
            "目标蛋白越长，翻译相关氨基酸与核糖体占用越高。",
        ),
        ProteinCostItem(
            "leader_signal_processing",
            "信号肽/leader 处理",
            "leader length + signal peptide length",
            leader_len + signal_len,
            0.0,
            "分泌 leader 和信号肽会增加进入分泌途径前后的处理负担。",
        ),
        ProteinCostItem(
            "er_translocation",
            "ER 转运",
            "through ER route and translocation reactions",
            float(stage_counts.get("translocation", 0) * 20 + (10 if target.through_er else 0)),
            0.0,
            "经 ER 分泌的目标蛋白需要转运与 Sec 通路相关步骤。",
        ),
        ProteinCostItem(
            "folding_dsb",
            "二硫键折叠 DSB",
            "declared DSB count",
            float(dsb * 18),
            0.0,
            "声明的二硫键数量越多，PDI/氧化折叠相关负担越高。",
        ),
        ProteinCostItem(
            "n_glycosylation",
            "N-糖基化 NG",
            "declared NG count and motif count",
            float(ng * 25 + (motif_count * 2 if ng > 0 else 0)),
            0.0,
            "N-糖基化位点和可见 motif 会增加 ER/Golgi 加工解释成本。",
        ),
        ProteinCostItem(
            "o_glycosylation",
            "O-糖基化 OG",
            "declared OG count and Ser/Thr support",
            float(og * 12 + min(ser_thr_count, og)),
            0.0,
            "O-糖基化位点会增加 ER/Golgi O-linked processing 负担。",
        ),
        ProteinCostItem(
            "misfolding_erad",
            "错误折叠/ERAD 风险",
            "length, PTM burden, and misfolding plan rows",
            float(mature_len / 20 + ptm_burden * 3 + stage_counts.get("misfolding", 0) * 10),
            0.0,
            "长序列和复杂 PTM 更容易带来折叠质量控制压力。",
        ),
        ProteinCostItem(
            "transport_secretion",
            "囊泡运输与分泌",
            "ER-Golgi, Golgi, final transport, exchange rows",
            float(
                (
                    stage_counts.get("er_to_golgi", 0)
                    + stage_counts.get("golgi_processing", 0)
                    + stage_counts.get("final_transport", 0)
                    + stage_counts.get("exchange", 0)
                )
                * 10
            ),
            0.0,
            "ER-Golgi 与最终分泌步骤越多，分泌路径解释成本越高。",
        ),
    )


def _normalise_cost_scores(items: tuple[ProteinCostItem, ...]) -> tuple[ProteinCostItem, ...]:
    total = sum(max(0.0, item.raw_value) for item in items)
    if total <= 0:
        return items
    return tuple(
        ProteinCostItem(
            category=item.category,
            label=item.label,
            basis=item.basis,
            raw_value=round(item.raw_value, 6),
            relative_score=round(max(0.0, item.raw_value) / total * 100.0, 3),
            interpretation=item.interpretation,
        )
        for item in items
    )


def _dominant_categories(items: tuple[ProteinCostItem, ...], limit: int = 3) -> tuple[str, ...]:
    ordered = sorted(items, key=lambda item: item.relative_score, reverse=True)
    return tuple(item.category for item in ordered[: max(0, int(limit))] if item.relative_score > 0)


def _cost_warnings(target: TargetSpec) -> list[str]:
    warnings = [
        "当前成本分析是 Python draft explanatory score，不代表真实发酵产量或湿实验成本。",
        "该分析不使用 LP shadow price，也不改变 corrected pipeline 的求解目标或约束。",
    ]
    if target.target_id == "hLF":
        warnings.append("hLF 使用项目定义 710aa 序列；结果不是 old MATLAB hLF fully aligned。")
    if target.source.startswith("request") or target.source.startswith("json:") or target.target_id.upper().startswith("CUSTOM"):
        if target.disulfide_sites == 0 and target.n_glycosylation_sites == 0 and target.o_glycosylation_sites == 0:
            warnings.append("custom target 的 DSB/NG/OG 均为 0；未声明 PTM 不会被自动推断进成本项。")
    return warnings


def _objective_evidence(simulation_result: SecretionSimulationResult) -> dict[str, object]:
    return {
        "success": simulation_result.success,
        "objective_reaction": simulation_result.exchange_reaction_id,
        "objective_value": simulation_result.objective_value,
        "secretion_flux": simulation_result.secretion_flux,
        "growth_rate": simulation_result.growth_rate,
        "status": simulation_result.status,
        "message": simulation_result.message,
    }


def _constraint_blocks(
    constraint_counts: dict[str, int],
    eq_length: int,
    ub_length: int,
) -> tuple[dict[str, object], ...]:
    blocks: list[dict[str, object]] = []
    start = 0
    for name in (
        "stoichiometric",
        "metabolic_coupling",
        "secretory_coupling",
        "protein_mass",
        "proteasome",
        "ribosome_assembly",
        "ribosome_translation",
        "misfolding",
    ):
        count = max(0, int(constraint_counts.get(name, 0) or 0))
        if count <= 0:
            continue
        end = min(eq_length, start + count)
        blocks.append({"constraint_type": "eq", "block": name, "start": start, "end": end})
        start += count
    if start < eq_length:
        blocks.append({"constraint_type": "eq", "block": "unknown", "start": start, "end": eq_length})
    if ub_length > 0:
        mito_count = max(0, int(constraint_counts.get("mitochondrial", 0) or 0))
        mito_end = min(ub_length, mito_count)
        if mito_end > 0:
            blocks.append({"constraint_type": "ub", "block": "mitochondrial", "start": 0, "end": mito_end})
        if mito_end < ub_length:
            blocks.append({"constraint_type": "ub", "block": "unknown", "start": mito_end, "end": ub_length})
    return tuple(blocks)


def _dominant_constraint_blocks(
    blocks: tuple[dict[str, object], ...],
    eq_marginals: tuple[float, ...],
    ub_marginals: tuple[float, ...],
    top_n: int,
) -> tuple[dict[str, object], ...]:
    rows: list[dict[str, object]] = []
    for block in blocks:
        values = _block_values(block, eq_marginals, ub_marginals)
        if not values:
            continue
        abs_values = [abs(value) for value in values]
        rows.append(
            {
                "constraint_type": block["constraint_type"],
                "block": block["block"],
                "row_count": len(values),
                "nonzero_marginal_count": _nonzero_count(values),
                "sum_abs_marginal": round(sum(abs_values), 12),
                "max_abs_marginal": round(max(abs_values), 12),
            }
        )
    return tuple(sorted(rows, key=lambda row: float(row["sum_abs_marginal"]), reverse=True)[:top_n])


def _top_constraint_marginals(
    blocks: tuple[dict[str, object], ...],
    eq_marginals: tuple[float, ...],
    ub_marginals: tuple[float, ...],
    top_n: int,
) -> tuple[dict[str, object], ...]:
    rows: list[dict[str, object]] = []
    for block in blocks:
        values = _block_values(block, eq_marginals, ub_marginals)
        start = int(block["start"])
        for offset, marginal in enumerate(values):
            if abs(marginal) <= 1e-12:
                continue
            row_index = start + offset
            rows.append(
                {
                    "constraint_type": block["constraint_type"],
                    "block": block["block"],
                    "row_index_0based": row_index,
                    "row_index_1based": row_index + 1,
                    "marginal": round(float(marginal), 12),
                    "abs_marginal": round(abs(float(marginal)), 12),
                }
            )
    return tuple(sorted(rows, key=lambda row: float(row["abs_marginal"]), reverse=True)[:top_n])


def _top_bound_marginals(
    lower_marginals: tuple[float, ...],
    upper_marginals: tuple[float, ...],
    reaction_ids: tuple[str, ...],
    target: TargetSpec,
    plan: SecretionPlanResult,
    top_n: int,
) -> tuple[dict[str, object], ...]:
    rows: list[dict[str, object]] = []
    for bound_type, values in (("lower", lower_marginals), ("upper", upper_marginals)):
        for index, marginal in enumerate(values):
            if abs(marginal) <= 1e-12:
                continue
            reaction_id = reaction_ids[index] if index < len(reaction_ids) else f"X{index + 1}"
            rows.append(
                {
                    "bound_type": bound_type,
                    "variable_index_0based": index,
                    "variable_id": f"X{index + 1}",
                    "reaction_id": reaction_id,
                    "secretory_process": _lp_reaction_process(reaction_id, target, plan),
                    "marginal": round(float(marginal), 12),
                    "abs_marginal": round(abs(float(marginal)), 12),
                }
            )
    return tuple(sorted(rows, key=lambda row: float(row["abs_marginal"]), reverse=True)[:top_n])


def _target_related_fluxes(
    target: TargetSpec,
    plan: SecretionPlanResult,
    simulation_result: SecretionSimulationResult,
    reaction_ids: tuple[str, ...],
) -> tuple[dict[str, object], ...]:
    fluxes = simulation_result.key_fluxes or {}
    rows: list[dict[str, object]] = []
    for reaction_id, flux in sorted(fluxes.items()):
        rows.append(
            {
                "reaction_id": reaction_id,
                "flux": flux,
                "secretory_process": _lp_reaction_process(reaction_id, target, plan),
                "is_target_exchange": reaction_id == simulation_result.exchange_reaction_id,
            }
        )
    if simulation_result.exchange_reaction_id and simulation_result.exchange_reaction_id not in fluxes:
        rows.append(
            {
                "reaction_id": simulation_result.exchange_reaction_id,
                "flux": simulation_result.secretion_flux,
                "secretory_process": "target_exchange",
                "is_target_exchange": True,
            }
        )
    return tuple(rows)


def _block_values(
    block: dict[str, object],
    eq_marginals: tuple[float, ...],
    ub_marginals: tuple[float, ...],
) -> tuple[float, ...]:
    values = eq_marginals if block["constraint_type"] == "eq" else ub_marginals
    start = int(block["start"])
    end = int(block["end"])
    return values[start:end]


def _lp_reaction_process(reaction_id: str, target: TargetSpec, plan: SecretionPlanResult) -> str:
    if reaction_id == "":
        return "unknown"
    if reaction_id in plan.reaction_ids:
        return "target_secretory_reaction"
    if reaction_id == f"r_{target.protein_id}_exchange" or reaction_id == f"{target.protein_id}_exchange":
        return "target_exchange"
    if target.protein_id and target.protein_id in reaction_id:
        return "target_related"
    lowered = reaction_id.lower()
    if "ribosome" in lowered:
        return "ribosome"
    if "misfold" in lowered or "erad" in lowered:
        return "misfolding_erad"
    if "glyco" in lowered or "_ng" in lowered or "_og" in lowered:
        return "glycosylation"
    if "biomass" in lowered:
        return "growth"
    if reaction_id in {"Ex_glc_D", "Ex_o2"}:
        return "medium_exchange"
    return "unknown"


def _nonzero_count(values: tuple[float, ...] | list[float]) -> int:
    return sum(1 for value in values if abs(float(value)) > 1e-12)


def _growth_tradeoff_points(tradeoff: GrowthTradeoffResult) -> list[GrowthTradeoffPoint]:
    points: list[GrowthTradeoffPoint] = []
    for row in tradeoff.tradeoff_rows:
        success = _optional_bool(row.get("success"))
        status = str(row.get("status") or "")
        flux = _optional_float(row.get("secretion_flux"))
        per_biomass = _optional_float(row.get("secretion_per_biomass"))
        interpretation = "可比较生长点" if success and flux is not None else "求解失败或无可比较分泌值"
        points.append(
            GrowthTradeoffPoint(
                mu=float(row.get("mu") or 0.0),
                success=success,
                secretion_flux=flux,
                secretion_per_biomass=per_biomass,
                status=status,
                interpretation=interpretation,
            )
        )
    return sorted(points, key=lambda point: point.mu)


def _with_growth_point_interpretation(
    point: GrowthTradeoffPoint,
    best_flux: GrowthTradeoffPoint | None,
    best_per_biomass: GrowthTradeoffPoint | None,
    baseline_growth_rate: float | None,
) -> GrowthTradeoffPoint:
    if not point.success or point.secretion_flux is None:
        interpretation = "求解失败或无可比较分泌值"
    else:
        labels: list[str] = []
        if best_flux is not None and point.mu == best_flux.mu:
            labels.append("最高分泌通量")
        if best_per_biomass is not None and point.mu == best_per_biomass.mu:
            labels.append("最高单位生物量分泌")
        if baseline_growth_rate is not None and abs(point.mu - float(baseline_growth_rate)) <= 1e-9:
            labels.append("基准生长点")
        interpretation = "；".join(labels) if labels else "可比较生长点"
    return GrowthTradeoffPoint(
        mu=point.mu,
        success=point.success,
        secretion_flux=point.secretion_flux,
        secretion_per_biomass=point.secretion_per_biomass,
        status=point.status,
        interpretation=interpretation,
    )


def _best_point(points: tuple[GrowthTradeoffPoint, ...], key: str) -> GrowthTradeoffPoint | None:
    comparable = [point for point in points if getattr(point, key) is not None]
    if not comparable:
        return None
    return max(comparable, key=lambda point: float(getattr(point, key) or 0.0))


def _growth_sensitivity(points: tuple[GrowthTradeoffPoint, ...]) -> tuple[str, str]:
    if len(points) < 2:
        return "insufficient_points", "insufficient_tradeoff_rows"
    if any((not point.success) or point.secretion_flux is None for point in points):
        return "mixed", "contains_failed_or_missing_points"
    valid_points = tuple(point for point in points if point.success and point.secretion_flux is not None)
    if len(valid_points) < 2:
        return "insufficient_points", "insufficient_comparable_points"
    ordered = sorted(valid_points, key=lambda point: point.mu)
    values = [float(point.secretion_flux or 0.0) for point in ordered]
    deltas = [right - left for left, right in zip(values, values[1:])]
    eps = 1e-12
    if all(delta > eps for delta in deltas):
        return "increasing", "monotonic_increasing_successful_grid"
    if all(delta < -eps for delta in deltas):
        return "decreasing", "monotonic_decreasing_successful_grid"
    return "mixed", "non_monotonic_successful_grid"


def _optional_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y", "success", "ok"}:
        return True
    if text in {"false", "0", "no", "n", "failed", "failure", "none", ""}:
        return False
    return False


def _optional_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


__all__ = [
    "GrowthTradeoffPoint",
    "ProteinLpAttributionResult",
    "ProteinCostAnalysisResult",
    "ProteinCostItem",
    "TargetGrowthAnalysisResult",
    "YieldImprovementCandidateRecommendation",
    "YieldImprovementRecommendationResult",
    "analyze_target_protein_cost",
    "analyze_target_protein_lp_attribution",
    "analyze_target_growth_impact",
    "analyze_yield_improvement_candidates",
    "build_cost_item_table",
    "build_growth_tradeoff_item_table",
    "build_yield_recommendation_table",
    "summarize_protein_lp_attribution",
    "summarize_protein_cost_slope_compatibility",
    "summarize_protein_cost_analysis",
    "summarize_target_growth_analysis",
    "summarize_yield_improvement_recommendations",
]
