from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pcsec_pichia.screens._prototype_adapter import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    TargetSpec,
    build_supported_target_model,
    build_target_enzymedata,
    classify_candidate_effect,
    classify_secretory_process,
    default_ko_genes,
    default_oe_reactions,
    run_ko_screen,
    run_oe_screen,
    run_pcsec_growth_tradeoff,
    run_pcsec_ko_screen,
    run_pcsec_oe_screen,
    run_pcsec_reaction_ko_screen,
    solve_pcsec_maximize,
)
from pcsec_pichia.screens.candidate_resolution import (
    reactions_for_gene,
    resolve_oe_gene_reactions,
    split_existing_genes,
    split_existing_reactions,
)
from pcsec_pichia.screens.gene_perturbation_map import (
    GenePerturbationMapping,
    GenePerturbationMapResult,
    GeneReactionMapping,
    build_gene_perturbation_map,
)


@dataclass(frozen=True)
class ScreenResult:
    target_id: str
    screen_type: str
    success: bool
    candidate_count: int
    rows: tuple[dict[str, Any], ...]
    constraint_counts: dict[str, int]
    baseline_objective_value: float | None
    result_status: str
    matlab_alignment_status: str


def run_knockout_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    genes: list[str],
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> ScreenResult:
    prepared = _prepare_screen_inputs(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    if not prepared["baseline_success"]:
        return _empty_screen_result(target.target_id, "knockout", prepared)

    raw_rows = run_pcsec_ko_screen(
        prepared["fixed_model"],
        prepared["baseline"],
        genes,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    rows = tuple(
        _normalize_screen_row(
            row,
            target_id=target.target_id,
            screen_type="knockout",
            intervention_type="KO",
            baseline_objective_value=prepared["baseline"].objective_value,
            complex_subunits=prepared["secretory"].complex_subunits,
            input_gene_id=str(row.get("gene")) if row.get("gene") is not None else None,
        )
        for row in raw_rows
    )
    return _screen_result(target.target_id, "knockout", rows, prepared)


def run_overexpression_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    reactions: list[str],
    growth_rate: float = 0.10,
    factor: float = 2.0,
    intervention_type: str = "OE_reaction",
    input_gene_ids_by_reaction: dict[str, str] | None = None,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> ScreenResult:
    prepared = _prepare_screen_inputs(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    if not prepared["baseline_success"]:
        return _empty_screen_result(target.target_id, "overexpression", prepared)

    raw_rows = run_pcsec_oe_screen(
        prepared["fixed_model"],
        prepared["baseline"],
        reactions,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        factor=factor,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    rows = tuple(
        _normalize_screen_row(
            row,
            target_id=target.target_id,
            screen_type="overexpression",
            intervention_type=intervention_type,
            baseline_objective_value=prepared["baseline"].objective_value,
            complex_subunits=prepared["secretory"].complex_subunits,
            input_gene_id=(input_gene_ids_by_reaction or {}).get(str(row.get("reaction"))),
        )
        for row in raw_rows
    )
    return _screen_result(target.target_id, "overexpression", rows, prepared)


def run_reaction_knockout_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    reactions: list[str],
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> ScreenResult:
    prepared = _prepare_screen_inputs(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    if not prepared["baseline_success"]:
        return _empty_screen_result(target.target_id, "knockout", prepared)

    raw_rows = run_pcsec_reaction_ko_screen(
        prepared["fixed_model"],
        prepared["baseline"],
        reactions,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    rows = tuple(
        _normalize_screen_row(
            row,
            target_id=target.target_id,
            screen_type="knockout",
            intervention_type="KO_reaction",
            baseline_objective_value=prepared["baseline"].objective_value,
            complex_subunits=prepared["secretory"].complex_subunits,
        )
        for row in raw_rows
    )
    return _screen_result(target.target_id, "knockout", rows, prepared)


def summarize_screen_result(result: ScreenResult) -> dict[str, Any]:
    return {
        "target_id": result.target_id,
        "screen_type": result.screen_type,
        "success": result.success,
        "candidate_count": result.candidate_count,
        "rows": result.rows,
        "constraint_counts": result.constraint_counts,
        "baseline_objective_value": result.baseline_objective_value,
        "result_status": result.result_status,
        "matlab_alignment_status": result.matlab_alignment_status,
    }


def _prepare_screen_inputs(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    growth_rate: float,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> dict[str, Any]:
    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        return {
            "baseline_success": False,
            "baseline": None,
            "constraint_counts": {},
            "fixed_model": None,
            "secretory": None,
            "combined": None,
            "exchange_reaction_id": build.exchange_reaction_id,
        }

    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = combined.with_target(target_enzymedata)
    fixed_model = build.model.with_bounds({"BIOMASS": (growth_rate, growth_rate)})
    baseline, counts = solve_pcsec_maximize(
        fixed_model,
        build.exchange_reaction_id,
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu=growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    return {
        "baseline_success": baseline.success,
        "baseline": baseline,
        "constraint_counts": {str(key): int(value) for key, value in counts.items()},
        "fixed_model": fixed_model,
        "secretory": target_secretory,
        "combined": target_combined,
        "exchange_reaction_id": build.exchange_reaction_id,
    }


def _normalize_screen_row(
    row: dict[str, Any],
    target_id: str,
    screen_type: str,
    intervention_type: str,
    baseline_objective_value: float | None,
    complex_subunits: dict[str, list[dict[str, object]]] | None,
    input_gene_id: str | None = None,
    resolved_reaction_id: str | None = None,
) -> dict[str, Any]:
    gene_id = row.get("gene")
    reaction_id = row.get("reaction")
    if reaction_id is None:
        preview = row.get("inactive_reactions_preview") or []
        reaction_id = preview[0] if preview else None
    complex_id = str(reaction_id).replace("_formation", "") if reaction_id else ""
    subunits = (complex_subunits or {}).get(complex_id, [])
    objective_value = row.get("objective_value")
    delta = row.get("delta_vs_baseline")
    normalized_delta = float(delta) if delta is not None else None
    normalized_baseline = float(baseline_objective_value) if baseline_objective_value is not None else None
    relative_delta = _relative_delta(normalized_delta, normalized_baseline)
    resolved_reaction = resolved_reaction_id or (str(reaction_id) if reaction_id is not None else None)
    process_code = classify_secretory_process(resolved_reaction)
    effect_code = classify_candidate_effect(bool(row.get("success")), relative_delta)
    solver_status_label = _solver_status_label(row.get("status"), bool(row.get("success")))
    return {
        **row,
        "target_id": target_id,
        "screen_type": screen_type,
        "intervention_type": intervention_type,
        "candidate_id": str(gene_id or reaction_id or ""),
        "gene_id": str(gene_id) if gene_id is not None else None,
        "reaction_id": str(reaction_id) if reaction_id is not None else None,
        "input_gene_id": input_gene_id,
        "resolved_reaction_id": resolved_reaction,
        "objective_value": float(objective_value) if objective_value is not None else None,
        "baseline_objective_value": baseline_objective_value,
        "delta_objective": normalized_delta,
        "effect_label": _effect_label(effect_code, row.get("status")),
        "solver_status_label": solver_status_label,
        "failure_reason": None if row.get("success") else solver_status_label,
        "secretory_process": _secretory_process_label(process_code),
        "complex_subunit_ids": [str(item["subunit_id"]) for item in subunits],
        "complex_subunit_stoichiometry": [float(item["stoichiometry"]) for item in subunits],
    }


def _relative_delta(delta: float | None, baseline_value: float | None) -> float | None:
    if delta is None or baseline_value is None or abs(baseline_value) < 1e-15:
        return None
    return float(delta) / abs(float(baseline_value))


def _effect_label(effect_code: str, raw_status: object | None = None) -> str:
    if effect_code in {"strong_improvement", "weak_improvement"}:
        return "提升分泌"
    if effect_code in {"strong_decrease", "weak_decrease"}:
        return "降低分泌"
    if effect_code == "neutral":
        return "无明显变化"
    if effect_code == "infeasible_at_fixed_mu":
        if str(raw_status) == "2":
            return "约束不可行"
        if str(raw_status) in {"missing_reaction", "unresolved_gene", "unresolved_reaction"}:
            return "未解析"
        return "求解失败"
    return "未解析"


def _solver_status_label(raw_status: object | None, success: bool) -> str:
    if success:
        return "求解成功"
    status = "" if raw_status is None else str(raw_status)
    labels = {
        "2": "约束不可行",
        "3": "目标无界",
        "4": "求解器数值错误",
        "missing_reaction": "反应未找到",
        "unresolved_gene": "基因未解析",
        "unresolved_reaction": "反应未解析",
        "missing_objective": "目标反应未找到",
    }
    return labels.get(status, "求解失败")


def _secretory_process_label(process_code: str) -> str:
    labels = {
        "ribosome": "翻译",
        "proteasome_degradation": "蛋白降解",
        "disulfide_folding": "ER 折叠 / DSB",
        "n_glycan_processing": "N-糖基化 NG",
        "o_glycan_processing": "O-糖基化 OG",
        "chaperone_folding": "ER 折叠 / 分子伴侣",
        "erad_misfolding": "错误折叠 / ERAD",
        "er_translocation": "ER 转运",
        "er_to_golgi_transport": "ER 到 Golgi 转运",
        "golgi_surface_transport": "Golgi 到胞外运输",
        "secretory_capacity": "分泌容量",
        "metabolic_or_other": "代谢或其它反应",
        "unknown": "未解析",
    }
    return labels.get(process_code, process_code)


def _screen_result(
    target_id: str,
    screen_type: str,
    rows: tuple[dict[str, Any], ...],
    prepared: dict[str, Any],
) -> ScreenResult:
    return ScreenResult(
        target_id=target_id,
        screen_type=screen_type,
        success=bool(rows) and all(bool(row.get("success")) for row in rows),
        candidate_count=len(rows),
        rows=rows,
        constraint_counts=prepared["constraint_counts"],
        baseline_objective_value=prepared["baseline"].objective_value,
        result_status="draft",
        matlab_alignment_status="pending",
    )


def _empty_screen_result(target_id: str, screen_type: str, prepared: dict[str, Any]) -> ScreenResult:
    baseline = prepared.get("baseline")
    return ScreenResult(
        target_id=target_id,
        screen_type=screen_type,
        success=False,
        candidate_count=0,
        rows=(),
        constraint_counts=prepared.get("constraint_counts", {}),
        baseline_objective_value=getattr(baseline, "objective_value", None),
        result_status="draft",
        matlab_alignment_status="pending",
    )


__all__ = [
    "ScreenResult",
    "GenePerturbationMapping",
    "GenePerturbationMapResult",
    "GeneReactionMapping",
    "build_gene_perturbation_map",
    "default_ko_genes",
    "default_oe_reactions",
    "reactions_for_gene",
    "resolve_oe_gene_reactions",
    "split_existing_genes",
    "run_knockout_screen",
    "run_ko_screen",
    "run_overexpression_screen",
    "run_oe_screen",
    "run_reaction_knockout_screen",
    "run_pcsec_growth_tradeoff",
    "run_pcsec_ko_screen",
    "run_pcsec_oe_screen",
    "run_pcsec_reaction_ko_screen",
    "split_existing_reactions",
    "summarize_screen_result",
]
