from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pcsec_pichia.analysis import analyze_target_protein_cost, summarize_protein_cost_analysis
from pcsec_pichia.alignment import (
    build_alignment_summary,
    hlf_project_710_known_matlab_compatibility_exceptions,
    opn_known_matlab_compatibility_exceptions,
    summarize_alignment,
)
from pcsec_pichia.constraints import build_pcsec_constraints
from pcsec_pichia.engines.base import PichiaSimulationRequest, PichiaSimulationRunResult
from pcsec_pichia.loading import CobraModel, load_pcsec_pichia_inputs
from pcsec_pichia.reports import write_simulation_outputs
from pcsec_pichia.screens import (
    ScreenResult,
    default_ko_genes,
    default_oe_reactions,
    resolve_oe_gene_reactions,
    run_knockout_screen,
    run_overexpression_screen,
    run_reaction_knockout_screen,
    split_existing_genes,
    split_existing_reactions,
)
from pcsec_pichia.secretion_plan import build_secretion_plan
from pcsec_pichia.simulation import run_growth_tradeoff, solve_secretion_capacity
from pcsec_pichia.targets import (
    TargetSpec,
    load_builtin_targets,
    load_custom_targets_json,
    target_spec_from_input,
    target_spec_from_mapping,
)


DEFAULT_ALIGNMENT_ARTIFACT = Path("local_runs") / "pichia_hlf_opn_probe" / "matlab_stage3_alignment" / "matlab_stage3_alignment_summary.json"
HLF_PROJECT_710_ALIGNMENT_ARTIFACT = (
    Path("local_runs")
    / "pichia_hlf_opn_probe"
    / "hlf_project_sequence_matlab_harness_2026-06-26"
    / "hlf_project_sequence_matlab_harness_summary.json"
)


def run_pichia_secretion_simulation(
    request: PichiaSimulationRequest,
    output_dir: Path,
) -> PichiaSimulationRunResult:
    root = Path(__file__).resolve().parents[3]
    inputs = load_pcsec_pichia_inputs(root, media_type=request.media_type, compatibility_mode=request.compatibility_mode)
    target = _resolve_target(request, root)

    plan = build_secretion_plan(target)
    protein_cost = analyze_target_protein_cost(target, plan)
    constraint_result = build_pcsec_constraints(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        mu=request.mu,
        write_ribosome_translation_constraint=request.enable_ribosome_translation_constraint,
        write_misfolding_constraints=request.enable_misfolding_constraint,
    )
    simulation = solve_secretion_capacity(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_rate=request.mu,
        write_ribosome_translation_constraint=request.enable_ribosome_translation_constraint,
        write_misfolding_constraints=request.enable_misfolding_constraint,
    )
    tradeoff = run_growth_tradeoff(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_points=_growth_points(request.growth_points, request.mu),
        write_ribosome_translation_constraint=request.enable_ribosome_translation_constraint,
        write_misfolding_constraints=request.enable_misfolding_constraint,
    )
    screen_plan = _build_screen_plan(inputs.prepared_model, request)

    # Keep smoke screens sequential so shared model/enzyme inputs stay easy to
    # reason about while the pipeline is still draft-alignment code.
    model = inputs.prepared_model
    aa, met, sec, comb = inputs.amino_acids, inputs.metabolic, inputs.secretory, inputs.combined
    mu = request.mu
    rib = request.enable_ribosome_translation_constraint
    mis = request.enable_misfolding_constraint

    ko = run_knockout_screen(
        model,
        target,
        aa,
        met,
        sec,
        comb,
        genes=screen_plan["ko_gene_ids"],
        growth_rate=mu,
        write_ribosome_translation_constraint=rib,
        write_misfolding_constraints=mis,
    )
    ko_reaction = run_reaction_knockout_screen(
        model,
        target,
        aa,
        met,
        sec,
        comb,
        reactions=screen_plan["ko_reaction_ids"],
        growth_rate=mu,
        write_ribosome_translation_constraint=rib,
        write_misfolding_constraints=mis,
    )
    oe_gene = run_overexpression_screen(
        model,
        target,
        aa,
        met,
        sec,
        comb,
        reactions=screen_plan["oe_gene_reaction_ids"],
        growth_rate=mu,
        intervention_type="OE_gene_proxy",
        input_gene_ids_by_reaction=screen_plan["oe_gene_by_reaction"],
        write_ribosome_translation_constraint=rib,
        write_misfolding_constraints=mis,
    )
    oe_reaction = run_overexpression_screen(
        model,
        target,
        aa,
        met,
        sec,
        comb,
        reactions=screen_plan["oe_reaction_ids"],
        growth_rate=mu,
        intervention_type="OE_reaction",
        write_ribosome_translation_constraint=rib,
        write_misfolding_constraints=mis,
    )

    ko = _append_unresolved_rows(ko, tuple(
        _unresolved_row(target.target_id, "knockout", "KO", gid, "unresolved_gene", input_gene_id=gid,
                        baseline_objective_value=ko.baseline_objective_value)
        for gid in screen_plan["unresolved_ko_gene_ids"]))
    ko_reaction = _append_unresolved_rows(ko_reaction, tuple(
        _unresolved_row(target.target_id, "knockout", "KO_reaction", rid, "unresolved_reaction", reaction_id=rid,
                        baseline_objective_value=ko_reaction.baseline_objective_value)
        for rid in screen_plan["unresolved_ko_reaction_ids"]))
    oe_gene = _append_unresolved_rows(oe_gene, tuple(
        _unresolved_row(target.target_id, "overexpression", "OE_gene_proxy", gid, "unresolved_gene", input_gene_id=gid,
                        baseline_objective_value=oe_gene.baseline_objective_value)
        for gid in screen_plan["unresolved_oe_gene_ids"]))
    oe_reaction = _append_unresolved_rows(oe_reaction, tuple(
        _unresolved_row(target.target_id, "overexpression", "OE_reaction", rid, "unresolved_reaction", reaction_id=rid,
                        baseline_objective_value=oe_reaction.baseline_objective_value)
        for rid in screen_plan["unresolved_oe_reaction_ids"]))
    screens_list = [ko, ko_reaction, oe_gene, oe_reaction]
    report = write_simulation_outputs(
        simulation,
        tradeoff,
        tuple(screens_list),
        output_dir=output_dir,
        output_prefix=f"{target.target_id}_draft",
    )

    alignment_request = _alignment_request_for_target(target.target_id, request.compatibility_mode, root)
    alignment = build_alignment_summary(
        alignment_request["target_id"],
        python_result_status="corrected_condition" if request.compatibility_mode == "corrected" else simulation.result_status,
        artifact_path=alignment_request["artifact_path"],
        rows_python=(
            constraint_result.eq_shape[0] + constraint_result.ub_shape[0]
            if constraint_result.eq_shape and constraint_result.ub_shape
            else None
        ),
        cols_python=constraint_result.eq_shape[1] if constraint_result.eq_shape else None,
        objective_python=simulation.objective_value,
        constraint_diff_status=alignment_request["constraint_diff_status"],
        compatibility_exceptions=alignment_request["compatibility_exceptions"],
    )
    alignment_payload = summarize_alignment(alignment)
    alignment_payload["python_target_id"] = target.target_id
    alignment_payload["alignment_artifact_target_id"] = alignment_request["target_id"]
    summary_payload = _attach_pipeline_metadata(
        report.summary_path,
        {
            "secretion_plan": {
                "supported": plan.supported,
                "route_kind": plan.route_kind,
                "reaction_count": plan.reaction_count,
                "ptm_counts": plan.ptm_counts,
            },
            "protein_cost_analysis": summarize_protein_cost_analysis(protein_cost),
            "alignment_summary": alignment_payload,
            "target_metadata": _target_metadata(target, request),
            "target_warnings": _target_warnings(target),
            "compatibility_mode": request.compatibility_mode,
            "glycosylation_mode": request.glycosylation_mode,
            "screen_warnings": screen_plan["warnings"],
            "screen_request": {
                "ko_gene_ids": list(screen_plan["requested_ko_gene_ids"]),
                "ko_reaction_ids": list(screen_plan["requested_ko_rxn_ids"]),
                "oe_gene_ids": list(screen_plan["requested_oe_gene_ids"]),
                "oe_reaction_ids": list(screen_plan["requested_oe_reaction_ids"]),
                "screen_candidate_limit": screen_plan["candidate_limit"],
            },
        },
    )
    report.report_path.write_text(_build_pipeline_report(summary_payload), encoding="utf-8")

    return PichiaSimulationRunResult(
        success=simulation.success and report.summary_path.exists(),
        target_id=target.target_id,
        candidate_id=request.candidate_id,
        mu=request.mu,
        production_ratio=simulation.objective_value,
        media_type=request.media_type,
        message=simulation.message,
        output_file=report.report_path,
        objective_value=str(simulation.objective_value) if simulation.objective_value is not None else None,
        result_status="corrected_condition" if request.compatibility_mode == "corrected" else simulation.result_status,
        summary_path=report.summary_path,
        report_path=report.report_path,
        matlab_alignment_status=alignment.matlab_alignment_status,
        constraint_counts=simulation.constraint_counts,
        candidate_table_path=report.candidate_table_path,
        tradeoff_path=report.tradeoff_path,
        alignment_summary=alignment_payload,
    )


def _build_screen_plan(model: CobraModel, request: PichiaSimulationRequest) -> dict[str, Any]:
    limit = _candidate_limit(request.screen_candidate_limit)
    requested_ko = _dedupe_nonempty(request.ko_gene_ids)[:limit]
    requested_ko_rxns = _dedupe_nonempty(request.ko_reaction_ids)[:limit]
    requested_oe_genes = _dedupe_nonempty(request.oe_gene_ids)[:limit]
    requested_oe_reactions = _dedupe_nonempty(request.oe_reaction_ids)[:limit]
    has_user_ko = bool(requested_ko or requested_ko_rxns)
    has_user_oe = bool(requested_oe_genes or requested_oe_reactions)

    ko_gene_ids, unresolved_ko = split_existing_genes(model, tuple(requested_ko))
    if not has_user_ko and limit > 0:
        ko_gene_ids = default_ko_genes(model, 1)

    # Reaction-level KO IDs (for complex KO, no gene fallback needed)
    existing_ko_rxns, unresolved_ko_rxns = split_existing_reactions(model, requested_ko_rxns)

    resolved_oe_from_genes, oe_gene_by_reaction, unresolved_oe_genes, gene_warnings = resolve_oe_gene_reactions(
        model, requested_oe_genes, limit
    )
    existing_oe_reactions, unresolved_oe_reactions = split_existing_reactions(model, requested_oe_reactions)
    if not has_user_oe and limit > 0:
        existing_oe_reactions = default_oe_reactions(model, 1)

    warnings = list(gene_warnings)
    warnings.extend(f"敲除基因未在模型中找到：{gene_id}" for gene_id in unresolved_ko)
    warnings.extend(f"敲除反应未在模型中找到：{reaction_id}" for reaction_id in unresolved_ko_rxns)
    warnings.extend(f"过表达基因无法解析到模型反应：{gene_id}" for gene_id in unresolved_oe_genes)
    warnings.extend(f"过表达反应未在模型中找到：{reaction_id}" for reaction_id in unresolved_oe_reactions)
    if requested_oe_genes:
        warnings.append("过表达基因当前按 reaction-level OE proxy 运行：先解析到该基因参与的反应，再逐个模拟反应容量上调。")

    return {
        "candidate_limit": limit,
        "requested_ko_gene_ids": tuple(requested_ko),
        "requested_ko_rxn_ids": tuple(requested_ko_rxns),
        "requested_oe_gene_ids": tuple(requested_oe_genes),
        "requested_oe_reaction_ids": tuple(requested_oe_reactions),
        "ko_gene_ids": ko_gene_ids[:limit],
        "unresolved_ko_gene_ids": unresolved_ko,
        "ko_reaction_ids": existing_ko_rxns[:limit],
        "unresolved_ko_reaction_ids": unresolved_ko_rxns,
        "oe_gene_reaction_ids": resolved_oe_from_genes[:limit],
        "oe_gene_by_reaction": oe_gene_by_reaction,
        "unresolved_oe_gene_ids": unresolved_oe_genes,
        "oe_reaction_ids": existing_oe_reactions[:limit],
        "unresolved_oe_reaction_ids": unresolved_oe_reactions,
        "warnings": warnings,
    }


def _append_unresolved_rows(result: ScreenResult, rows: tuple[dict[str, Any], ...]) -> ScreenResult:
    if not rows:
        return result
    merged = (*result.rows, *rows)
    return ScreenResult(
        target_id=result.target_id,
        screen_type=result.screen_type,
        success=bool(merged) and all(bool(row.get("success")) for row in merged),
        candidate_count=len(merged),
        rows=merged,
        constraint_counts=result.constraint_counts,
        baseline_objective_value=result.baseline_objective_value,
        result_status=result.result_status,
        matlab_alignment_status=result.matlab_alignment_status,
    )


def _unresolved_row(
    target_id: str,
    screen_type: str,
    intervention_type: str,
    candidate_id: str,
    status: str,
    baseline_objective_value: float | None,
    input_gene_id: str | None = None,
    reaction_id: str | None = None,
) -> dict[str, Any]:
    return {
        "target_id": target_id,
        "screen_type": screen_type,
        "candidate_id": candidate_id,
        "gene_id": input_gene_id if intervention_type == "KO" else None,
        "reaction_id": reaction_id,
        "input_gene_id": input_gene_id,
        "resolved_reaction_id": None,
        "intervention_type": intervention_type,
        "success": False,
        "status": status,
        "objective_value": None,
        "baseline_objective_value": baseline_objective_value,
        "delta_objective": None,
        "effect_label": "未解析",
        "secretory_process": "未解析",
        "mapping_level": "unresolved",
        "mapping_confidence": "unresolved",
        "mapping_interpretation": "未解析到可解释的模型反应。",
        "complex_id": None,
        "complex_subunit_ids": [],
        "complex_subunit_stoichiometry": [],
    }


def _candidate_limit(value: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return 20
    return max(0, min(parsed, 20))


def _dedupe_nonempty(values: tuple[str, ...]) -> list[str]:
    deduped: list[str] = []
    for value in values:
        item = str(value).strip()
        if item and item not in deduped:
            deduped.append(item)
    return deduped


def _growth_points(values: tuple[float, ...], fallback_mu: float) -> tuple[float, ...]:
    points = tuple(float(value) for value in values if float(value) > 0)
    return points or (fallback_mu,)


def _resolve_target(request: PichiaSimulationRequest, root: Path) -> TargetSpec:
    if request.target_input is not None and request.leader_candidate is not None:
        return target_spec_from_input(request.target_input, request.leader_candidate)
    if isinstance(request.target_input, (str, Path)):
        return _select_target(load_custom_targets_json(Path(request.target_input)), request.target_id)
    if isinstance(request.target_input, dict):
        if "targets" in request.target_input:
            targets = [target_spec_from_mapping(item, source="request.target_input") for item in request.target_input["targets"]]
            return _select_target(targets, request.target_id)
        return target_spec_from_mapping(request.target_input, source="request.target_input")
    return _select_target(load_builtin_targets(root), request.target_id)


def _select_target(targets: list[TargetSpec], target_id: str) -> TargetSpec:
    for target in targets:
        if target.target_id == target_id or target.protein_id == target_id:
            return target
    available = ", ".join(target.target_id for target in targets)
    raise ValueError(f"Target {target_id!r} was not found. Available targets: {available}")


def _attach_pipeline_metadata(summary_path: Path, metadata: dict[str, Any]) -> dict[str, Any]:
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    payload.update(metadata)
    summary_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return payload


def _uses_opn_known_alignment_exceptions(target_id: str, compatibility_mode: str) -> bool:
    return compatibility_mode == "corrected" and target_id == "OPN_ALPHA_FULL_PROJECT"


def _uses_hlf_project_710_alignment_exceptions(target_id: str, compatibility_mode: str) -> bool:
    return compatibility_mode == "corrected" and target_id == "hLF"


def _alignment_request_for_target(target_id: str, compatibility_mode: str, root: Path) -> dict[str, Any]:
    if _uses_opn_known_alignment_exceptions(target_id, compatibility_mode):
        return {
            "target_id": target_id,
            "artifact_path": root / DEFAULT_ALIGNMENT_ARTIFACT,
            "constraint_diff_status": "known_matlab_compatibility_differences",
            "compatibility_exceptions": opn_known_matlab_compatibility_exceptions(),
        }
    if _uses_hlf_project_710_alignment_exceptions(target_id, compatibility_mode):
        return {
            "target_id": "hLF_PROJECT_710",
            "artifact_path": root / HLF_PROJECT_710_ALIGNMENT_ARTIFACT,
            "constraint_diff_status": "known_matlab_compatibility_differences",
            "compatibility_exceptions": hlf_project_710_known_matlab_compatibility_exceptions(),
        }
    return {
        "target_id": target_id,
        "artifact_path": root / DEFAULT_ALIGNMENT_ARTIFACT,
        "constraint_diff_status": None,
        "compatibility_exceptions": None,
    }


def _target_metadata(target: TargetSpec, request: PichiaSimulationRequest | None = None) -> dict[str, Any]:
    metadata = {
        "target_id": target.target_id,
        "protein_id": target.protein_id,
        "source": target.source,
        "sequence_role": _target_sequence_role(target),
        "normalization_mode": _target_normalization_mode(target),
        "alignment_target_kind": _alignment_target_kind(target),
        "full_sequence_length": len(target.full_sequence),
        "mature_sequence_length": len(target.mature_sequence),
        "leader_sequence_length": len(target.leader_sequence),
        "signal_peptide_length": len(target.signal_peptide_sequence),
        "contains_leader": bool(target.leader_sequence),
        "contains_signal_peptide": bool(target.signal_peptide_sequence),
        "disulfide_sites": target.disulfide_sites,
        "n_glycosylation_sites": target.n_glycosylation_sites,
        "o_glycosylation_sites": target.o_glycosylation_sites,
        "target_parameter_status": _target_parameter_status_for_report(target),
    }
    metadata.update(_request_sequence_contract_metadata(target, request))
    return metadata


def _request_sequence_contract_metadata(target: TargetSpec, request: PichiaSimulationRequest | None) -> dict[str, Any]:
    if request is None:
        return {}
    has_contract = (
        request.target_input is not None
        or request.sequence_role != "unknown"
        or request.normalization_mode != "as_provided"
        or request.terminal_stop_policy != "allow_for_record_only"
        or request.contains_signal_peptide is not None
        or request.contains_leader is not None
        or request.original_sequence_length is not None
    )
    if not has_contract:
        return {}
    return {
        "sequence_role": (
            _target_sequence_role(target) if request.sequence_role in {"", "unknown"} else request.sequence_role
        ),
        "normalization_mode": request.normalization_mode or _target_normalization_mode(target),
        "contains_signal_peptide": (
            bool(target.signal_peptide_sequence)
            if request.contains_signal_peptide is None
            else bool(request.contains_signal_peptide)
        ),
        "contains_leader": bool(target.leader_sequence) if request.contains_leader is None else bool(request.contains_leader),
        "terminal_stop_policy": request.terminal_stop_policy,
        "original_sequence_length": request.original_sequence_length,
        "normalized_sequence_length": request.normalized_sequence_length or len(target.mature_sequence),
        "original_full_sequence_length": request.original_full_sequence_length,
        "normalized_full_sequence_length": request.normalized_full_sequence_length or len(target.full_sequence),
        "original_leader_sequence_length": request.original_leader_sequence_length,
        "normalized_leader_sequence_length": request.normalized_leader_sequence_length or len(target.leader_sequence),
        "original_signal_peptide_length": request.original_signal_peptide_length,
        "normalized_signal_peptide_length": request.normalized_signal_peptide_length or len(target.signal_peptide_sequence),
        "terminal_stop_present": request.terminal_stop_present,
        "terminal_stop_removed": request.terminal_stop_removed,
    }


def _target_sequence_role(target: TargetSpec) -> str:
    if target.target_id == "hLF":
        return "native_signal_plus_mature_hLF"
    if target.target_id.startswith("OPN_"):
        return "mature_secreted_with_leader_candidate"
    return "custom_user_sequence"


def _target_normalization_mode(target: TargetSpec) -> str:
    if target.target_id == "hLF":
        return "user_provided_as_provided"
    return "as_provided"


def _alignment_target_kind(target: TargetSpec) -> str:
    if target.target_id == "hLF":
        return "project_defined_hLF"
    if target.target_id.startswith("OPN_"):
        return "opn_leader_candidate"
    return "custom"


def _target_parameter_status_for_report(target: TargetSpec) -> str:
    if target.target_id == "hLF":
        return "draft_matlab_alignment_pending"
    return "draft"


def _target_warnings(target: TargetSpec) -> list[str]:
    if target.target_id == "hLF":
        return [
            "hLF 使用用户提供的 710aa 目标序列：人源天然信号肽 19aa + mature hLF 691aa。",
            "当前 Python target `hLF` 对应 MATLAB artifact target `hLF_PROJECT_710`，可报告为 `aligned_except_known_matlab_compatibility_differences`，但不是 fully aligned。",
            "旧 MATLAB hLF baseline 仍保持 matlab_failed；Python corrected 结果不能声明为旧 MATLAB fully aligned。",
        ]
    return []


def _build_pipeline_report(summary: dict[str, Any]) -> str:
    alignment = summary.get("alignment_summary") or {}
    exceptions = alignment.get("compatibility_exceptions") or []
    target_metadata = summary.get("target_metadata") or {}
    target_warnings = summary.get("target_warnings") or []
    protein_cost = summary.get("protein_cost_analysis") or {}
    lines = [
        f"# pcSecPichia Python 分泌仿真报告: {summary.get('target_id')}",
        "",
        "## 状态",
        "",
        "- 当前结果是 Python 草稿结果。",
        f"- 结果状态: `{summary.get('result_status')}`.",
        f"- MATLAB 对齐状态: `{summary.get('matlab_alignment_status')}`.",
        f"- 目标参数状态: `{summary.get('target_parameter_status')}`.",
    ]
    if alignment.get("alignment_artifact_target_id"):
        lines.append(f"- MATLAB artifact target: `{alignment.get('alignment_artifact_target_id')}`.")
    if target_metadata:
        lines.extend(
            [
                "",
                "## 目标蛋白元数据",
                "",
                f"- alignment target kind: `{target_metadata.get('alignment_target_kind')}`.",
                f"- sequence role: `{target_metadata.get('sequence_role')}`.",
                f"- normalization mode: `{target_metadata.get('normalization_mode')}`.",
                f"- full sequence length: `{target_metadata.get('full_sequence_length')}` aa.",
                f"- mature sequence length: `{target_metadata.get('mature_sequence_length')}` aa.",
                f"- leader sequence length: `{target_metadata.get('leader_sequence_length')}` aa.",
                f"- signal peptide length: `{target_metadata.get('signal_peptide_length')}` aa.",
                f"- DSB/NG/OG: `{target_metadata.get('disulfide_sites')}` / `{target_metadata.get('n_glycosylation_sites')}` / `{target_metadata.get('o_glycosylation_sites')}`.",
            ]
        )
        if "terminal_stop_policy" in target_metadata:
            lines.extend(
                [
                    f"- terminal stop policy: `{target_metadata.get('terminal_stop_policy')}`.",
                    f"- raw/normalized mature length: `{target_metadata.get('original_sequence_length')}` / `{target_metadata.get('normalized_sequence_length')}` aa.",
                    f"- raw/normalized full length: `{target_metadata.get('original_full_sequence_length')}` / `{target_metadata.get('normalized_full_sequence_length')}` aa.",
                    f"- terminal stop present/removed: `{target_metadata.get('terminal_stop_present')}` / `{target_metadata.get('terminal_stop_removed')}`.",
                ]
            )
    if target_warnings:
        lines.extend(["", "## 目标输入边界", "", *[f"- {item}" for item in target_warnings]])
    if protein_cost:
        lines.extend(_protein_cost_report_lines(protein_cost))
    screen_warnings = summary.get("screen_warnings") or []
    if screen_warnings:
        lines.extend(["", "## 基因扰动提示", "", *[f"- {item}" for item in screen_warnings]])
    if exceptions:
        lines.extend(
            [
                "- 修正条件下的 Python 结果仍标记为 `corrected_condition`，不等同于旧 MATLAB 基线完全对齐。",
                "- 已知 MATLAB 兼容差异：",
                *[f"  - `{item.get('id')}`: `{item.get('count')}` {item.get('category')}" for item in exceptions],
            ]
        )
    lines.extend(
        [
            "",
            "## 仿真结果",
            "",
            f"- 是否成功: `{summary.get('success')}`.",
            f"- 目标函数值: `{summary.get('objective_value')}`.",
            f"- 生长速率: `{summary.get('growth_rate')}`.",
            f"- 约束计数: `{summary.get('constraint_counts') or {}}`.",
            "",
            "## 输出",
            "",
            f"- 候选行数: `{summary.get('candidate_count')}`.",
            f"- 生长权衡行数: `{len((summary.get('tradeoff') or {}).get('tradeoff_rows') or [])}`.",
            "",
        ]
    )
    candidate_interpretation = summary.get("candidate_interpretation") or {}
    if candidate_interpretation:
        lines.extend(
            [
                "## 候选解释",
                "",
                f"- 分类汇总: `{candidate_interpretation.get('effect_counts') or {}}`.",
                "- `约束不可行` 表示在当前固定生长速率和约束组合下，该扰动没有可行解；这不等同于真实生物系统必然不可行。",
                "",
            ]
        )
        for row in candidate_interpretation.get("top_explanations") or []:
            lines.append(f"- {row.get('summary')}")
        lines.append("")
    return "\n".join(lines)


def _protein_cost_report_lines(protein_cost: dict[str, Any]) -> list[str]:
    items = protein_cost.get("cost_items") or []
    dominant = protein_cost.get("dominant_cost_categories") or []
    lines = [
        "",
        "## 目标蛋白成本分析",
        "",
        "- 当前结果是 Python draft explanatory score，不代表真实发酵产量或湿实验成本。",
        f"- 成本分析状态: `{protein_cost.get('result_status')}`.",
        f"- 总相对成本分: `{protein_cost.get('total_relative_score')}`.",
        f"- 主要成本类别: `{', '.join(str(item) for item in dominant)}`.",
        "",
        "| 类别 | 成本项 | 相对分 | 依据 |",
        "| --- | --- | ---: | --- |",
    ]
    for item in items:
        lines.append(
            f"| `{item.get('category')}` | {item.get('label')} | "
            f"{item.get('relative_score')} | {item.get('basis')} |"
        )
    warnings = protein_cost.get("warnings") or []
    if warnings:
        lines.extend(["", "提示:"])
        lines.extend(f"- {warning}" for warning in warnings)
    return lines


__all__ = ["run_pichia_secretion_simulation"]
