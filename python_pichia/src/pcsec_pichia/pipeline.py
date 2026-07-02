from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pcsec_pichia.analysis import (
    analyze_target_growth_impact,
    analyze_target_protein_cost,
    analyze_target_protein_lp_attribution,
    analyze_yield_improvement_candidates,
    summarize_protein_cost_analysis,
    summarize_protein_cost_slope_compatibility,
    summarize_protein_lp_attribution,
    summarize_target_growth_analysis,
    summarize_yield_improvement_recommendations,
)
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
from pcsec_pichia.screens.planning import build_screen_plan
from pcsec_pichia.services.gene_evidence import (
    DEFAULT_GENE_EVIDENCE_CACHE,
    evidence_for_gene,
    load_gene_evidence_cache,
    recommendation_tier_for_candidate,
)
from pcsec_pichia.services.gene_catalog import build_lightweight_secretion_gene_evidence
from pcsec_pichia.simulation import (
    run_growth_tradeoff,
    run_protein_cost_slope_compatibility,
    solve_secretion_capacity,
)
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
    lp_attribution = analyze_target_protein_lp_attribution(
        target,
        plan,
        constraint_result.constraint_counts,
        simulation,
        reaction_ids=tuple(inputs.prepared_model.rxns),
    )
    cost_slope_compatibility = None
    if request.enable_cost_slope_compatibility:
        cost_slope_ratios, cost_slope_policy = _cost_slope_secretion_ratio_policy(request, simulation)
        cost_slope_compatibility = run_protein_cost_slope_compatibility(
            inputs.prepared_model,
            target,
            inputs.amino_acids,
            inputs.metabolic,
            inputs.secretory,
            inputs.combined,
            growth_rates=request.cost_slope_growth_rates,
            secretion_ratios=cost_slope_ratios,
            write_ribosome_translation_constraint=request.enable_ribosome_translation_constraint,
            write_misfolding_constraints=request.enable_misfolding_constraint,
            medium_compatibility_mode=request.cost_slope_medium_compatibility_mode,
        )
        cost_slope_compatibility = _with_cost_slope_ratio_policy(
            cost_slope_compatibility,
            cost_slope_policy,
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
    growth_analysis = analyze_target_growth_impact(tradeoff, baseline_growth_rate=request.mu)

    evidence_by_gene = load_gene_evidence_cache(root / DEFAULT_GENE_EVIDENCE_CACHE)
    screen_plan = build_screen_plan(
        inputs.prepared_model,
        request,
        inputs.secretory.complex_subunits,
        evidence_by_gene=evidence_by_gene,
    )

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
    screens_list = [_annotate_gene_evidence(result) for result in screens_list]
    candidate_rows = tuple(row for result in screens_list for row in result.rows)
    yield_recommendations = analyze_yield_improvement_candidates(
        candidate_rows,
        target_id=target.target_id,
        baseline_objective=simulation.objective_value,
        baseline_secretion_flux=simulation.secretion_flux,
        curated_gene_evidence=tuple(build_lightweight_secretion_gene_evidence(model)),
    )
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
    protein_cost_payload = summarize_protein_cost_analysis(protein_cost)
    protein_cost_payload["lp_attribution"] = summarize_protein_lp_attribution(lp_attribution)
    protein_cost_payload["cost_slope_compatibility"] = summarize_protein_cost_slope_compatibility(
        cost_slope_compatibility
    )
    summary_payload = _attach_pipeline_metadata(
        report.summary_path,
        {
            "secretion_plan": {
                "supported": plan.supported,
                "route_kind": plan.route_kind,
                "reaction_count": plan.reaction_count,
                "ptm_counts": plan.ptm_counts,
            },
            "protein_cost_analysis": protein_cost_payload,
            "target_growth_analysis": summarize_target_growth_analysis(growth_analysis),
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


def _annotate_gene_evidence(
    result: ScreenResult,
    evidence_by_gene: dict[str, Any] | None = None,
) -> ScreenResult:
    evidence_by_gene = load_gene_evidence_cache() if evidence_by_gene is None else evidence_by_gene
    rows: list[dict[str, Any]] = []
    for row in result.rows:
        gene_id = str(row.get("canonical_gene_id") or row.get("gene_id") or row.get("input_gene_id") or "").strip()
        record = evidence_for_gene(gene_id, evidence_by_gene)
        if record is None:
            rows.append(
                {
                    **row,
                    "standard_gene_symbol": "",
                    "display_name": row.get("input_gene_id") or row.get("gene_id") or row.get("candidate_id"),
                    "protein_name": "",
                    "function_annotation": "",
                    "external_ids": {},
                    "ec_numbers": [],
                    "go_terms": [],
                    "ortholog_symbol": "",
                    "evidence_sources": [],
                    "evidence_confidence": "low_model_only",
                    "wet_lab_readiness": "model_only_not_experiment_ready",
                    **_candidate_evidence_tier_fields(row, None, result.target_id),
                }
            )
            continue
        rows.append(
            {
                **row,
                "standard_gene_symbol": record.standard_gene_symbol,
                "display_name": record.display_name,
                "protein_name": record.protein_name,
                "function_annotation": record.function_annotation,
                "external_ids": record.external_ids or {},
                "ec_numbers": list(record.ec_numbers),
                "go_terms": list(record.go_terms),
                "ortholog_symbol": record.ortholog_symbol,
                "evidence_sources": list(record.evidence_sources),
                "evidence_confidence": record.evidence_confidence,
                "wet_lab_readiness": record.wet_lab_readiness,
                **_candidate_evidence_tier_fields(row, record, result.target_id),
            }
        )
    return ScreenResult(
        target_id=result.target_id,
        screen_type=result.screen_type,
        success=result.success,
        candidate_count=result.candidate_count,
        rows=tuple(rows),
        constraint_counts=result.constraint_counts,
        baseline_objective_value=result.baseline_objective_value,
        result_status=result.result_status,
        matlab_alignment_status=result.matlab_alignment_status,
    )


def _candidate_evidence_tier_fields(
    row: dict[str, Any],
    record: Any | None,
    target_id: str,
) -> dict[str, Any]:
    intervention_type = str(row.get("intervention_type") or row.get("screen_type") or "")
    gene_id = str(row.get("canonical_gene_id") or row.get("gene_id") or row.get("input_gene_id") or "")
    status = str(row.get("status") or "")
    resolved = not status.startswith("unresolved")
    model_gpr_executable = (
        intervention_type == "KO"
        and str(row.get("ko_support_status") or "") == "ko_runnable_gpr_gene_deletion"
    )
    oe_reaction_proxy = (
        intervention_type in {"OE_gene_proxy", "OE_reaction"}
        and resolved
        and status
        not in {"not_run_complex_subunit_limited", "not_run_gene_oe_proxy", "not_run_no_gpr_effect"}
    )
    aliases = tuple(str(row.get(key) or "") for key in ("input_gene_id", "display_name", "standard_gene_symbol"))
    tier, reason, phenotype = recommendation_tier_for_candidate(
        gene_id=gene_id,
        intervention_type=intervention_type,
        target_protein_context=target_id,
        model_gpr_executable=model_gpr_executable,
        oe_reaction_proxy=oe_reaction_proxy,
        resolved=resolved,
        database_annotation_available=record is not None,
        aliases=aliases,
    )
    return {
        "database_annotation_sources": list(record.evidence_sources) if record else [],
        "database_annotation_confidence": record.evidence_confidence if record else "",
        "model_gpr_executable": model_gpr_executable,
        "oe_reaction_proxy": oe_reaction_proxy,
        "phenotype_evidence": phenotype.to_dict() if phenotype else {},
        "recommendation_tier": tier,
        "recommendation_tier_reason": reason,
    }

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


def _cost_slope_secretion_ratio_policy(
    request: PichiaSimulationRequest,
    simulation: Any,
) -> tuple[tuple[float, ...], dict[str, Any]]:
    explicit = tuple(float(value) for value in request.cost_slope_secretion_ratios if float(value) > 0)
    if explicit:
        return explicit, {
            "secretion_ratio_policy": "explicit_absolute_ratios",
            "capacity_reference": None,
            "capacity_fractions": (),
            "warnings": (),
        }
    capacity = simulation.secretion_flux if getattr(simulation, "success", False) else None
    fractions = tuple(float(value) for value in request.cost_slope_capacity_fractions if 0 < float(value) < 1)
    if capacity is not None and capacity > 0 and fractions:
        ratios = tuple(float(capacity) * fraction for fraction in fractions)
        return ratios, {
            "secretion_ratio_policy": "capacity_fraction_ratios",
            "capacity_reference": float(capacity),
            "capacity_fractions": fractions,
            "warnings": (
                "No explicit target secretion ratios were provided; cost slope ratios were generated from current secretion capacity fractions.",
            ),
        }
    fallback = (5e-7, 1e-6, 5e-6, 1e-5, 2e-5)
    return fallback, {
        "secretion_ratio_policy": "fallback_absolute_ratios_capacity_unavailable",
        "capacity_reference": None if capacity is None else float(capacity),
        "capacity_fractions": fractions,
        "warnings": (
            "Secretion capacity was unavailable, so historical absolute protein-cost ratios were used as fallback.",
        ),
    }


def _cost_slope_ratio_policy_lines(cost_slope: dict[str, Any]) -> list[str]:
    policy = str(cost_slope.get("secretion_ratio_policy") or "explicit_absolute_ratios")
    capacity = cost_slope.get("capacity_reference")
    fractions = tuple(cost_slope.get("capacity_fractions") or ())
    lines = [f"- secretion ratio policy: `{policy}`."]
    if policy == "capacity_fraction_ratios":
        fraction_text = ", ".join(f"{float(value):.0%}" for value in fractions)
        lines.append(
            "- target secretion ratios: generated from current corrected secretion capacity "
            f"`{capacity}` using capacity fractions `{fraction_text}`."
        )
        lines.append(
            "- interpretation: this is the default when experimental target secretion ratios are unknown; "
            "explicit user ratios should override it when available."
        )
    elif policy == "explicit_absolute_ratios":
        lines.append(
            "- target secretion ratios: explicit absolute ratios supplied by the request; "
            "these are treated as the historical MATLAB-style fixed secretion requirements."
        )
    else:
        lines.append(
            "- target secretion ratios: fallback historical absolute ratios because current secretion capacity "
            "was unavailable; treat this as a diagnostic fallback, not calibrated biology."
        )
    return lines


def _with_cost_slope_ratio_policy(result: Any, policy: dict[str, Any]) -> Any:
    return type(result)(
        **{
            **result.__dict__,
            "warnings": (*result.warnings, *tuple(policy.get("warnings") or ())),
            "secretion_ratio_policy": str(policy["secretion_ratio_policy"]),
            "capacity_reference": policy.get("capacity_reference"),
            "capacity_fractions": tuple(policy.get("capacity_fractions") or ()),
        }
    )


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
    target_growth = summary.get("target_growth_analysis") or {}
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
    if target_growth:
        lines.extend(_target_growth_report_lines(target_growth))
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
    lp_attribution = protein_cost.get("lp_attribution") or {}
    if lp_attribution:
        lines.extend(_lp_attribution_report_lines(lp_attribution))
    cost_slope = protein_cost.get("cost_slope_compatibility") or {}
    if cost_slope:
        lines.extend(_cost_slope_report_lines(cost_slope))
    return lines


def _lp_attribution_report_lines(lp_attribution: dict[str, Any]) -> list[str]:
    objective = lp_attribution.get("objective_evidence") or {}
    lines = [
        "",
        "### LP 级归因证据",
        "",
        "- 当前结果是 Python draft LP sensitivity，基于 SciPy HiGHS marginals；不是 MATLAB/SoPlex fully aligned shadow price。",
        f"- LP 归因状态: `{lp_attribution.get('result_status')}`.",
        f"- 目标反应: `{objective.get('objective_reaction')}`.",
        f"- 分泌通量: `{objective.get('secretion_flux')}`.",
        "",
    ]
    lines.extend(_markdown_table("主导约束块", lp_attribution.get("dominant_constraint_blocks") or ()))
    lines.extend(_markdown_table("Top constraint marginals", lp_attribution.get("top_constraint_marginals") or ()))
    lines.extend(_markdown_table("Top bound marginals", lp_attribution.get("top_bound_marginals") or ()))
    lines.extend(_markdown_table("目标相关 flux", lp_attribution.get("target_related_fluxes") or ()))
    warnings = lp_attribution.get("warnings") or []
    if warnings:
        lines.extend(["", "LP 归因提示:"])
        lines.extend(f"- {warning}" for warning in warnings)
    return lines


def _cost_slope_report_lines(cost_slope: dict[str, Any]) -> list[str]:
    lines = [
        "",
        "### MATLAB-compatible 成本 slope（可选）",
        "",
        f"- 开启状态: `{cost_slope.get('enabled')}`.",
        f"- 结果状态: `{cost_slope.get('result_status')}`.",
        "- 当前默认路线: 固定生长率、corrected medium、最大化目标蛋白分泌通量。",
        "- 历史成本路线: 固定目标蛋白 exchange ratio 和生长率，优化 `Ex_glc_D`，再估计 glucose/ribosome cost slope。",
        "- 该模式用于历史 MATLAB `Protein_cost_TP` 定义对比，不替代默认 Python corrected pipeline。",
    ]
    lines.extend(_cost_slope_ratio_policy_lines(cost_slope))
    lines.append(f"- medium compatibility mode: `{cost_slope.get('medium_compatibility_mode', 'corrected')}`.")
    overrides = cost_slope.get("medium_bound_overrides") or []
    if overrides:
        lines.extend(_markdown_table("medium bound overrides", overrides, limit=12))
    lines.extend(_markdown_table("glucose cost slopes", cost_slope.get("glucose_cost_slopes") or ()))
    lines.extend(_markdown_table("ribosome cost slopes", cost_slope.get("ribosome_cost_slopes") or ()))
    lines.extend(_markdown_table("cost slope rows", cost_slope.get("rows") or (), limit=10))
    warnings = cost_slope.get("warnings") or []
    if warnings:
        lines.extend(["", "cost slope 提示:"])
        lines.extend(f"- {warning}" for warning in warnings)
    return lines


def _markdown_table(title: str, rows: Any, limit: int = 8) -> list[str]:
    rows = list(rows or [])[:limit]
    if not rows:
        return []
    keys = list(rows[0].keys())[:6]
    lines = ["", f"#### {title}", "", "| " + " | ".join(keys) + " |", "| " + " | ".join("---" for _ in keys) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(key, "")) for key in keys) + " |")
    return lines


def _target_growth_report_lines(target_growth: dict[str, Any]) -> list[str]:
    items = target_growth.get("tradeoff_points") or []
    best_flux = target_growth.get("best_secretion_point") or {}
    best_per_biomass = target_growth.get("best_secretion_per_biomass_point") or {}
    lines = [
        "",
        "## 目标蛋白生长分析",
        "",
        "- 当前结果是 Python draft explanatory tradeoff，不代表真实发酵生长预测。",
        f"- 生长分析状态: `{target_growth.get('result_status')}`.",
        f"- 趋势标签: `{target_growth.get('growth_sensitivity_label')}`.",
        f"- 趋势原因: `{target_growth.get('growth_sensitivity_reason')}`.",
        f"- 可比较生长点数量: `{target_growth.get('valid_point_count')}`.",
        f"- 最高分泌通量生长点: `{best_flux.get('mu')}`.",
        f"- 最高单位生物量分泌生长点: `{best_per_biomass.get('mu')}`.",
        "",
        "| mu | success | secretion flux | secretion / biomass | interpretation |",
        "| ---: | --- | ---: | ---: | --- |",
    ]
    for item in items:
        lines.append(
            f"| {item.get('mu')} | {item.get('success')} | {item.get('secretion_flux')} | "
            f"{item.get('secretion_per_biomass')} | {item.get('interpretation')} |"
        )
    warnings = target_growth.get("warnings") or []
    if warnings:
        lines.extend(["", "提示:"])
        lines.extend(f"- {warning}" for warning in warnings)
    return lines


__all__ = ["run_pichia_secretion_simulation"]
