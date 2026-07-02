from __future__ import annotations

from pathlib import Path
from typing import Any

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.core.paths import ProjectPaths

from app.services.pichia_secretion_schema import SecretionRunRequest


def preview_screen_inputs(
    request: SecretionRunRequest,
    paths: ProjectPaths | None = None,
) -> dict[str, Any]:
    """Preview KO/OE input resolution without running pcSec simulation.

    This is a UI/API explanation helper. It only reads model indices and GPR
    rules, then reports which manual candidates are resolvable. It must not
    change engine behavior or perform optimization.
    """
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    ensure_python_pichia_on_path()

    from pcsec_pichia.loading import load_pcsec_pichia_inputs

    inputs = load_pcsec_pichia_inputs(resolved_paths.repo_root)
    model = inputs.prepared_model
    return _preview_screen_inputs_for_model(
        model,
        request,
        complex_subunits=inputs.secretory.complex_subunits,
        repo_root=resolved_paths.repo_root,
    )


def _preview_screen_inputs_for_model(
    model: Any,
    request: SecretionRunRequest,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    repo_root: Path | None = None,
) -> dict[str, Any]:
    overlay_summary: dict[str, Any] = {}
    overlay_warnings: list[str] = []
    gene_rule_evidence_by_name = {}
    gene_rule_overlay = None
    if request.enable_gene_rule_overlay:
        from pcsec_pichia.services.gene_rule_overlay import (
            DEFAULT_GENE_RULE_EVIDENCE_CACHE,
            apply_gpr_overlay_for_analysis,
            build_gpr_overlay,
            load_gene_rule_evidence_cache,
        )

        cache_path = (repo_root / DEFAULT_GENE_RULE_EVIDENCE_CACHE) if repo_root is not None else None
        gene_rule_evidence_by_name = load_gene_rule_evidence_cache(cache_path)
        gene_rule_overlay = build_gpr_overlay(model, gene_rule_evidence_by_name)
        model = apply_gpr_overlay_for_analysis(model, gene_rule_overlay)
        overlay_summary = gene_rule_overlay.to_dict()
        overlay_warnings.extend(str(item) for item in gene_rule_overlay.warnings)
        if not gene_rule_overlay.entries:
            overlay_warnings.append(
                "External GPR overlay was enabled, but no high-confidence executable rules were available."
            )
    from pcsec_pichia.services.gene_evidence import DEFAULT_GENE_EVIDENCE_CACHE, load_gene_evidence_cache
    from pcsec_pichia.screens.planning import build_screen_plan

    evidence_cache_path = (repo_root / DEFAULT_GENE_EVIDENCE_CACHE) if repo_root is not None else None
    evidence_by_gene = load_gene_evidence_cache(evidence_cache_path)

    screen_plan = build_screen_plan(
        model,
        request,
        complex_subunits=complex_subunits,
        gene_rule_evidence_by_name=gene_rule_evidence_by_name,
        gene_rule_overlay=gene_rule_overlay,
        evidence_by_gene=evidence_by_gene,
    )
    limit = screen_plan.candidate_limit
    ko_gene_ids = screen_plan.requested_ko_gene_ids
    ko_reaction_ids = screen_plan.requested_ko_rxn_ids
    oe_gene_ids = screen_plan.requested_oe_gene_ids
    oe_reaction_ids = screen_plan.requested_oe_reaction_ids

    from pcsec_pichia.screens import (
        build_gene_capability_profile,
        build_gene_perturbation_map,
    )

    target_context = str(getattr(request, "target_id", "") or "")
    ko_canonical_by_input = _canonical_by_input_from_plan(ko_gene_ids, screen_plan.ko_input_by_gene)
    oe_canonical_by_input = _canonical_by_input_from_plan(oe_gene_ids, screen_plan.oe_input_by_gene)
    ko_canonical_ids = _dedupe_text_tuple(tuple(ko_canonical_by_input.values()))
    oe_canonical_ids = _dedupe_text_tuple(tuple(oe_canonical_by_input.values()))

    gene_capabilities = {
        gene_id: build_gene_capability_profile(
            model,
            gene_id,
            complex_subunits=complex_subunits,
            evidence_by_gene=evidence_by_gene,
            target_protein_context=target_context,
        ).to_dict()
        for gene_id in _dedupe_text_tuple((*ko_gene_ids, *oe_gene_ids))
    }
    canonical_gene_capabilities = {
        canonical_gene_id: build_gene_capability_profile(
            model,
            canonical_gene_id,
            complex_subunits=complex_subunits,
            aliases=tuple(
                input_id
                for input_id, mapped_id in {**ko_canonical_by_input, **oe_canonical_by_input}.items()
                if mapped_id == canonical_gene_id and input_id != canonical_gene_id
            ),
            evidence_by_gene=evidence_by_gene,
            target_protein_context=target_context,
        ).to_dict()
        for canonical_gene_id in _dedupe_text_tuple((*ko_canonical_ids, *oe_canonical_ids))
    }

    ko_genes = [
        _gene_row(
            gene_id,
            "KO",
            resolved=ko_canonical_by_input[gene_id] not in screen_plan.unresolved_ko_gene_ids,
            capability=_capability_for_preview_row(
                gene_capabilities.get(gene_id),
                canonical_gene_capabilities.get(ko_canonical_by_input[gene_id]),
            ),
            canonical_gene_id=ko_canonical_by_input[gene_id],
        )
        for gene_id in ko_gene_ids
    ]
    ko_reactions = [
        _reaction_row(reaction_id, "KO_reaction", resolved=reaction_id in screen_plan.ko_reaction_ids)
        for reaction_id in ko_reaction_ids
    ]
    oe_genes = [
        _oe_gene_row_from_plan(
            screen_plan.oe_gene_plans_by_gene[oe_canonical_by_input[gene_id]],
            all_reactions=list(screen_plan.oe_gene_plans_by_gene[oe_canonical_by_input[gene_id]].affected_reactions),
            limit=limit,
            input_id=gene_id,
            capability=_capability_for_preview_row(
                gene_capabilities.get(gene_id),
                canonical_gene_capabilities.get(oe_canonical_by_input[gene_id]),
            ),
        )
        for gene_id in oe_gene_ids
    ]
    oe_reactions = [
        _reaction_row(reaction_id, "OE_reaction", resolved=reaction_id in screen_plan.oe_reaction_ids)
        for reaction_id in oe_reaction_ids
    ]
    gene_mapping = build_gene_perturbation_map(
        model,
        _dedupe_text_tuple((*ko_canonical_ids, *oe_canonical_ids)),
        complex_subunits=complex_subunits,
    ).to_dict()
    _attach_gene_mapping_input_ids(
        gene_mapping,
        _input_ids_by_canonical({**ko_canonical_by_input, **oe_canonical_by_input}),
    )
    warnings = _dedupe_text_list([
        *overlay_warnings,
        *screen_plan.warnings,
        *_screen_preview_warnings(ko_genes, ko_reactions, oe_genes, oe_reactions, request),
    ])

    return {
        "candidate_limit": limit,
        "semantics": {
            "KO": "按模型基因 ID 敲除；未在模型 gene_index 中找到的 ID 会写入 unresolved_gene。",
            "KO_reaction": "按模型反应 ID 直接敲除，主要用于复合体级扰动。",
            "OE_gene_proxy": "过表达基因会先做 GPR-aware 规划；只有单基因/同工酶等可解释场景才运行 reaction-level OE proxy。",
            "OE_reaction": "按模型反应 ID 直接做反应级过表达代理，作为高级诊断入口保留。",
        },
        "ko_genes": ko_genes,
        "ko_reactions": ko_reactions,
        "oe_genes": oe_genes,
        "oe_reactions": oe_reactions,
        "gene_mapping": gene_mapping,
        "gene_mapping_rows": gene_mapping["rows"],
        "gene_capabilities": list(gene_capabilities.values()),
        "gene_rule_overlay": overlay_summary,
        "warnings": warnings,
    }


def _gene_row(
    gene_id: str,
    intervention_type: str,
    resolved: bool,
    capability: dict[str, Any] | None = None,
    canonical_gene_id: str | None = None,
) -> dict[str, Any]:
    return {
        "input_id": gene_id,
        "canonical_gene_id": canonical_gene_id or (capability or {}).get("canonical_gene_id", gene_id),
        "intervention_type": intervention_type,
        "resolved": resolved,
        "status": "resolved" if resolved else "unresolved_gene",
        "resolved_reaction_count": None,
        "resolved_reactions_preview": [],
        "ko_support_status": (capability or {}).get("ko_support_status", ""),
        "oe_support_status": (capability or {}).get("oe_support_status", ""),
        "gpr_role": (capability or {}).get("gpr_role", ""),
        "support_reason": (capability or {}).get("support_reason", ""),
        "missing_information": list((capability or {}).get("missing_information") or []),
        "confidence": (capability or {}).get("confidence", ""),
        **_capability_evidence_preview_fields(capability, intervention_type),
    }


def _reaction_row(reaction_id: str, intervention_type: str, resolved: bool) -> dict[str, Any]:
    return {
        "input_id": reaction_id,
        "intervention_type": intervention_type,
        "resolved": resolved,
        "status": "resolved" if resolved else "unresolved_reaction",
        "resolved_reaction_count": 1 if resolved else 0,
        "resolved_reactions_preview": [reaction_id] if resolved else [],
    }


def _oe_gene_row_from_plan(
    plan: Any,
    all_reactions: list[str],
    limit: int,
    input_id: str | None = None,
    capability: dict[str, Any] | None = None,
) -> dict[str, Any]:
    executable = list(getattr(plan, "executable_reactions", ()))
    resolved = bool(getattr(plan, "resolved", False))
    status = "resolved" if resolved and executable else "not_run_complex_subunit_limited"
    if not resolved:
        status = "unresolved_gene"
    elif getattr(plan, "oe_support_status", "") == "oe_no_gpr_effect":
        status = "not_run_no_gpr_effect"
    elif getattr(plan, "oe_support_status", "") == "oe_explain_only_no_capacity_model":
        status = "not_run_gene_oe_proxy"
    return {
        "input_id": input_id or plan.gene_id,
        "canonical_gene_id": plan.gene_id,
        "intervention_type": "OE_gene_proxy",
        "resolved": resolved,
        "status": status,
        "resolved_reaction_count": len(executable),
        "resolved_reactions_preview": executable[: max(0, min(limit, 5))],
        "affected_reaction_count": len(all_reactions),
        "affected_reactions_preview": all_reactions[: max(0, min(limit, 5))],
        "gpr_role": plan.gpr_role,
        "capacity_effect": plan.capacity_effect,
        "simulation_basis": plan.simulation_basis,
        "mapping_confidence": plan.mapping_confidence,
        "ko_support_status": (capability or {}).get("ko_support_status", plan.ko_support_status),
        "oe_support_status": (capability or {}).get("oe_support_status", plan.oe_support_status),
        "support_reason": (capability or {}).get("support_reason", plan.support_reason),
        "missing_information": list((capability or {}).get("missing_information") or plan.missing_information),
        "warnings": list(plan.warnings),
        "truncated": len(executable) > limit,
        **_capability_evidence_preview_fields(capability, "OE"),
    }


def _capability_evidence_preview_fields(
    capability: dict[str, Any] | None,
    intervention_type: str,
) -> dict[str, Any]:
    data = capability or {}
    key = "KO" if str(intervention_type).upper().startswith("KO") else "OE"
    tier = data.get("recommendation_tier") if isinstance(data.get("recommendation_tier"), dict) else {}
    reason = (
        data.get("recommendation_tier_reason")
        if isinstance(data.get("recommendation_tier_reason"), dict)
        else {}
    )
    phenotype = data.get("phenotype_evidence") if isinstance(data.get("phenotype_evidence"), dict) else {}
    return {
        "database_annotation_sources": list(data.get("database_annotation_sources") or []),
        "database_annotation_confidence": data.get("database_annotation_confidence", ""),
        "model_gpr_executable": bool(data.get("model_gpr_executable")),
        "oe_reaction_proxy": bool(data.get("oe_reaction_proxy")),
        "phenotype_evidence": phenotype.get(key, {}),
        "recommendation_tier": tier.get(key, "manual_review_required"),
        "recommendation_tier_reason": reason.get(key, ""),
    }


def _screen_preview_warnings(
    ko_genes: list[dict[str, Any]],
    ko_reactions: list[dict[str, Any]],
    oe_genes: list[dict[str, Any]],
    oe_reactions: list[dict[str, Any]],
    request: SecretionRunRequest,
) -> list[str]:
    warnings: list[str] = []
    if not any((ko_genes, ko_reactions, oe_genes, oe_reactions)):
        warnings.append("未填写 KO/OE 候选；正式运行会使用小规模默认 smoke 候选。")
    for row in ko_genes:
        if not row["resolved"]:
            warnings.append(f"敲除基因未在模型中找到：{row['input_id']}")
    for row in ko_reactions:
        if not row["resolved"]:
            warnings.append(f"敲除反应未在模型中找到：{row['input_id']}")
    for row in oe_genes:
        if not row["resolved"]:
            warnings.append(f"过表达基因无法解析到模型反应：{row['input_id']}")
        elif row.get("status") == "not_run_no_gpr_effect":
            warnings.append(f"过表达基因 {row['input_id']} 只做解释：当前模型没有 GPR 反应引用它。")
        elif row.get("status") == "not_run_gene_oe_proxy":
            warnings.append(f"过表达基因 {row['input_id']} 只做解释：当前缺少可解释的容量增强模型。")
        elif row.get("status") == "not_run_complex_subunit_limited":
            warnings.append(f"过表达基因 {row['input_id']} 只做解释：复合体亚基或混合 GPR 不直接运行 capacity proxy。")
        elif row.get("truncated"):
            warnings.append(
                f"过表达基因 {row['input_id']} 解析到的反应超过候选上限 "
                f"{request.screen_candidate_limit}，正式运行会截断。"
            )
        for item in row.get("warnings") or []:
            warnings.append(str(item))
    for row in oe_reactions:
        if not row["resolved"]:
            warnings.append(f"过表达反应未在模型中找到：{row['input_id']}")
    if oe_genes:
        warnings.append("过表达基因当前是 GPR-aware planning + reaction-level proxy；这不是完整的基因表达调控模拟。")
    return warnings


def _dedupe_text_tuple(values: tuple[str, ...]) -> tuple[str, ...]:
    deduped: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in deduped:
            deduped.append(text)
    return tuple(deduped)


def _dedupe_text_list(values: list[str]) -> list[str]:
    deduped: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in deduped:
            deduped.append(text)
    return deduped


def _capability_for_preview_row(
    input_capability: dict[str, Any] | None,
    canonical_capability: dict[str, Any] | None,
) -> dict[str, Any] | None:
    if input_capability and input_capability.get("resolved"):
        return input_capability
    return canonical_capability or input_capability


def _canonical_by_input_from_plan(input_ids: tuple[str, ...], input_by_gene: dict[str, str]) -> dict[str, str]:
    canonical_by_input = {input_id: input_id for input_id in input_ids}
    for canonical_id, input_id in input_by_gene.items():
        if input_id in canonical_by_input:
            canonical_by_input[input_id] = canonical_id
    return canonical_by_input


def _input_ids_by_canonical(canonical_by_input: dict[str, str]) -> dict[str, list[str]]:
    input_ids: dict[str, list[str]] = {}
    for input_id, canonical_id in canonical_by_input.items():
        input_ids.setdefault(canonical_id, []).append(input_id)
    return input_ids


def _attach_gene_mapping_input_ids(
    gene_mapping: dict[str, Any],
    input_ids_by_canonical: dict[str, list[str]],
) -> None:
    for gene in gene_mapping.get("genes") or []:
        canonical_id = str(gene.get("gene_id") or "")
        input_ids = input_ids_by_canonical.get(canonical_id, [canonical_id])
        gene["input_gene_ids"] = input_ids
        gene["input_gene_id"] = input_ids[0] if len(input_ids) == 1 else ", ".join(input_ids)
        gene["canonical_gene_id"] = canonical_id
    for row in gene_mapping.get("rows") or []:
        canonical_id = str(row.get("gene_id") or "")
        input_ids = input_ids_by_canonical.get(canonical_id, [canonical_id])
        row["input_gene_ids"] = input_ids
        row["input_gene_id"] = input_ids[0] if len(input_ids) == 1 else ", ".join(input_ids)
        row["canonical_gene_id"] = canonical_id


__all__ = ["preview_screen_inputs", "_preview_screen_inputs_for_model"]
