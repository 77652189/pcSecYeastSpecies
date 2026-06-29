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
    return _preview_screen_inputs_for_model(
        inputs.prepared_model,
        request,
        complex_subunits=inputs.secretory.complex_subunits,
    )


def _preview_screen_inputs_for_model(
    model: Any,
    request: SecretionRunRequest,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> dict[str, Any]:
    limit = max(0, min(int(request.screen_candidate_limit), 20))
    ko_gene_ids = _dedupe_text_tuple(request.ko_gene_ids or request.ko_candidates)[:limit]
    ko_reaction_ids = _dedupe_text_tuple(request.ko_reaction_ids)[:limit]
    oe_gene_ids = _dedupe_text_tuple(request.oe_gene_ids)[:limit]
    oe_reaction_ids = _dedupe_text_tuple(request.oe_reaction_ids or request.oe_candidates)[:limit]

    from pcsec_pichia.screens import (
        build_gene_perturbation_map,
        reactions_for_gene,
        resolve_oe_gene_reactions,
        split_existing_genes,
        split_existing_reactions,
    )

    existing_ko_genes, _ = split_existing_genes(model, ko_gene_ids)
    existing_ko_reactions, _ = split_existing_reactions(model, ko_reaction_ids)
    _, _, unresolved_oe_genes, oe_gene_warnings = resolve_oe_gene_reactions(
        model,
        oe_gene_ids,
        limit,
    )
    existing_oe_reactions, _ = split_existing_reactions(model, oe_reaction_ids)

    ko_genes = [
        _gene_row(gene_id, "KO", resolved=gene_id in existing_ko_genes)
        for gene_id in ko_gene_ids
    ]
    ko_reactions = [
        _reaction_row(reaction_id, "KO_reaction", resolved=reaction_id in existing_ko_reactions)
        for reaction_id in ko_reaction_ids
    ]
    oe_genes = [
        _oe_gene_row(
            gene_id,
            reactions_for_gene(model, gene_id),
            resolved=gene_id not in unresolved_oe_genes,
            limit=limit,
            truncated=any(gene_id in warning for warning in oe_gene_warnings),
        )
        for gene_id in oe_gene_ids
    ]
    oe_reactions = [
        _reaction_row(reaction_id, "OE_reaction", resolved=reaction_id in existing_oe_reactions)
        for reaction_id in oe_reaction_ids
    ]
    gene_mapping = build_gene_perturbation_map(
        model,
        _dedupe_text_tuple((*ko_gene_ids, *oe_gene_ids)),
        complex_subunits=complex_subunits,
    ).to_dict()
    warnings = _screen_preview_warnings(ko_genes, ko_reactions, oe_genes, oe_reactions, request)

    return {
        "candidate_limit": limit,
        "semantics": {
            "KO": "按模型基因 ID 敲除；未在模型 gene_index 中找到的 ID 会写入 unresolved_gene。",
            "KO_reaction": "按模型反应 ID 直接敲除，主要用于复合体级扰动。",
            "OE_gene_proxy": "过表达基因会先解析到该基因参与的反应，再按 reaction-level OE proxy 逐个模拟。",
            "OE_reaction": "按模型反应 ID 直接做反应级过表达代理。",
        },
        "ko_genes": ko_genes,
        "ko_reactions": ko_reactions,
        "oe_genes": oe_genes,
        "oe_reactions": oe_reactions,
        "gene_mapping": gene_mapping,
        "gene_mapping_rows": gene_mapping["rows"],
        "warnings": warnings,
    }


def _gene_row(gene_id: str, intervention_type: str, resolved: bool) -> dict[str, Any]:
    return {
        "input_id": gene_id,
        "intervention_type": intervention_type,
        "resolved": resolved,
        "status": "resolved" if resolved else "unresolved_gene",
        "resolved_reaction_count": None,
        "resolved_reactions_preview": [],
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


def _oe_gene_row(
    gene_id: str,
    reactions: list[str],
    resolved: bool,
    limit: int,
    truncated: bool,
) -> dict[str, Any]:
    return {
        "input_id": gene_id,
        "intervention_type": "OE_gene_proxy",
        "resolved": resolved,
        "status": "resolved" if resolved else "unresolved_gene",
        "resolved_reaction_count": len(reactions),
        "resolved_reactions_preview": reactions[: max(0, min(limit, 5))],
        "truncated": truncated,
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
        elif row.get("truncated"):
            warnings.append(
                f"过表达基因 {row['input_id']} 解析到的反应超过候选上限 "
                f"{request.screen_candidate_limit}，正式运行会截断。"
            )
    for row in oe_reactions:
        if not row["resolved"]:
            warnings.append(f"过表达反应未在模型中找到：{row['input_id']}")
    if oe_genes:
        warnings.append("过表达基因当前按 reaction-level OE proxy 运行；这不是完整的基因表达调控模拟。")
    return warnings


def _dedupe_text_tuple(values: tuple[str, ...]) -> tuple[str, ...]:
    deduped: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in deduped:
            deduped.append(text)
    return tuple(deduped)


__all__ = ["preview_screen_inputs", "_preview_screen_inputs_for_model"]
