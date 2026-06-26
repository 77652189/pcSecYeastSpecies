from __future__ import annotations

import re
from typing import Any


def reactions_for_gene(model: Any, gene_id: str) -> list[str]:
    gene_index = getattr(model, "gene_index", {})
    if gene_id not in gene_index:
        return []
    gene_number = int(gene_index[gene_id]) + 1
    gene_pattern = re.compile(rf"(?<![A-Za-z0-9_.-]){re.escape(gene_id)}(?![A-Za-z0-9_.-])")
    reactions: list[str] = []
    for reaction_id, rule, gr_rule in zip(
        getattr(model, "rxns", []),
        getattr(model, "rules", []),
        getattr(model, "gr_rules", []),
    ):
        rule_text = str(rule or "")
        gr_rule_text = str(gr_rule or "")
        rule_matches = any(int(match.group(1)) == gene_number for match in re.finditer(r"x\((\d+)\)", rule_text))
        gr_rule_matches = bool(gene_pattern.search(gr_rule_text))
        if rule_matches or gr_rule_matches:
            reactions.append(str(reaction_id))
    return reactions


def resolve_oe_gene_reactions(
    model: Any,
    gene_ids: tuple[str, ...],
    limit: int,
) -> tuple[list[str], dict[str, str], tuple[str, ...], list[str]]:
    reactions: list[str] = []
    gene_by_reaction: dict[str, str] = {}
    unresolved: list[str] = []
    warnings: list[str] = []
    for gene_id in gene_ids:
        if gene_id not in getattr(model, "gene_index", {}):
            unresolved.append(gene_id)
            continue
        matched = reactions_for_gene(model, gene_id)
        if not matched:
            unresolved.append(gene_id)
            continue
        for reaction_id in matched:
            if len(reactions) >= limit:
                warnings.append(f"过表达基因 {gene_id} 解析到的反应超过候选上限 {limit}，已截断。")
                break
            if reaction_id not in gene_by_reaction:
                reactions.append(reaction_id)
                gene_by_reaction[reaction_id] = gene_id
            elif gene_id not in gene_by_reaction[reaction_id].split(","):
                gene_by_reaction[reaction_id] = f"{gene_by_reaction[reaction_id]},{gene_id}"
    return reactions, gene_by_reaction, tuple(unresolved), warnings


def split_existing_genes(model: Any, gene_ids: tuple[str, ...]) -> tuple[list[str], tuple[str, ...]]:
    gene_index = getattr(model, "gene_index", {})
    existing: list[str] = []
    unresolved: list[str] = []
    for gene_id in gene_ids:
        if gene_id in gene_index:
            existing.append(gene_id)
        else:
            unresolved.append(gene_id)
    return existing, tuple(unresolved)


def split_existing_reactions(model: Any, reaction_ids: tuple[str, ...]) -> tuple[list[str], tuple[str, ...]]:
    reaction_index = getattr(model, "reaction_index", {})
    existing: list[str] = []
    unresolved: list[str] = []
    for reaction_id in reaction_ids:
        if reaction_id in reaction_index:
            existing.append(reaction_id)
        else:
            unresolved.append(reaction_id)
    return existing, tuple(unresolved)


__all__ = ["reactions_for_gene", "resolve_oe_gene_reactions", "split_existing_genes", "split_existing_reactions"]
