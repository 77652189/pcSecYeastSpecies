from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from pcsec_pichia.screens._prototype_adapter import default_ko_genes, default_oe_reactions
from pcsec_pichia.screens.candidate_resolution import split_existing_genes, split_existing_reactions
from pcsec_pichia.screens.gene_interventions import GeneInterventionPlan, plan_gene_overexpression
from pcsec_pichia.services import gene_evidence
from pcsec_pichia.services.gene_evidence import GeneExternalEvidence
from pcsec_pichia.services.gene_rule_overlay import GeneRuleEvidence, overlay_aliases_for_executable_rules


@dataclass(frozen=True)
class ScreenPlanResult:
    candidate_limit: int
    requested_ko_gene_ids: tuple[str, ...]
    requested_ko_rxn_ids: tuple[str, ...]
    requested_oe_gene_ids: tuple[str, ...]
    requested_oe_reaction_ids: tuple[str, ...]
    ko_gene_ids: list[str]
    ko_input_by_gene: dict[str, str]
    unresolved_ko_gene_ids: tuple[str, ...]
    ko_reaction_ids: list[str]
    unresolved_ko_reaction_ids: tuple[str, ...]
    oe_gene_reaction_ids: list[str]
    oe_gene_by_reaction: dict[str, str]
    oe_gene_plans_by_gene: dict[str, GeneInterventionPlan]
    oe_input_by_gene: dict[str, str]
    oe_gene_explain_only_plans: tuple[GeneInterventionPlan, ...]
    unresolved_oe_gene_ids: tuple[str, ...]
    oe_reaction_ids: list[str]
    unresolved_oe_reaction_ids: tuple[str, ...]
    warnings: list[str]
    has_user_ko: bool
    has_user_oe: bool

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_screen_plan(
    model: Any,
    request: Any,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    gene_rule_evidence_by_name: dict[str, GeneRuleEvidence] | None = None,
    gene_rule_overlay: Any | None = None,
    evidence_by_gene: dict[str, GeneExternalEvidence] | None = None,
) -> ScreenPlanResult:
    limit = _candidate_limit(getattr(request, "screen_candidate_limit", 20))
    requested_ko = _dedupe_nonempty(_request_tuple(request, "ko_gene_ids", "ko_candidates"))[:limit]
    requested_ko_rxns = _dedupe_nonempty(_request_tuple(request, "ko_reaction_ids"))[:limit]
    requested_oe_genes = _dedupe_nonempty(_request_tuple(request, "oe_gene_ids"))[:limit]
    requested_oe_reactions = _dedupe_nonempty(_request_tuple(request, "oe_reaction_ids", "oe_candidates"))[:limit]
    has_user_ko = bool(requested_ko or requested_ko_rxns)
    has_user_oe = bool(requested_oe_genes or requested_oe_reactions)

    evidence_by_gene = gene_evidence.load_gene_evidence_cache() if evidence_by_gene is None else evidence_by_gene
    overlay_aliases = (
        overlay_aliases_for_executable_rules(gene_rule_evidence_by_name or {}, gene_rule_overlay)
        if getattr(request, "enable_gene_rule_overlay", False)
        else {}
    )
    ko_gene_candidates, ko_alias_warnings, ko_input_by_gene = _canonicalize_gene_ids(
        model,
        requested_ko,
        evidence_by_gene,
        overlay_aliases,
    )
    oe_gene_candidates, oe_alias_warnings, oe_input_by_gene = _canonicalize_gene_ids(
        model,
        requested_oe_genes,
        evidence_by_gene,
        overlay_aliases,
    )

    ko_gene_ids, unresolved_ko = split_existing_genes(model, tuple(ko_gene_candidates))
    if not has_user_ko and not has_user_oe and limit > 0:
        ko_gene_ids = default_ko_genes(model, 1)

    existing_ko_rxns, unresolved_ko_rxns = split_existing_reactions(model, requested_ko_rxns)

    oe_gene_plans = tuple(
        plan_gene_overexpression(model, gene_id, complex_subunits=complex_subunits)
        for gene_id in oe_gene_candidates
    )
    (
        resolved_oe_from_genes,
        oe_gene_by_reaction,
        unresolved_oe_genes,
        gene_warnings,
        oe_gene_explain_only_plans,
    ) = _resolve_oe_gene_plans(oe_gene_plans, limit)
    existing_oe_reactions, unresolved_oe_reactions = split_existing_reactions(model, requested_oe_reactions)
    if not has_user_ko and not has_user_oe and limit > 0:
        existing_oe_reactions = default_oe_reactions(model, 1)

    warnings = [*ko_alias_warnings, *oe_alias_warnings, *gene_warnings]
    warnings.extend(f"敲除基因未在模型中找到：{gene_id}" for gene_id in unresolved_ko)
    warnings.extend(f"敲除反应未在模型中找到：{reaction_id}" for reaction_id in unresolved_ko_rxns)
    warnings.extend(f"过表达基因无法解析到模型反应：{gene_id}" for gene_id in unresolved_oe_genes)
    warnings.extend(f"过表达反应未在模型中找到：{reaction_id}" for reaction_id in unresolved_oe_reactions)
    if requested_oe_genes:
        warnings.append("过表达基因先进行 GPR-aware 规划；可解释的单基因/同工酶反应才运行 reaction-level OE proxy。")
        warnings.append("复合体亚基单基因 OE 默认只解释不增强 capacity，避免虚构分泌提升。")

    return ScreenPlanResult(
        candidate_limit=limit,
        requested_ko_gene_ids=tuple(requested_ko),
        requested_ko_rxn_ids=tuple(requested_ko_rxns),
        requested_oe_gene_ids=tuple(requested_oe_genes),
        requested_oe_reaction_ids=tuple(requested_oe_reactions),
        ko_gene_ids=ko_gene_ids[:limit],
        ko_input_by_gene=ko_input_by_gene,
        unresolved_ko_gene_ids=unresolved_ko,
        ko_reaction_ids=existing_ko_rxns[:limit],
        unresolved_ko_reaction_ids=unresolved_ko_rxns,
        oe_gene_reaction_ids=resolved_oe_from_genes[:limit],
        oe_gene_by_reaction=oe_gene_by_reaction,
        oe_gene_plans_by_gene={plan.gene_id: plan for plan in oe_gene_plans},
        oe_input_by_gene=oe_input_by_gene,
        oe_gene_explain_only_plans=oe_gene_explain_only_plans,
        unresolved_oe_gene_ids=unresolved_oe_genes,
        oe_reaction_ids=existing_oe_reactions[:limit],
        unresolved_oe_reaction_ids=unresolved_oe_reactions,
        warnings=warnings,
        has_user_ko=has_user_ko,
        has_user_oe=has_user_oe,
    )


def _request_tuple(request: Any, primary_name: str, fallback_name: str | None = None) -> tuple[str, ...]:
    primary = getattr(request, primary_name, ()) or ()
    if primary or fallback_name is None:
        return tuple(primary)
    return tuple(getattr(request, fallback_name, ()) or ())


def _resolve_oe_gene_plans(
    plans: tuple[GeneInterventionPlan, ...],
    limit: int,
) -> tuple[list[str], dict[str, str], tuple[str, ...], list[str], tuple[GeneInterventionPlan, ...]]:
    reactions: list[str] = []
    gene_by_reaction: dict[str, str] = {}
    unresolved: list[str] = []
    warnings: list[str] = []
    explain_only: list[GeneInterventionPlan] = []
    for plan in plans:
        warnings.extend(str(item) for item in plan.warnings)
        if not plan.resolved:
            unresolved.append(plan.gene_id)
            continue
        if not plan.executable_reactions:
            explain_only.append(plan)
            continue
        if plan.explain_only_reactions:
            explain_only.append(plan)
        for reaction_id in plan.executable_reactions:
            if len(reactions) >= limit:
                warnings.append(f"过表达基因 {plan.gene_id} 解析到的可运行反应超过候选上限 {limit}，已截断。")
                break
            if reaction_id not in gene_by_reaction:
                reactions.append(reaction_id)
                gene_by_reaction[reaction_id] = plan.gene_id
            elif plan.gene_id not in gene_by_reaction[reaction_id].split(","):
                gene_by_reaction[reaction_id] = f"{gene_by_reaction[reaction_id]},{plan.gene_id}"
    return reactions, gene_by_reaction, tuple(unresolved), warnings, tuple(explain_only)


def _candidate_limit(value: int) -> int:
    try:
        return max(0, min(int(value), 20))
    except (TypeError, ValueError):
        return 20


def _dedupe_nonempty(values: tuple[str, ...]) -> list[str]:
    deduped: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in deduped:
            deduped.append(text)
    return deduped


def _canonicalize_gene_ids(
    model: Any,
    gene_ids: list[str],
    evidence_by_gene: dict[str, GeneExternalEvidence],
    overlay_aliases: dict[str, str] | None = None,
) -> tuple[list[str], list[str], dict[str, str]]:
    canonical_ids: list[str] = []
    warnings: list[str] = []
    input_by_gene: dict[str, str] = {}
    for gene_id in gene_ids:
        canonical = _resolve_overlay_gene_identifier(gene_id, overlay_aliases or {})
        if canonical is None:
            canonical = gene_evidence.resolve_gene_identifier(gene_id, evidence_by_gene, model=model)
        canonical = str(canonical or gene_id).strip()
        canonical = canonical or gene_id
        if canonical not in canonical_ids:
            canonical_ids.append(canonical)
        input_by_gene.setdefault(canonical, gene_id)
        if canonical != gene_id:
            warnings.append(f"基因别名 `{gene_id}` 已解析为模型基因 ID `{canonical}`。")
    return canonical_ids, warnings, input_by_gene


def _resolve_overlay_gene_identifier(gene_id: str, overlay_aliases: dict[str, str]) -> str | None:
    query = str(gene_id or "").strip()
    if not query:
        return None
    for alias, locus_tag in overlay_aliases.items():
        if query.lower() == str(alias).lower() and locus_tag:
            return str(locus_tag)
    return None


__all__ = ["ScreenPlanResult", "build_screen_plan"]
