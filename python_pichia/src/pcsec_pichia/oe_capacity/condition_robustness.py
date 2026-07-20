from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class ConditionRankingCase:
    """One (carbon_source, growth_rate) context's relative-OE ranking.

    Built from a single submit_oe_capacity_screen() result; carries no
    solver call of its own.
    """

    context_id: str
    carbon_source_id: str
    growth_rate: float
    ranked_gene_ids: tuple[str, ...]
    gene_deltas: tuple[tuple[str, float], ...]
    missing_gene_ids: tuple[str, ...]


@dataclass(frozen=True)
class ConditionRobustnessResult:
    target_id: str
    gene_ids: tuple[str, ...]
    baseline_case: ConditionRankingCase
    cases: tuple[ConditionRankingCase, ...]
    stable_cases: tuple[ConditionRankingCase, ...]
    unstable_cases: tuple[ConditionRankingCase, ...]
    top1_stable_cases: tuple[ConditionRankingCase, ...]
    top1_unstable_cases: tuple[ConditionRankingCase, ...]

    @property
    def full_order_is_stable(self) -> bool:
        return not self.unstable_cases

    @property
    def top1_is_stable(self) -> bool:
        return not self.top1_unstable_cases


def ranking_case_from_screen_result(
    *,
    context_id: str,
    carbon_source_id: str,
    growth_rate: float,
    screen_result: Mapping[str, Any],
) -> ConditionRankingCase:
    """Rank genes by relative_vs_baseline_delta from one screen's rows.

    This is the same comparison value _render_row_comparison() surfaces in
    the UI (app/ui/views/oe_capacity.py): the relative-uncalibrated OE
    objective minus this gene's own no-OE baseline, under this one context.
    Genes without a usable delta (not executable, failed, or the requested
    scenario set didn't include "nominal") are reported separately rather
    than silently dropped, so an unstable ranking can't be an artifact of
    quietly shrinking the compared set.
    """

    scored: list[tuple[str, float]] = []
    missing: list[str] = []
    for row in (*(screen_result.get("rows") or []), *(screen_result.get("failures") or [])):
        gene_id = str(row.get("gene_id"))
        delta = row.get("relative_vs_baseline_delta")
        if delta is None:
            missing.append(gene_id)
        else:
            scored.append((gene_id, float(delta)))
    scored.sort(key=lambda item: item[1], reverse=True)
    return ConditionRankingCase(
        context_id=context_id,
        carbon_source_id=carbon_source_id,
        growth_rate=float(growth_rate),
        ranked_gene_ids=tuple(gene_id for gene_id, _ in scored),
        gene_deltas=tuple(scored),
        missing_gene_ids=tuple(missing),
    )


def evaluate_condition_robustness(
    *,
    target_id: str,
    gene_ids: tuple[str, ...],
    cases: tuple[ConditionRankingCase, ...],
) -> ConditionRobustnessResult:
    """Partition condition cases into stable/unstable against the first case.

    Mirrors analysis.shadow_lp.validation.run_shadow_validation_matrix's
    "compute every case, then partition into within_tolerance vs
    review_required" shape, with the ranking's own first case standing in
    for shadow_lp's external reference. This function never solves
    anything; the caller (a tools/ script, not this module) is responsible
    for actually calling submit_oe_capacity_screen once per context.
    """

    if not cases:
        raise ValueError("evaluate_condition_robustness requires at least one case.")
    baseline = cases[0]
    stable = tuple(case for case in cases[1:] if case.ranked_gene_ids == baseline.ranked_gene_ids)
    unstable = tuple(case for case in cases[1:] if case.ranked_gene_ids != baseline.ranked_gene_ids)
    baseline_top1 = baseline.ranked_gene_ids[0] if baseline.ranked_gene_ids else None
    top1_stable = tuple(
        case
        for case in cases[1:]
        if case.ranked_gene_ids and case.ranked_gene_ids[0] == baseline_top1
    )
    top1_unstable = tuple(
        case
        for case in cases[1:]
        if not case.ranked_gene_ids or case.ranked_gene_ids[0] != baseline_top1
    )
    return ConditionRobustnessResult(
        target_id=target_id,
        gene_ids=tuple(gene_ids),
        baseline_case=baseline,
        cases=tuple(cases),
        stable_cases=stable,
        unstable_cases=unstable,
        top1_stable_cases=top1_stable,
        top1_unstable_cases=top1_unstable,
    )


__all__ = [
    "ConditionRankingCase",
    "ConditionRobustnessResult",
    "evaluate_condition_robustness",
    "ranking_case_from_screen_result",
]
