from __future__ import annotations

import pytest

from pcsec_pichia.oe_capacity import (
    ConditionRankingCase,
    evaluate_condition_robustness,
    ranking_case_from_screen_result,
)


def _screen_result(deltas: dict[str, float | None]) -> dict:
    rows = [{"gene_id": gene_id, "relative_vs_baseline_delta": delta} for gene_id, delta in deltas.items()]
    return {"rows": rows, "failures": []}


def test_ranking_case_orders_genes_by_delta_descending() -> None:
    case = ranking_case_from_screen_result(
        context_id="glucose_mu_0.1",
        carbon_source_id="glucose",
        growth_rate=0.1,
        screen_result=_screen_result({"G1": 0.01, "G2": 0.05, "G3": -0.02}),
    )

    assert case.ranked_gene_ids == ("G2", "G1", "G3")
    assert case.gene_deltas == (("G2", 0.05), ("G1", 0.01), ("G3", -0.02))
    assert case.missing_gene_ids == ()


def test_ranking_case_reports_missing_deltas_separately_instead_of_dropping_them() -> None:
    case = ranking_case_from_screen_result(
        context_id="glucose_mu_0.1",
        carbon_source_id="glucose",
        growth_rate=0.1,
        screen_result=_screen_result({"G1": 0.01, "G2": None}),
    )

    assert case.ranked_gene_ids == ("G1",)
    assert case.missing_gene_ids == ("G2",)


def test_ranking_case_reads_failures_too() -> None:
    result = {
        "rows": [{"gene_id": "G1", "relative_vs_baseline_delta": 0.01}],
        "failures": [{"gene_id": "G2", "relative_vs_baseline_delta": None}],
    }
    case = ranking_case_from_screen_result(
        context_id="glucose_mu_0.1", carbon_source_id="glucose", growth_rate=0.1, screen_result=result
    )

    assert case.ranked_gene_ids == ("G1",)
    assert case.missing_gene_ids == ("G2",)


def _case(context_id: str, ranked_gene_ids: tuple[str, ...]) -> ConditionRankingCase:
    return ConditionRankingCase(
        context_id=context_id,
        carbon_source_id=context_id.split("_mu_")[0],
        growth_rate=0.1,
        ranked_gene_ids=ranked_gene_ids,
        gene_deltas=tuple((gene_id, 0.0) for gene_id in ranked_gene_ids),
        missing_gene_ids=(),
    )


def test_identical_ranking_across_contexts_is_reported_stable() -> None:
    baseline = _case("glucose_mu_0.1", ("G2", "G1", "G3"))
    same_order = _case("glucose_mu_0.15", ("G2", "G1", "G3"))

    result = evaluate_condition_robustness(target_id="hLF", gene_ids=("G1", "G2", "G3"), cases=(baseline, same_order))

    assert result.full_order_is_stable
    assert result.top1_is_stable
    assert result.stable_cases == (same_order,)
    assert result.unstable_cases == ()


def test_reordered_ranking_is_reported_unstable_even_if_top1_unchanged() -> None:
    baseline = _case("glucose_mu_0.1", ("G2", "G1", "G3"))
    reordered_tail = _case("glucose_mu_0.15", ("G2", "G3", "G1"))

    result = evaluate_condition_robustness(
        target_id="hLF", gene_ids=("G1", "G2", "G3"), cases=(baseline, reordered_tail)
    )

    assert not result.full_order_is_stable
    assert result.top1_is_stable
    assert result.unstable_cases == (reordered_tail,)
    assert result.top1_unstable_cases == ()


def test_top1_change_is_reported_unstable() -> None:
    baseline = _case("glucose_mu_0.1", ("G2", "G1", "G3"))
    new_leader = _case("glycerol_mu_0.1", ("G1", "G2", "G3"))

    result = evaluate_condition_robustness(target_id="hLF", gene_ids=("G1", "G2", "G3"), cases=(baseline, new_leader))

    assert not result.full_order_is_stable
    assert not result.top1_is_stable
    assert result.top1_unstable_cases == (new_leader,)


def test_evaluate_condition_robustness_requires_at_least_one_case() -> None:
    with pytest.raises(ValueError, match="at least one case"):
        evaluate_condition_robustness(target_id="hLF", gene_ids=("G1",), cases=())
