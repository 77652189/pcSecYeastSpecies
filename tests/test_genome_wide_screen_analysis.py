from __future__ import annotations

import pandas as pd

from app.services.genome_wide_screen_analysis import (
    analyze_single_target,
    complex_subunit_oe_hypothesis_candidates,
    load_gene_tradeoff_csv,
)


def _row(
    gene_id: str,
    intervention_type: str,
    *,
    secretion_ratio: float | None,
    growth_retention: float | None = 1.0,
    max_feasible_mu: float | None = 0.10,
    skipped_reason: str | None = None,
    gpr_role: str = "single_gene",
    secretory_process: str = "代谢或其它反应",
    candidate_kind: str = "gene",
    hypothesis_note: str = "",
    feasibility_interpretation: str = "definitive",
    has_timeout: bool = False,
    timeout_mu_points: str = "",
    proven_infeasible_mu_points: str = "",
    other_solver_failure_mu_points: str = "",
    solver_retry_count: int = 0,
    timeout_retry_mu_points: str = "",
) -> dict[str, object]:
    return {
        "target_id": "hLF",
        "gene_id": gene_id,
        "common_name": "",
        "candidate_kind": candidate_kind,
        "intervention_type": intervention_type,
        "secretion_ratio_vs_wildtype": secretion_ratio,
        "growth_retention_ratio": growth_retention,
        "max_feasible_mu": max_feasible_mu,
        "skipped_reason": skipped_reason,
        "gpr_role": gpr_role,
        "secretory_process": secretory_process,
        "affected_reactions": "R1",
        "hypothesis_note": hypothesis_note,
        "feasibility_interpretation": feasibility_interpretation,
        "has_timeout": has_timeout,
        "timeout_mu_points": timeout_mu_points,
        "proven_infeasible_mu_points": proven_infeasible_mu_points,
        "other_solver_failure_mu_points": other_solver_failure_mu_points,
        "solver_retry_count": solver_retry_count,
        "timeout_retry_mu_points": timeout_retry_mu_points,
    }


def _frame(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_essential_gene_is_ko_infeasible_and_excluded_from_yield_down() -> None:
    frame = _frame(
        [
            _row("G_ESSENTIAL", "KO", secretion_ratio=None, growth_retention=None, max_feasible_mu=None),
            _row("G_ESSENTIAL", "OE", secretion_ratio=1.0),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert list(result.essential_genes["gene_id"]) == ["G_ESSENTIAL"]
    assert result.ko_yield_down.empty
    assert result.ko_yield_up_growth_cost.empty
    assert result.ko_clean_wins.empty


def test_timeout_ko_is_solver_inconclusive_not_essential() -> None:
    frame = _frame(
        [
            _row(
                "G_TIMEOUT",
                "KO",
                secretion_ratio=None,
                growth_retention=None,
                max_feasible_mu=None,
                feasibility_interpretation="inconclusive_due_to_timeout",
                has_timeout=True,
                timeout_mu_points="[0.1]",
            ),
            _row(
                "G_PROVEN_ESSENTIAL",
                "KO",
                secretion_ratio=None,
                growth_retention=None,
                max_feasible_mu=None,
                feasibility_interpretation="definitive",
                proven_infeasible_mu_points="[0.1]",
            ),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert list(result.essential_genes["gene_id"]) == ["G_PROVEN_ESSENTIAL"]
    assert list(result.solver_inconclusive_ko["gene_id"]) == ["G_TIMEOUT"]
    assert result.solver_inconclusive_ko.iloc[0]["feasibility_interpretation"] == "inconclusive_due_to_timeout"
    assert bool(result.solver_inconclusive_ko.iloc[0]["has_timeout"]) is True


def test_retry_resolved_success_is_retry_evidence_not_solver_inconclusive() -> None:
    frame = _frame(
        [
            _row(
                "G_RETRY_SUCCESS",
                "KO",
                secretion_ratio=1.05,
                growth_retention=1.0,
                max_feasible_mu=0.1,
                solver_retry_count=1,
                timeout_retry_mu_points="[0.1]",
                feasibility_interpretation="definitive",
                has_timeout=False,
            ),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert result.solver_inconclusive_ko.empty
    assert list(result.solver_retry_evidence["gene_id"]) == ["G_RETRY_SUCCESS"]
    assert int(result.solver_retry_evidence.iloc[0]["solver_retry_count"]) == 1
    assert result.solver_retry_evidence.iloc[0]["timeout_retry_mu_points"] == "[0.1]"


def test_retry_still_timeout_remains_solver_inconclusive_with_retry_evidence() -> None:
    frame = _frame(
        [
            _row(
                "G_RETRY_TIMEOUT",
                "KO",
                secretion_ratio=None,
                growth_retention=None,
                max_feasible_mu=None,
                solver_retry_count=1,
                timeout_retry_mu_points="[0.1]",
                feasibility_interpretation="inconclusive_due_to_timeout",
                has_timeout=True,
                timeout_mu_points="[0.1]",
            ),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert list(result.solver_inconclusive_ko["gene_id"]) == ["G_RETRY_TIMEOUT"]
    assert list(result.solver_inconclusive_rows["gene_id"]) == ["G_RETRY_TIMEOUT"]
    assert list(result.solver_retry_evidence["gene_id"]) == ["G_RETRY_TIMEOUT"]
    assert result.solver_inconclusive_ko.iloc[0]["timeout_retry_mu_points"] == "[0.1]"


def test_retry_still_timeout_oe_row_enters_solver_inconclusive_rows() -> None:
    frame = _frame(
        [
            _row(
                "G_RETRY_OE_TIMEOUT",
                "OE",
                secretion_ratio=None,
                growth_retention=None,
                max_feasible_mu=None,
                solver_retry_count=1,
                timeout_retry_mu_points="[0.1]",
                feasibility_interpretation="inconclusive_due_to_timeout",
                has_timeout=True,
                timeout_mu_points="[0.1]",
            ),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert result.solver_inconclusive_ko.empty
    assert list(result.solver_inconclusive_rows["gene_id"]) == ["G_RETRY_OE_TIMEOUT"]
    assert result.solver_inconclusive_rows.iloc[0]["intervention_type"] == "OE"


def test_load_gene_tradeoff_csv_preserves_retry_columns(tmp_path) -> None:
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    pd.DataFrame(
        [
            _row(
                "G_RETRY_SUCCESS",
                "KO",
                secretion_ratio=1.05,
                solver_retry_count=1,
                timeout_retry_mu_points="[0.1]",
            )
        ]
    ).to_csv(csv_path, index=False)

    frame = load_gene_tradeoff_csv(str(csv_path))

    assert "solver_retry_count" in frame.columns
    assert "timeout_retry_mu_points" in frame.columns
    assert int(frame.iloc[0]["solver_retry_count"]) == 1
    assert frame.iloc[0]["timeout_retry_mu_points"] == "[0.1]"


def test_ko_yield_down_catches_feasible_but_worse_than_wildtype_knockouts() -> None:
    frame = _frame(
        [
            # A viable complex-subunit KO that hurts secretion - this used to be invisible:
            # not essential (feasible), not yield-up (ratio < 1), so it fell through every
            # existing bucket.
            _row(
                "G_COMPLEX_HURTS",
                "KO",
                secretion_ratio=0.75,
                growth_retention=1.0,
                gpr_role="complex_subunit",
            ),
            _row("G_COMPLEX_HURTS", "OE", secretion_ratio=None, skipped_reason="no_structural_effect"),
            # A neutral gene: ratio essentially 1.0, should not appear in yield_down.
            _row("G_NEUTRAL", "KO", secretion_ratio=1.0005),
            _row("G_NEUTRAL", "OE", secretion_ratio=1.0),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert list(result.ko_yield_down["gene_id"]) == ["G_COMPLEX_HURTS"]
    assert result.ko_yield_down.iloc[0]["gpr_role"] == "complex_subunit"
    assert "G_NEUTRAL" not in set(result.ko_yield_down["gene_id"])
    assert result.essential_genes.empty


def test_ko_yield_down_sorted_worst_first() -> None:
    frame = _frame(
        [
            _row("G_MILD", "KO", secretion_ratio=0.98),
            _row("G_MILD", "OE", secretion_ratio=1.0),
            _row("G_SEVERE", "KO", secretion_ratio=0.5),
            _row("G_SEVERE", "OE", secretion_ratio=1.0),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert list(result.ko_yield_down["gene_id"]) == ["G_SEVERE", "G_MILD"]


def test_yield_up_dimensions_still_split_by_growth_cost() -> None:
    frame = _frame(
        [
            _row("G_COSTLY_WIN", "KO", secretion_ratio=1.05, growth_retention=0.8),
            _row("G_COSTLY_WIN", "OE", secretion_ratio=1.0),
            _row("G_CLEAN_WIN", "KO", secretion_ratio=1.08, growth_retention=1.0),
            _row("G_CLEAN_WIN", "OE", secretion_ratio=1.0),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert list(result.ko_yield_up_growth_cost["gene_id"]) == ["G_COSTLY_WIN"]
    assert list(result.ko_clean_wins["gene_id"]) == ["G_CLEAN_WIN"]
    assert result.ko_yield_down.empty


def test_to_summary_dict_includes_yield_down() -> None:
    frame = _frame(
        [
            _row("G_HURTS", "KO", secretion_ratio=0.6, gpr_role="complex_subunit"),
            _row("G_HURTS", "OE", secretion_ratio=None, skipped_reason="no_structural_effect"),
        ]
    )

    summary = analyze_single_target(frame, "hLF").to_summary_dict()

    assert summary["ko_yield_down_count"] == 1
    assert summary["ko_yield_down_sample"][0]["gene_id"] == "G_HURTS"


def test_to_summary_dict_includes_solver_inconclusive_ko() -> None:
    frame = _frame(
        [
            _row(
                "G_TIMEOUT",
                "KO",
                secretion_ratio=None,
                growth_retention=None,
                max_feasible_mu=None,
                feasibility_interpretation="inconclusive_due_to_timeout",
                has_timeout=True,
            ),
        ]
    )

    summary = analyze_single_target(frame, "hLF").to_summary_dict()

    assert summary["solver_inconclusive_ko_count"] == 1
    assert summary["solver_inconclusive_ko_sample"][0]["gene_id"] == "G_TIMEOUT"
    assert summary["solver_inconclusive_row_count"] == 1
    assert summary["solver_inconclusive_rows_sample"][0]["gene_id"] == "G_TIMEOUT"


def test_to_summary_dict_includes_solver_retry_evidence() -> None:
    frame = _frame(
        [
            _row(
                "G_RETRY_SUCCESS",
                "OE",
                secretion_ratio=1.05,
                solver_retry_count=1,
                timeout_retry_mu_points="[0.1]",
            ),
        ]
    )

    summary = analyze_single_target(frame, "hLF").to_summary_dict()

    assert summary["solver_retry_evidence_count"] == 1
    assert summary["solver_retry_evidence_sample"][0]["gene_id"] == "G_RETRY_SUCCESS"
    assert summary["solver_retry_evidence_sample"][0]["timeout_retry_mu_points"] == "[0.1]"


def test_complex_subunit_oe_hypothesis_candidates_requires_ko_decrease_and_complex_role() -> None:
    frame = _frame(
        [
            # Qualifies: feasible KO decrease + complex_subunit + no OE data.
            _row("G_HYPOTHESIS_TARGET", "KO", secretion_ratio=0.7, gpr_role="complex_subunit"),
            _row("G_HYPOTHESIS_TARGET", "OE", secretion_ratio=None, skipped_reason="no_structural_effect"),
            # Not qualifying: single_gene role already has real OE data, not what this test is for.
            _row("G_SINGLE_GENE", "KO", secretion_ratio=0.7, gpr_role="single_gene"),
            _row("G_SINGLE_GENE", "OE", secretion_ratio=1.02),
            # Not qualifying: complex_subunit but KO ratio is neutral (no real effect to rescue).
            _row("G_NEUTRAL_COMPLEX", "KO", secretion_ratio=1.0, gpr_role="complex_subunit"),
            _row("G_NEUTRAL_COMPLEX", "OE", secretion_ratio=None, skipped_reason="no_structural_effect"),
        ]
    )

    candidates = complex_subunit_oe_hypothesis_candidates(frame, "hLF")

    assert candidates == ["G_HYPOTHESIS_TARGET"]


def test_complex_oe_hypothesis_dimension_shows_neutral_rows_not_just_wins() -> None:
    """Unlike oe_yield_up, a hypothesis row with ratio ~1.0 IS the reportable result
    ("no rescue") - it must not be filtered out the way a neutral ordinary OE row would be.
    """
    frame = _frame(
        [
            _row(
                "ATPS3m_no_1_fwd",
                "OE",
                secretion_ratio=0.999998,
                candidate_kind="complex_oe_hypothesis",
                hypothesis_note="assumes proportional whole-complex co-overexpression",
            ),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert list(result.complex_oe_hypothesis["gene_id"]) == ["ATPS3m_no_1_fwd"]
    assert result.complex_oe_hypothesis.iloc[0]["hypothesis_note"] == "assumes proportional whole-complex co-overexpression"
    assert result.oe_yield_up.empty  # ratio ~1.0 correctly does NOT count as an ordinary OE win


def test_complex_oe_hypothesis_dimension_excludes_ordinary_oe_rows() -> None:
    frame = _frame(
        [
            _row("G_ORDINARY", "OE", secretion_ratio=1.05, candidate_kind="gene"),
        ]
    )

    result = analyze_single_target(frame, "hLF")

    assert result.complex_oe_hypothesis.empty
    assert list(result.oe_yield_up["gene_id"]) == ["G_ORDINARY"]


def test_to_summary_dict_includes_complex_oe_hypothesis() -> None:
    frame = _frame(
        [
            _row("R1", "OE", secretion_ratio=1.0, candidate_kind="complex_oe_hypothesis", hypothesis_note="note"),
        ]
    )

    summary = analyze_single_target(frame, "hLF").to_summary_dict()

    assert summary["complex_oe_hypothesis_count"] == 1
    assert summary["complex_oe_hypothesis_sample"][0]["gene_id"] == "R1"


def test_complex_subunit_oe_hypothesis_candidates_scoped_to_target() -> None:
    frame = _frame(
        [
            _row("G_HYPOTHESIS_TARGET", "KO", secretion_ratio=0.7, gpr_role="complex_subunit"),
            _row("G_HYPOTHESIS_TARGET", "OE", secretion_ratio=None, skipped_reason="no_structural_effect"),
        ]
    )
    frame.loc[frame["gene_id"] == "G_HYPOTHESIS_TARGET", "target_id"] = "OPN"

    assert complex_subunit_oe_hypothesis_candidates(frame, "hLF") == []
    assert complex_subunit_oe_hypothesis_candidates(frame, "OPN") == ["G_HYPOTHESIS_TARGET"]
