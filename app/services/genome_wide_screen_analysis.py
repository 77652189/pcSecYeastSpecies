"""Turns raw gene_tradeoff_rows.csv output into the analysis dimensions used
to review a genome-wide KO/OE screen:

1. Essential genes (KO infeasible even at the lowest growth rate tested)
2. Solver-inconclusive KO rows (timeout/solver failure; not proven essential)
3. KO candidates that raise secretion but cost growth (need a rescue strategy)
4. KO candidates that raise secretion with growth fully retained (clean wins)
5. KO candidates that lower secretion (still viable, just worse - not "essential" and not
   "no effect"; without this dimension these rows are invisible in every other bucket)
6. OE candidates that raise secretion
7. Hypothetical whole-complex OE test results (see candidate_kind == "complex_oe_hypothesis";
   these are ~always near ratio 1.0 by nature of what they're testing, so without this
   dimension they too would fall through every ratio-threshold bucket above)
8. Target-specific divergence (same gene, different effect per target)

Each function returns a plain pandas DataFrame so the UI can render it
directly and the report generator can serialize it to text/JSON.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

SECRETION_UP_THRESHOLD = 1.01
SECRETION_DOWN_THRESHOLD = 0.99
GROWTH_COST_THRESHOLD = 0.99
GROWTH_FULLY_RETAINED_THRESHOLD = 0.999
DIVERGENCE_TOP_N = 20


@dataclass(frozen=True)
class DimensionalResults:
    target_id: str
    essential_genes: pd.DataFrame
    solver_inconclusive_ko: pd.DataFrame
    solver_inconclusive_rows: pd.DataFrame
    solver_retry_evidence: pd.DataFrame
    ko_yield_up_growth_cost: pd.DataFrame
    ko_clean_wins: pd.DataFrame
    ko_yield_down: pd.DataFrame
    oe_yield_up: pd.DataFrame
    complex_oe_hypothesis: pd.DataFrame
    row_count: int
    skipped_count: int

    def to_summary_dict(self, max_rows_per_dimension: int = 15) -> dict[str, object]:
        """Compact, LLM-prompt-friendly summary: counts + top rows per dimension."""
        return {
            "target_id": self.target_id,
            "row_count": self.row_count,
            "skipped_count": self.skipped_count,
            "essential_gene_count": len(self.essential_genes),
            "essential_genes_sample": self.essential_genes.head(max_rows_per_dimension).to_dict("records"),
            "solver_inconclusive_ko_count": len(self.solver_inconclusive_ko),
            "solver_inconclusive_ko_sample": self.solver_inconclusive_ko.head(max_rows_per_dimension).to_dict("records"),
            "solver_inconclusive_row_count": len(self.solver_inconclusive_rows),
            "solver_inconclusive_rows_sample": self.solver_inconclusive_rows.head(max_rows_per_dimension).to_dict("records"),
            "solver_retry_evidence_count": len(self.solver_retry_evidence),
            "solver_retry_evidence_sample": self.solver_retry_evidence.head(max_rows_per_dimension).to_dict("records"),
            "ko_yield_up_growth_cost_count": len(self.ko_yield_up_growth_cost),
            "ko_yield_up_growth_cost": self.ko_yield_up_growth_cost.head(max_rows_per_dimension).to_dict("records"),
            "ko_clean_win_count": len(self.ko_clean_wins),
            "ko_clean_wins": self.ko_clean_wins.head(max_rows_per_dimension).to_dict("records"),
            "ko_yield_down_count": len(self.ko_yield_down),
            "ko_yield_down_sample": self.ko_yield_down.head(max_rows_per_dimension).to_dict("records"),
            "oe_yield_up_count": len(self.oe_yield_up),
            "oe_yield_up": self.oe_yield_up.head(max_rows_per_dimension).to_dict("records"),
            "complex_oe_hypothesis_count": len(self.complex_oe_hypothesis),
            "complex_oe_hypothesis_sample": self.complex_oe_hypothesis.head(max_rows_per_dimension).to_dict("records"),
        }


def load_gene_tradeoff_csv(csv_path: str) -> pd.DataFrame:
    frame = pd.read_csv(csv_path)
    # common_name/candidate_kind were added when the curated-catalog reaction screen shipped;
    # older gene-only CSVs won't have them, so backfill defaults instead of KeyError-ing downstream.
    if "candidate_kind" not in frame.columns:
        frame["candidate_kind"] = "gene"
    if "common_name" not in frame.columns:
        frame["common_name"] = ""
    if "hypothesis_note" not in frame.columns:
        frame["hypothesis_note"] = ""
    _ensure_solver_outcome_columns(frame)
    return frame


def _ensure_solver_outcome_columns(frame: pd.DataFrame) -> None:
    defaults = {
        "feasibility_interpretation": "definitive",
        "has_timeout": False,
        "timeout_mu_points": "",
        "proven_infeasible_mu_points": "",
        "other_solver_failure_mu_points": "",
        "solver_retry_count": 0,
        "timeout_retry_mu_points": "",
    }
    for column, default in defaults.items():
        if column not in frame.columns:
            frame[column] = default


def analyze_single_target(frame: pd.DataFrame, target_id: str) -> DimensionalResults:
    """Compute the single-target review dimensions for one target's rows."""
    target_rows = frame[frame.target_id == target_id].copy()
    _ensure_solver_outcome_columns(target_rows)
    ko = target_rows[target_rows.intervention_type == "KO"].dropna(subset=["secretion_ratio_vs_wildtype"])
    oe = target_rows[target_rows.intervention_type == "OE"].dropna(subset=["secretion_ratio_vs_wildtype"])

    display_cols = ["gene_id", "common_name", "candidate_kind"]
    solver_cols = [
        "secretory_process",
        "affected_reactions",
        "feasibility_interpretation",
        "has_timeout",
        "timeout_mu_points",
        "proven_infeasible_mu_points",
        "other_solver_failure_mu_points",
        "solver_retry_count",
        "timeout_retry_mu_points",
    ]
    infeasible_ko_mask = (
        (target_rows.intervention_type == "KO")
        & target_rows.max_feasible_mu.isna()
        & target_rows.skipped_reason.isna()
    )
    inconclusive_mask = target_rows.feasibility_interpretation.isin(
        {"inconclusive_due_to_timeout", "inconclusive_due_to_solver_failure"}
    )

    essential = target_rows[infeasible_ko_mask & ~inconclusive_mask][
        display_cols + ["secretory_process", "affected_reactions"]
    ].reset_index(drop=True)
    solver_inconclusive_rows = target_rows[target_rows.skipped_reason.isna() & inconclusive_mask][
        display_cols + ["intervention_type", "max_feasible_mu", "secretion_ratio_vs_wildtype"] + solver_cols
    ].reset_index(drop=True)
    solver_inconclusive = target_rows[infeasible_ko_mask & inconclusive_mask][display_cols + solver_cols].reset_index(drop=True)
    retry_rows = target_rows[pd.to_numeric(target_rows.solver_retry_count, errors="coerce").fillna(0).astype(int) > 0][
        display_cols + ["intervention_type", "max_feasible_mu", "secretion_ratio_vs_wildtype"] + solver_cols
    ].reset_index(drop=True)

    yield_up_growth_cost = ko[
        (ko.secretion_ratio_vs_wildtype > SECRETION_UP_THRESHOLD) & (ko.growth_retention_ratio < GROWTH_COST_THRESHOLD)
    ][display_cols + ["secretion_ratio_vs_wildtype", "growth_retention_ratio", "secretory_process", "affected_reactions"]]
    yield_up_growth_cost = yield_up_growth_cost.sort_values("secretion_ratio_vs_wildtype", ascending=False).reset_index(drop=True)

    clean_wins = ko[
        (ko.secretion_ratio_vs_wildtype > SECRETION_UP_THRESHOLD) & (ko.growth_retention_ratio >= GROWTH_FULLY_RETAINED_THRESHOLD)
    ][display_cols + ["secretion_ratio_vs_wildtype", "secretory_process", "affected_reactions"]]
    clean_wins = clean_wins.sort_values("secretion_ratio_vs_wildtype", ascending=False).reset_index(drop=True)

    # Feasible (not essential - that's its own dimension) KO rows that make secretion worse.
    # These fall outside both "yield up" buckets above, so without this dimension they are
    # invisible anywhere in the productized results even though the raw CSV has them -
    # e.g. every complex-subunit gene whose knockout hurts a shared metabolic complex.
    yield_down = ko[ko.secretion_ratio_vs_wildtype < SECRETION_DOWN_THRESHOLD][
        display_cols + ["secretion_ratio_vs_wildtype", "growth_retention_ratio", "gpr_role", "secretory_process", "affected_reactions"]
    ]
    yield_down = yield_down.sort_values("secretion_ratio_vs_wildtype", ascending=True).reset_index(drop=True)

    oe_yield_up = oe[oe.secretion_ratio_vs_wildtype > SECRETION_UP_THRESHOLD][
        display_cols + ["secretion_ratio_vs_wildtype", "growth_retention_ratio", "secretory_process", "affected_reactions"]
    ]
    oe_yield_up = oe_yield_up.sort_values("secretion_ratio_vs_wildtype", ascending=False).reset_index(drop=True)

    # Hypothetical whole-complex OE test rows (see COMPLEX_OE_HYPOTHESIS_ASSUMPTION). Shown
    # in full regardless of ratio - unlike oe_yield_up, a ratio near 1.0 here *is* the result
    # ("no rescue"), not a row to discard, so filtering by SECRETION_UP_THRESHOLD would hide
    # the finding instead of surfacing it.
    hypothesis_rows = target_rows[target_rows.candidate_kind == "complex_oe_hypothesis"].dropna(
        subset=["secretion_ratio_vs_wildtype"]
    )
    complex_oe_hypothesis = hypothesis_rows[
        display_cols
        + ["secretion_ratio_vs_wildtype", "growth_retention_ratio", "secretory_process", "affected_reactions", "hypothesis_note"]
    ]
    complex_oe_hypothesis = complex_oe_hypothesis.sort_values("secretion_ratio_vs_wildtype", ascending=False).reset_index(drop=True)

    return DimensionalResults(
        target_id=target_id,
        essential_genes=essential,
        solver_inconclusive_ko=solver_inconclusive,
        solver_inconclusive_rows=solver_inconclusive_rows,
        solver_retry_evidence=retry_rows,
        ko_yield_up_growth_cost=yield_up_growth_cost,
        ko_clean_wins=clean_wins,
        ko_yield_down=yield_down,
        oe_yield_up=oe_yield_up,
        complex_oe_hypothesis=complex_oe_hypothesis,
        row_count=len(target_rows),
        skipped_count=int(target_rows.skipped_reason.notna().sum()),
    )


def complex_subunit_oe_hypothesis_candidates(frame: pd.DataFrame, target_id: str) -> list[str]:
    """Gene ids worth a "hypothetical whole-complex OE" test: KO is feasible but lowers
    secretion (same criterion as ko_yield_down), and the gene's GPR role is complex_subunit
    - meaning the gene-level screen already correctly skipped single-gene OE for it (see
    resolve_complex_subunit_oe_hypothesis_candidates for what "hypothetical" means here).
    """
    ko = frame[(frame.target_id == target_id) & (frame.intervention_type == "KO")].dropna(
        subset=["secretion_ratio_vs_wildtype"]
    )
    candidates = ko[(ko.secretion_ratio_vs_wildtype < SECRETION_DOWN_THRESHOLD) & (ko.gpr_role == "complex_subunit")]
    return candidates["gene_id"].astype(str).tolist()


def analyze_target_divergence(frame: pd.DataFrame, target_ids: list[str], top_n: int = DIVERGENCE_TOP_N) -> pd.DataFrame:
    """Same gene KO, compared pairwise across the first two targets given. Requires >= 2 targets."""
    if len(target_ids) < 2:
        return pd.DataFrame()
    a_id, b_id = target_ids[0], target_ids[1]
    ko = frame[frame.intervention_type == "KO"]
    a = ko[ko.target_id == a_id][["gene_id", "common_name", "secretion_ratio_vs_wildtype", "growth_retention_ratio"]]
    b = ko[ko.target_id == b_id][["gene_id", "secretion_ratio_vs_wildtype", "growth_retention_ratio"]]
    merged = a.merge(b, on="gene_id", suffixes=(f"_{a_id}", f"_{b_id}"))
    ratio_a_col, ratio_b_col = f"secretion_ratio_vs_wildtype_{a_id}", f"secretion_ratio_vs_wildtype_{b_id}"
    merged = merged.dropna(subset=[ratio_a_col, ratio_b_col])
    merged["divergence"] = (merged[ratio_a_col] - merged[ratio_b_col]).abs()
    return merged.sort_values("divergence", ascending=False).head(top_n).reset_index(drop=True)


__all__ = [
    "DIVERGENCE_TOP_N",
    "GROWTH_COST_THRESHOLD",
    "GROWTH_FULLY_RETAINED_THRESHOLD",
    "SECRETION_DOWN_THRESHOLD",
    "SECRETION_UP_THRESHOLD",
    "DimensionalResults",
    "analyze_single_target",
    "analyze_target_divergence",
    "complex_subunit_oe_hypothesis_candidates",
    "load_gene_tradeoff_csv",
]
