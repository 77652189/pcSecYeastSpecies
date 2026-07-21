from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

from pcsec_pichia.core.modes import (
    CompatibilityMode,
    DEFAULT_COMPATIBILITY_MODE,
    GlycosylationMode,
    ResultStatus,
)


PICHIA_MU_MIN = 0.01
PICHIA_MU_MAX = 0.44


@dataclass(frozen=True)
class PichiaSimulationRequest:
    target_id: str
    candidate_id: str
    mu: float = 0.10
    production_ratio: float = 1e-8
    media_type: int = 4
    carbon_source_id: str = "glucose"
    timeout_seconds: int = 300
    target_input: Any | None = None
    leader_candidate: Any | None = None
    compatibility_mode: CompatibilityMode = DEFAULT_COMPATIBILITY_MODE
    glycosylation_mode: GlycosylationMode = "native"
    enable_ribosome_translation_constraint: bool = False
    # direction_3_erad_constraint_activation (2026-07-17): stays optional, not
    # a new default. Real kcat backs these 1418 rows and hLF/OPN both stay
    # feasible with it on (test_pipeline_runs_builtin_{hlf,opn}_with_optional_
    # constraints); a small ERAD/proteasome candidate set showed 5%-14% real
    # secretion_ratio_vs_wildtype swings when toggled. Flipping the default
    # anyway would (a) silently change every existing direction-2 screen's
    # candidate-ranking semantics without an explicit opt-in, which
    # EXECUTION_PLAN.md's ERAD acceptance line forbids, and (b) add 1418
    # constraint rows to every solve genome-wide, including the vast majority
    # of genes that never touch the ERAD/misfolding pathway and gain nothing
    # from it. Turn this on explicitly per-run when a candidate is already
    # known/suspected ERAD/proteasome-pathway-adjacent (see
    # services/gene_catalog.py's CAT_ERAD/CAT_PROTEASOME), not globally.
    enable_misfolding_constraint: bool = False
    growth_points: tuple[float, ...] = ()  # empty = auto grid around mu; see pipeline._growth_points
    ko_gene_ids: tuple[str, ...] = ()
    ko_reaction_ids: tuple[str, ...] = ()
    oe_gene_ids: tuple[str, ...] = ()
    oe_reaction_ids: tuple[str, ...] = ()
    screen_candidate_limit: int = 20
    enable_gene_rule_overlay: bool = False
    enable_cost_slope_compatibility: bool = False
    cost_slope_growth_rates: tuple[float, ...] = (0.05, 0.10)
    cost_slope_secretion_ratios: tuple[float, ...] = ()
    cost_slope_capacity_fractions: tuple[float, ...] = (0.10, 0.25, 0.50, 0.75, 0.90)
    cost_slope_medium_compatibility_mode: str = "corrected"
    # R1 (ADR-004): opt-in solver-robustness re-solve of the LP-attribution bottleneck.
    # Off by default because each extra method re-solves the full pcSec LP (slow); turn on
    # per-run to check whether the top bottleneck is stable across HiGHS algorithms or is a
    # degenerate/numerical artifact. Does not change the default solver or the objective.
    enable_solver_robustness_check: bool = False
    # Methods to re-solve with and compare against the deterministic default (highs-ds).
    solver_robustness_methods: tuple[str, ...] = ("highs", "highs-ipm")
    sequence_role: str = "unknown"
    normalization_mode: str = "as_provided"
    contains_signal_peptide: bool | None = None
    contains_leader: bool | None = None
    terminal_stop_policy: str = "allow_for_record_only"
    original_sequence_length: int | None = None
    normalized_sequence_length: int | None = None
    original_full_sequence_length: int | None = None
    normalized_full_sequence_length: int | None = None
    original_leader_sequence_length: int | None = None
    normalized_leader_sequence_length: int | None = None
    original_signal_peptide_length: int | None = None
    normalized_signal_peptide_length: int | None = None
    terminal_stop_present: bool | None = None
    terminal_stop_removed: bool | None = None


@dataclass(frozen=True)
class PichiaSimulationRunResult:
    success: bool
    target_id: str
    candidate_id: str
    mu: float
    production_ratio: float | None
    media_type: int
    message: str
    lp_file: Path | None = None
    output_file: Path | None = None
    objective_value: str | None = None
    command_output: str = ""
    result_status: ResultStatus = "draft"
    summary_path: Path | None = None
    report_path: Path | None = None
    matlab_alignment_status: str = "pending"
    constraint_counts: dict[str, int] = field(default_factory=dict)
    candidate_table_path: Path | None = None
    tradeoff_path: Path | None = None
    alignment_summary: dict[str, Any] = field(default_factory=dict)


class PichiaEngine(Protocol):
    engine_name: str

    def run_target_smoke(self, request: PichiaSimulationRequest) -> PichiaSimulationRunResult:
        """Run a small target-protein secretion simulation."""
