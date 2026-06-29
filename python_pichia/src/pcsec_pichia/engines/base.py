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
    timeout_seconds: int = 300
    target_input: Any | None = None
    leader_candidate: Any | None = None
    compatibility_mode: CompatibilityMode = DEFAULT_COMPATIBILITY_MODE
    glycosylation_mode: GlycosylationMode = "native"
    enable_ribosome_translation_constraint: bool = False
    enable_misfolding_constraint: bool = False
    growth_points: tuple[float, ...] = (0.10,)
    ko_gene_ids: tuple[str, ...] = ()
    ko_reaction_ids: tuple[str, ...] = ()
    oe_gene_ids: tuple[str, ...] = ()
    oe_reaction_ids: tuple[str, ...] = ()
    screen_candidate_limit: int = 20
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
