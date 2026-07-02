from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal


TargetSource = Literal["builtin", "custom_json", "custom_sequence"]
SequenceRole = Literal["preprotein", "mature_secreted", "custom_user_sequence", "unknown"]
NormalizationMode = Literal["as_provided", "remove_terminal_stop", "none"]
TerminalStopPolicy = Literal["reject", "strip", "allow_for_record_only"]


@dataclass(frozen=True)
class SecretionRunRequest:
    target_source: TargetSource
    target_id: str
    target_name: str = ""
    sequence: str = ""
    leader_sequence: str = ""
    signal_peptide_sequence: str = ""
    sequence_role: SequenceRole = "unknown"
    normalization_mode: NormalizationMode = "as_provided"
    contains_signal_peptide: bool | None = None
    contains_leader: bool | None = None
    terminal_stop_policy: TerminalStopPolicy = "allow_for_record_only"
    disulfide_sites: int = 0
    n_glycosylation_sites: int = 0
    o_glycosylation_sites: int = 0
    custom_json_path: Path | None = None
    enable_ribosome_translation_constraint: bool = False
    enable_misfolding_constraint: bool = False
    ko_gene_ids: tuple[str, ...] = ()
    oe_gene_ids: tuple[str, ...] = ()
    oe_reaction_ids: tuple[str, ...] = ()
    ko_reaction_ids: tuple[str, ...] = ()
    screen_candidate_limit: int = 20
    enable_gene_rule_overlay: bool = False
    ko_candidates: tuple[str, ...] = ()
    oe_candidates: tuple[str, ...] = ()
    growth_points: tuple[float, ...] = (0.10,)
    enable_cost_slope_compatibility: bool = False
    cost_slope_growth_rates: tuple[float, ...] = (0.05, 0.10)
    cost_slope_secretion_ratios: tuple[float, ...] = ()
    cost_slope_capacity_fractions: tuple[float, ...] = (0.10, 0.25, 0.50, 0.75, 0.90)
    cost_slope_medium_compatibility_mode: str = "corrected"
    mu: float = 0.10
    media_type: int = 4
    carbon_source_id: str = "glucose"
    output_dir: Path | None = None


@dataclass(frozen=True)
class SecretionRunResponse:
    success: bool
    target_id: str
    result_status: str
    matlab_alignment_status: str
    objective_value: str | None = None
    secretion_flux: str | None = None
    constraint_counts: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    output_dir: Path | None = None
    summary_path: Path | None = None
    report_path: Path | None = None
    candidate_table_path: Path | None = None
    tradeoff_path: Path | None = None
    alignment_summary: dict[str, Any] = field(default_factory=dict)
    target_metadata: dict[str, Any] = field(default_factory=dict)
    target_warnings: list[str] = field(default_factory=list)
    protein_cost_analysis: dict[str, Any] = field(default_factory=dict)
    target_growth_analysis: dict[str, Any] = field(default_factory=dict)
    yield_improvement_recommendations: dict[str, Any] = field(default_factory=dict)
    medium_condition: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class BuiltinTargetTemplate:
    target_id: str
    label: str
    parameter_status: str
    source: str
    leader_length: int
    signal_peptide_length: int
    full_sequence_length: int
    alignment_target_kind: str
    sequence_role: str
    normalization_mode: str
    mature_sequence_length: int = 0
    target_warning: str = ""
    note: str = ""


__all__ = [
    "BuiltinTargetTemplate",
    "NormalizationMode",
    "SecretionRunRequest",
    "SecretionRunResponse",
    "SequenceRole",
    "TargetSource",
    "TerminalStopPolicy",
]
