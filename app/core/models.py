from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field


SpeciesCode = Literal["SCE", "PPA", "KMX", "Unknown"]


class DatasetInfo(BaseModel):
    id: str
    name: str
    path: Path
    category: str
    suffix: str
    species: SpeciesCode = "Unknown"
    size_bytes: int
    modified_at: str


class LoadedDataset(BaseModel):
    info: DatasetInfo
    kind: Literal["table", "mat"]
    tables: dict[str, list[dict[str, Any]]] = Field(default_factory=dict)
    variable_summary: list[dict[str, Any]] = Field(default_factory=list)


class HealthItem(BaseModel):
    name: str
    status: Literal["ok", "warning", "missing", "error"]
    detail: str = ""


class HealthReport(BaseModel):
    items: list[HealthItem]
    preflight_output: str = ""


class SoplexSummary(BaseModel):
    optimal: bool
    objective_value: str | None = None
    status_line: str | None = None


class SimulationResult(BaseModel):
    success: bool
    mu: float
    message: str
    lp_file: Path | None = None
    output_file: Path | None = None
    objective_value: str | None = None
    command_output: str = ""


class OpnCandidate(BaseModel):
    candidate_id: str
    category: str
    leader_sequence: str
    signal_peptide_sequence: str
    processing_route: str
    source_note: str
    rationale: str
    caution: str
    leader_length: int
    construct_length: int
    mature_opn_internal_kex2_like_sites: str = ""
    construct_nxs_t_motifs: str = ""


class OpnSimulationResult(SimulationResult):
    candidate_id: str
    production_ratio: float
    media_type: int


class OpnCandidateRank(BaseModel):
    candidate_id: str
    category: str
    category_label: str
    recommendation: str
    experimental_role: str
    rank: int
    model_rank: int | None = None
    objective_value: float | None = None
    objective_text: str | None = None
    objective_delta_percent: float | None = None
    optimal: bool = False
    leader_length: int
    construct_length: int
    processing_route: str
    evidence_level: str
    risk_level: str
    reason: str
    output_file: Path | None = None


class OpnConstructDesign(BaseModel):
    candidate_id: str
    experimental_role: str
    recommendation: str
    leader_sequence: str
    signal_peptide_sequence: str
    mature_opn_sequence: str
    full_protein_sequence: str
    leader_length: int
    signal_peptide_length: int
    mature_opn_length: int
    full_protein_length: int
    contains_alpha_pro_region: bool
    processing_route: str
    kex2_risk: str
    codon_optimization_next: str
    note: str


class CdsDesignRecord(BaseModel):
    construct_id: str
    experimental_role: str = ""
    recommendation: str = ""
    cds_candidate_rank: int
    recommended_subset: bool = False
    leader_sequence: str = ""
    mature_opn_sequence: str = ""
    full_protein_sequence: str = ""
    aa_length: int
    cds_length: int
    gc_percent: float | None = None
    gc_status: str = ""
    cai_training: float | None = None
    cai_public: float | None = None
    quality_status: str = ""
    warnings: int = 0
    restriction_sites: int = 0
    motif_hits: int = 0
    length_multiple_of_three: bool = False
    translation_matches_input: bool = False
    internal_stop_codons: list[int] = Field(default_factory=list)
    kex2_risk: str = ""
    risk_note: str = ""
    source: str = ""
    cds: str


class CdsDesignResult(BaseModel):
    available: bool
    message: str
    records: list[CdsDesignRecord] = Field(default_factory=list)
    csv_file: Path | None = None
    xlsx_file: Path | None = None
    fasta_file: Path | None = None
