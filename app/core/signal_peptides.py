from __future__ import annotations

import re
from dataclasses import dataclass

from pydantic import BaseModel


AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


class SignalPeptideCandidate(BaseModel):
    candidate_id: str
    leader_sequence: str
    signal_peptide_sequence: str
    category: str
    processing_route: str
    source_note: str
    rationale: str
    caution: str
    category_label: str = ""
    library_stage: str = "候选库"
    source_type: str = "待补充"
    leader_length: int | None = None
    signal_peptide_length: int | None = None
    construct_length: int | None = None
    accession: str = ""
    uniprot_id: str = ""
    protein_name: str = ""
    organism_name: str = ""
    protein_sequence: str = ""
    protein_length: int | None = None
    uniprot_signal_start: int | None = None
    uniprot_signal_end: int | None = None
    already_in_formal_library: bool = False

    def as_row(self) -> dict[str, object]:
        row = self.model_dump()
        row["leader_length"] = self.leader_length if self.leader_length is not None else len(self.leader_sequence)
        row["signal_peptide_length"] = (
            self.signal_peptide_length
            if self.signal_peptide_length is not None
            else len(self.signal_peptide_sequence)
        )
        if self.protein_sequence and row["protein_length"] is None:
            row["protein_length"] = len(self.protein_sequence)
        return row


@dataclass(frozen=True)
class CandidateDiscoveryResult:
    rows: list[dict[str, object]]
    source_url: str
    errors: list[str]


@dataclass(frozen=True)
class UniProtCandidateLibraryResult:
    rows: list[dict[str, object]]
    source_url: str
    errors: list[str]
    initial_hit_count: int
    fetched_record_count: int
    extracted_signal_count: int
    deduplicated_count: int
