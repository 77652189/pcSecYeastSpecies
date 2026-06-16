from __future__ import annotations

from dataclasses import dataclass

from app.core.signal_peptides import SignalPeptideCandidate
from app.services.opn import OPN_SHORTLIST, OpnCandidateCatalog, opn_category_label


@dataclass(frozen=True)
class OpnSignalPeptideCandidateSource:
    catalog: OpnCandidateCatalog

    def list_candidates(self) -> list[SignalPeptideCandidate]:
        candidates = []
        for candidate in self.catalog.list_candidates():
            candidates.append(
                SignalPeptideCandidate(
                    candidate_id=candidate.candidate_id,
                    category=candidate.category,
                    category_label=opn_category_label(candidate.category),
                    library_stage=_library_stage(candidate.candidate_id, candidate.category),
                    leader_sequence=candidate.leader_sequence,
                    signal_peptide_sequence=candidate.signal_peptide_sequence,
                    leader_length=candidate.leader_length,
                    signal_peptide_length=len(candidate.signal_peptide_sequence),
                    construct_length=candidate.construct_length,
                    processing_route=candidate.processing_route,
                    source_type=_source_type(candidate.source_note),
                    source_note=candidate.source_note,
                    rationale=candidate.rationale,
                    caution=candidate.caution,
                )
            )
        return candidates


def _library_stage(candidate_id: str, category: str) -> str:
    if candidate_id in OPN_SHORTLIST:
        return "首轮推荐"
    if category == "project_baseline":
        return "对照基线"
    return "候选库"


def _source_type(source_note: str) -> str:
    lowered = source_note.lower()
    if "uniprot" in lowered:
        return "UniProt"
    if "reported" in lowered or "paper" in lowered or "pmcid" in lowered:
        return "文献"
    if "project" in lowered:
        return "项目基线"
    return "待补充"
