from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Literal, Protocol


AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
TargetParameterStatus = Literal["draft", "ready_for_model", "validated"]

PCSEC_TARGET_MODEL_COLUMNS = [
    "Protein name",
    "abbreviation",
    "ThroughER",
    "Signal peptide ",
    "Disulfide site",
    "N-glycosylation site",
    "O-linked glycisylation ",
    "Transmembrane",
    "GPI site",
    "Localization",
    "sequence",
    "Length",
    "sp sequence",
    "Signal peptide length",
    "Cotranslation",
]


def clean_amino_acid_sequence(value: object) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


@dataclass(frozen=True)
class TargetProteinInput:
    target_id: str
    protein_name: str
    abbreviation: str
    mature_sequence: str
    through_er: int = 1
    signal_peptide: int = 1
    disulfide_sites: int = 0
    n_glycosylation_sites: int = 0
    o_glycosylation_sites: int = 0
    transmembrane: int = 0
    gpi_sites: int = 0
    localization: str = "e"
    cotranslation: int = 0
    parameter_status: TargetParameterStatus = "draft"
    missing_parameters: tuple[str, ...] = field(default_factory=tuple)
    evidence_note: str = ""

    def normalized(self) -> "TargetProteinInput":
        return TargetProteinInput(
            target_id=self.target_id.strip(),
            protein_name=self.protein_name.strip(),
            abbreviation=self.abbreviation.strip(),
            mature_sequence=clean_amino_acid_sequence(self.mature_sequence),
            through_er=self.through_er,
            signal_peptide=self.signal_peptide,
            disulfide_sites=self.disulfide_sites,
            n_glycosylation_sites=self.n_glycosylation_sites,
            o_glycosylation_sites=self.o_glycosylation_sites,
            transmembrane=self.transmembrane,
            gpi_sites=self.gpi_sites,
            localization=self.localization.strip() or "e",
            cotranslation=self.cotranslation,
            parameter_status=self.parameter_status,
            missing_parameters=self.missing_parameters,
            evidence_note=self.evidence_note,
        )

    def validation_errors(self) -> list[str]:
        target = self.normalized()
        errors: list[str] = []
        if not target.target_id:
            errors.append("target_id is required.")
        if not target.protein_name:
            errors.append("protein_name is required.")
        if not target.abbreviation:
            errors.append("abbreviation is required.")
        if not target.mature_sequence:
            errors.append("mature_sequence is required before model simulation.")
        elif not AA_PATTERN.fullmatch(target.mature_sequence):
            errors.append("mature_sequence must use standard amino-acid letters.")
        if target.missing_parameters:
            missing = ", ".join(target.missing_parameters)
            errors.append(f"{target.target_id} 参数待确认：{missing}")
        elif target.parameter_status == "draft":
            missing = ", ".join(target.missing_parameters) if target.missing_parameters else "model-ready parameters"
            errors.append(f"{target.target_id} 参数待确认：{missing}")
        return errors

    @property
    def ready_for_model(self) -> bool:
        return not self.validation_errors()

    def readiness_message(self) -> str:
        errors = self.validation_errors()
        if not errors:
            return "参数完整，可进入 pcSecPichia 模型。"
        return "；".join(errors)


@dataclass(frozen=True)
class LeaderCandidateInput:
    candidate_id: str
    leader_sequence: str
    signal_peptide_sequence: str
    category: str = ""
    processing_route: str = ""
    source_note: str = ""
    rationale: str = ""
    caution: str = ""

    def normalized(self) -> "LeaderCandidateInput":
        return LeaderCandidateInput(
            candidate_id=self.candidate_id.strip(),
            leader_sequence=clean_amino_acid_sequence(self.leader_sequence),
            signal_peptide_sequence=clean_amino_acid_sequence(self.signal_peptide_sequence),
            category=self.category.strip(),
            processing_route=self.processing_route.strip(),
            source_note=self.source_note.strip(),
            rationale=self.rationale.strip(),
            caution=self.caution.strip(),
        )

    def validation_errors(self) -> list[str]:
        candidate = self.normalized()
        errors: list[str] = []
        if not candidate.candidate_id:
            errors.append("candidate_id is required.")
        if not AA_PATTERN.fullmatch(candidate.leader_sequence):
            errors.append(f"{candidate.candidate_id or 'candidate'} leader_sequence must use standard amino-acid letters.")
        if not AA_PATTERN.fullmatch(candidate.signal_peptide_sequence):
            errors.append(
                f"{candidate.candidate_id or 'candidate'} signal_peptide_sequence must use standard amino-acid letters."
            )
        if (
            candidate.leader_sequence
            and candidate.signal_peptide_sequence
            and candidate.signal_peptide_sequence not in candidate.leader_sequence
        ):
            errors.append(
                f"{candidate.candidate_id or 'candidate'} signal_peptide_sequence must be contained in leader_sequence."
            )
        return errors


@dataclass(frozen=True)
class TargetInputSet:
    target: TargetProteinInput
    leaders: list[LeaderCandidateInput]
    source_name: str = "target input set"
    errors: list[str] = field(default_factory=list)

    def validation_errors(self) -> list[str]:
        errors = list(self.errors)
        errors.extend(self.target.validation_errors())
        seen: set[str] = set()
        for leader in self.leaders:
            normalized = leader.normalized()
            if normalized.candidate_id in seen:
                errors.append(f"duplicate leader candidate_id: {normalized.candidate_id}")
            seen.add(normalized.candidate_id)
            errors.extend(normalized.validation_errors())
        return errors

    @property
    def ready_for_model(self) -> bool:
        return not self.validation_errors()


class TargetInputProvider(Protocol):
    source_name: str

    def load_input_set(self) -> TargetInputSet:
        """Return target-protein and leader-candidate inputs."""


@dataclass(frozen=True)
class StaticTargetInputProvider:
    input_set: TargetInputSet
    source_name: str = "static target inputs"

    def load_input_set(self) -> TargetInputSet:
        return self.input_set


@dataclass
class TargetRegistry:
    _input_sets: dict[str, TargetInputSet] = field(default_factory=dict)

    def register(self, input_set: TargetInputSet) -> None:
        target_id = input_set.target.normalized().target_id
        if not target_id:
            raise ValueError("target_id is required.")
        self._input_sets[target_id] = input_set

    def list_targets(self) -> list[TargetProteinInput]:
        return [input_set.target for input_set in self._input_sets.values()]

    def get(self, target_id: str) -> TargetInputSet:
        try:
            return self._input_sets[target_id]
        except KeyError as exc:
            raise KeyError(f"未找到目标蛋白：{target_id}") from exc


def build_pcsec_target_row(target: TargetProteinInput, leader: LeaderCandidateInput) -> dict[str, object]:
    normalized_target = target.normalized()
    normalized_leader = leader.normalized()
    sequence = normalized_leader.leader_sequence + normalized_target.mature_sequence
    return {
        "Protein name": normalized_leader.candidate_id,
        "abbreviation": normalized_leader.candidate_id,
        "ThroughER": normalized_target.through_er,
        "Signal peptide ": normalized_target.signal_peptide,
        "Disulfide site": normalized_target.disulfide_sites,
        "N-glycosylation site": normalized_target.n_glycosylation_sites,
        "O-linked glycisylation ": normalized_target.o_glycosylation_sites,
        "Transmembrane": normalized_target.transmembrane,
        "GPI site": normalized_target.gpi_sites,
        "Localization": normalized_target.localization,
        "sequence": sequence,
        "Length": len(sequence),
        "sp sequence": normalized_leader.signal_peptide_sequence,
        "Signal peptide length": len(normalized_leader.signal_peptide_sequence),
        "Cotranslation": normalized_target.cotranslation,
    }


def target_input_set_from_opn_input_set(opn_input_set: object) -> TargetInputSet:
    target = getattr(opn_input_set, "target")
    leaders = getattr(opn_input_set, "candidates")
    generic_target = TargetProteinInput(
        target_id=str(getattr(target, "abbreviation", "OPN")),
        protein_name=str(getattr(target, "protein_name", "OPN")),
        abbreviation=str(getattr(target, "abbreviation", "OPN")),
        mature_sequence=str(getattr(target, "mature_sequence", "")),
        through_er=int(getattr(target, "through_er", 1)),
        signal_peptide=int(getattr(target, "signal_peptide", 1)),
        disulfide_sites=int(getattr(target, "disulfide_sites", 0)),
        n_glycosylation_sites=int(getattr(target, "n_glycosylation_sites", 0)),
        o_glycosylation_sites=int(getattr(target, "o_glycosylation_sites", 0)),
        transmembrane=int(getattr(target, "transmembrane", 0)),
        gpi_sites=int(getattr(target, "gpi_sites", 0)),
        localization=str(getattr(target, "localization", "e")),
        cotranslation=int(getattr(target, "cotranslation", 0)),
        parameter_status="ready_for_model",
        evidence_note="Converted from the existing OPN input provider.",
    )
    generic_leaders = [
        LeaderCandidateInput(
            candidate_id=str(getattr(leader, "candidate_id")),
            leader_sequence=str(getattr(leader, "leader_sequence")),
            signal_peptide_sequence=str(getattr(leader, "signal_peptide_sequence")),
            category=str(getattr(leader, "category", "")),
            processing_route=str(getattr(leader, "processing_route", "")),
            source_note=str(getattr(leader, "source_note", "")),
            rationale=str(getattr(leader, "rationale", "")),
            caution=str(getattr(leader, "caution", "")),
        )
        for leader in leaders
    ]
    return TargetInputSet(
        target=generic_target,
        leaders=generic_leaders,
        source_name=str(getattr(opn_input_set, "source_name", "converted OPN inputs")),
        errors=list(getattr(opn_input_set, "errors", [])),
    )


def draft_hlf_input_set() -> TargetInputSet:
    return TargetInputSet(
        target=TargetProteinInput(
            target_id="hLF",
            protein_name="human lactoferrin",
            abbreviation="hLF",
            mature_sequence="",
            through_er=1,
            signal_peptide=1,
            parameter_status="draft",
            missing_parameters=(
                "mature_sequence",
                "disulfide_sites",
                "n_glycosylation_sites",
                "o_glycosylation_sites",
                "PTM evidence source",
            ),
            evidence_note="hLF is reserved as a target draft until sequence and PTM parameters are confirmed.",
        ),
        leaders=[],
        source_name="builtin hLF draft",
    )
