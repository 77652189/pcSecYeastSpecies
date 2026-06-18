from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol


AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

DEFAULT_MATURE_OPN = (
    "IPVKQADSGSSEEKQLYNKYPDAVATWLNPDPSQKQNLLAPQNAVSSEETNDFKQETLPSKSNES"
    "HDHMDDMDDEDDDDHVDSQDSIDSNDSDDVDDTDDSHQSDESHHSDESDELVTDFPTDLPATEVFT"
    "PVVPTVDTYDGRGDSVVYGLRSKSKKFRRPDIQYPDATDEDITSHMESEELNGAYKAIPVAQDLN"
    "APSDWDSRGKDSYETSQLDDQSAETHSHKQSRLYKRKANDESNEHSDVIDSQELSKVSREFHSHEF"
    "HSHEDMLVVDPKSKEEDKHLKFRISHELDSASSEVN"
)

PROJECT_ALPHA_FACTOR_LEADER = (
    "MRFPSIFTAVLFAASSALAAPVNTTTEDETAQIPAEAVIGYSDLEGDFDVAVLPFSNSTNNGLLFI"
    "NTTIASIAAKEEGVSLEKREAEA"
)
PROJECT_ALPHA_FACTOR_SP = "MRFPSIFTAVLFAASSALA"
PROJECT_ALPHA_FACTOR_PRO = PROJECT_ALPHA_FACTOR_LEADER[len(PROJECT_ALPHA_FACTOR_SP) :]

HUMAN_SPP1_NATIVE_SP = "MRIAVICFCLLGITCA"
SC_OST1_N23 = "MRQVWFSWIVGLFLCFFNVSSAA"
PPA_DDDK18_SP = "MFNLKTILISTLASIAVA"
PPA_PAS_CHR3_0030_SP = "MKFAISTLLIILQAAAVFAA"
PPA_EPX1_SA_SP = "MKLSTNLILAIAAASAVVSA"

CANDIDATE_INPUT_COLUMNS = [
    "candidate_id",
    "leader_sequence",
    "signal_peptide_sequence",
    "category",
    "processing_route",
    "source_note",
    "rationale",
    "caution",
]


@dataclass(frozen=True)
class OpnTargetProteinInput:
    protein_name: str = "OPN"
    abbreviation: str = "OPN"
    mature_sequence: str = DEFAULT_MATURE_OPN
    through_er: int = 1
    signal_peptide: int = 1
    disulfide_sites: int = 0
    n_glycosylation_sites: int = 0
    o_glycosylation_sites: int = 7
    transmembrane: int = 0
    gpi_sites: int = 0
    localization: str = "e"
    cotranslation: int = 0


@dataclass(frozen=True)
class OpnLeaderCandidateInput:
    candidate_id: str
    leader_sequence: str
    signal_peptide_sequence: str
    category: str
    processing_route: str
    source_note: str
    rationale: str
    caution: str


@dataclass(frozen=True)
class OpnInputSet:
    target: OpnTargetProteinInput
    candidates: list[OpnLeaderCandidateInput]
    source_name: str = "builtin OPN inputs"
    errors: list[str] = field(default_factory=list)


class OpnInputProvider(Protocol):
    source_name: str

    def load_input_set(self) -> OpnInputSet:
        """Return target-protein and leader-candidate inputs for pcSecPichia."""


@dataclass(frozen=True)
class StaticOpnInputProvider:
    input_set: OpnInputSet
    source_name: str = "static OPN inputs"

    def load_input_set(self) -> OpnInputSet:
        return self.input_set


@dataclass(frozen=True)
class BuiltinOpnInputProvider:
    source_name: str = "builtin OPN inputs"

    def load_input_set(self) -> OpnInputSet:
        return default_opn_input_set(source_name=self.source_name)


@dataclass(frozen=True)
class CsvOpnCandidateInputProvider:
    path: Path | None = None
    content: bytes | str | None = None
    target: OpnTargetProteinInput = field(default_factory=OpnTargetProteinInput)
    source_name: str = "CSV OPN leader candidates"

    def load_input_set(self) -> OpnInputSet:
        text, read_errors = self._read_text()
        if read_errors:
            return OpnInputSet(self.target, [], source_name=self.source_name, errors=read_errors)
        reader = csv.DictReader(io.StringIO(text))
        missing = [column for column in CANDIDATE_INPUT_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            return OpnInputSet(
                self.target,
                [],
                source_name=self.source_name,
                errors=[f"Missing required OPN candidate columns: {', '.join(missing)}"],
            )
        candidates: list[OpnLeaderCandidateInput] = []
        errors: list[str] = []
        seen_ids: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            candidate, row_errors = _candidate_from_row(row, row_number, seen_ids)
            if row_errors:
                errors.extend(row_errors)
                continue
            if candidate is not None:
                seen_ids.add(candidate.candidate_id)
                candidates.append(candidate)
        return OpnInputSet(self.target, candidates, source_name=self.source_name, errors=errors)

    def _read_text(self) -> tuple[str, list[str]]:
        if self.content is not None:
            if isinstance(self.content, bytes):
                return self.content.decode("utf-8-sig"), []
            return str(self.content), []
        if self.path is None:
            return "", ["No OPN candidate input file or content was provided."]
        try:
            return self.path.read_text(encoding="utf-8-sig"), []
        except OSError as exc:
            return "", [f"Failed to read OPN candidate input file: {exc}"]


def default_opn_input_set(source_name: str = "builtin OPN inputs") -> OpnInputSet:
    return OpnInputSet(
        target=OpnTargetProteinInput(),
        candidates=[
            OpnLeaderCandidateInput(
                candidate_id="OPN_ALPHA_FULL_PROJECT",
                leader_sequence=PROJECT_ALPHA_FACTOR_LEADER,
                signal_peptide_sequence=PROJECT_ALPHA_FACTOR_SP,
                category="project_baseline",
                processing_route="alpha-factor prepro; signal peptidase plus Kex2/Ste13-like pro-leader processing",
                source_note="Project pcSecPichia alpha-factor leader sequence used by the existing target-protein workflow.",
                rationale="Best baseline because it matches the current project modeling convention and common Pichia secretion practice.",
                caution="Mature OPN contains internal dibasic Kex2-like motifs, so wet-lab constructs should check proteolytic clipping risk.",
            ),
            OpnLeaderCandidateInput(
                candidate_id="OPN_ALPHA_PRE_ONLY",
                leader_sequence=PROJECT_ALPHA_FACTOR_SP,
                signal_peptide_sequence=PROJECT_ALPHA_FACTOR_SP,
                category="yeast_signal_only",
                processing_route="signal peptidase only",
                source_note="Signal peptide portion of the project alpha-factor leader, without the pro region.",
                rationale="Separates the cost of the alpha pre signal from the cost and processing risk of the alpha pro region.",
                caution="May secrete less efficiently than the full alpha-factor prepro leader; this is a comparison arm, not a recommended final construct.",
            ),
            OpnLeaderCandidateInput(
                candidate_id="OPN_NATIVE_SPP1",
                leader_sequence=HUMAN_SPP1_NATIVE_SP,
                signal_peptide_sequence=HUMAN_SPP1_NATIVE_SP,
                category="target_native_signal",
                processing_route="signal peptidase only",
                source_note="Human SPP1/osteopontin native N-terminal signal peptide, UniProt P10451 residues 1-16.",
                rationale="Tests whether the target's native mammalian signal peptide is a viable low-length secretion signal in the model.",
                caution="Mammalian native signal peptides are not automatically optimal in Pichia; use this mainly as a biological reference.",
            ),
            OpnLeaderCandidateInput(
                candidate_id="OPN_OST1N23_ALPHA_PRO",
                leader_sequence=SC_OST1_N23 + PROJECT_ALPHA_FACTOR_PRO,
                signal_peptide_sequence=SC_OST1_N23,
                category="hybrid_yeast_leader",
                processing_route="Ost1 N-terminal signal peptide plus alpha-factor pro region",
                source_note="S. cerevisiae OST1 N-terminal pre sequence from UniProt P41543 residues 1-23 combined with the project alpha pro region.",
                rationale="Hybrid Ost1-alpha leaders are a common secretion-engineering comparison against standard alpha-factor leaders.",
                caution="Still uses the alpha pro processing route, so it does not remove Kex2-like clipping concerns.",
            ),
            OpnLeaderCandidateInput(
                candidate_id="OPN_PPA_DDDK18",
                leader_sequence=PPA_DDDK18_SP,
                signal_peptide_sequence=PPA_DDDK18_SP,
                category="pichia_native_signal",
                processing_route="signal peptidase only",
                source_note="Reported Pichia DDDK 18-aa signal peptide candidate.",
                rationale="Useful for OPN because it avoids alpha pro/Kex2 processing while staying in a Pichia-derived signal family.",
                caution="Candidate sequence should be experimentally confirmed in the final strain and vector context.",
            ),
            OpnLeaderCandidateInput(
                candidate_id="OPN_PPA_PASCHR3_0030",
                leader_sequence=PPA_PAS_CHR3_0030_SP,
                signal_peptide_sequence=PPA_PAS_CHR3_0030_SP,
                category="pichia_native_signal",
                processing_route="signal peptidase only",
                source_note="Pichia pastoris PAS_chr3_0030 signal peptide reported in secretion-leader screening resources.",
                rationale="A short Pichia-native leader candidate for comparing secretion cost without a pro peptide.",
                caution="Source protein and industrial performance may be product-dependent; treat as a screen candidate.",
            ),
            OpnLeaderCandidateInput(
                candidate_id="OPN_PPA_EPX1_SA",
                leader_sequence=PPA_EPX1_SA_SP,
                signal_peptide_sequence=PPA_EPX1_SA_SP,
                category="pichia_native_signal",
                processing_route="signal peptidase only",
                source_note="Pichia EPX1 signal-anchor/signal-peptide fragment reported in secretion-leader screening resources.",
                rationale="Adds a second Pichia-native short leader with different hydrophobic-core composition.",
                caution="Needs wet-lab confirmation because signal-anchor behavior can depend strongly on the downstream protein.",
            ),
        ],
        source_name=source_name,
    )


def clean_amino_acid_sequence(value: object) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


def _candidate_from_row(
    row: dict[str, str],
    row_number: int,
    seen_ids: set[str],
) -> tuple[OpnLeaderCandidateInput | None, list[str]]:
    candidate_id = str(row.get("candidate_id", "")).strip()
    leader = clean_amino_acid_sequence(row.get("leader_sequence", ""))
    signal = clean_amino_acid_sequence(row.get("signal_peptide_sequence", ""))
    errors: list[str] = []
    if not candidate_id:
        errors.append(f"Row {row_number}: candidate_id is required.")
    if candidate_id in seen_ids:
        errors.append(f"Row {row_number}: duplicate candidate_id in input file: {candidate_id}")
    if not AA_PATTERN.fullmatch(leader):
        errors.append(f"Row {row_number}: leader_sequence must use standard amino-acid letters.")
    if not AA_PATTERN.fullmatch(signal):
        errors.append(f"Row {row_number}: signal_peptide_sequence must use standard amino-acid letters.")
    if leader and signal and signal not in leader:
        errors.append(f"Row {row_number}: signal_peptide_sequence must be contained in leader_sequence.")
    if errors:
        return None, errors
    return (
        OpnLeaderCandidateInput(
            candidate_id=candidate_id,
            leader_sequence=leader,
            signal_peptide_sequence=signal,
            category=str(row.get("category", "")).strip(),
            processing_route=str(row.get("processing_route", "")).strip(),
            source_note=str(row.get("source_note", "")).strip(),
            rationale=str(row.get("rationale", "")).strip(),
            caution=str(row.get("caution", "")).strip(),
        ),
        [],
    )
