from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from pcsec_pichia.core.target_inputs import (
    AA_PATTERN,
    LeaderCandidateInput,
    StaticTargetInputProvider,
    TargetInputProvider,
    TargetInputSet,
    TargetProteinInput,
    TargetRegistry,
    build_pcsec_target_row,
    clean_amino_acid_sequence,
)
from pcsec_pichia.probe import (
    TargetSpec,
    clean_sequence,
    load_hlf_default,
    load_opn_default,
    load_targets,
    target_features,
    target_from_mapping,
)


OPN_BASELINE_TARGET_ID = "OPN_ALPHA_FULL_PROJECT"
OPN_DEFAULT_O_GLYCOSYLATION_SITES = 7


@dataclass(frozen=True)
class BuiltinTargetSummary:
    target_id: str
    protein_id: str
    source: str
    parameter_status: str


def list_supported_builtin_targets() -> tuple[BuiltinTargetSummary, ...]:
    return (
        *tuple(
            BuiltinTargetSummary(
                target_id=candidate_id,
                protein_id=candidate_id,
                source="Data/pcSecPichia/TargetProtein_OPN_candidates.csv",
                parameter_status="ready_for_model",
            )
            for candidate_id in _opn_candidate_ids(_repo_root())
        ),
        BuiltinTargetSummary(
            target_id="hLF",
            protein_id="hLF",
            source="用户提供: hLF native signal (19aa) + mature hLF (691aa)",
            parameter_status="draft_matlab_alignment_pending",
        ),
    )


def load_builtin_targets(root: Path | None = None) -> list[TargetSpec]:
    resolved_root = root or _repo_root()
    return [
        *load_opn_candidate_targets(resolved_root),
        load_hlf_default(resolved_root),
    ]


def load_opn_candidate_targets(root: Path | None = None) -> list[TargetSpec]:
    resolved_root = root or _repo_root()
    return [load_opn_candidate_target(candidate_id, resolved_root) for candidate_id in _opn_candidate_ids(resolved_root)]


def load_opn_candidate_target(candidate_id: str, root: Path | None = None) -> TargetSpec:
    resolved_root = root or _repo_root()
    if candidate_id == OPN_BASELINE_TARGET_ID:
        return load_opn_default(resolved_root)

    rows = _opn_candidate_rows(resolved_root)
    meta_rows = _opn_candidate_meta_rows(resolved_root)
    try:
        row = rows[candidate_id]
    except KeyError as exc:
        raise KeyError(f"Unknown OPN candidate target: {candidate_id}") from exc
    try:
        meta = meta_rows[candidate_id]
    except KeyError as exc:
        raise KeyError(f"Missing OPN candidate metadata: {candidate_id}") from exc

    full_sequence = clean_sequence(row.get("sequence", ""))
    leader_sequence = clean_sequence(meta.get("leader_sequence", ""))
    signal_peptide_sequence = clean_sequence(meta.get("signal_peptide_sequence", "") or row.get("sp sequence", ""))
    if not full_sequence.startswith(leader_sequence):
        raise ValueError(f"{candidate_id} sequence does not start with its OPN leader sequence.")

    return TargetSpec(
        target_id=candidate_id,
        protein_id=candidate_id,
        mature_sequence=full_sequence[len(leader_sequence) :],
        leader_sequence=leader_sequence,
        signal_peptide_sequence=signal_peptide_sequence,
        through_er=_parse_bool_int(row.get("ThroughER", "1")),
        localization=str(row.get("Localization", "e") or "e"),
        disulfide_sites=_parse_int(row.get("Disulfide site"), default=0),
        n_glycosylation_sites=_parse_int(row.get("N-glycosylation site"), default=0),
        o_glycosylation_sites=_parse_int(
            row.get("O-linked glycisylation "),
            default=OPN_DEFAULT_O_GLYCOSYLATION_SITES,
        ),
        transmembrane=_parse_int(row.get("Transmembrane"), default=0),
        gpi_sites=_parse_int(row.get("GPI site"), default=0),
        cotranslation=_parse_int(row.get("Cotranslation"), default=0),
        source="Data/pcSecPichia/TargetProtein_OPN_candidates.csv",
    )


def load_custom_targets_json(path: Path) -> list[TargetSpec]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    targets = payload.get("targets", [])
    if not isinstance(targets, list):
        raise ValueError("targets JSON must contain a list field named 'targets'.")
    return [target_spec_from_mapping(item, source=f"json:{path}") for item in targets]


def target_spec_from_input(target: TargetProteinInput, leader: LeaderCandidateInput) -> TargetSpec:
    target_errors = target.validation_errors()
    leader_errors = leader.validation_errors()
    if target_errors or leader_errors:
        raise ValueError("; ".join([*target_errors, *leader_errors]))

    normalized_target = target.normalized()
    normalized_leader = leader.normalized()
    return TargetSpec(
        target_id=normalized_target.target_id,
        protein_id=normalized_target.abbreviation,
        mature_sequence=normalized_target.mature_sequence,
        leader_sequence=normalized_leader.leader_sequence,
        signal_peptide_sequence=normalized_leader.signal_peptide_sequence,
        through_er=bool(normalized_target.through_er),
        localization=normalized_target.localization,
        disulfide_sites=normalized_target.disulfide_sites,
        n_glycosylation_sites=normalized_target.n_glycosylation_sites,
        o_glycosylation_sites=normalized_target.o_glycosylation_sites,
        transmembrane=normalized_target.transmembrane,
        gpi_sites=normalized_target.gpi_sites,
        cotranslation=normalized_target.cotranslation,
        source="TargetProteinInput",
    )


def target_spec_from_mapping(item: dict[str, object], source: str) -> TargetSpec:
    _validate_mapping_sequence(item, "mature_sequence")
    _validate_mapping_sequence(item, "leader_sequence", required=False)
    _validate_mapping_sequence(item, "signal_peptide_sequence", required=False)
    return target_from_mapping(item, source)


def _validate_mapping_sequence(item: dict[str, object], field_name: str, required: bool = True) -> None:
    value = item.get(field_name, "")
    cleaned = clean_amino_acid_sequence(value)
    if not cleaned:
        if required:
            raise ValueError(f"{field_name} is required.")
        return
    if not AA_PATTERN.fullmatch(cleaned):
        target_id = item.get("target_id") or item.get("protein_id") or "target"
        raise ValueError(f"{target_id} {field_name} must use standard amino-acid letters.")


def _opn_candidate_ids(root: Path) -> tuple[str, ...]:
    return tuple(_opn_candidate_rows(root).keys())


def _opn_candidate_rows(root: Path) -> dict[str, dict[str, str]]:
    csv_path = root / "Data" / "pcSecPichia" / "TargetProtein_OPN_candidates.csv"
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        return {str(row["abbreviation"]): row for row in csv.DictReader(handle)}


def _opn_candidate_meta_rows(root: Path) -> dict[str, dict[str, str]]:
    csv_path = root / "Data" / "pcSecPichia" / "TargetProtein_OPN_candidates_meta.csv"
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        return {str(row["candidate_id"]): row for row in csv.DictReader(handle)}


def _parse_int(value: object, default: int) -> int:
    if value is None:
        return default
    text = str(value).strip()
    if not text:
        return default
    return int(float(text))


def _parse_bool_int(value: object) -> bool:
    return bool(_parse_int(value, default=0))


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


__all__ = [
    "BuiltinTargetSummary",
    "LeaderCandidateInput",
    "StaticTargetInputProvider",
    "TargetInputProvider",
    "TargetInputSet",
    "TargetProteinInput",
    "TargetRegistry",
    "TargetSpec",
    "build_pcsec_target_row",
    "clean_sequence",
    "list_supported_builtin_targets",
    "load_builtin_targets",
    "load_custom_targets_json",
    "load_hlf_default",
    "load_opn_candidate_target",
    "load_opn_candidate_targets",
    "load_opn_default",
    "load_targets",
    "target_features",
    "target_spec_from_input",
    "target_spec_from_mapping",
    "target_from_mapping",
]
