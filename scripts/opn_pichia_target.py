"""Prepare OPN target protein inputs for pcSecPichia.

The user supplied the mature human osteopontin (OPN/SPP1) sequence with the
native signal peptide removed. This module keeps the mature sequence fixed and
builds model-ready secretory leader candidates that can be read by
``addTargetProtein`` through the existing fakeProteinInfo path.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from dataclasses import asdict
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from app.core.opn_inputs import (
    BuiltinOpnInputProvider,
    DEFAULT_MATURE_OPN,
    HUMAN_SPP1_NATIVE_SP,
    OpnInputProvider,
    OpnInputSet,
    OpnLeaderCandidateInput,
    OpnTargetProteinInput,
    PPA_DDDK18_SP,
    PPA_EPX1_SA_SP,
    PPA_PAS_CHR3_0030_SP,
    PROJECT_ALPHA_FACTOR_LEADER,
    PROJECT_ALPHA_FACTOR_PRO,
    PROJECT_ALPHA_FACTOR_SP,
    SC_OST1_N23,
)


MATURE_OPN = DEFAULT_MATURE_OPN
SignalPeptideCandidate = OpnLeaderCandidateInput

MODEL_COLUMNS = [
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

AA_MW = {
    "A": 71.08,
    "B": 114.60,
    "C": 103.14,
    "D": 115.09,
    "E": 129.11,
    "F": 147.17,
    "G": 57.05,
    "H": 137.14,
    "I": 113.16,
    "J": 113.16,
    "K": 128.17,
    "L": 113.16,
    "M": 131.20,
    "N": 114.10,
    "O": 255.31,
    "P": 97.12,
    "Q": 128.13,
    "R": 156.19,
    "S": 87.08,
    "T": 101.10,
    "U": 150.04,
    "V": 99.13,
    "W": 186.21,
    "X": 126.50,
    "Y": 163.17,
    "Z": 128.62,
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def target_csv_path(root: Path | None = None) -> Path:
    base = root or repo_root()
    return base / "Data" / "pcSecPichia" / "TargetProtein_OPN.csv"


def candidate_csv_path(root: Path | None = None) -> Path:
    base = root or repo_root()
    return base / "Data" / "pcSecPichia" / "TargetProtein_OPN_candidates.csv"


def candidate_metadata_path(root: Path | None = None) -> Path:
    base = root or repo_root()
    return base / "Data" / "pcSecPichia" / "TargetProtein_OPN_candidates_meta.csv"


def opn_construct_sequence(input_provider: OpnInputProvider | None = None) -> str:
    input_set = resolve_input_set(input_provider)
    return PROJECT_ALPHA_FACTOR_LEADER + input_set.target.mature_sequence


def candidate_definitions() -> list[SignalPeptideCandidate]:
    return resolve_input_set().candidates


def resolve_input_set(input_provider: OpnInputProvider | None = None) -> OpnInputSet:
    provider = input_provider or BuiltinOpnInputProvider()
    input_set = provider.load_input_set()
    if input_set.errors:
        raise ValueError("; ".join(input_set.errors))
    return input_set


def build_model_row(
    candidate: SignalPeptideCandidate,
    target: OpnTargetProteinInput | None = None,
) -> dict[str, object]:
    target = target or OpnTargetProteinInput()
    sequence = candidate.leader_sequence + target.mature_sequence
    return {
        "Protein name": candidate.candidate_id,
        "abbreviation": candidate.candidate_id,
        "ThroughER": target.through_er,
        "Signal peptide ": target.signal_peptide,
        "Disulfide site": target.disulfide_sites,
        "N-glycosylation site": target.n_glycosylation_sites,
        "O-linked glycisylation ": target.o_glycosylation_sites,
        "Transmembrane": target.transmembrane,
        "GPI site": target.gpi_sites,
        "Localization": target.localization,
        "sequence": sequence,
        "Length": len(sequence),
        "sp sequence": candidate.signal_peptide_sequence,
        "Signal peptide length": len(candidate.signal_peptide_sequence),
        "Cotranslation": target.cotranslation,
    }


def build_candidate_rows(input_provider: OpnInputProvider | None = None) -> list[dict[str, object]]:
    input_set = resolve_input_set(input_provider)
    return [build_model_row(candidate, input_set.target) for candidate in input_set.candidates]


def build_candidate_metadata(input_provider: OpnInputProvider | None = None) -> list[dict[str, object]]:
    input_set = resolve_input_set(input_provider)
    rows = []
    for candidate in input_set.candidates:
        model_row = build_model_row(candidate, input_set.target)
        rows.append(
            {
                **asdict(candidate),
                "leader_length": len(candidate.leader_sequence),
                "construct_length": model_row["Length"],
                "mature_opn_internal_kex2_like_sites": dibasic_kex2_like_sites(input_set.target.mature_sequence),
                "construct_nxs_t_motifs": n_glycosylation_motifs(str(model_row["sequence"])),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_candidate_files(
    root: Path | None = None,
    input_provider: OpnInputProvider | None = None,
) -> None:
    base = root or repo_root()
    write_csv(candidate_csv_path(base), build_candidate_rows(input_provider), MODEL_COLUMNS)
    metadata = build_candidate_metadata(input_provider)
    write_csv(candidate_metadata_path(base), metadata, list(metadata[0]))


def n_glycosylation_motifs(seq: str) -> list[tuple[int, str]]:
    return [(match.start() + 1, match.group(0)) for match in re.finditer(r"N[^P][ST]", seq)]


def dibasic_kex2_like_sites(seq: str) -> list[tuple[int, str]]:
    return [(match.start() + 1, match.group(0)) for match in re.finditer(r"KR|RR", seq)]


def protein_mw(seq: str) -> float:
    counts = Counter(seq)
    return 18.0 + sum(counts[aa] * mw for aa, mw in AA_MW.items())


def load_opn_row(path: Path | None = None) -> dict[str, str]:
    csv_path = path or target_csv_path()
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one OPN row in {csv_path}, found {len(rows)}.")
    return rows[0]


def load_candidate_rows(path: Path | None = None) -> list[dict[str, str]]:
    csv_path = path or candidate_csv_path()
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No OPN candidate rows found in {csv_path}.")
    return rows


def build_summary(
    row: dict[str, str],
    target: OpnTargetProteinInput | None = None,
) -> dict[str, object]:
    target = target or OpnTargetProteinInput()
    seq = row["sequence"]
    mature_motifs = n_glycosylation_motifs(target.mature_sequence)
    construct_motifs = n_glycosylation_motifs(seq)
    return {
        "target_id": row["Protein name"],
        "protein_id": row["abbreviation"],
        "mature_opn_length": len(target.mature_sequence),
        "leader_length": len(seq) - len(target.mature_sequence),
        "construct_length": len(seq),
        "model_length_field": int(row["Length"]),
        "signal_peptide_flag": int(row["Signal peptide "]),
        "signal_peptide_sequence": row["sp sequence"],
        "signal_peptide_length": int(row["Signal peptide length"]),
        "through_er": int(row["ThroughER"]),
        "localization": row["Localization"],
        "disulfide_sites": int(row["Disulfide site"]),
        "n_glycosylation_sites_for_model": int(row["N-glycosylation site"]),
        "o_glycosylation_sites_for_model": int(row["O-linked glycisylation "]),
        "transmembrane": int(row["Transmembrane"]),
        "gpi_sites": int(row["GPI site"]),
        "cotranslation": int(row["Cotranslation"]),
        "mature_cysteine_count": target.mature_sequence.count("C"),
        "mature_ser_thr_count": target.mature_sequence.count("S") + target.mature_sequence.count("T"),
        "mature_nxs_t_motifs": mature_motifs,
        "construct_nxs_t_motifs": construct_motifs,
        "mature_opn_internal_kex2_like_sites": dibasic_kex2_like_sites(target.mature_sequence),
        "construct_mw_da_model_formula": round(protein_mw(seq), 2),
        "extra_mw_da_default_pp": int(row["N-glycosylation site"]) * 3346
        + int(row["O-linked glycisylation "]) * 1080
        + int(row["GPI site"]) * 2009,
    }


def validate_model_row(
    row: dict[str, object | str],
    target: OpnTargetProteinInput | None = None,
) -> None:
    target = target or OpnTargetProteinInput()
    seq = str(row["sequence"])
    if not seq.endswith(target.mature_sequence):
        raise ValueError("OPN construct does not end with the mature OPN sequence.")
    if int(row["Length"]) != len(seq):
        raise ValueError("OPN Length field does not match sequence length.")
    if int(row["Signal peptide "]) != 1:
        raise ValueError("OPN candidate must be marked as a signal-peptide secretory protein.")
    if not str(row["sp sequence"]):
        raise ValueError("OPN candidate is missing the signal peptide sequence.")
    if int(row["Signal peptide length"]) != len(str(row["sp sequence"])):
        raise ValueError("Signal peptide length does not match the signal-peptide sequence.")
    if int(row["Disulfide site"]) != target.mature_sequence.count("C") // 2:
        raise ValueError("OPN disulfide count is inconsistent with mature cysteine count.")


def validate_row(
    row: dict[str, str],
    target: OpnTargetProteinInput | None = None,
) -> None:
    target = target or OpnTargetProteinInput()
    expected = PROJECT_ALPHA_FACTOR_LEADER + target.mature_sequence
    validate_model_row(row, target)
    if row["sequence"] != expected:
        raise ValueError("OPN construct sequence does not match the expected project alpha-factor leader + mature OPN.")


def validate_candidate_rows(
    rows: list[dict[str, object | str]],
    target: OpnTargetProteinInput | None = None,
) -> None:
    target = target or OpnTargetProteinInput()
    seen: set[str] = set()
    for row in rows:
        candidate_id = str(row["Protein name"])
        if candidate_id in seen:
            raise ValueError(f"Duplicate OPN candidate id: {candidate_id}")
        seen.add(candidate_id)
        if row["Protein name"] != row["abbreviation"]:
            raise ValueError(f"Protein name and abbreviation must match for fakeProteinInfo: {candidate_id}")
        validate_model_row(row, target)


def write_summary(summary: dict[str, object], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "OPN_target_summary.json"
    md_path = output_dir / "OPN_target_summary.md"
    json_path.write_text(json.dumps(summary, indent=2, ensure_ascii=True), encoding="utf-8")
    lines = [
        "# OPN pcSecPichia target summary",
        "",
        f"- Target id: `{summary['target_id']}`",
        f"- Protein id: `{summary['protein_id']}`",
        f"- Mature OPN length: {summary['mature_opn_length']} aa",
        f"- Leader length: {summary['leader_length']} aa",
        f"- Full construct length used by the model: {summary['construct_length']} aa",
        f"- Signal peptide flag: {summary['signal_peptide_flag']}",
        f"- Secretory localization: `{summary['localization']}`",
        f"- Disulfide sites: {summary['disulfide_sites']}",
        f"- N-glycosylation sites used by model: {summary['n_glycosylation_sites_for_model']}",
        f"- O-glycosylation sites used by model: {summary['o_glycosylation_sites_for_model']}",
        f"- Mature N-X-S/T motifs observed but not modeled as N-glycosylation: {summary['mature_nxs_t_motifs']}",
        f"- Mature Kex2-like dibasic motifs: {summary['mature_opn_internal_kex2_like_sites']}",
        f"- Approximate construct MW by project formula: {summary['construct_mw_da_model_formula']} Da",
        f"- Default Pichia extra PTM mass used by model: {summary['extra_mw_da_default_pp']} Da",
    ]
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=target_csv_path())
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=repo_root() / "local_runs" / "OPN_PPA_glc_smoke",
    )
    parser.add_argument("--write-candidates", action="store_true")
    args = parser.parse_args()

    if args.write_candidates:
        rows = build_candidate_rows()
        validate_candidate_rows(rows)
        write_candidate_files()
        print(f"Wrote {candidate_csv_path()}")
        print(f"Wrote {candidate_metadata_path()}")

    row = load_opn_row(args.csv)
    validate_row(row)
    summary = build_summary(row)
    write_summary(summary, args.summary_dir)
    print(json.dumps(summary, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
