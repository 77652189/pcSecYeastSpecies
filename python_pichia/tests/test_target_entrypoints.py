from __future__ import annotations

import json
from pathlib import Path

import pytest

from pcsec_pichia.core.target_inputs import LeaderCandidateInput, TargetProteinInput
from pcsec_pichia.targets import (
    list_supported_builtin_targets,
    load_builtin_targets,
    load_custom_targets_json,
    load_opn_candidate_target,
    load_opn_candidate_targets,
    target_spec_from_input,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
OPN_CANDIDATE_IDS = {
    "OPN_ALPHA_FULL_PROJECT",
    "OPN_ALPHA_PRE_ONLY",
    "OPN_NATIVE_SPP1",
    "OPN_OST1N23_ALPHA_PRO",
    "OPN_PPA_DDDK18",
    "OPN_PPA_PASCHR3_0030",
    "OPN_PPA_EPX1_SA",
}


def test_builtin_target_catalog_lists_opn_and_hlf() -> None:
    catalog = list_supported_builtin_targets()
    by_id = {item.target_id: item for item in catalog}

    assert OPN_CANDIDATE_IDS.issubset(by_id)
    assert by_id["OPN_ALPHA_FULL_PROJECT"].parameter_status == "ready_for_model"
    assert by_id["OPN_PPA_DDDK18"].parameter_status == "ready_for_model"
    assert by_id["OPN_PPA_PASCHR3_0030"].parameter_status == "ready_for_model"
    assert "hLF" in by_id
    assert by_id["hLF"].parameter_status == "draft_matlab_alignment_pending"
    assert "hLF_CLEAN" not in by_id
    assert "hLF_MATURE_SECRETED" not in by_id


def test_load_builtin_targets_returns_opn_and_hlf_specs() -> None:
    targets = load_builtin_targets(REPO_ROOT)
    by_id = {target.target_id: target for target in targets}

    assert OPN_CANDIDATE_IDS.issubset(by_id)
    assert by_id["OPN_ALPHA_FULL_PROJECT"].o_glycosylation_sites == 7
    assert by_id["OPN_ALPHA_FULL_PROJECT"].localization == "e"
    assert len(by_id["OPN_ALPHA_FULL_PROJECT"].full_sequence) == 383
    assert by_id["OPN_PPA_DDDK18"].o_glycosylation_sites == 7
    assert len(by_id["OPN_PPA_DDDK18"].leader_sequence) == 18
    assert len(by_id["OPN_PPA_DDDK18"].full_sequence) == 316
    assert by_id["OPN_PPA_PASCHR3_0030"].o_glycosylation_sites == 7
    assert len(by_id["OPN_PPA_PASCHR3_0030"].leader_sequence) == 20
    assert len(by_id["OPN_PPA_PASCHR3_0030"].full_sequence) == 318
    assert by_id["hLF"].disulfide_sites == 21
    assert by_id["hLF"].n_glycosylation_sites == 4
    assert by_id["hLF"].o_glycosylation_sites == 0
    assert by_id["hLF"].leader_sequence == "MKLVFLVLLFLGALGLCLA"
    assert by_id["hLF"].signal_peptide_sequence == "MKLVFLVLLFLGALGLCLA"
    assert len(by_id["hLF"].mature_sequence) == 691
    assert len(by_id["hLF"].full_sequence) == 710
    assert "用户提供" in by_id["hLF"].source


def test_load_opn_candidate_targets_returns_all_candidate_specs() -> None:
    targets = load_opn_candidate_targets(REPO_ROOT)
    by_id = {target.target_id: target for target in targets}

    assert set(by_id) == OPN_CANDIDATE_IDS
    for target_id, target in by_id.items():
        assert target.protein_id == target_id
        assert target.through_er is True
        assert target.localization == "e"
        assert target.disulfide_sites == 0
        assert target.n_glycosylation_sites == 0
        assert target.o_glycosylation_sites == 7
        assert target.transmembrane == 0
        assert target.gpi_sites == 0
        assert target.cotranslation == 0
        assert target.full_sequence == target.leader_sequence + target.mature_sequence


def test_load_opn_candidate_target_preserves_baseline_compatibility() -> None:
    baseline = load_opn_candidate_target("OPN_ALPHA_FULL_PROJECT", REPO_ROOT)

    assert len(baseline.leader_sequence) == 85
    assert len(baseline.mature_sequence) == 298
    assert len(baseline.full_sequence) == 383
    assert baseline.full_sequence.startswith("MRFPSIFTAVLFAASSALA")


def test_load_custom_targets_json_returns_target_specs() -> None:
    targets = load_custom_targets_json(REPO_ROOT / "local_runs" / "pichia_hlf_opn_probe" / "targets.example.json")
    by_id = {target.target_id: target for target in targets}

    assert by_id["OPN_CUSTOM"].o_glycosylation_sites == 7
    assert by_id["HLF_CUSTOM"].disulfide_sites == 21
    assert by_id["HLF_CUSTOM"].n_glycosylation_sites == 4


def test_custom_targets_json_rejects_nonstandard_amino_acids(tmp_path: Path) -> None:
    path = tmp_path / "bad_targets.json"
    path.write_text(
        json.dumps(
            {
                "targets": [
                    {
                        "target_id": "BAD",
                        "protein_id": "BAD",
                        "mature_sequence": "ACDUZ",
                        "leader_sequence": "MAGA",
                        "signal_peptide_sequence": "MA",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="standard amino-acid"):
        load_custom_targets_json(path)


def test_target_spec_from_input_builds_custom_sequence_entry() -> None:
    target = TargetProteinInput(
        target_id="CUSTOM",
        protein_name="custom protein",
        abbreviation="CUSTOM",
        mature_sequence="acdefg",
        disulfide_sites=1,
        n_glycosylation_sites=2,
        o_glycosylation_sites=3,
        parameter_status="ready_for_model",
    )
    leader = LeaderCandidateInput(
        candidate_id="CUSTOM_LEADER",
        leader_sequence="mmmaaa",
        signal_peptide_sequence="mmm",
    )

    spec = target_spec_from_input(target, leader)

    assert spec.target_id == "CUSTOM"
    assert spec.protein_id == "CUSTOM"
    assert spec.mature_sequence == "ACDEFG"
    assert spec.leader_sequence == "MMMAAA"
    assert spec.signal_peptide_sequence == "MMM"
    assert spec.full_sequence == "MMMAAAACDEFG"
    assert spec.disulfide_sites == 1
    assert spec.n_glycosylation_sites == 2
    assert spec.o_glycosylation_sites == 3


def test_target_spec_from_input_rejects_draft_parameters() -> None:
    target = TargetProteinInput(
        target_id="DRAFT",
        protein_name="draft protein",
        abbreviation="DRAFT",
        mature_sequence="ACDE",
    )
    leader = LeaderCandidateInput(
        candidate_id="LEADER",
        leader_sequence="MMAA",
        signal_peptide_sequence="MM",
    )

    with pytest.raises(ValueError, match="参数待确认"):
        target_spec_from_input(target, leader)
