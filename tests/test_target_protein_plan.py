from __future__ import annotations

import pytest

from app.core.opn_inputs import BuiltinOpnInputProvider
from pcsec_pichia.core.target_inputs import (
    LeaderCandidateInput,
    TargetProteinInput,
    draft_hlf_input_set,
    target_input_set_from_opn_input_set,
)
from pcsec_pichia.core.target_protein_plan import build_target_protein_plan, protein_extra_mw, protein_mw


def test_target_protein_plan_builds_opn_alpha_reference_reaction_plan() -> None:
    generic = target_input_set_from_opn_input_set(BuiltinOpnInputProvider().load_input_set())
    leader = next(item for item in generic.leaders if item.candidate_id == "OPN_ALPHA_FULL_PROJECT")

    plan = build_target_protein_plan(generic.target, leader)

    assert plan.target_id == "OPN"
    assert plan.protein_id == "OPN_ALPHA_FULL_PROJECT"
    assert plan.full_length == 387
    assert plan.signal_peptide_length == 19
    assert plan.localization == "e"
    assert plan.through_er is True
    assert plan.protein_extra_mw == pytest.approx(7560.0)
    assert plan.protein_mw == pytest.approx(protein_mw(plan.full_sequence))

    reaction_ids = set(plan.reaction_ids)
    assert "r_OPN_ALPHA_FULL_PROJECT_peptide_translation" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT_Post_translation_PSTA_sec_RAC_complex" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT_ERAD7B_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT_Mature" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT_HDSVII_sec_Vps1p_Chc1p_Clc1p_complex" in reaction_ids
    assert "r_OPN_ALPHA_FULL_PROJECT_SP_degradation" in reaction_ids
    assert "r_OPN_ALPHA_FULL_PROJECT_subunit_degradation" in reaction_ids
    assert "OPN_ALPHA_FULL_PROJECT exchange" in reaction_ids

    assert len(plan.reaction_ids_by_stage("translation")) == 1
    assert len(plan.reaction_ids_by_stage("translocation")) == 6
    assert len(plan.reaction_ids_by_stage("misfolding")) == 12
    assert len(plan.reaction_ids_by_stage("er_to_golgi")) == 4
    assert len(plan.reaction_ids_by_stage("final_transport")) == 2


def test_target_protein_plan_rejects_draft_hlf_until_parameters_are_confirmed() -> None:
    hlf = draft_hlf_input_set()
    leader = LeaderCandidateInput(
        candidate_id="HLF_TEST_LEADER",
        leader_sequence="MKFAISTLLIILQAAAVFAA",
        signal_peptide_sequence="MKFAISTLLIILQAAAVFAA",
    )

    with pytest.raises(ValueError, match="参数待确认"):
        build_target_protein_plan(hlf.target, leader)


def test_target_protein_plan_supports_future_ready_target_without_engine_changes() -> None:
    target = TargetProteinInput(
        target_id="TEST_PROTEIN",
        protein_name="test protein",
        abbreviation="TEST",
        mature_sequence="MSTNPKPQR",
        through_er=0,
        signal_peptide=0,
        localization="c",
        parameter_status="ready_for_model",
    )
    leader = LeaderCandidateInput(
        candidate_id="TEST_NO_LEADER",
        leader_sequence="M",
        signal_peptide_sequence="M",
    )

    plan = build_target_protein_plan(target, leader)

    assert plan.protein_id == "TEST_NO_LEADER"
    assert plan.full_sequence == "MMSTNPKPQR"
    assert plan.protein_extra_mw == pytest.approx(protein_extra_mw(0, 0, 0, 0))
    assert "TEST_NO_LEADER_folding_c" in plan.reaction_ids
    assert "r_TEST_NO_LEADER_subunit_degradation" in plan.reaction_ids
