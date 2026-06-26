from __future__ import annotations

import pytest

from pcsec_pichia.adapters.mat_loader import MatStructLoader
from pcsec_pichia.adapters.pichia_target_enzymedata import (
    simulate_target_reaction_coefficients,
    simulate_target_reaction_coefficients_matlab_compatible,
    target_enzymedata_from_plan,
)
from pcsec_pichia.adapters.aa_stoichiometry import AminoAcidStoichiometry
from app.core.opn_inputs import BuiltinOpnInputProvider
from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.core.target_inputs import target_input_set_from_opn_input_set
from pcsec_pichia.core.target_protein_plan import build_target_protein_plan, protein_mw
from pcsec_pichia.services._pichia_opn_builder import add_opn_target_reactions


def _opn_plan(candidate_id: str):
    generic = target_input_set_from_opn_input_set(BuiltinOpnInputProvider().load_input_set())
    leader = next(item for item in generic.leaders if item.candidate_id == candidate_id)
    return build_target_protein_plan(generic.target, leader)


def _opn_model(plan):
    paths = ProjectPaths.discover()
    loader = MatStructLoader(paths)
    model = loader.load_pcsec_pichia_model()
    amino_acids = AminoAcidStoichiometry.from_workbook(paths.pichia_aa_id_xlsx)
    return add_opn_target_reactions(model, plan, amino_acids)


def test_target_enzymedata_from_plan_matches_calculate_mw_fields_for_opn() -> None:
    plan = _opn_plan("OPN_ALPHA_FULL_PROJECT")

    enzymedata = target_enzymedata_from_plan(plan)

    assert enzymedata.proteins == ["OPN_ALPHA_FULL_PROJECT"]
    assert enzymedata.protein_count == 1
    assert enzymedata.kdeg[0] == pytest.approx(0.0)
    assert enzymedata.protein_length[0] == pytest.approx(387.0)
    assert enzymedata.protein_mw[0] == pytest.approx(protein_mw(plan.full_sequence))
    assert enzymedata.protein_pst_info == ("DSB", "NG", "OG", "GPI")
    assert enzymedata.protein_pst.tolist() == [[0.0, 0.0, 7.0, 0.0]]
    assert enzymedata.protein_extra_mw_specific.tolist() == [[0.0, 0.0, 7560.0, 0.0]]
    assert enzymedata.protein_extra_mw[0] == pytest.approx(7560.0)
    assert enzymedata.protein_loc == ["e"]
    assert enzymedata.pst_value("OPN_ALPHA_FULL_PROJECT", "OG") == pytest.approx(7.0)


def test_simulate_target_reaction_coefficients_matches_sec_coefref_for_current_opn_branch() -> None:
    paths = ProjectPaths.discover()
    plan = _opn_plan("OPN_ALPHA_FULL_PROJECT")
    model_result = _opn_model(plan)
    secretory_enzymedata = MatStructLoader(paths).load_secretory_enzymedata()
    target_enzymedata = target_enzymedata_from_plan(plan)

    with_coefficients = simulate_target_reaction_coefficients(
        model_result.model,
        secretory_enzymedata,
        target_enzymedata,
    )

    assert len(with_coefficients.reaction_coefficients) == 18
    assert with_coefficients.reaction_coefficients[
        "OPN_ALPHA_FULL_PROJECT_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex"
    ] == pytest.approx(7.0)
    assert with_coefficients.reaction_coefficients[
        "OPN_ALPHA_FULL_PROJECT_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex"
    ] == pytest.approx(7.0)
    assert with_coefficients.reaction_coefficients[
        "OPN_ALPHA_FULL_PROJECT_Post_translation_PSTA_sec_BIP_NEFS_complex"
    ] == pytest.approx(387.0 / 40.0)
    assert with_coefficients.reaction_coefficients[
        "OPN_ALPHA_FULL_PROJECT_cycle_accumulation_sec_acc_Kar2p_complex"
    ] == pytest.approx(10.0 * 387.0 / 40.0)
    assert with_coefficients.reaction_coefficients[
        "OPN_ALPHA_FULL_PROJECT_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex"
    ] == pytest.approx(plan.protein_mw / 54580.0)


def test_matlab_compatible_target_reaction_coefficients_match_third_underscore_behavior() -> None:
    paths = ProjectPaths.discover()
    loader = MatStructLoader(paths)
    secretory_enzymedata = loader.load_secretory_enzymedata()

    alpha_plan = _opn_plan("OPN_ALPHA_FULL_PROJECT")
    alpha_model = _opn_model(alpha_plan)
    alpha_full = simulate_target_reaction_coefficients(
        alpha_model.model,
        secretory_enzymedata,
        target_enzymedata_from_plan(alpha_plan),
    )
    alpha_matlab = simulate_target_reaction_coefficients_matlab_compatible(
        alpha_model.model,
        secretory_enzymedata,
        target_enzymedata_from_plan(alpha_plan),
    )

    dddk_plan = _opn_plan("OPN_PPA_DDDK18")
    dddk_model = _opn_model(dddk_plan)
    dddk_matlab = simulate_target_reaction_coefficients_matlab_compatible(
        dddk_model.model,
        secretory_enzymedata,
        target_enzymedata_from_plan(dddk_plan),
    )

    assert len(alpha_full.reaction_coefficients) == 18
    assert len(alpha_matlab.reaction_coefficients) == 0
    assert len(dddk_matlab.reaction_coefficients) == 18
    assert "OPN_PPA_DDDK18_misfold_ERAD_sec_Kar2p_complex" in dddk_matlab.reaction_coefficients


def test_target_enzymedata_merges_into_secretory_and_combined_enzymedata_for_lp_constraints() -> None:
    paths = ProjectPaths.discover()
    loader = MatStructLoader(paths)
    plan = _opn_plan("OPN_ALPHA_FULL_PROJECT")
    model_result = _opn_model(plan)
    secretory_enzymedata = loader.load_secretory_enzymedata()
    combined_enzymedata = loader.load_combined_enzymedata()
    target_enzymedata = simulate_target_reaction_coefficients(
        model_result.model,
        secretory_enzymedata,
        target_enzymedata_from_plan(plan),
    )

    secretory_with_target = secretory_enzymedata.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    combined_with_target = combined_enzymedata.with_target_proteins(target_enzymedata)

    assert len(secretory_with_target.reaction_coefficients) == 8242
    assert "OPN_ALPHA_FULL_PROJECT_misfold_ERAD_sec_Kar2p_complex" in secretory_with_target.reaction_coefficients
    assert combined_with_target.protein_count == 1418
    assert combined_with_target.exact_protein_length("OPN_ALPHA_FULL_PROJECT") == pytest.approx(387.0)
    assert combined_with_target.molecular_weight_for_dilution_reaction(
        "OPN_ALPHA_FULL_PROJECT_dilution_misfolding_er"
    ) == pytest.approx(plan.protein_mw)
