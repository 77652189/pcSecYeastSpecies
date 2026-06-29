"""Minimal OPN target reaction builder for ``PichiaModel``.

This module contains the bare minimum needed to add OPN target-protein
reactions to a ``PichiaModel``. It replaces the legacy 7 000-line
``pichia_target_extender.py`` which had 25+ route-specific wrappers.

Only the OPN-like (OG-only extracellular) path is kept — all other
synthetic PTM routes (DSB-only, NG-only, GPI-only, transmembrane, etc.)
have been removed. New target-planning work should use the validated
probe prototype (``pcsec_pichia.probe``) with its ``CobraModel``-based
workflow instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from pcsec_pichia.core.pichia_model import PichiaModel
from pcsec_pichia.core.target_protein_plan import TargetProteinBuildPlan


@dataclass(frozen=True)
class OpnBuildResult:
    """Result of adding OPN target reactions to a PichiaModel."""

    model: PichiaModel
    protein_id: str
    added_reaction_ids: tuple[str, ...]
    added_reaction_count: int
    added_metabolite_count: int
    exchange_reaction_id: str
    status: str = "opn_built"


def add_opn_target_reactions(
    model: PichiaModel,
    plan: TargetProteinBuildPlan,
    amino_acids: Any,
) -> OpnBuildResult:
    """Add OPN target-protein reactions to a *PichiaModel*.

    This is the single remaining entry point for PichiaModel-based target
    extension. It chains five sub-steps: translocation, OG-misfolding,
    ER-Golgi transport, translation-degradation, and exchange.
    """
    original_rxn_count = len(model.rxns)
    original_met_count = len(model.mets)
    added_ids: list[str] = []

    # 1. Post-translation translocation
    tloc = _add_translocation(model, plan)
    added_ids.extend(tloc["added_reaction_ids"])

    # 2. OG glycosylation + misfolding
    misfold = _add_og_misfolding(tloc["model"], plan)
    added_ids.extend(misfold["added_reaction_ids"])

    # 3. ER -> Golgi transport, maturation, final transport
    transport = _add_er_golji_transport(misfold["model"], plan)
    added_ids.extend(transport["added_reaction_ids"])

    # 4. Translation + degradation (requires amino acid stoichiometry)
    core = _add_translation_degradation(transport["model"], plan, amino_acids)
    added_ids.extend(core["added_reaction_ids"])

    # 5. Exchange reaction
    exchange = _add_exchange(core["model"], plan)
    added_ids.append(exchange["exchange_reaction_id"])

    planned = plan.reaction_ids
    if set(added_ids) != set(planned) or len(added_ids) != len(planned):
        raise ValueError(
            f"OPN target extension produced {len(added_ids)} reactions but plan expects {len(planned)}."
        )

    return OpnBuildResult(
        model=exchange["model"],
        protein_id=plan.protein_id,
        added_reaction_ids=planned,
        added_reaction_count=len(exchange["model"].rxns) - original_rxn_count,
        added_metabolite_count=len(exchange["model"].mets) - original_met_count,
        exchange_reaction_id=exchange["exchange_reaction_id"],
    )


# ---------------------------------------------------------------------------
# Sub-step helpers (portable stoichiometry knowledge)
# ---------------------------------------------------------------------------

def _translocation_reaction_ids(protein_id: str) -> tuple[str, ...]:
    return (
        f"{protein_id}_Post_translation_PSTA_sec_RAC_complex",
        f"{protein_id}_Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex",
        f"{protein_id}_Post_translation_PSTA_sec_SEC61SEC63C_complex",
        f"{protein_id}_Post_translation_PSTA_sec_BIP_NEFS_complex",
        f"{protein_id}_Post_translation_TC_sec_SPC_complex",
        f"{protein_id}_export_sp_to_c",
    )


def _add_translocation(model: PichiaModel, plan: TargetProteinBuildPlan) -> dict[str, Any]:
    pid = plan.protein_id
    rids = _translocation_reaction_ids(pid)
    length_factor = _mround(len(plan.full_sequence) / 40.0)
    stoichs = (
        {f"{pid}_peptide[c]": -1.0, f"{pid}_translocate_1[c]": 1.0},
        {f"{pid}_translocate_1[c]": -1.0, "atp[c]": -1.0, "h2o[c]": -1.0,
         f"{pid}_translocate_2[c]": 1.0, "adp[c]": 1.0, "h[c]": 1.0, "pi[c]": 1.0},
        {f"{pid}_translocate_2[c]": -1.0, f"{pid}_translocate_3[c]": 1.0},
        {f"{pid}_translocate_3[c]": -1.0, "atp[c]": -length_factor, "h2o[c]": -length_factor,
         f"{pid}[er]": 1.0, "adp[c]": length_factor, "h[c]": length_factor, "pi[c]": length_factor},
        {f"{pid}_translocate_3[c]": -1.0, "h2o[c]": -1.0, f"{pid}[er]": 1.0, f"{pid}_sp[er]": 1.0},
        {f"{pid}_sp[er]": -1.0, f"{pid}_sp[c]": 1.0},
    )
    for rid, stoich in zip(rids, stoichs):
        model = model.add_reaction(rid, stoich, lower_bound=0.0, upper_bound=1000.0)
    return {"model": model, "added_reaction_ids": list(rids)}

def _add_og_misfolding(model: PichiaModel, plan: TargetProteinBuildPlan) -> dict[str, Any]:
    pid = plan.protein_id
    length_factor = _mround(len(plan.full_sequence) / 40.0)
    acc_factor = 10.0 * length_factor
    og = float(plan.o_glycosylation_sites)

    base = f"{pid}"
    rids = [
        f"{pid}_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex",
        f"{pid}_misfold_ERAD_sec_Kar2p_complex",
        f"{pid}_ERAD2B",
        f"{pid}_ERAD3B",
        f"{pid}_ERAD4B",
        f"{pid}_ERAD5B",
        f"{pid}_ERADL_sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex",
        f"{pid}_ERADL_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex",
        f"{pid}_ERAD7B_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex",
        f"{pid}_degradation_misfolding_c",
        f"{pid}_cycle_accumulation",
        f"{pid}_cycle_accumulation_sec_acc_Kar2p_complex",
        f"{pid}_dilution_misfolding_er",
    ]
    stoichs = [
        {f"{base}[er]": -1.0, "dolmanp[er]": -og, f"{base}_OG_M1[er]": 1.0, "dolp[er]": og},
        {f"{base}_OG_M1[er]": -1.0, "atp[er]": -length_factor, "h2o[er]": -length_factor,
         f"{base}_OG_M1_misf[er]": 1.0, "adp[er]": length_factor, "h[er]": length_factor, "pi[er]": length_factor},
        {f"{base}_OG_M1_misf[er]": -1.0, f"{base}_OG_M1_misf_G1[er]": 1.0},
        {f"{base}_OG_M1_misf_G1[er]": -1.0, f"{base}_OG_M1_misf_G2[er]": 1.0},
        {f"{base}_OG_M1_misf_G2[er]": -1.0, f"{base}_OG_M1_misf_G3[er]": 1.0},
        {f"{base}_OG_M1_misf_G3[er]": -1.0, f"{base}_OG_M1_misf_G4[er]": 1.0},
        {f"{base}_OG_M1_misf_G4[er]": -1.0, f"{base}_OG_M1_misf_G5[er]": 1.0},
        {f"{base}_OG_M1_misf_G5[er]": -1.0, "Ubiquitin_for_Transfer[c]": -8.0,
         f"{base}_OG_M1_misf_G6[c]": 1.0, "Ubiquitin[c]": 8.0},
        {f"{base}_OG_M1_misf_G6[c]": -1.0, f"{pid}_misfolding[c]": 1.0, "man[er]": og},
        {f"{pid}_misfolding[c]": -1.0, f"{pid}_subunit[c]": 1.0},
        {f"{base}_OG_M1_misf[er]": -1.0, f"{base}_OG_M1_misf2[er]": 1.0},
        {f"{base}_OG_M1_misf2[er]": -1.0, "atp[er]": -acc_factor, "h2o[er]": -acc_factor,
         f"{base}_OG_M1_misfolding_acc[er]": 1.0, "adp[er]": acc_factor, "h[er]": acc_factor, "pi[er]": acc_factor},
        {f"{base}_OG_M1_misfolding_acc[er]": -1.0},
    ]
    for rid, stoich in zip(rids, stoichs):
        model = model.add_reaction(rid, stoich, lower_bound=0.0, upper_bound=1000.0)
    return {"model": model, "added_reaction_ids": rids}


def _add_er_golji_transport(model: PichiaModel, plan: TargetProteinBuildPlan) -> dict[str, Any]:
    pid = plan.protein_id
    og = float(plan.o_glycosylation_sites)
    mw_ratio = float(plan.protein_mw) / 54580.0

    rids = [
        f"{pid}_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex",
        f"{pid}_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Shl23p_Lst1p_Erv29p_Bet1p_Bos1p_complex",
        f"{pid}_COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex",
        f"{pid}_COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex",
        f"{pid}_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex",
        f"{pid}_Mature",
        f"{pid}_HDSVI_sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex",
        f"{pid}_HDSVII_sec_Vps1p_Chc1p_Clc1p_complex",
    ]
    stoichs = [
        # COPII coat complexes (1-3 use mw_ratio, 4 is plain dilution)
        {f"{pid}_OG_M1[er]": -1.0, f"{pid}_OG_M1[Golgi]": mw_ratio},
        {f"{pid}_OG_M1[er]": -1.0, f"{pid}_OG_M1[Golgi]": mw_ratio},
        {f"{pid}_OG_M1[er]": -1.0, f"{pid}[Golgi]": 1.0},
        {f"{pid}_OG_M1[er]": -1.0, f"{pid}[Golgi]": 1.0},
        # Golgi O-mannosylation
        {f"{pid}_OG_M1[Golgi]": -1.0, "dolmanp[Golgi]": -og, f"{pid}_OG[Golgi]": 1.0, "dolp[Golgi]": og},
        # Mature
        {f"{pid}_OG[Golgi]": -1.0, f"{pid}[e]": 1.0},
        # Final transport
        {f"{pid}_OG[Golgi]": -1.0, f"{pid}_secreted_vesicle[Golgi]": 1.0},
        {f"{pid}_secreted_vesicle[Golgi]": -1.0, f"{pid}[e]": 1.0},
    ]
    for rid, stoich in zip(rids, stoichs):
        model = model.add_reaction(rid, stoich, lower_bound=0.0, upper_bound=1000.0)
    return {"model": model, "added_reaction_ids": rids}


def _add_translation_degradation(
    model: PichiaModel,
    plan: TargetProteinBuildPlan,
    amino_acids: Any,
) -> dict[str, Any]:
    pid = plan.protein_id
    full_seq = str(plan.full_sequence)
    sp = str(plan.signal_peptide_sequence) or ""
    sp_len = len(sp)

    rids = [
        f"r_{pid}_peptide_translation",
        f"r_{pid}_SP_degradation",
        f"r_{pid}_subunit_degradation",
    ]
    stoichs = [
        amino_acids.translation_stoichiometry(pid, full_seq),
        amino_acids.signal_peptide_degradation_stoichiometry(pid, sp) if sp_len else {f"{pid}_sp[c]": -1.0, f"{pid}_subunit[c]": 1.0},
        amino_acids.subunit_degradation_stoichiometry(pid, full_seq[sp_len:]) if sp_len else amino_acids.subunit_degradation_stoichiometry(pid, full_seq),
    ]
    for rid, stoich in zip(rids, stoichs):
        model = model.add_reaction(rid, stoich, lower_bound=0.0, upper_bound=1000.0)
    return {"model": model, "added_reaction_ids": rids}


def _add_exchange(model: PichiaModel, plan: TargetProteinBuildPlan) -> dict[str, Any]:
    pid = plan.protein_id
    rid = f"{pid} exchange"
    model = model.add_reaction(rid, {f"{pid}[e]": -1.0}, lower_bound=0.0, upper_bound=1000.0)
    return {"model": model, "exchange_reaction_id": rid}


def _mround(value: float) -> float:
    """MATLAB's round() behaviour for positive values (round-half-up)."""
    return float(int(value + 0.5))


__all__ = [
    "OpnBuildResult",
    "add_opn_target_reactions",
]
