from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Literal

from pcsec_pichia.core.target_inputs import LeaderCandidateInput, TargetProteinInput, build_pcsec_target_row


PlanStage = Literal[
    "translation",
    "translocation",
    "folding",
    "misfolding",
    "er_to_golgi",
    "golgi_processing",
    "maturation",
    "final_transport",
    "degradation",
    "exchange",
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


@dataclass(frozen=True)
class PlannedTargetReaction:
    reaction_id: str
    stage: PlanStage
    source_function: str


@dataclass(frozen=True)
class TargetProteinBuildPlan:
    target_id: str
    protein_id: str
    full_sequence: str
    mature_sequence: str
    leader_sequence: str
    signal_peptide_sequence: str
    signal_peptide_length: int
    localization: str
    through_er: bool
    disulfide_sites: int
    n_glycosylation_sites: int
    o_glycosylation_sites: int
    transmembrane: int
    gpi_sites: int
    cotranslation: int
    protein_mw: float
    protein_extra_mw: float
    reactions: tuple[PlannedTargetReaction, ...]

    @property
    def full_length(self) -> int:
        return len(self.full_sequence)

    @property
    def reaction_ids(self) -> tuple[str, ...]:
        return tuple(reaction.reaction_id for reaction in self.reactions)

    def reaction_ids_by_stage(self, stage: PlanStage) -> tuple[str, ...]:
        return tuple(reaction.reaction_id for reaction in self.reactions if reaction.stage == stage)


def build_target_protein_plan(
    target: TargetProteinInput,
    leader: LeaderCandidateInput,
) -> TargetProteinBuildPlan:
    target = target.normalized()
    leader = leader.normalized()
    errors = target.validation_errors() + leader.validation_errors()
    if errors:
        raise ValueError("; ".join(errors))

    row = build_pcsec_target_row(target, leader)
    protein_id = str(row["abbreviation"])
    full_sequence = str(row["sequence"])
    localization = str(row["Localization"]).split(",", maxsplit=1)[0] or "e"
    reactions: list[PlannedTargetReaction] = []

    if int(row["ThroughER"]) == 1:
        reactions.extend(_secretory_reactions(row))
    else:
        reactions.extend(_non_secretory_folding_reactions(row))
    _add(reactions, f"r_{protein_id}_peptide_translation", "translation", "addTranslationRxns")
    reactions.extend(_degradation_reactions(row))
    _add(reactions, f"{protein_id} exchange", "exchange", "addTargetProtein")

    return TargetProteinBuildPlan(
        target_id=target.target_id,
        protein_id=protein_id,
        full_sequence=full_sequence,
        mature_sequence=target.mature_sequence,
        leader_sequence=leader.leader_sequence,
        signal_peptide_sequence=str(row["sp sequence"]),
        signal_peptide_length=int(row["Signal peptide length"]),
        localization=localization,
        through_er=int(row["ThroughER"]) == 1,
        disulfide_sites=int(row["Disulfide site"]),
        n_glycosylation_sites=int(row["N-glycosylation site"]),
        o_glycosylation_sites=int(row["O-linked glycisylation "]),
        transmembrane=int(row["Transmembrane"]),
        gpi_sites=int(row["GPI site"]),
        cotranslation=int(row["Cotranslation"]),
        protein_mw=protein_mw(full_sequence),
        protein_extra_mw=protein_extra_mw(
            disulfides=int(row["Disulfide site"]),
            n_glycans=int(row["N-glycosylation site"]),
            o_glycans=int(row["O-linked glycisylation "]),
            gpi_sites=int(row["GPI site"]),
        ),
        reactions=tuple(reactions),
    )


def protein_mw(sequence: str) -> float:
    counts = Counter(sequence)
    return 18.0 + sum(counts[aa] * mw for aa, mw in AA_MW.items())


def protein_extra_mw(disulfides: int, n_glycans: int, o_glycans: int, gpi_sites: int) -> float:
    return disulfides * 0.0 + n_glycans * 3346.0 + o_glycans * 1080.0 + gpi_sites * 2009.0


def _secretory_reactions(row: dict[str, object]) -> list[PlannedTargetReaction]:
    protein_id = str(row["abbreviation"])
    full_length = int(row["Length"])
    signal_peptide = int(row["Signal peptide "])
    disulfides = int(row["Disulfide site"])
    n_glycans = int(row["N-glycosylation site"])
    o_glycans = int(row["O-linked glycisylation "])
    transmembrane = int(row["Transmembrane"])
    gpi_sites = int(row["GPI site"])
    localization = str(row["Localization"]).split(",", maxsplit=1)[0] or "e"
    cotranslation = int(row["Cotranslation"])

    reactions: list[PlannedTargetReaction] = []
    reactions.extend(_translocation_reactions(protein_id, cotranslation, signal_peptide, gpi_sites))
    reactions.extend(_dsb_reactions(protein_id, disulfides))
    reactions.extend(_gpi_reactions(protein_id, gpi_sites))
    reactions.extend(_og_reactions(protein_id, o_glycans))
    reactions.extend(_ng_reactions(protein_id, n_glycans))
    reactions.extend(_misfolding_reactions(protein_id, full_length, n_glycans, o_glycans, disulfides, gpi_sites, transmembrane, localization))
    reactions.extend(_er_to_golgi_reactions(protein_id, gpi_sites, transmembrane))
    reactions.extend(_golgi_n_reactions(protein_id, n_glycans))
    reactions.extend(_golgi_o_reactions(protein_id, o_glycans))
    _add(reactions, f"{protein_id}_Mature", "maturation", "mature")
    reactions.extend(_final_transport_reactions(protein_id, localization))
    return reactions


def _non_secretory_folding_reactions(row: dict[str, object]) -> list[PlannedTargetReaction]:
    protein_id = str(row["abbreviation"])
    localization = str(row["Localization"]).split(",", maxsplit=1)[0] or "c"
    reactions: list[PlannedTargetReaction] = []
    if localization != "c":
        _add(reactions, f"{protein_id}_importing_{localization}", "translocation", "addfolding")
        _add(reactions, f"{protein_id}_subunit_export_{localization}", "final_transport", "addfolding")
    _add(reactions, f"{protein_id}_folding_{localization}", "folding", "addfolding")
    _add(reactions, f"{protein_id}_misfold_{localization}", "misfolding", "addfolding")
    _add(reactions, f"{protein_id}_degradation_misfolding_{localization}", "misfolding", "addfolding")
    _add(reactions, f"{protein_id}_dilution_misfolding_{localization}", "misfolding", "addfolding")
    return reactions


def _translocation_reactions(
    protein_id: str,
    cotranslation: int,
    signal_peptide: int,
    gpi_sites: int,
) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if cotranslation == 1:
        suffixes = (
            "co_translation_TC_sec_SRPC_complex",
            "co_translation_TC_sec_SRC_complex",
            "co_translation_TC_sec_SEC61C_complex",
            "co_translation_TC_sec_SSH1C_complex",
            "co_translation_TC_sec_SPC_complex",
            "export_sp_to_c",
        )
        source = "co_Translation_translocation"
    elif gpi_sites > 0 and signal_peptide == 0:
        suffixes = (
            "Post_translation_PSTB_sec_Sgt2_Get4_Get5_complex",
            "Post_translation_PSTB_sec_Get3_complex",
            "Post_translation_PSTB_sec_Get1_Get2_complex",
        )
        source = "Post_Translation_translocation_tail"
    else:
        suffixes = (
            "Post_translation_PSTA_sec_RAC_complex",
            "Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex",
            "Post_translation_PSTA_sec_SEC61SEC63C_complex",
            "Post_translation_PSTA_sec_BIP_NEFS_complex",
            "Post_translation_TC_sec_SPC_complex",
            "export_sp_to_c",
        )
        source = "Post_Translation_translocation"
    for suffix in suffixes:
        _add(reactions, f"{protein_id}_{suffix}", "translocation", source)
    return reactions


def _dsb_reactions(protein_id: str, disulfides: int) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if disulfides > 0:
        _add(reactions, f"{protein_id}_DSB_sec_BIP_NEFS_complex", "folding", "addDSB")
        _add(reactions, f"{protein_id}_DSB_PDI_II_sec_PDI1_ERV2_Ero1p_complex", "folding", "addDSB")
    return reactions


def _gpi_reactions(protein_id: str, gpi_sites: int) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if gpi_sites > 0:
        for suffix in (
            "GPIRI_sec_GPIR_complex",
            "GPIRII_sec_Bst1p_complex",
            "GPIRIII_sec_Per1p_complex",
            "GPIRIV_sec_Gup1p_complex",
            "GPIRV_sec_Cwh43p_Las21p_Mcd4p_complex",
            "GPIRVI_sec_Ted1p_complex",
            "GPIRIB_sec_GPIR_complex",
        ):
            _add(reactions, f"{protein_id}_{suffix}", "folding", "addGPI")
    return reactions


def _og_reactions(protein_id: str, o_glycans: int) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if o_glycans > 0:
        _add(reactions, f"{protein_id}_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex", "folding", "addOG")
    return reactions


def _ng_reactions(protein_id: str, n_glycans: int) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if n_glycans > 0:
        for suffix in (
            "ERNG_NG_sec_OSTC_complex",
            "ERNG_FLI_NG_sec_Cwh41p_complex",
            "ERNG_FLII_NG_sec_Rot2p_complex",
            "ERNG_FLIII_NG_sec_Rot2p_complex",
            "ERNG_FLIV_NG_sec_Mns1p_complex",
        ):
            _add(reactions, f"{protein_id}_{suffix}", "folding", "addNG")
    return reactions


def _misfolding_reactions(
    protein_id: str,
    full_length: int,
    n_glycans: int,
    o_glycans: int,
    disulfides: int,
    gpi_sites: int,
    transmembrane: int,
    localization: str,
) -> list[PlannedTargetReaction]:
    del full_length
    reactions: list[PlannedTargetReaction] = []
    _add(reactions, f"{protein_id}_misfold_ERAD_sec_Kar2p_complex", "misfolding", "addMisfold")
    if disulfides > 0:
        _add(reactions, f"{protein_id}_ERAD2A_sec_Pdi1p_complex", "misfolding", "addMisfold")
    else:
        _add(reactions, f"{protein_id}_ERAD2B", "misfolding", "addMisfold")
    if n_glycans > 0:
        _add(reactions, f"{protein_id}_ERAD3A_sec_Mns1p_complex", "misfolding", "addMisfold")
        _add(reactions, f"{protein_id}_ERAD4A_sec_Mnl1p_Pdi1p_complex", "misfolding", "addMisfold")
    else:
        _add(reactions, f"{protein_id}_ERAD3B", "misfolding", "addMisfold")
        _add(reactions, f"{protein_id}_ERAD4B", "misfolding", "addMisfold")
    _add(reactions, f"{protein_id}_ERAD5A" if gpi_sites > 0 else f"{protein_id}_ERAD5B", "misfolding", "addMisfold")
    if transmembrane > 0 and localization == "c" and disulfides == 0 and n_glycans == 0:
        _add(reactions, f"{protein_id}_ERADC_sec_Ubc6p_Ubc7p_Doa10p_complex", "misfolding", "addMisfold")
        _add(reactions, f"{protein_id}_ERADC_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex", "misfolding", "addMisfold")
    elif transmembrane > 0 and localization in {"er", "erm"}:
        _add(reactions, f"{protein_id}_ERADM_sec_Ubc6p_Ubc7p_Hrd1p_Hrd3p_Der1p_complex", "misfolding", "addMisfold")
        _add(reactions, f"{protein_id}_ERADM_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex", "misfolding", "addMisfold")
        _add(reactions, f"{protein_id}_ERADM2_sec_Ubc6p_Ubc7p_Doa10p_complex", "misfolding", "addMisfold")
    else:
        _add(reactions, f"{protein_id}_ERADL_sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex", "misfolding", "addMisfold")
        _add(reactions, f"{protein_id}_ERADL_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex", "misfolding", "addMisfold")
    if n_glycans > 0:
        _add(reactions, f"{protein_id}_ERAD7A_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex", "misfolding", "addMisfold")
    elif o_glycans > 0:
        _add(reactions, f"{protein_id}_ERAD7B_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex", "misfolding", "addMisfold")
    else:
        _add(reactions, f"{protein_id}_ERAD7C_sec_Dsk2p_Rad23p_Uba1p_complex", "misfolding", "addMisfold")
    _add(reactions, f"{protein_id}_degradation_misfolding_c", "misfolding", "addMisfold")
    if disulfides > 0:
        _add(reactions, f"{protein_id}_cycle_accumulation_sec_pdi1p_ero1p_complex", "misfolding", "addMisfold")
    else:
        _add(reactions, f"{protein_id}_cycle_accumulation", "misfolding", "addMisfold")
    _add(reactions, f"{protein_id}_cycle_accumulation_sec_acc_Kar2p_complex", "misfolding", "addMisfold")
    _add(reactions, f"{protein_id}_dilution_misfolding_er", "misfolding", "addMisfold")
    return reactions


def _er_to_golgi_reactions(protein_id: str, gpi_sites: int, transmembrane: int) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if gpi_sites > 0:
        suffixes = (
            "COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex",
            "COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Shl23p_Lst1p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex",
            "COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex",
            "COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex",
        )
        source = "coat_GPI"
    elif transmembrane > 0:
        suffixes = (
            "COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Sec23p_Sec24p_Bet1p_Bos1p_complex",
            "COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Shl23p_Lst1p_Bet1p_Bos1p_complex",
            "COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex",
            "COPII_ERGL_sec_Ypt1p_Uso1p_bug1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex",
        )
        source = "coat_trans_membrane"
    else:
        suffixes = (
            "COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex",
            "COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Shl23p_Lst1p_Erv29p_Bet1p_Bos1p_complex",
            "COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex",
            "COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex",
        )
        source = "coat_other"
    for suffix in suffixes:
        _add(reactions, f"{protein_id}_{suffix}", "er_to_golgi", source)
    return reactions


def _golgi_n_reactions(protein_id: str, n_glycans: int) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if n_glycans > 0:
        for suffix in (
            "GLNG_Golgi_N_linked_glycosylation_I_sec_Och1p_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex",
            "GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pA_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pB_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pC_complex",
        ):
            _add(reactions, f"{protein_id}_{suffix}", "golgi_processing", "golgiProcessing_N_PP")
    return reactions


def _golgi_o_reactions(protein_id: str, o_glycans: int) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if o_glycans > 0:
        _add(reactions, f"{protein_id}_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex", "golgi_processing", "golgiProcessing_O")
    return reactions


def _final_transport_reactions(protein_id: str, localization: str) -> list[PlannedTargetReaction]:
    reactions: list[PlannedTargetReaction] = []
    if localization == "e":
        _add(
            reactions,
            f"{protein_id}_HDSVI_sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex",
            "final_transport",
            "transportFromGolgiToS",
        )
        _add(reactions, f"{protein_id}_HDSVII_sec_Vps1p_Chc1p_Clc1p_complex", "final_transport", "transportFromGolgiToS")
    elif localization == "ce":
        _add(
            reactions,
            f"{protein_id}_LDSV_sec_Arf1p_Sec3p_Sec5p_Sec6p_Sec8p_Sec10p_Sec15p_Exo70p_Exo84p_Sec4p_Chc1p_Clc1p_complex",
            "final_transport",
            "transportFromGolgiToCe",
        )
    elif localization == "vm":
        _add(
            reactions,
            f"{protein_id}_ALPtransport_sec_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p_complex",
            "final_transport",
            "transportFromGolgiToVM",
        )
    elif localization == "v":
        _add(
            reactions,
            f"{protein_id}_CPYI_sec_Gga1p_Gga2p_Arf1p_Apl4p_Apl2p_Apm1p_Aps1p_Chc1p_Clc1p_Pep12p_Vps45p_Vps5p_Swa2p_complex",
            "final_transport",
            "transportFromGolgiToV",
        )
        _add(
            reactions,
            f"{protein_id}_CPYII_sec_Vps4p_Vps27p_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_complex",
            "final_transport",
            "transportFromGolgiToV",
        )
    elif localization in {"er", "erm"}:
        _add(
            reactions,
            f"{protein_id}_GLER_COPI_formation_sec_Arf1p_Gea2p_Rer1p_Erd2p_Cop1p_Sec26p_Sec27p_Sec21p_Ret2p_Sec28p_Ret3p_complex",
            "final_transport",
            "transportFromGolgiToER",
        )
        _add(
            reactions,
            f"{protein_id}_GLER_COPI_uncoating_and_fission_sec_Rer1p_Ret2p_Cop1p_Sec27p_Sec21p_Bet1p_complex",
            "final_transport",
            "transportFromGolgiToER",
        )
        _add(reactions, f"{protein_id}_GLER3_Final_demand", "final_transport", "transportFromGolgiToER")
    else:
        _add(reactions, f"{protein_id}_transportFromGolgiToOthercompartment", "final_transport", "transportFromGolgiToOther")
    return reactions


def _degradation_reactions(row: dict[str, object]) -> list[PlannedTargetReaction]:
    protein_id = str(row["abbreviation"])
    through_er = int(row["ThroughER"]) == 1
    has_signal_peptide = int(row["Signal peptide "]) == 1
    reactions: list[PlannedTargetReaction] = []
    if through_er and has_signal_peptide:
        _add(reactions, f"r_{protein_id}_SP_degradation", "degradation", "addDegradationRxns")
    _add(reactions, f"r_{protein_id}_subunit_degradation", "degradation", "addDegradationRxns")
    return reactions


def _add(
    reactions: list[PlannedTargetReaction],
    reaction_id: str,
    stage: PlanStage,
    source_function: str,
) -> None:
    reactions.append(PlannedTargetReaction(reaction_id=reaction_id, stage=stage, source_function=source_function))
