from __future__ import annotations

import numpy as np

from pcsec_pichia.core.pichia_enzymes import SecretoryEnzymeData, TargetProteinEnzymeData
from pcsec_pichia.core.pichia_model import PichiaModel
from pcsec_pichia.core.target_protein_plan import TargetProteinBuildPlan


PST_INFO = ("DSB", "NG", "OG", "GPI")
PICHIA_PTM_EXTRA_MW = np.array([0.0, 3346.0, 1080.0, 2009.0], dtype=float)


def target_enzymedata_from_plan(plan: TargetProteinBuildPlan) -> TargetProteinEnzymeData:
    """Build the target-protein enzymedata record matching calculateMW.m.

    MATLAB creates enzymedata_TP with proteins/kdeg in addTargetProtein(), then
    calculateMW() fills proteinMWs, proteinLength, proteinPST, proteinExtraMW
    and proteinLoc. This function ports that narrow target-protein part.
    """

    protein_pst = np.array(
        [[plan.disulfide_sites, plan.n_glycosylation_sites, plan.o_glycosylation_sites, plan.gpi_sites]],
        dtype=float,
    )
    extra_specific = protein_pst * PICHIA_PTM_EXTRA_MW
    return TargetProteinEnzymeData(
        proteins=[plan.protein_id],
        protein_mw=np.array([plan.protein_mw], dtype=float),
        protein_length=np.array([plan.full_length], dtype=float),
        protein_pst=protein_pst,
        protein_pst_info=PST_INFO,
        protein_extra_mw_specific=extra_specific,
        protein_extra_mw=np.sum(extra_specific, axis=1),
        protein_loc=[plan.localization],
        kdeg=np.zeros(1, dtype=float),
    )


def simulate_target_reaction_coefficients(
    model: PichiaModel,
    secretory_enzymedata: SecretoryEnzymeData,
    target_enzymedata: TargetProteinEnzymeData,
) -> TargetProteinEnzymeData:
    return _simulate_target_reaction_coefficients(
        model,
        secretory_enzymedata,
        target_enzymedata,
        reaction_to_protein=_target_protein_for_reaction,
    )


def simulate_target_reaction_coefficients_matlab_compatible(
    model: PichiaModel,
    secretory_enzymedata: SecretoryEnzymeData,
    target_enzymedata: TargetProteinEnzymeData,
) -> TargetProteinEnzymeData:
    """Mimic SimulateRxnKcatCoef.m target-protein matching.

    The MATLAB implementation extracts the reaction prefix before the third
    underscore and matches that against enzymedata.proteins. This is not the
    cleanest long-term behavior for target ids containing more than three
    underscore-delimited fields, but preserving it is essential when comparing
    Python-generated smoke LPs against MATLAB baseline LP files.
    """

    return _simulate_target_reaction_coefficients(
        model,
        secretory_enzymedata,
        target_enzymedata,
        reaction_to_protein=_matlab_prefix_target_protein_for_reaction,
    )


def _simulate_target_reaction_coefficients(
    model: PichiaModel,
    secretory_enzymedata: SecretoryEnzymeData,
    target_enzymedata: TargetProteinEnzymeData,
    reaction_to_protein,
) -> TargetProteinEnzymeData:
    coefficients: dict[str, float] = {}
    for complex_id, coefficient_ref in zip(secretory_enzymedata.complexes, secretory_enzymedata.coefficient_refs):
        suffix = f"_{complex_id}"
        for reaction_id in model.rxns:
            if reaction_id.startswith("dummyER") or not reaction_id.endswith(suffix):
                continue
            protein_id = reaction_to_protein(reaction_id, target_enzymedata.proteins)
            if protein_id is None:
                continue
            coefficients[reaction_id] = _evaluate_coefficient_ref(coefficient_ref, target_enzymedata, protein_id)
    return TargetProteinEnzymeData(
        proteins=target_enzymedata.proteins,
        protein_mw=target_enzymedata.protein_mw,
        protein_length=target_enzymedata.protein_length,
        protein_pst=target_enzymedata.protein_pst,
        protein_pst_info=target_enzymedata.protein_pst_info,
        protein_extra_mw_specific=target_enzymedata.protein_extra_mw_specific,
        protein_extra_mw=target_enzymedata.protein_extra_mw,
        protein_loc=target_enzymedata.protein_loc,
        kdeg=target_enzymedata.kdeg,
        reaction_coefficients=coefficients,
    )


def _target_protein_for_reaction(reaction_id: str, protein_ids: list[str]) -> str | None:
    for protein_id in sorted(protein_ids, key=len, reverse=True):
        if reaction_id.startswith(f"{protein_id}_"):
            return protein_id
    return None


def _matlab_prefix_target_protein_for_reaction(reaction_id: str, protein_ids: list[str]) -> str | None:
    prefix = _matlab_prefix_before_third_underscore(reaction_id)
    if prefix is None:
        return None
    return prefix if prefix in protein_ids else None


def _matlab_prefix_before_third_underscore(reaction_id: str) -> str | None:
    positions = [index for index, character in enumerate(reaction_id) if character == "_"]
    if len(positions) < 3:
        return None
    return reaction_id[: positions[2]]


def _evaluate_coefficient_ref(
    coefficient_ref: str,
    target_enzymedata: TargetProteinEnzymeData,
    protein_id: str,
) -> float:
    coefficient_ref = coefficient_ref.strip()
    index = target_enzymedata.index_for_protein(protein_id)
    protein_length = float(target_enzymedata.protein_length[index])
    protein_mw = float(target_enzymedata.protein_mw[index])
    variables = {
        "DSB": target_enzymedata.pst_value(protein_id, "DSB"),
        "NG": target_enzymedata.pst_value(protein_id, "NG"),
        "OG": target_enzymedata.pst_value(protein_id, "OG"),
        "GPI": target_enzymedata.pst_value(protein_id, "GPI"),
        "proteinLength": protein_length,
        "ProteinMW": protein_mw,
    }
    if coefficient_ref == "proteinLength":
        return protein_length / 467.0
    if coefficient_ref == "ProteinMW":
        return protein_mw / 54580.0
    try:
        return float(eval(coefficient_ref, {"__builtins__": {}}, variables))
    except Exception as exc:  # noqa: BLE001
        raise ValueError(f"Unsupported target enzymedata coefficient reference: {coefficient_ref}") from exc
