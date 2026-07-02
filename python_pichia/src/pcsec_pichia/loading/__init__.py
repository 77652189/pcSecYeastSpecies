"""Formal pcSecPichia loading entrypoints.

This module is still a prototype-backed adapter: it is the reviewed boundary
where formal workflows may reach the original probe loaders while migration is
in progress. Downstream app/service code should import loading helpers from
here rather than importing ``pcsec_pichia.probe`` directly.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from pcsec_pichia.core.modes import DEFAULT_COMPATIBILITY_MODE, CompatibilityMode, validate_compatibility_mode
from pcsec_pichia.media import (
    CarbonSourceFormulation,
    compose_medium_condition_bounds,
    load_carbon_source_formulation,
    load_medium_condition_spec,
    summarize_medium_condition,
)
from pcsec_pichia.probe import (
    AminoAcidStoichiometry,
    CombinedEnzymeData,
    CobraModel,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    load_aa_stoichiometry,
    load_combined_enzymedata,
    load_metabolic_enzymedata,
    load_pcsec_pichia_model,
    load_secretory_enzymedata,
    repo_root,
    set_media_pp_like_matlab,
)


MINIMAL_EXCHANGE_REACTIONS: tuple[str, ...] = (
    "Ex_nh4",
    "Ex_o2",
    "Ex_pi",
    "Ex_so4",
    "Ex_fe2",
    "Ex_h",
    "Ex_h2o",
    "Ex_na1",
    "Ex_k",
    "Ex_co2",
)


@dataclass(frozen=True)
class PcSecPichiaInputs:
    root: Path
    model: CobraModel
    prepared_model: CobraModel
    amino_acids: AminoAcidStoichiometry
    metabolic: MetabolicEnzymeData
    secretory: SecretoryEnzymeData
    combined: CombinedEnzymeData
    compatibility_mode: CompatibilityMode
    media_type: int
    carbon_source_id: str
    carbon_source_formulation: CarbonSourceFormulation
    medium_condition_id: str


CARBON_SOURCE_IDS: tuple[str, ...] = (
    "glucose",
    "glycerol",
    "methanol",
    "glucose_glycerol",
    "glycerol_methanol",
)


def set_media_pp(
    model: CobraModel,
    media_type: int = 4,
    compatibility_mode: CompatibilityMode = DEFAULT_COMPATIBILITY_MODE,
) -> CobraModel:
    """Configure pcSecPichia medium with an explicit MATLAB compatibility mode."""

    mode = validate_compatibility_mode(compatibility_mode)
    configured = set_media_pp_like_matlab(model, media_type=media_type)
    if mode == "matlab_compat":
        return configured
    return configured.with_bounds({reaction_id: (-1000.0, None) for reaction_id in MINIMAL_EXCHANGE_REACTIONS})


def prepare_glucose_model(
    model: CobraModel,
    media_type: int = 4,
    compatibility_mode: CompatibilityMode = DEFAULT_COMPATIBILITY_MODE,
) -> CobraModel:
    return prepare_carbon_source_model(
        model,
        media_type=media_type,
        compatibility_mode=compatibility_mode,
        carbon_source_id="glucose",
    )


def prepare_carbon_source_model(
    model: CobraModel,
    media_type: int = 4,
    compatibility_mode: CompatibilityMode = DEFAULT_COMPATIBILITY_MODE,
    carbon_source_id: str = "glucose",
) -> CobraModel:
    carbon_source = validate_carbon_source_id(carbon_source_id)
    mode = validate_compatibility_mode(compatibility_mode)
    configured = set_media_pp(model, media_type=media_type, compatibility_mode=mode)
    condition_id = medium_condition_id_for(media_type, mode, carbon_source)
    try:
        return configured.with_bounds(_medium_condition_bounds(configured, condition_id))
    except KeyError:
        pass

    return configured.with_bounds(_legacy_carbon_source_bounds(carbon_source))


def _legacy_carbon_source_bounds(carbon_source: str) -> dict[str, tuple[float | None, float | None]]:
    changes: dict[str, tuple[float | None, float | None]] = {
        "BIOMASS": (None, 1000.0),
        "Ex_o2": (-1000.0, None),
    }
    if carbon_source == "glucose":
        changes.update(
            {
                "Ex_glc_D": (-1000.0, None),
                "Ex_glyc": (0.0, None),
                "BIOMASS_glyc": (0.0, 0.0),
                "Ex_meoh": (0.0, None),
                "BIOMASS_meoh": (0.0, 0.0),
            }
        )
    elif carbon_source == "glycerol":
        changes.update(
            {
                "Ex_glc_D": (0.0, None),
                "Ex_glyc": (-1000.0, None),
                "BIOMASS_glyc": (0.0, 1000.0),
                "Ex_meoh": (0.0, None),
                "BIOMASS_meoh": (0.0, 0.0),
            }
        )
    elif carbon_source == "methanol":
        changes.update(
            {
                "Ex_glc_D": (0.0, None),
                "Ex_glyc": (0.0, None),
                "BIOMASS_glyc": (0.0, 0.0),
                "Ex_meoh": (-1000.0, None),
                "BIOMASS_meoh": (0.0, 1000.0),
            }
        )
    elif carbon_source == "glucose_glycerol":
        changes.update(
            {
                "Ex_glc_D": (-1000.0, None),
                "Ex_glyc": (-1000.0, None),
                "BIOMASS_glyc": (0.0, 1000.0),
                "Ex_meoh": (0.0, None),
                "BIOMASS_meoh": (0.0, 0.0),
            }
        )
    else:
        changes.update(
            {
                "Ex_glc_D": (0.0, None),
                "Ex_glyc": (-1000.0, None),
                "BIOMASS_glyc": (0.0, 1000.0),
                "Ex_meoh": (-1000.0, None),
                "BIOMASS_meoh": (0.0, 1000.0),
            }
        )
    return changes


def load_pcsec_pichia_inputs(
    root: Path | None = None,
    media_type: int = 4,
    compatibility_mode: CompatibilityMode = DEFAULT_COMPATIBILITY_MODE,
    carbon_source_id: str = "glucose",
) -> PcSecPichiaInputs:
    """Load the prototype-backed pcSecPichia input bundle for formal workflows."""

    mode = validate_compatibility_mode(compatibility_mode)
    carbon_source = validate_carbon_source_id(carbon_source_id)
    resolved_root = root or repo_root()
    model = load_pcsec_pichia_model(resolved_root)
    medium_condition_id = medium_condition_id_for(media_type, mode, carbon_source)
    carbon_source_formulation = load_carbon_source_formulation(carbon_source)
    return PcSecPichiaInputs(
        root=resolved_root,
        model=model,
        prepared_model=prepare_carbon_source_model(
            model,
            media_type=media_type,
            compatibility_mode=mode,
            carbon_source_id=carbon_source,
        ),
        amino_acids=load_aa_stoichiometry(resolved_root),
        metabolic=load_metabolic_enzymedata(resolved_root),
        secretory=load_secretory_enzymedata(resolved_root),
        combined=load_combined_enzymedata(resolved_root),
        compatibility_mode=mode,
        media_type=media_type,
        carbon_source_id=carbon_source,
        carbon_source_formulation=carbon_source_formulation,
        medium_condition_id=medium_condition_id,
    )


def validate_carbon_source_id(carbon_source_id: str) -> str:
    normalized = str(carbon_source_id or "glucose").strip().lower()
    aliases = {
        "glc": "glucose",
        "meoh": "methanol",
        "glucose+glycerol": "glucose_glycerol",
        "glycerol_glucose": "glucose_glycerol",
        "methanol_glycerol": "glycerol_methanol",
        "glycerol+methanol": "glycerol_methanol",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in CARBON_SOURCE_IDS:
        raise ValueError(f"carbon_source_id must be one of {', '.join(CARBON_SOURCE_IDS)}.")
    return normalized


def medium_condition_id_for(
    media_type: int,
    compatibility_mode: CompatibilityMode,
    carbon_source_id: str,
) -> str:
    carbon_source = validate_carbon_source_id(carbon_source_id)
    base = {
        2: "ynb_minimal",
        3: "ynb_core_aa",
        4: "ynb_core_aa",
        5: "ynb_all_aa",
    }.get(int(media_type), "ynb_core_aa")
    suffix = "matlab_compat" if compatibility_mode == "matlab_compat" else "corrected"
    candidate = f"{carbon_source}_{base}_{suffix}"
    try:
        load_medium_condition_spec(candidate)
        return candidate
    except KeyError:
        return candidate


def medium_condition_summary_for_inputs(inputs: PcSecPichiaInputs) -> dict[str, object]:
    try:
        summary = summarize_medium_condition(inputs.medium_condition_id)
    except KeyError:
        summary = {
            "condition_id": inputs.medium_condition_id,
            "carbon_source_id": inputs.carbon_source_id,
            "carbon_source_formulation": _carbon_source_formulation_payload(inputs.carbon_source_formulation),
            "legacy_media_type": inputs.media_type,
            "status": "active_runtime_condition",
        }
    summary["runtime_carbon_source_formulation"] = _carbon_source_formulation_payload(inputs.carbon_source_formulation)
    summary["runtime_key_bounds"] = _runtime_key_bounds(inputs.prepared_model)
    return summary


def _carbon_source_formulation_payload(formulation: CarbonSourceFormulation) -> dict[str, object]:
    return {
        "carbon_source_id": formulation.carbon_source_id,
        "active_uptake_reaction_ids": formulation.active_uptake_reaction_ids,
        "blocked_uptake_reaction_ids": formulation.blocked_uptake_reaction_ids,
        "candidate_growth_reaction_ids": formulation.candidate_growth_reaction_ids,
        "selected_growth_reaction_id": formulation.selected_growth_reaction_id,
        "carbon_objective_weights": formulation.carbon_objective_weights,
        "formulation_status": formulation.formulation_status,
        "warnings": formulation.warnings,
        "matlab_alignment_note": formulation.matlab_alignment_note,
    }


def _medium_condition_bounds(
    model: CobraModel,
    condition_id: str,
) -> dict[str, tuple[float | None, float | None]]:
    return _existing_bounds(
        model,
        {
            bound.reaction_id: (bound.lower_bound, bound.upper_bound)
            for bound in compose_medium_condition_bounds(condition_id)
        },
    )


def _runtime_key_bounds(model: CobraModel) -> dict[str, dict[str, float]]:
    key_reactions = (
        "BIOMASS",
        "BIOMASS_glyc",
        "BIOMASS_meoh",
        "Ex_glc_D",
        "Ex_glyc",
        "Ex_meoh",
        "Ex_o2",
    )
    payload: dict[str, dict[str, float]] = {}
    for reaction_id in key_reactions:
        index = model.reaction_index.get(reaction_id)
        if index is None:
            continue
        payload[reaction_id] = {
            "lower_bound": float(model.lb[index]),
            "upper_bound": float(model.ub[index]),
        }
    return payload


def _existing_bounds(
    model: CobraModel,
    changes: dict[str, tuple[float | None, float | None]],
) -> dict[str, tuple[float | None, float | None]]:
    return {reaction_id: bounds for reaction_id, bounds in changes.items() if reaction_id in model.reaction_index}


__all__ = [
    "CobraModel",
    "CARBON_SOURCE_IDS",
    "MINIMAL_EXCHANGE_REACTIONS",
    "PcSecPichiaInputs",
    "load_aa_stoichiometry",
    "load_combined_enzymedata",
    "load_carbon_source_formulation",
    "load_metabolic_enzymedata",
    "load_pcsec_pichia_inputs",
    "load_pcsec_pichia_model",
    "load_secretory_enzymedata",
    "medium_condition_id_for",
    "medium_condition_summary_for_inputs",
    "prepare_carbon_source_model",
    "prepare_glucose_model",
    "repo_root",
    "set_media_pp",
    "validate_carbon_source_id",
]
