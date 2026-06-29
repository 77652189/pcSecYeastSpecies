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
    configured = set_media_pp(model, media_type=media_type, compatibility_mode=compatibility_mode)
    return configured.with_bounds(
        {
            "Ex_glc_D": (-1000.0, None),
            "BIOMASS": (None, 1000.0),
            "Ex_glyc": (0.0, None),
            "BIOMASS_glyc": (0.0, 0.0),
            "Ex_meoh": (0.0, None),
            "BIOMASS_meoh": (0.0, 0.0),
            "Ex_o2": (-1000.0, None),
        }
    )


def load_pcsec_pichia_inputs(
    root: Path | None = None,
    media_type: int = 4,
    compatibility_mode: CompatibilityMode = DEFAULT_COMPATIBILITY_MODE,
) -> PcSecPichiaInputs:
    """Load the prototype-backed pcSecPichia input bundle for formal workflows."""

    mode = validate_compatibility_mode(compatibility_mode)
    resolved_root = root or repo_root()
    model = load_pcsec_pichia_model(resolved_root)
    return PcSecPichiaInputs(
        root=resolved_root,
        model=model,
        prepared_model=prepare_glucose_model(model, media_type=media_type, compatibility_mode=mode),
        amino_acids=load_aa_stoichiometry(resolved_root),
        metabolic=load_metabolic_enzymedata(resolved_root),
        secretory=load_secretory_enzymedata(resolved_root),
        combined=load_combined_enzymedata(resolved_root),
        compatibility_mode=mode,
        media_type=media_type,
    )


__all__ = [
    "CobraModel",
    "MINIMAL_EXCHANGE_REACTIONS",
    "PcSecPichiaInputs",
    "load_aa_stoichiometry",
    "load_combined_enzymedata",
    "load_metabolic_enzymedata",
    "load_pcsec_pichia_inputs",
    "load_pcsec_pichia_model",
    "load_secretory_enzymedata",
    "prepare_glucose_model",
    "repo_root",
    "set_media_pp",
]
