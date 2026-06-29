from __future__ import annotations

from dataclasses import dataclass, field

from pcsec_pichia.core.pichia_model import PichiaModel


MINIMAL_EXCHANGES = (
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

YNB_VITAMINS = (
    "Ex_btn",
    "Ex_thm",
    "Ex_4abz",
    "Ex_pnto_R",
    "Ex_inost",
    "Ex_nac",
    "Ex_ribflv",
)

CORE_AMINO_ACIDS = (
    "Ex_arg_L",
    "Ex_asp_L",
    "Ex_glu_L",
    "Ex_gly",
    "Ex_his_L",
    "Ex_ile_L",
    "Ex_leu_L",
    "Ex_lys_L",
    "Ex_met_L",
    "Ex_phe_L",
    "Ex_thr_L",
    "Ex_trp_L",
    "Ex_tyr_L",
    "Ex_val_L",
    "Ex_ura",
)

ALL_AMINO_ACIDS = CORE_AMINO_ACIDS + (
    "Ex_ala_L",
    "Ex_asn_L",
    "Ex_cys_L",
    "Ex_gln_L",
    "Ex_pro_L",
    "Ex_ser_L",
)

OPTIONAL_GLUCOSE_UPPER_REACTIONS = (
    "LIPIDS",
    "PROTEINS",
    "STEROLS",
)

GLYCEROL_BLOCK_REACTIONS = (
    "BIOMASS_glyc",
    "LIPIDS_glyc",
    "PROTEINS_glyc",
    "STEROLS_glyc",
)

METHANOL_BLOCK_REACTIONS = (
    "BIOMASS_meoh",
    "LIPIDS_meoh",
    "PROTEINS_meoh",
    "STEROLS_meoh",
)


@dataclass(frozen=True)
class MediaConfiguration:
    model: PichiaModel
    media_type: int
    missing_reactions: tuple[str, ...] = field(default_factory=tuple)
    configured_reactions: tuple[str, ...] = field(default_factory=tuple)


def exchange_reaction_ids(model: PichiaModel) -> tuple[str, ...]:
    return tuple(reaction_id for reaction_id in model.rxns if reaction_id.startswith("Ex_"))


def set_media_pp(
    model: PichiaModel,
    media_type: int = 1,
    open_minimal_exchanges: bool = False,
) -> MediaConfiguration:
    if media_type not in {1, 2, 3, 4, 5}:
        raise ValueError("media_type must be one of 1, 2, 3, 4, 5.")

    changes: dict[str, tuple[float | None, float | None]] = {
        reaction_id: (0.0, 1000.0) for reaction_id in exchange_reaction_ids(model)
    }
    missing: list[str] = []
    configured: list[str] = list(changes)

    # Match the current MATLAB setMediaPP.m behavior: its minimal-media helper
    # call does not assign the returned model, so these exchanges remain closed
    # until individual simulation scripts reopen the required exchanges.
    if open_minimal_exchanges:
        _set_lower(changes, missing, configured, model, MINIMAL_EXCHANGES, -1000.0)
    if media_type == 2:
        _set_lower(changes, missing, configured, model, YNB_VITAMINS, -2.0)
    elif media_type in {3, 4}:
        _set_lower(changes, missing, configured, model, YNB_VITAMINS, -2.0)
        _set_lower(changes, missing, configured, model, CORE_AMINO_ACIDS, -0.08)
    elif media_type == 5:
        _set_lower(changes, missing, configured, model, YNB_VITAMINS, -2.0)
        _set_lower(changes, missing, configured, model, ALL_AMINO_ACIDS, -0.08)

    return MediaConfiguration(
        model=model.with_reaction_bounds(changes),
        media_type=media_type,
        missing_reactions=tuple(dict.fromkeys(missing)),
        configured_reactions=tuple(dict.fromkeys(configured)),
    )


def apply_glucose_reference_conditions(
    model: PichiaModel,
    media_type: int = 2,
    block_misfold_dilution: bool = True,
    open_minimal_exchanges: bool = False,
) -> MediaConfiguration:
    media = set_media_pp(model, media_type, open_minimal_exchanges=open_minimal_exchanges)
    configured_model = media.model
    changes: dict[str, tuple[float | None, float | None]] = {
        "Ex_glc_D": (-1000.0, None),
        "BIOMASS": (None, 1000.0),
        "Ex_glyc": (0.0, None),
        "Ex_meoh": (0.0, None),
        "Ex_o2": (-1000.0, None),
    }
    missing = list(media.missing_reactions)
    configured = list(media.configured_reactions)
    for reaction_id in OPTIONAL_GLUCOSE_UPPER_REACTIONS:
        _set_optional_upper(changes, configured_model, reaction_id, 1000.0)
    for reaction_id in GLYCEROL_BLOCK_REACTIONS + METHANOL_BLOCK_REACTIONS:
        _set_optional_both(changes, configured_model, reaction_id, 0.0)
    for reaction_id in list(changes):
        if reaction_id not in configured_model.reaction_index:
            missing.append(reaction_id)
            changes.pop(reaction_id)
            continue
        configured.append(reaction_id)
    configured_model = configured_model.with_reaction_bounds(changes)
    if block_misfold_dilution:
        configured_model = _block_reactions_containing(configured_model, "dilution_misfolding")
    return MediaConfiguration(
        model=configured_model,
        media_type=media_type,
        missing_reactions=tuple(dict.fromkeys(missing)),
        configured_reactions=tuple(dict.fromkeys(configured)),
    )


def set_fixed_growth_rate(model: PichiaModel, mu: float) -> PichiaModel:
    return model.change_rxn_bounds("BIOMASS", lower=mu, upper=mu)


def _set_lower(
    changes: dict[str, tuple[float | None, float | None]],
    missing: list[str],
    configured: list[str],
    model: PichiaModel,
    reaction_ids: tuple[str, ...],
    lower_bound: float,
) -> None:
    for reaction_id in reaction_ids:
        if reaction_id not in model.reaction_index:
            missing.append(reaction_id)
            continue
        _lower, upper = changes.get(reaction_id, (None, None))
        changes[reaction_id] = (lower_bound, upper)
        configured.append(reaction_id)


def _set_optional_upper(
    changes: dict[str, tuple[float | None, float | None]],
    model: PichiaModel,
    reaction_id: str,
    upper_bound: float,
) -> None:
    if reaction_id in model.reaction_index:
        lower, _upper = changes.get(reaction_id, (None, None))
        changes[reaction_id] = (lower, upper_bound)


def _set_optional_both(
    changes: dict[str, tuple[float | None, float | None]],
    model: PichiaModel,
    reaction_id: str,
    value: float,
) -> None:
    if reaction_id in model.reaction_index:
        changes[reaction_id] = (value, value)


def _block_reactions_containing(model: PichiaModel, text: str) -> PichiaModel:
    changes = {reaction_id: (None, 0.0) for reaction_id in model.rxns if text in reaction_id}
    if not changes:
        return model
    return model.with_reaction_bounds(changes)
