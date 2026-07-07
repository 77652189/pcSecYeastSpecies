from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from pcsec_pichia.analysis.shadow_lp.constraint_spec import ShadowConstraintConfig
from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.probe import build_supported_target_model, build_target_enzymedata
from pcsec_pichia.targets import TargetSpec, load_builtin_targets


@dataclass(frozen=True)
class ShadowTargetPreparation:
    """Target-extended pcSec inputs for symbolic shadow LP constraint builders."""

    target_id: str
    target: Any
    model: Any
    fixed_model: Any
    exchange_reaction_id: str
    metabolic: Any
    secretory: Any
    combined: Any
    added_reaction_count: int
    added_metabolite_count: int
    warnings: tuple[str, ...] = ()


def fixed_growth_bounds(
    model: Any,
    growth_rate: float,
) -> dict[str, tuple[float | None, float | None]]:
    """Return the pcSec fixed-growth bound overrides without applying a solve."""

    reaction_index = _reaction_index(model)
    bounds: dict[str, tuple[float | None, float | None]] = {"BIOMASS": (float(growth_rate), float(growth_rate))}
    for reaction_id in ("BIOMASS_glyc", "BIOMASS_meoh"):
        if reaction_id in reaction_index:
            bounds[reaction_id] = (0.0, 0.0)
    return bounds


def prepare_shadow_target(
    inputs: PcSecPichiaInputs,
    target: TargetSpec,
    config: ShadowConstraintConfig | None = None,
) -> ShadowTargetPreparation:
    """Prepare target-extended model/enzyme inputs for symbolic constraint builders."""

    resolved_config = config or ShadowConstraintConfig()
    build = build_supported_target_model(inputs.prepared_model, target, inputs.amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        raise ValueError(f"Target build failed for {target.target_id}: {build.reason}")

    target_enzymedata = build_target_enzymedata(target, build.model, inputs.secretory)
    target_secretory = inputs.secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = _with_target_enzymedata(inputs.combined, target_enzymedata)
    fixed_model = _apply_bounds(build.model, fixed_growth_bounds(build.model, resolved_config.growth_rate))
    return ShadowTargetPreparation(
        target_id=target.target_id,
        target=target,
        model=build.model,
        fixed_model=fixed_model,
        exchange_reaction_id=build.exchange_reaction_id,
        metabolic=inputs.metabolic,
        secretory=target_secretory,
        combined=target_combined,
        added_reaction_count=int(build.added_reaction_count),
        added_metabolite_count=int(build.added_metabolite_count),
    )


def prepare_builtin_shadow_target(
    target_id: str,
    root: Path | None = None,
    config: ShadowConstraintConfig | None = None,
) -> ShadowTargetPreparation:
    """Load pcSec inputs and prepare one built-in target without running a solver."""

    resolved_root = root or _repo_root()
    targets = {target.target_id: target for target in load_builtin_targets(resolved_root)}
    try:
        target = targets[target_id]
    except KeyError as exc:
        raise KeyError(f"Unknown built-in target: {target_id}") from exc
    inputs = load_pcsec_pichia_inputs(resolved_root)
    return prepare_shadow_target(inputs, target, config=config)


def _with_target_enzymedata(combined: Any, target_enzymedata: Any) -> Any:
    if hasattr(combined, "with_target"):
        return combined.with_target(target_enzymedata)
    if hasattr(combined, "with_target_proteins"):
        return combined.with_target_proteins(target_enzymedata)
    raise TypeError("combined enzyme data must provide with_target() or with_target_proteins().")


def _apply_bounds(model: Any, bounds: Mapping[str, tuple[float | None, float | None]]) -> Any:
    if hasattr(model, "with_bounds"):
        return model.with_bounds(dict(bounds))
    if hasattr(model, "with_reaction_bounds"):
        return model.with_reaction_bounds(dict(bounds))
    raise TypeError("pcSec model must provide with_bounds() or with_reaction_bounds().")


def _reaction_index(model: Any) -> Mapping[str, int]:
    reaction_index = getattr(model, "reaction_index")
    if callable(reaction_index):
        reaction_index = reaction_index()
    return reaction_index


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[5]
