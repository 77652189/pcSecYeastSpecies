from __future__ import annotations

from typing import Any, Mapping

from pcsec_pichia.analysis.shadow_lp.backends import ScipyHighsBackend, SolverBackend
from pcsec_pichia.analysis.shadow_lp.constraint_spec import ShadowConstraintConfig
from pcsec_pichia.analysis.shadow_lp.ladder import ShadowLadderResult, run_shadow_ladder_for_prepared_target
from pcsec_pichia.analysis.shadow_lp.model_adapter import ShadowTargetPreparation, fixed_growth_bounds
from pcsec_pichia.probe import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    TargetSpec,
    build_supported_target_model,
    build_target_enzymedata,
)
from pcsec_pichia.simulation import SecretionSimulationResult


SHADOW_CONSTRAINT_COUNT_KEYS: tuple[str, ...] = (
    "stoichiometric",
    "metabolic_coupling",
    "secretory_coupling",
    "protein_mass",
    "proteasome",
    "ribosome_assembly",
    "ribosome_translation",
    "misfolding",
    "mitochondrial",
    "eq_total",
    "ub_total",
)


def solve_shadow_secretion_capacity(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    backend: SolverBackend | None = None,
    solver_options: Mapping[str, Any] | None = None,
) -> SecretionSimulationResult:
    """Solve fixed-growth secretion capacity through the formal shadow LP path."""

    config = ShadowConstraintConfig(
        growth_rate=float(growth_rate),
        enable_ribosome_translation=bool(write_ribosome_translation_constraint),
        enable_misfolding=bool(write_misfolding_constraints),
    )
    try:
        prep = _prepare_shadow_target_from_components(
            model,
            target,
            amino_acids,
            metabolic,
            secretory,
            combined,
            config=config,
        )
    except _ShadowTargetBuildError as exc:
        return _build_failed_result(target, growth_rate, exc.status, exc.message)
    except Exception as exc:
        return _build_failed_result(target, growth_rate, "exception", str(exc))

    resolved_backend = backend or ScipyHighsBackend()
    ladder = run_shadow_ladder_for_prepared_target(
        prep,
        config=config,
        backend=resolved_backend,
        solver_options=solver_options,
    )
    final = ladder.final_layer
    success = bool(final.success)
    objective = final.objective if success else None
    return SecretionSimulationResult(
        success=success,
        target_id=target.target_id,
        objective_value=objective,
        growth_rate=float(growth_rate),
        secretion_flux=objective if success else None,
        status=str(final.status),
        message=str(final.message),
        constraint_counts=_shadow_constraint_counts(final.backend_metadata, final.eq_constraint_count, final.ub_constraint_count),
        result_status="shadow_lp_capacity",
        target_parameter_status=_target_parameter_status(target),
        matlab_alignment_status="shadow_validation_pending",
        exchange_reaction_id=prep.exchange_reaction_id,
        build_status="supported",
        lp_sensitivity=None,
        key_fluxes={reaction_id: float(value) for reaction_id, value in final.key_fluxes.items() if value is not None},
        growth_reaction_id="BIOMASS",
        open_growth_reaction_ids=("BIOMASS",),
        growth_reaction_status="single_growth_reaction",
        warnings=_shadow_warnings(ladder),
        solver_mode="shadow",
        shadow_metadata=_shadow_metadata(ladder, final, resolved_backend),
    )


class _ShadowTargetBuildError(Exception):
    def __init__(self, status: str, message: str) -> None:
        super().__init__(message)
        self.status = status
        self.message = message


def _prepare_shadow_target_from_components(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    config: ShadowConstraintConfig,
) -> ShadowTargetPreparation:
    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        raise _ShadowTargetBuildError(str(build.status), str(build.reason))
    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = _with_target_enzymedata(combined, target_enzymedata)
    fixed_model = _apply_bounds(build.model, fixed_growth_bounds(build.model, config.growth_rate))
    return ShadowTargetPreparation(
        target_id=target.target_id,
        target=target,
        model=build.model,
        fixed_model=fixed_model,
        exchange_reaction_id=build.exchange_reaction_id,
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        added_reaction_count=int(build.added_reaction_count),
        added_metabolite_count=int(build.added_metabolite_count),
    )


def _shadow_constraint_counts(
    backend_metadata: Mapping[str, Any],
    eq_constraint_count: int,
    ub_constraint_count: int,
) -> dict[str, int]:
    layer_counts = backend_metadata.get("layer_counts", {})
    counts = {
        "stoichiometric": int(backend_metadata.get("stoichiometric_constraint_count", 0)),
        "metabolic_coupling": int(layer_counts.get("metabolic_coupling", 0)),
        "secretory_coupling": int(layer_counts.get("secretory_coupling", 0)),
        "protein_mass": int(layer_counts.get("protein_mass", 0)),
        "proteasome": int(layer_counts.get("proteasome", 0)),
        "ribosome_assembly": int(layer_counts.get("ribosome_assembly", 0)),
        "ribosome_translation": int(layer_counts.get("ribosome_translation", 0)),
        "misfolding": int(layer_counts.get("misfolding", 0)),
        "mitochondrial": int(layer_counts.get("mitochondrial", 0)),
        "eq_total": int(eq_constraint_count),
        "ub_total": int(ub_constraint_count),
    }
    return {key: counts[key] for key in SHADOW_CONSTRAINT_COUNT_KEYS}


def _shadow_metadata(
    ladder: ShadowLadderResult,
    final: Any,
    backend: SolverBackend,
) -> dict[str, Any]:
    return {
        "solver_mode": "shadow",
        "backend": backend.name,
        "default_large_model_backend": "ScipyHighsBackend",
        "canonical_final_layer": ladder.final_layer_id,
        "final_layer": final.to_dict(),
        "layer_order": tuple(layer.layer_id for layer in ladder.layers),
        "skipped_layers": {
            layer.layer_id: dict(layer.skipped_layers)
            for layer in ladder.layers
            if layer.skipped_layers
        },
        "shadow_ladder_warning_count": len(_shadow_warnings(ladder)),
        "reference_solver_used": False,
    }


def _shadow_warnings(ladder: ShadowLadderResult) -> tuple[str, ...]:
    warnings: list[str] = list(ladder.warnings)
    for layer in ladder.layers:
        warnings.extend(str(warning) for warning in layer.warnings)
    return tuple(dict.fromkeys(warnings))


def _build_failed_result(
    target: TargetSpec,
    growth_rate: float,
    status: str,
    message: str,
) -> SecretionSimulationResult:
    return SecretionSimulationResult(
        success=False,
        target_id=target.target_id,
        objective_value=None,
        growth_rate=float(growth_rate),
        secretion_flux=None,
        status=str(status),
        message=str(message),
        constraint_counts={key: 0 for key in SHADOW_CONSTRAINT_COUNT_KEYS},
        result_status="shadow_lp_capacity",
        target_parameter_status=_target_parameter_status(target),
        matlab_alignment_status="shadow_validation_pending",
        exchange_reaction_id=None,
        build_status=str(status),
        lp_sensitivity=None,
        key_fluxes={},
        warnings=(str(message),),
        solver_mode="shadow",
        shadow_metadata={
            "solver_mode": "shadow",
            "reference_solver_used": False,
            "build_status": str(status),
        },
    )


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


def _target_parameter_status(target: TargetSpec) -> str:
    pending_targets = {"hlf"}
    if target.target_id.lower() in pending_targets or target.protein_id.lower() in pending_targets:
        return "draft_matlab_alignment_pending"
    return "draft"


__all__ = [
    "SHADOW_CONSTRAINT_COUNT_KEYS",
    "solve_shadow_secretion_capacity",
]
