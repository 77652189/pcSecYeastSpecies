from __future__ import annotations

import importlib
import importlib.util
from dataclasses import asdict, dataclass, field
from typing import Any, Literal

import numpy as np
from scipy import sparse

ShadowSense = Literal["maximize", "minimize"]


UNAVAILABLE_MESSAGE = "COBRApy is not installed; optional shadow FBA is unavailable."


@dataclass(frozen=True)
class CobraPyShadowFlux:
    reaction_id: str
    flux: float | None
    lower_bound: float | None = None
    upper_bound: float | None = None


@dataclass(frozen=True)
class CobraPyShadowModelBuildResult:
    available: bool
    status: str
    message: str
    model: Any | None = None
    model_summary: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "available": self.available,
            "status": self.status,
            "message": self.message,
            "model_summary": dict(self.model_summary),
        }


@dataclass(frozen=True)
class CobraPyShadowFBAResult:
    available: bool
    success: bool
    status: str
    message: str
    objective_reaction: str
    sense: ShadowSense
    objective_value: float | None = None
    fluxes: dict[str, float] = field(default_factory=dict)
    key_fluxes: tuple[CobraPyShadowFlux, ...] = ()
    model_summary: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "available": self.available,
            "success": self.success,
            "status": self.status,
            "message": self.message,
            "objective_reaction": self.objective_reaction,
            "sense": self.sense,
            "objective_value": self.objective_value,
            "fluxes": dict(self.fluxes),
            "key_fluxes": [asdict(item) for item in self.key_fluxes],
            "model_summary": dict(self.model_summary),
        }


@dataclass(frozen=True)
class CobraPyShadowComparison:
    comparable: bool
    status: str
    message: str
    objective_abs_diff: float | None = None
    objective_rel_diff: float | None = None
    within_tolerance: bool = False
    key_flux_diffs: dict[str, float | None] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def cobrapy_available() -> bool:
    return importlib.util.find_spec("cobra") is not None


def convert_to_cobrapy_model(
    model: Any,
    objective_reaction: str | None = None,
    sense: ShadowSense = "maximize",
) -> CobraPyShadowModelBuildResult:
    cobra = _import_cobra()
    if cobra is None:
        return CobraPyShadowModelBuildResult(
            available=False,
            status="unavailable",
            message=UNAVAILABLE_MESSAGE,
            model_summary=_model_summary(model),
        )
    if sense not in {"maximize", "minimize"}:
        return CobraPyShadowModelBuildResult(
            available=True,
            status="invalid_sense",
            message=f"Unsupported objective sense: {sense!r}.",
            model_summary=_model_summary(model),
        )
    if objective_reaction is not None and objective_reaction not in _reaction_index(model):
        return CobraPyShadowModelBuildResult(
            available=True,
            status="missing_objective",
            message=f"Objective reaction not found: {objective_reaction}",
            model_summary=_model_summary(model),
        )
    if _has_nonzero_rhs(model):
        return CobraPyShadowModelBuildResult(
            available=True,
            status="unsupported_nonzero_rhs",
            message="COBRApy shadow adapter only supports base GEM FBA with zero b vector; pcSec/nonzero RHS constraints are not converted.",
            model_summary=_model_summary(model),
        )

    cobra_model = cobra.Model(str(getattr(model, "model_id", None) or getattr(model, "source_file", "pcsec_shadow")))
    metabolites = {met_id: cobra.Metabolite(met_id) for met_id in model.mets}
    matrix = _s_matrix_csc(model)

    reactions = []
    for reaction_offset, reaction_id in enumerate(model.rxns):
        reaction = cobra.Reaction(reaction_id)
        reaction.lower_bound = float(model.lb[reaction_offset])
        reaction.upper_bound = float(model.ub[reaction_offset])
        stoichiometry = {}
        column = matrix.getcol(reaction_offset)
        for metabolite_offset, coefficient in zip(column.indices, column.data):
            if coefficient:
                stoichiometry[metabolites[model.mets[int(metabolite_offset)]]] = float(coefficient)
        if stoichiometry:
            reaction.add_metabolites(stoichiometry)
        reactions.append(reaction)
    cobra_model.add_reactions(reactions)

    if objective_reaction is not None:
        cobra_model.objective = cobra_model.reactions.get_by_id(objective_reaction)
        cobra_model.objective_direction = "max" if sense == "maximize" else "min"

    return CobraPyShadowModelBuildResult(
        available=True,
        status="converted",
        message="Converted base GEM stoichiometry, bounds, and objective for optional COBRApy shadow FBA.",
        model=cobra_model,
        model_summary=_model_summary(model),
    )


build_cobrapy_shadow_model = convert_to_cobrapy_model


def solve_cobrapy_shadow_fba(
    model: Any,
    objective_reaction: str,
    sense: ShadowSense = "maximize",
    key_reactions: tuple[str, ...] | list[str] = (),
) -> CobraPyShadowFBAResult:
    if objective_reaction not in _reaction_index(model):
        return CobraPyShadowFBAResult(
            available=cobrapy_available(),
            success=False,
            status="missing_objective",
            message=f"Objective reaction not found: {objective_reaction}",
            objective_reaction=objective_reaction,
            sense=sense,
            model_summary=_model_summary(model),
        )

    build = convert_to_cobrapy_model(model, objective_reaction=objective_reaction, sense=sense)
    if not build.available or build.model is None:
        return CobraPyShadowFBAResult(
            available=build.available,
            success=False,
            status=build.status,
            message=build.message,
            objective_reaction=objective_reaction,
            sense=sense,
            model_summary=build.model_summary,
        )
    if build.status != "converted":
        return CobraPyShadowFBAResult(
            available=True,
            success=False,
            status=build.status,
            message=build.message,
            objective_reaction=objective_reaction,
            sense=sense,
            model_summary=build.model_summary,
        )

    solution = build.model.optimize()
    success = str(getattr(solution, "status", "")).lower() == "optimal"
    fluxes: dict[str, float] = {}
    key_fluxes: list[CobraPyShadowFlux] = []
    if success:
        for reaction_id, value in getattr(solution, "fluxes", {}).items():
            value_float = float(value)
            if abs(value_float) > 1e-12:
                fluxes[str(reaction_id)] = value_float
        cobra_reaction_index = {reaction.id: reaction for reaction in build.model.reactions}
        for reaction_id in key_reactions:
            reaction = cobra_reaction_index.get(reaction_id)
            key_fluxes.append(
                CobraPyShadowFlux(
                    reaction_id=reaction_id,
                    flux=float(solution.fluxes[reaction_id]) if reaction_id in solution.fluxes.index else None,
                    lower_bound=float(reaction.lower_bound) if reaction is not None else None,
                    upper_bound=float(reaction.upper_bound) if reaction is not None else None,
                )
            )

    return CobraPyShadowFBAResult(
        available=True,
        success=success,
        status=str(getattr(solution, "status", "")),
        message="COBRApy shadow FBA solved base GEM model." if success else "COBRApy shadow FBA did not reach optimal status.",
        objective_reaction=objective_reaction,
        sense=sense,
        objective_value=float(getattr(solution, "objective_value", 0.0)) if success else None,
        fluxes=fluxes,
        key_fluxes=tuple(key_fluxes),
        model_summary=build.model_summary,
    )


def compare_shadow_fba(
    current_result: Any,
    shadow_result: CobraPyShadowFBAResult,
    objective_tolerance: float = 1e-7,
    relative_tolerance: float = 1e-6,
    key_reactions: tuple[str, ...] | list[str] = (),
) -> CobraPyShadowComparison:
    if not shadow_result.available:
        return CobraPyShadowComparison(
            comparable=False,
            status="unavailable",
            message=shadow_result.message,
        )
    if not shadow_result.success:
        return CobraPyShadowComparison(
            comparable=False,
            status=shadow_result.status,
            message="Shadow FBA did not solve successfully.",
        )
    current_objective = getattr(current_result, "objective_value", None)
    if current_objective is None or shadow_result.objective_value is None:
        return CobraPyShadowComparison(
            comparable=False,
            status="missing_objective_value",
            message="Both current and shadow results need objective values for comparison.",
        )

    objective_abs_diff = abs(float(current_objective) - float(shadow_result.objective_value))
    denominator = max(abs(float(current_objective)), abs(float(shadow_result.objective_value)), 1e-12)
    objective_rel_diff = objective_abs_diff / denominator
    current_fluxes = getattr(current_result, "fluxes", {}) or {}
    key_flux_diffs = {
        reaction_id: (
            abs(float(current_fluxes[reaction_id]) - float(shadow_result.fluxes[reaction_id]))
            if reaction_id in current_fluxes and reaction_id in shadow_result.fluxes
            else None
        )
        for reaction_id in key_reactions
    }
    within_tolerance = objective_abs_diff <= objective_tolerance or objective_rel_diff <= relative_tolerance
    return CobraPyShadowComparison(
        comparable=True,
        status="compared",
        message="Compared current SciPy/HiGHS baseline FBA with optional COBRApy shadow FBA.",
        objective_abs_diff=objective_abs_diff,
        objective_rel_diff=objective_rel_diff,
        within_tolerance=within_tolerance,
        key_flux_diffs=key_flux_diffs,
    )


def _import_cobra() -> Any | None:
    if not cobrapy_available():
        return None
    return importlib.import_module("cobra")


def _reaction_index(model: Any) -> dict[str, int]:
    reaction_index = getattr(model, "reaction_index")
    return dict(reaction_index() if callable(reaction_index) else reaction_index)


def _s_matrix_csc(model: Any) -> sparse.csc_matrix:
    matrix = getattr(model, "s_matrix")
    return matrix.tocsc() if sparse.issparse(matrix) else sparse.csc_matrix(matrix)


def _has_nonzero_rhs(model: Any) -> bool:
    b = getattr(model, "b", None)
    if b is None:
        return False
    return bool(np.any(np.abs(np.asarray(b, dtype=float).reshape(-1)) > 1e-12))


def _model_summary(model: Any) -> dict[str, object]:
    matrix = getattr(model, "s_matrix", None)
    shape = tuple(matrix.shape) if matrix is not None else None
    nnz = int(matrix.nnz) if hasattr(matrix, "nnz") else None
    return {
        "reaction_count": len(getattr(model, "rxns", ())),
        "metabolite_count": len(getattr(model, "mets", ())),
        "gene_count": len(getattr(model, "genes", ())),
        "stoichiometric_shape": shape,
        "stoichiometric_nnz": nnz,
        "supports_pcsec_constraints": False,
        "gene_reaction_rule_semantics": "not_converted_in_phase1",
    }


__all__ = [
    "CobraPyShadowComparison",
    "CobraPyShadowFBAResult",
    "CobraPyShadowFlux",
    "CobraPyShadowModelBuildResult",
    "build_cobrapy_shadow_model",
    "cobrapy_available",
    "compare_shadow_fba",
    "convert_to_cobrapy_model",
    "solve_cobrapy_shadow_fba",
]
