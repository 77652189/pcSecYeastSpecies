from __future__ import annotations

import importlib
import time
from dataclasses import asdict, dataclass, field
from typing import Any, Mapping, Protocol, runtime_checkable

import numpy as np

from pcsec_pichia.analysis.shadow_lp.constraint_spec import LPProblem
from pcsec_pichia.analysis.shadow_lp.lp_problem import assemble_lp_problem


@dataclass(frozen=True)
class SolverResult:
    """Backend-neutral result for a shadow LP solve."""

    success: bool
    status: str
    objective: float | None
    fluxes: Mapping[str, float]
    message: str
    key_fluxes: Mapping[str, float | None] = field(default_factory=dict)
    timings: Mapping[str, float] = field(default_factory=dict)
    duals: Mapping[str, Any] | None = None
    backend_metadata: Mapping[str, Any] = field(default_factory=dict)


@runtime_checkable
class SolverBackend(Protocol):
    """Protocol implemented by replaceable shadow LP solver backends."""

    name: str
    supports_duals: bool
    supports_time_limit: bool

    def available(self) -> bool:
        """Return whether this backend can solve in the current environment."""

    def solve(self, problem: LPProblem, options: Mapping[str, Any] | None = None) -> SolverResult:
        """Solve a backend-neutral LP problem."""


@dataclass(frozen=True)
class ScipyHighsBackend:
    """Solve assembled shadow LP problems with scipy.optimize.linprog(method='highs')."""

    name: str = "scipy-highs"
    supports_duals: bool = True
    supports_time_limit: bool = True

    def available(self) -> bool:
        try:
            _import_linprog()
        except ImportError:
            return False
        return True

    def solve(self, problem: LPProblem, options: Mapping[str, Any] | None = None) -> SolverResult:
        resolved_options = dict(options or {})
        total_started = time.perf_counter()
        assemble_started = time.perf_counter()
        assembled = assemble_lp_problem(problem)
        assemble_seconds = time.perf_counter() - assemble_started

        linprog = _import_linprog()
        scipy_options: dict[str, Any] = {
            "presolve": bool(resolved_options.get("presolve", True)),
            "disp": bool(resolved_options.get("disp", False)),
        }
        time_limit = resolved_options.get("time_limit", resolved_options.get("time_limit_seconds"))
        if time_limit is not None:
            scipy_options["time_limit"] = float(time_limit)

        scipy_objective = assembled.objective_vector
        if assembled.objective_sense == "maximize":
            scipy_objective = -scipy_objective
        elif assembled.objective_sense != "minimize":
            raise ValueError(f"Unsupported objective sense: {assembled.objective_sense}")

        solve_started = time.perf_counter()
        result = linprog(
            c=scipy_objective,
            A_eq=assembled.A_eq,
            b_eq=assembled.b_eq,
            A_ub=assembled.A_ub,
            b_ub=assembled.b_ub,
            bounds=list(assembled.bounds),
            method="highs",
            options=scipy_options,
        )
        solve_seconds = time.perf_counter() - solve_started

        objective_value = None
        fluxes: dict[str, float] = {}
        key_fluxes = {reaction_id: None for reaction_id in assembled.key_reaction_ids}
        if bool(result.success) and result.x is not None:
            x = np.asarray(result.x, dtype=float)
            objective_value = float(np.dot(assembled.objective_vector, x))
            threshold = float(resolved_options.get("flux_threshold", 1e-12))
            fluxes = {
                reaction_id: float(x[index])
                for reaction_id, index in assembled.reaction_index.items()
                if abs(float(x[index])) > threshold
            }
            key_fluxes = {
                reaction_id: (
                    float(x[assembled.reaction_index[reaction_id]])
                    if reaction_id in assembled.reaction_index
                    else None
                )
                for reaction_id in assembled.key_reaction_ids
            }

        timings = {
            "assemble_seconds": assemble_seconds,
            "solve_seconds": solve_seconds,
            "total_seconds": time.perf_counter() - total_started,
        }
        return SolverResult(
            success=bool(result.success),
            status=str(result.status),
            objective=objective_value,
            fluxes=fluxes,
            key_fluxes=key_fluxes,
            message=str(result.message),
            timings=timings,
            duals=_extract_duals(result) if bool(resolved_options.get("include_duals", False)) else None,
            backend_metadata={
                "backend": self.name,
                "scipy_method": "highs",
                "scipy_options": scipy_options,
                "status_code": int(result.status),
                "nit": _optional_int(getattr(result, "nit", None)),
                "crossover_nit": _optional_int(getattr(result, "crossover_nit", None)),
                "objective_sense": assembled.objective_sense,
                "variable_count": assembled.diagnostics.variable_count,
                "eq_constraint_count": assembled.diagnostics.eq_constraint_count,
                "ub_constraint_count": assembled.diagnostics.ub_constraint_count,
                "constraint_count": assembled.diagnostics.constraint_count,
                "stoichiometric_constraint_count": assembled.diagnostics.stoichiometric_constraint_count,
                "layer_counts": dict(assembled.diagnostics.layer_counts),
                "block_order": assembled.diagnostics.block_order,
                "constraint_order": tuple(asdict(item) for item in assembled.diagnostics.constraint_order[:16]),
            },
        )


def _import_linprog() -> Any:
    optimize = importlib.import_module("scipy.optimize")
    return optimize.linprog


def _extract_duals(result: Any) -> Mapping[str, Any] | None:
    duals: dict[str, Any] = {}
    for key in ("eqlin", "ineqlin", "lower", "upper"):
        value = getattr(result, key, None)
        if value is None:
            continue
        marginals = getattr(value, "marginals", None)
        residual = getattr(value, "residual", None)
        duals[key] = {
            "marginals": None if marginals is None else np.asarray(marginals, dtype=float),
            "residual": None if residual is None else np.asarray(residual, dtype=float),
        }
    return duals or None


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    return int(value)
