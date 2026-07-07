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


@dataclass(frozen=True)
class CobraOptlangBackend:
    """Optional tiny-LP semantic backend backed by COBRApy/optlang.

    This backend is intentionally limited to small validation fixtures. The
    full pcSec shadow ladder should keep using ScipyHighsBackend.
    """

    name: str = "cobra-optlang"
    supports_duals: bool = False
    supports_time_limit: bool = False
    default_max_variables: int = 128
    default_max_constraints: int = 512

    def available(self) -> bool:
        try:
            _import_cobra()
        except ImportError:
            return False
        return True

    def solve(self, problem: LPProblem, options: Mapping[str, Any] | None = None) -> SolverResult:
        resolved_options = dict(options or {})
        total_started = time.perf_counter()
        availability = self.available()
        if not availability:
            return SolverResult(
                success=False,
                status="unavailable",
                objective=None,
                fluxes={},
                message="COBRApy/optlang is not available in this environment.",
                timings={"total_seconds": time.perf_counter() - total_started},
                backend_metadata={"backend": self.name, "available": False},
            )

        assemble_started = time.perf_counter()
        assembled = assemble_lp_problem(problem)
        assemble_seconds = time.perf_counter() - assemble_started
        max_variables = int(resolved_options.get("max_variables", self.default_max_variables))
        max_constraints = int(resolved_options.get("max_constraints", self.default_max_constraints))
        if assembled.diagnostics.variable_count > max_variables or assembled.diagnostics.constraint_count > max_constraints:
            return SolverResult(
                success=False,
                status="too_large",
                objective=None,
                fluxes={},
                key_fluxes={reaction_id: None for reaction_id in assembled.key_reaction_ids},
                message="CobraOptlangBackend is limited to tiny validation LPs.",
                timings={
                    "assemble_seconds": assemble_seconds,
                    "total_seconds": time.perf_counter() - total_started,
                },
                backend_metadata={
                    "backend": self.name,
                    "available": True,
                    "max_variables": max_variables,
                    "max_constraints": max_constraints,
                    "variable_count": assembled.diagnostics.variable_count,
                    "constraint_count": assembled.diagnostics.constraint_count,
                },
            )

        cobra = _import_cobra()
        model = cobra.Model("shadow_lp_tiny_validation")
        reaction_ids_by_original = _unique_cobra_ids(assembled.reaction_ids)
        reactions = []
        for reaction_id, (lower, upper) in zip(assembled.reaction_ids, assembled.bounds):
            reaction = cobra.Reaction(reaction_ids_by_original[reaction_id])
            reaction.name = reaction_id
            reaction.lower_bound = -1_000_000.0 if lower is None else float(lower)
            reaction.upper_bound = 1_000_000.0 if upper is None else float(upper)
            reactions.append(reaction)
        model.add_reactions(reactions)
        reaction_by_id = {original: reaction for original, reaction in zip(assembled.reaction_ids, reactions)}

        constraints = []
        for row_index in range(assembled.A_eq.shape[0]):
            expr = _optlang_row_expression(assembled.A_eq.getrow(row_index), assembled.reaction_ids, reaction_by_id)
            rhs = float(assembled.b_eq[row_index])
            constraints.append(model.problem.Constraint(expr, lb=rhs, ub=rhs, name=f"eq_{row_index}"))
        if assembled.A_ub is not None and assembled.b_ub is not None:
            for row_index in range(assembled.A_ub.shape[0]):
                expr = _optlang_row_expression(assembled.A_ub.getrow(row_index), assembled.reaction_ids, reaction_by_id)
                constraints.append(
                    model.problem.Constraint(
                        expr,
                        lb=None,
                        ub=float(assembled.b_ub[row_index]),
                        name=f"ub_{row_index}",
                    )
                )
        model.add_cons_vars(constraints)

        objective_expr = 0
        for reaction_id, coefficient in zip(assembled.reaction_ids, assembled.objective_vector):
            if coefficient:
                objective_expr += float(coefficient) * reaction_by_id[reaction_id].flux_expression
        direction = "max" if assembled.objective_sense == "maximize" else "min"
        model.objective = model.problem.Objective(objective_expr, direction=direction)

        solve_started = time.perf_counter()
        solution = model.optimize()
        solve_seconds = time.perf_counter() - solve_started
        success = str(solution.status).lower() == "optimal"
        fluxes: dict[str, float] = {}
        key_fluxes = {reaction_id: None for reaction_id in assembled.key_reaction_ids}
        objective_value = None
        if success:
            objective_value = float(solution.objective_value)
            threshold = float(resolved_options.get("flux_threshold", 1e-12))
            for reaction_id, reaction in reaction_by_id.items():
                value = float(solution.fluxes[reaction.id])
                if abs(value) > threshold:
                    fluxes[reaction_id] = value
            key_fluxes = {
                reaction_id: (
                    float(solution.fluxes[reaction_by_id[reaction_id].id])
                    if reaction_id in reaction_by_id
                    else None
                )
                for reaction_id in assembled.key_reaction_ids
            }
        return SolverResult(
            success=success,
            status=str(solution.status),
            objective=objective_value,
            fluxes=fluxes,
            key_fluxes=key_fluxes,
            message=str(solution.status),
            timings={
                "assemble_seconds": assemble_seconds,
                "solve_seconds": solve_seconds,
                "total_seconds": time.perf_counter() - total_started,
            },
            backend_metadata={
                "backend": self.name,
                "available": True,
                "status": str(solution.status),
                "variable_count": assembled.diagnostics.variable_count,
                "constraint_count": assembled.diagnostics.constraint_count,
                "limited_to_tiny_lp": True,
            },
        )


def _import_cobra() -> Any:
    try:
        cobra = importlib.import_module("cobra")
        importlib.import_module("optlang")
    except Exception as exc:
        raise ImportError("COBRApy/optlang is not available.") from exc
    return cobra


def _cobra_id(reaction_id: str) -> str:
    cleaned = [character if character.isalnum() or character == "_" else "_" for character in reaction_id]
    text = "".join(cleaned).strip("_")
    return f"R_{text or 'reaction'}"


def _unique_cobra_ids(reaction_ids: tuple[str, ...]) -> Mapping[str, str]:
    used: set[str] = set()
    result: dict[str, str] = {}
    for index, reaction_id in enumerate(reaction_ids):
        base_id = _cobra_id(reaction_id)
        cobra_id = base_id
        if cobra_id in used:
            cobra_id = f"{base_id}_{index}"
        while cobra_id in used:
            cobra_id = f"{base_id}_{index}_{len(used)}"
        used.add(cobra_id)
        result[reaction_id] = cobra_id
    return result


def _optlang_row_expression(row: Any, reaction_ids: tuple[str, ...], reaction_by_id: Mapping[str, Any]) -> Any:
    expr = 0
    coo = row.tocoo()
    for column, value in zip(coo.col, coo.data):
        if value:
            reaction_id = reaction_ids[int(column)]
            expr += float(value) * reaction_by_id[reaction_id].flux_expression
    return expr


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
