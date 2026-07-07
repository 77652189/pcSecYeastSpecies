from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from scipy import sparse

from pcsec_pichia.analysis.shadow_lp.constraint_builders import (
    REFERENCE_LAYER_ORDER,
    build_shadow_constraint_blocks,
)
from pcsec_pichia.analysis.shadow_lp.constraint_spec import ConstraintBlock, ConstraintSpec, LPProblem, OptimizationSense
from pcsec_pichia.analysis.shadow_lp.model_adapter import ShadowTargetPreparation


RESOURCE_LAYER_BLOCK_IDS: tuple[str, ...] = (
    "metabolic_coupling",
    "secretory_coupling",
    "protein_mass",
    "proteasome",
    "ribosome_assembly",
    "mitochondrial",
)

LADDER_LAYER_BLOCK_IDS: Mapping[str, tuple[str, ...]] = {
    "base_cobrapy_fba": (),
    "fixed_growth": (),
    "metabolic_coupling": ("metabolic_coupling",),
    "secretory_coupling": ("metabolic_coupling", "secretory_coupling"),
    "protein_mass": ("metabolic_coupling", "secretory_coupling", "protein_mass"),
    "resource": RESOURCE_LAYER_BLOCK_IDS,
    "mitochondrial": RESOURCE_LAYER_BLOCK_IDS,
}


@dataclass(frozen=True)
class ConstraintOrderEntry:
    """One assembled symbolic constraint row and its final matrix location."""

    row_kind: str
    row_index: int
    layer_id: str
    constraint_name: str
    original_sense: str
    term_count: int
    rhs: float


@dataclass(frozen=True)
class LPAssemblyDiagnostics:
    """Auditable metadata produced while assembling symbolic constraints."""

    variable_count: int
    stoichiometric_constraint_count: int
    eq_constraint_count: int
    ub_constraint_count: int
    constraint_count: int
    objective_reaction_ids: tuple[str, ...]
    layer_counts: Mapping[str, int]
    block_order: tuple[str, ...]
    constraint_order: tuple[ConstraintOrderEntry, ...]


@dataclass(frozen=True)
class AssembledLPProblem:
    """Sparse matrix form consumed by solver backends."""

    reaction_ids: tuple[str, ...]
    reaction_index: Mapping[str, int]
    bounds: tuple[tuple[float | None, float | None], ...]
    objective_vector: np.ndarray
    objective_sense: OptimizationSense
    A_eq: sparse.csr_matrix
    b_eq: np.ndarray
    A_ub: sparse.csr_matrix | None
    b_ub: np.ndarray | None
    key_reaction_ids: tuple[str, ...]
    diagnostics: LPAssemblyDiagnostics
    metadata: Mapping[str, Any] = field(default_factory=dict)


def lp_problem_from_model(
    model: Any,
    objective_reaction_id: str,
    constraint_blocks: tuple[ConstraintBlock, ...] = (),
    objective_sense: OptimizationSense = "maximize",
    key_reaction_ids: tuple[str, ...] = (),
    metadata: Mapping[str, Any] | None = None,
) -> LPProblem:
    """Create a backend-neutral LPProblem from a pcSec model and symbolic blocks."""

    reaction_ids = tuple(str(reaction_id) for reaction_id in getattr(model, "rxns"))
    reaction_index = {reaction_id: index for index, reaction_id in enumerate(reaction_ids)}
    if objective_reaction_id not in reaction_index:
        raise KeyError(f"Objective reaction not found: {objective_reaction_id}")

    bounds = tuple(
        (_normalize_bound(lower, lower=True), _normalize_bound(upper, lower=False))
        for lower, upper in zip(getattr(model, "lb"), getattr(model, "ub"))
    )
    if len(bounds) != len(reaction_ids):
        raise ValueError(f"Bounds length does not match reaction count: {len(bounds)} != {len(reaction_ids)}")
    return LPProblem(
        reaction_ids=reaction_ids,
        bounds=bounds,
        objective={objective_reaction_id: 1.0},
        stoichiometric_matrix=getattr(model, "s_matrix"),
        rhs=np.array(getattr(model, "b"), dtype=float).reshape(-1),
        constraint_blocks=tuple(constraint_blocks),
        objective_sense=objective_sense,
        key_reaction_ids=_unique_reaction_ids((*key_reaction_ids, objective_reaction_id)),
        metadata=dict(metadata or {}),
    )


def build_shadow_ladder_lp_problems(
    prep: ShadowTargetPreparation,
    blocks: tuple[ConstraintBlock, ...] | None = None,
    layer_ids: tuple[str, ...] = (
        "base_cobrapy_fba",
        "fixed_growth",
        "metabolic_coupling",
        "secretory_coupling",
        "protein_mass",
        "resource",
    ),
    key_reaction_ids: tuple[str, ...] = ("BIOMASS", "Ex_glc_D", "Ex_o2"),
) -> dict[str, LPProblem]:
    """Build pcSec shadow ladder LP problems without solving them."""

    resolved_blocks = blocks or build_shadow_constraint_blocks(prep)
    block_by_id = {block.layer_id: block for block in resolved_blocks}
    problems: dict[str, LPProblem] = {}
    for layer_id in layer_ids:
        try:
            selected_block_ids = LADDER_LAYER_BLOCK_IDS[layer_id]
        except KeyError as exc:
            raise KeyError(f"Unknown shadow LP ladder layer: {layer_id}") from exc
        model = prep.model if layer_id == "base_cobrapy_fba" else prep.fixed_model
        selected_blocks = tuple(block_by_id[block_id] for block_id in selected_block_ids)
        problems[layer_id] = lp_problem_from_model(
            model,
            prep.exchange_reaction_id,
            selected_blocks,
            key_reaction_ids=key_reaction_ids,
            metadata={
                "target_id": prep.target_id,
                "layer_id": layer_id,
                "exchange_reaction_id": prep.exchange_reaction_id,
                "selected_block_ids": selected_block_ids,
                "reference_layer_order": REFERENCE_LAYER_ORDER,
            },
        )
    return problems


def assemble_lp_problem(problem: LPProblem) -> AssembledLPProblem:
    """Assemble symbolic pcSec constraints into sparse equality and inequality matrices."""

    reaction_index = {reaction_id: index for index, reaction_id in enumerate(problem.reaction_ids)}
    objective_vector = _objective_vector(problem.objective, reaction_index, len(problem.reaction_ids))

    stoichiometric_matrix = problem.stoichiometric_matrix.tocsr()
    rhs = np.array(problem.rhs, dtype=float).reshape(-1)
    if stoichiometric_matrix.shape[1] != len(problem.reaction_ids):
        raise ValueError(
            "Stoichiometric matrix column count does not match reaction count: "
            f"{stoichiometric_matrix.shape[1]} != {len(problem.reaction_ids)}"
        )
    if stoichiometric_matrix.shape[0] != rhs.shape[0]:
        raise ValueError(
            "Stoichiometric matrix row count does not match RHS length: "
            f"{stoichiometric_matrix.shape[0]} != {rhs.shape[0]}"
        )

    eq_blocks: list[sparse.csr_matrix] = [stoichiometric_matrix]
    eq_rhs: list[float] = []
    ub_blocks: list[sparse.csr_matrix] = []
    ub_rhs: list[float] = []
    order: list[ConstraintOrderEntry] = []
    layer_counts: dict[str, int] = {}

    for block in problem.constraint_blocks:
        layer_counts[block.layer_id] = int(block.counts.get(block.layer_id, len(block.constraints)))
        for constraint in block.constraints:
            row = _constraint_row(constraint, reaction_index, len(problem.reaction_ids))
            if constraint.sense == "eq":
                row_index = stoichiometric_matrix.shape[0] + len(eq_rhs)
                eq_blocks.append(row)
                eq_rhs.append(float(constraint.rhs))
                order.append(_order_entry("eq", row_index, constraint))
            elif constraint.sense == "le":
                row_index = len(ub_rhs)
                ub_blocks.append(row)
                ub_rhs.append(float(constraint.rhs))
                order.append(_order_entry("ub", row_index, constraint))
            elif constraint.sense == "ge":
                row_index = len(ub_rhs)
                ub_blocks.append(-row)
                ub_rhs.append(-float(constraint.rhs))
                order.append(_order_entry("ub", row_index, constraint))
            else:
                raise ValueError(f"Unsupported constraint sense: {constraint.sense}")

    A_eq = sparse.vstack(eq_blocks, format="csr")
    b_eq = np.concatenate([rhs, np.array(eq_rhs, dtype=float)])
    if ub_blocks:
        A_ub = sparse.vstack(ub_blocks, format="csr")
        b_ub = np.array(ub_rhs, dtype=float)
    else:
        A_ub = None
        b_ub = None

    ub_count = 0 if A_ub is None else int(A_ub.shape[0])
    diagnostics = LPAssemblyDiagnostics(
        variable_count=len(problem.reaction_ids),
        stoichiometric_constraint_count=int(stoichiometric_matrix.shape[0]),
        eq_constraint_count=int(A_eq.shape[0]),
        ub_constraint_count=ub_count,
        constraint_count=int(A_eq.shape[0] + ub_count),
        objective_reaction_ids=tuple(problem.objective),
        layer_counts=layer_counts,
        block_order=tuple(block.layer_id for block in problem.constraint_blocks),
        constraint_order=tuple(order),
    )
    return AssembledLPProblem(
        reaction_ids=problem.reaction_ids,
        reaction_index=reaction_index,
        bounds=problem.bounds,
        objective_vector=objective_vector,
        objective_sense=problem.objective_sense,
        A_eq=A_eq,
        b_eq=b_eq,
        A_ub=A_ub,
        b_ub=b_ub,
        key_reaction_ids=problem.key_reaction_ids,
        diagnostics=diagnostics,
        metadata=problem.metadata,
    )


def _constraint_row(
    constraint: ConstraintSpec,
    reaction_index: Mapping[str, int],
    variable_count: int,
) -> sparse.csr_matrix:
    columns: list[int] = []
    values: list[float] = []
    for reaction_id, coefficient in constraint.terms.items():
        try:
            column = reaction_index[reaction_id]
        except KeyError as exc:
            raise KeyError(
                f"Constraint {constraint.name!r} references unknown reaction {reaction_id!r}."
            ) from exc
        value = float(coefficient)
        if not value:
            continue
        columns.append(column)
        values.append(value)
    return sparse.csr_matrix((values, ([0] * len(columns), columns)), shape=(1, variable_count))


def _objective_vector(
    objective: Mapping[str, float],
    reaction_index: Mapping[str, int],
    variable_count: int,
) -> np.ndarray:
    vector = np.zeros(variable_count, dtype=float)
    for reaction_id, coefficient in objective.items():
        try:
            vector[reaction_index[reaction_id]] = float(coefficient)
        except KeyError as exc:
            raise KeyError(f"Objective references unknown reaction {reaction_id!r}.") from exc
    return vector


def _order_entry(row_kind: str, row_index: int, constraint: ConstraintSpec) -> ConstraintOrderEntry:
    return ConstraintOrderEntry(
        row_kind=row_kind,
        row_index=int(row_index),
        layer_id=constraint.layer,
        constraint_name=constraint.name,
        original_sense=constraint.sense,
        term_count=len(constraint.terms),
        rhs=float(constraint.rhs),
    )


def _normalize_bound(value: Any, lower: bool) -> float | None:
    number = float(value)
    if np.isneginf(number) and lower:
        return None
    if np.isposinf(number) and not lower:
        return None
    return number


def _unique_reaction_ids(reaction_ids: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(reaction_id) for reaction_id in reaction_ids))
