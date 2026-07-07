from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import pytest
from scipy import sparse

from pcsec_pichia.analysis.shadow_lp import (
    ConstraintBlock,
    ConstraintSpec,
    LPProblem,
    ScipyHighsBackend,
    assemble_lp_problem,
    build_shadow_constraint_blocks,
    build_shadow_ladder_lp_problems,
    prepare_builtin_shadow_target,
)


LADDER_LAYERS = (
    "base_cobrapy_fba",
    "fixed_growth",
    "metabolic_coupling",
    "secretory_coupling",
    "protein_mass",
    "resource",
)

FINAL_CONSTRAINT_COUNTS = {
    "hLF": 23025,
    "OPN_ALPHA_FULL_PROJECT": 23016,
}

FINAL_OBJECTIVES = {
    "hLF": 0.0032850100270232106,
    "OPN_ALPHA_FULL_PROJECT": 0.006572021526431409,
}


def test_lp_assembler_uses_reaction_ids_and_records_constraint_order() -> None:
    block = ConstraintBlock(
        layer_id="toy_layer",
        constraints=(
            ConstraintSpec(
                name="toy_eq",
                layer="toy_layer",
                sense="eq",
                terms={"R2": 2.0, "R1": -1.0},
                rhs=3.0,
                source="unit test",
                enabled_by_default=True,
            ),
            ConstraintSpec(
                name="toy_le",
                layer="toy_layer",
                sense="le",
                terms={"R1": 1.5},
                rhs=4.0,
                source="unit test",
                enabled_by_default=True,
            ),
        ),
        counts={"toy_layer": 2},
    )
    problem = LPProblem(
        reaction_ids=("R1", "R2"),
        bounds=((0.0, 10.0), (None, None)),
        objective={"R2": 1.0},
        stoichiometric_matrix=sparse.csr_matrix([[1.0, 0.0]]),
        rhs=np.array([0.0]),
        constraint_blocks=(block,),
        key_reaction_ids=("R1", "R2"),
    )

    assembled = assemble_lp_problem(problem)

    assert assembled.reaction_index == {"R1": 0, "R2": 1}
    assert assembled.bounds == ((0.0, 10.0), (None, None))
    assert assembled.objective_vector.tolist() == [0.0, 1.0]
    assert assembled.A_eq.shape == (2, 2)
    assert assembled.b_eq.tolist() == [0.0, 3.0]
    assert assembled.A_eq.toarray().tolist() == [[1.0, 0.0], [-1.0, 2.0]]
    assert assembled.A_ub is not None
    assert assembled.A_ub.toarray().tolist() == [[1.5, 0.0]]
    assert assembled.b_ub is not None
    assert assembled.b_ub.tolist() == [4.0]
    assert assembled.diagnostics.constraint_count == 3
    assert [entry.constraint_name for entry in assembled.diagnostics.constraint_order] == ["toy_eq", "toy_le"]
    assert [entry.row_index for entry in assembled.diagnostics.constraint_order] == [1, 0]


@pytest.fixture(scope="module")
def solved_shadow_ladders() -> Mapping[str, Mapping[str, dict[str, Any]]]:
    backend = ScipyHighsBackend()
    assert backend.available()
    payload: dict[str, Mapping[str, dict[str, Any]]] = {}
    for target_id in ("hLF", "OPN_ALPHA_FULL_PROJECT"):
        prep = prepare_builtin_shadow_target(target_id)
        blocks = build_shadow_constraint_blocks(prep)
        problems = build_shadow_ladder_lp_problems(prep, blocks, layer_ids=LADDER_LAYERS)
        target_payload: dict[str, dict[str, Any]] = {}
        for layer_id, problem in problems.items():
            assembled = assemble_lp_problem(problem)
            result = backend.solve(problem, options={"time_limit": 600.0, "presolve": True})
            assert result.success, f"{target_id}/{layer_id}: {result.status} {result.message}"
            target_payload[layer_id] = {"assembled": assembled, "result": result}
        payload[target_id] = target_payload
    return payload


def test_scipy_highs_backend_solves_hlf_and_opn_shadow_ladder_layers(
    solved_shadow_ladders: Mapping[str, Mapping[str, dict[str, Any]]],
) -> None:
    for target_id, layers in solved_shadow_ladders.items():
        assert tuple(layers) == LADDER_LAYERS
        for layer_id, payload in layers.items():
            assembled = payload["assembled"]
            result = payload["result"]
            assert result.status == "0"
            assert result.objective is not None
            assert result.objective > 0.0
            assert result.key_fluxes["BIOMASS"] is not None
            assert result.backend_metadata["backend"] == "scipy-highs"
            assert result.backend_metadata["constraint_count"] == assembled.diagnostics.constraint_count
            assert result.timings["assemble_seconds"] >= 0.0
            assert result.timings["solve_seconds"] >= 0.0


def test_final_resource_layer_matrix_counts_and_objectives_match_reference_payload(
    solved_shadow_ladders: Mapping[str, Mapping[str, dict[str, Any]]],
) -> None:
    for target_id, expected_constraints in FINAL_CONSTRAINT_COUNTS.items():
        assembled = solved_shadow_ladders[target_id]["resource"]["assembled"]
        result = solved_shadow_ladders[target_id]["resource"]["result"]

        assert assembled.diagnostics.constraint_count == expected_constraints
        assert assembled.diagnostics.ub_constraint_count == 1
        assert result.backend_metadata["constraint_count"] == expected_constraints
        assert result.objective == pytest.approx(FINAL_OBJECTIVES[target_id], rel=1e-4)
