from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from pcsec_pichia.analysis.shadow_lp import CobraOptlangBackend, ConstraintBlock, ConstraintSpec, LPProblem, ScipyHighsBackend
from pcsec_pichia.analysis.shadow_lp.backends import _unique_cobra_ids


def _tiny_problem() -> LPProblem:
    limit_block = ConstraintBlock(
        layer_id="toy_limit",
        constraints=(
            ConstraintSpec(
                name="target_limit",
                layer="toy_limit",
                sense="le",
                terms={"TARGET": 1.0},
                rhs=4.0,
                source="tiny semantic test",
                enabled_by_default=True,
            ),
        ),
        counts={"toy_limit": 1},
    )
    return LPProblem(
        reaction_ids=("SOURCE", "TARGET"),
        bounds=((0.0, 10.0), (0.0, 10.0)),
        objective={"TARGET": 1.0},
        stoichiometric_matrix=sparse.csr_matrix([[1.0, -1.0]]),
        rhs=np.array([0.0]),
        constraint_blocks=(limit_block,),
        key_reaction_ids=("SOURCE", "TARGET"),
    )


def test_cobra_optlang_backend_skips_when_optional_dependencies_are_unavailable() -> None:
    backend = CobraOptlangBackend()
    if backend.available():
        pytest.skip("COBRApy/optlang is available; unavailable skip path is not active in this environment.")

    result = backend.solve(_tiny_problem())

    assert result.success is False
    assert result.status == "unavailable"
    assert result.backend_metadata["available"] is False


def test_cobra_optlang_backend_generates_unique_sanitized_reaction_ids() -> None:
    cobra_ids = _unique_cobra_ids(("TARGET-A", "TARGET_A", "TARGET A"))

    assert len(set(cobra_ids.values())) == 3
    assert all(cobra_id.startswith("R_TARGET") for cobra_id in cobra_ids.values())


def test_cobra_optlang_backend_matches_scipy_highs_on_tiny_lp() -> None:
    cobra_backend = CobraOptlangBackend()
    if not cobra_backend.available():
        pytest.skip("COBRApy/optlang optional dependency is unavailable.")
    problem = _tiny_problem()

    scipy_result = ScipyHighsBackend().solve(problem)
    cobra_result = cobra_backend.solve(problem)

    assert scipy_result.success is True
    assert cobra_result.success is True
    assert cobra_result.objective == pytest.approx(scipy_result.objective)
    assert cobra_result.backend_metadata["backend"] == "cobra-optlang"
    assert cobra_result.backend_metadata["limited_to_tiny_lp"] is True
