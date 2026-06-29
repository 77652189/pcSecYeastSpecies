from __future__ import annotations

import pytest

from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.probe import solve_maximize, prepare_glucose_model, load_pcsec_pichia_model


def test_baseline_growth_solve_with_prototype() -> None:
    """Solve base pcSecPichia model with prototype (glucose conditions)."""
    root = ProjectPaths.discover().repo_root
    model = load_pcsec_pichia_model(root)
    glucose = prepare_glucose_model(model, media_type=2)
    # Maximise glucose uptake (minimise Ex_glc_D) - the standard FBA objective
    result = solve_maximize(glucose, "Ex_glc_D", key_reactions=["BIOMASS", "Ex_glc_D"])

    assert result.success is True
    assert result.status == "0"
    assert result.objective == "Ex_glc_D"
    assert result.objective_value is not None
    # Glucose uptake should be negative (consumed)
    assert result.objective_value < 0


def test_lp_solver_rejects_unknown_objective_reaction() -> None:
    root = ProjectPaths.discover().repo_root
    model = load_pcsec_pichia_model(root)
    glucose = prepare_glucose_model(model, media_type=2)
    result = solve_maximize(glucose, "NonExistentReaction")

    assert result.success is False
    assert result.objective == "NonExistentReaction"
