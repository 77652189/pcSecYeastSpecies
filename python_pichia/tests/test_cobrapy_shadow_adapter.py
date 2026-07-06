from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from pcsec_pichia.adapters.cobrapy_shadow import (
    UNAVAILABLE_MESSAGE,
    cobrapy_available,
    compare_shadow_fba,
    convert_to_cobrapy_model,
    solve_cobrapy_shadow_fba,
)
from pcsec_pichia.adapters.lp_solver import ScipyHiGHSSolver
from pcsec_pichia.core.pichia_model import PichiaModel


REPO_ROOT = Path(__file__).resolve().parents[2]


def _tiny_base_gem_model() -> PichiaModel:
    return PichiaModel(
        model_id="tiny_shadow_fba",
        source_file=Path("tiny_shadow_fba.mat"),
        rxns=["EX_A", "R_AB", "EX_B"],
        mets=["A_c", "B_c"],
        genes=[],
        lb=np.array([0.0, 0.0, 0.0], dtype=float),
        ub=np.array([10.0, 1000.0, 1000.0], dtype=float),
        c=np.zeros(3, dtype=float),
        b=np.zeros(2, dtype=float),
        s_matrix=sparse.csc_matrix(
            [
                [1.0, -1.0, 0.0],
                [0.0, 1.0, -1.0],
            ]
        ),
        rules=["", "", ""],
        gr_rules=["", "", ""],
    )


def test_cobrapy_shadow_unavailable_is_stable_without_importerror(monkeypatch: pytest.MonkeyPatch) -> None:
    original_find_spec = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None if name == "cobra" else original_find_spec(name))
    model = _tiny_base_gem_model()

    build = convert_to_cobrapy_model(model, objective_reaction="EX_B")
    solved = solve_cobrapy_shadow_fba(model, "EX_B", key_reactions=("EX_A", "EX_B"))

    assert build.available is False
    assert build.status == "unavailable"
    assert UNAVAILABLE_MESSAGE in build.message
    assert build.to_dict()["model_summary"]["supports_pcsec_constraints"] is False
    assert solved.available is False
    assert solved.success is False
    assert solved.status == "unavailable"
    assert "COBRApy is not installed" in solved.message
    assert solved.to_dict()["key_fluxes"] == []


def test_missing_objective_is_reported_without_requiring_cobrapy(monkeypatch: pytest.MonkeyPatch) -> None:
    original_find_spec = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None if name == "cobra" else original_find_spec(name))

    solved = solve_cobrapy_shadow_fba(_tiny_base_gem_model(), "NO_SUCH_REACTION")

    assert solved.status == "missing_objective"
    assert solved.available is False
    assert solved.success is False
    assert "Objective reaction not found" in solved.message


def test_no_default_paths_import_or_call_cobrapy_shadow_layer() -> None:
    default_path_files = [
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "pipeline.py",
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "simulation" / "__init__.py",
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "screens" / "__init__.py",
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "reports" / "__init__.py",
        REPO_ROOT / "app" / "services" / "pichia_request_mapping_service.py",
        REPO_ROOT / "app" / "ui" / "streamlit_app.py",
    ]

    for path in default_path_files:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported_modules = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        imported_from_modules = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        assert "pcsec_pichia.adapters.cobrapy_shadow" not in imported_modules
        assert "pcsec_pichia.adapters.cobrapy_shadow" not in imported_from_modules
        assert "cobrapy_shadow" not in path.read_text(encoding="utf-8")


def test_shadow_result_comparison_reports_unavailable_without_failing(monkeypatch: pytest.MonkeyPatch) -> None:
    original_find_spec = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None if name == "cobra" else original_find_spec(name))
    current_result = ScipyHiGHSSolver().solve(_tiny_base_gem_model(), "EX_B", key_reactions=("EX_A", "EX_B"))
    shadow_result = solve_cobrapy_shadow_fba(_tiny_base_gem_model(), "EX_B", key_reactions=("EX_A", "EX_B"))

    comparison = compare_shadow_fba(current_result, shadow_result, key_reactions=("EX_A", "EX_B"))

    assert comparison.comparable is False
    assert comparison.status == "unavailable"
    assert "COBRApy is not installed" in comparison.message


def test_cobrapy_shadow_tiny_model_parity_when_optional_dependency_is_installed() -> None:
    if not cobrapy_available():
        pytest.skip("COBRApy is not installed; optional shadow parity test skipped.")

    model = _tiny_base_gem_model()
    current_result = ScipyHiGHSSolver().solve(model, "EX_B", key_reactions=("EX_A", "EX_B"))
    shadow_result = solve_cobrapy_shadow_fba(model, "EX_B", key_reactions=("EX_A", "EX_B"))
    comparison = compare_shadow_fba(current_result, shadow_result, key_reactions=("EX_A", "EX_B"))

    assert shadow_result.available is True
    assert shadow_result.success is True
    assert shadow_result.objective_value == pytest.approx(current_result.objective_value)
    assert comparison.comparable is True
    assert comparison.within_tolerance is True
