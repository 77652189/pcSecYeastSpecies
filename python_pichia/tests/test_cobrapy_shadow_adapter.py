from __future__ import annotations

import ast
import importlib
import importlib.machinery
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import sparse

from pcsec_pichia.adapters.cobrapy_shadow import (
    CobraPyShadowFBAResult,
    CobraPyShadowFlux,
    UNAVAILABLE_MESSAGE,
    _cobra_id_map,
    _cobra_safe_id,
    _s_matrix_csc,
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


def test_cobrapy_shadow_unavailable_when_cobra_import_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    original_find_spec = importlib.util.find_spec
    original_import_module = importlib.import_module

    def fake_find_spec(name: str, *args: object, **kwargs: object) -> importlib.machinery.ModuleSpec | None:
        if name == "cobra":
            return importlib.machinery.ModuleSpec(name="cobra", loader=None)
        return original_find_spec(name, *args, **kwargs)

    def fake_import_module(name: str, package: str | None = None) -> object:
        if name == "cobra":
            raise ImportError("broken cobra dependency")
        return original_import_module(name, package)

    monkeypatch.setattr(importlib.util, "find_spec", fake_find_spec)
    monkeypatch.setattr(importlib, "import_module", fake_import_module)
    model = _tiny_base_gem_model()

    build = convert_to_cobrapy_model(model, objective_reaction="EX_B")
    solved = solve_cobrapy_shadow_fba(model, "EX_B", key_reactions=("EX_A", "EX_B"))

    assert build.available is False
    assert build.status == "unavailable"
    assert "COBRApy is not installed" in build.message
    assert solved.available is False
    assert solved.success is False
    assert solved.status == "unavailable"
    assert "COBRApy is not installed" in solved.message


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


def test_shadow_comparison_uses_key_fluxes_for_zero_flux_reactions() -> None:
    current_result = SimpleNamespace(
        objective_value=10.0,
        fluxes={"EX_B": 10.0},
        key_fluxes=[SimpleNamespace(reaction_id="R_ZERO", flux=0.0)],
    )
    shadow_result = CobraPyShadowFBAResult(
        available=True,
        success=True,
        status="optimal",
        message="solved",
        objective_reaction="EX_B",
        sense="maximize",
        objective_value=10.0,
        fluxes={"EX_B": 10.0},
        key_fluxes=(CobraPyShadowFlux(reaction_id="R_ZERO", flux=0.0),),
    )

    comparison = compare_shadow_fba(current_result, shadow_result, key_reactions=("R_ZERO",))

    assert comparison.comparable is True
    assert comparison.key_flux_diffs["R_ZERO"] == 0.0


def test_shadow_adapter_normalizes_sparse_array_to_csc_matrix() -> None:
    if not hasattr(sparse, "csc_array"):
        pytest.skip("SciPy sparse arrays are not available in this environment.")

    model = _tiny_base_gem_model()
    model = PichiaModel(
        model_id=model.model_id,
        source_file=model.source_file,
        rxns=model.rxns,
        mets=model.mets,
        genes=model.genes,
        lb=model.lb,
        ub=model.ub,
        c=model.c,
        b=model.b,
        s_matrix=sparse.csc_array(model.s_matrix),
        rules=model.rules,
        gr_rules=model.gr_rules,
    )

    matrix = _s_matrix_csc(model)

    assert isinstance(matrix, sparse.csc_matrix)
    assert matrix.getcol(0).nnz == 1


def test_cobrapy_safe_id_replaces_whitespace_and_deduplicates() -> None:
    used_ids: set[str] = set()

    first = _cobra_safe_id("metabolite A [c]", used_ids)
    second = _cobra_safe_id("metabolite A [c]", used_ids)

    assert first == "metabolite__A__[c]"
    assert second == "metabolite__A__[c]__2"


def test_cobrapy_id_map_preserves_original_keys_for_spaced_reactions() -> None:
    id_map = _cobra_id_map(("GPI-anchor assembly, step 2", "GPI-anchor__assembly,__step__2"))

    assert id_map["GPI-anchor assembly, step 2"] == "GPI-anchor__assembly,__step__2"
    assert id_map["GPI-anchor__assembly,__step__2"] == "GPI-anchor__assembly,__step__2__2"


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
