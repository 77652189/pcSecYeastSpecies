from __future__ import annotations

import ast
import importlib
import importlib.machinery
import importlib.util
import json
import shutil
from pathlib import Path

import pytest

from pcsec_pichia.analysis.cobrapy_shadow_baseline import (
    CobraPyShadowBaselineCase,
    run_cobrapy_shadow_baseline,
)
from pcsec_pichia.analysis import cobrapy_shadow_baseline
from pcsec_pichia.core.paths import ProjectPaths


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_real_model_shadow_baseline_unavailable_writes_only_local_runs(monkeypatch: pytest.MonkeyPatch) -> None:
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

    paths = ProjectPaths.discover(REPO_ROOT)
    output_dir = paths.local_runs_dir / "cobrapy_shadow_baseline_pytest"
    shutil.rmtree(output_dir, ignore_errors=True)
    try:
        run = run_cobrapy_shadow_baseline(
            paths=paths,
            objective_reactions=("Ex_glc_D",),
            output_dir=output_dir,
        )

        assert run.result_status == "completed_shadow_unavailable"
        assert run.summary_json_path is not None
        assert run.report_markdown_path is not None
        assert Path(run.summary_json_path).is_file()
        assert Path(run.report_markdown_path).is_file()
        assert Path(run.output_dir).resolve().is_relative_to(paths.local_runs_dir.resolve())
        assert run.cases
        assert all(case.shadow_available is False for case in run.cases)
        assert all(case.shadow_success is False for case in run.cases)
        assert all(case.current_success in {True, False} for case in run.cases)

        payload = json.loads(Path(run.summary_json_path).read_text(encoding="utf-8"))
        assert payload["cases"][0]["case_id"] == "Ex_glc_D_maximize"
        assert payload["cases"][0]["shadow_available"] is False
        assert "pcSec protein/secretion constraints" in "\n".join(payload["warnings"])
        for artifact_path in (Path(run.summary_json_path), Path(run.report_markdown_path)):
            assert artifact_path.resolve().is_relative_to(paths.local_runs_dir.resolve())
            assert not any(part in {"Data", "Model", "Enzymedata", "Results"} for part in artifact_path.parts)
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)


def test_shadow_baseline_rejects_output_outside_local_runs() -> None:
    paths = ProjectPaths.discover(REPO_ROOT)

    with pytest.raises(ValueError, match="local_runs"):
        run_cobrapy_shadow_baseline(
            paths=paths,
            objective_reactions=("Ex_glc_D",),
            output_dir=REPO_ROOT / "Results" / "cobrapy_shadow_baseline",
            write_artifacts=False,
        )


def test_shadow_baseline_status_does_not_pass_when_shadow_solves_fail() -> None:
    case = CobraPyShadowBaselineCase(
        case_id="BIOMASS_maximize",
        objective_reaction="BIOMASS",
        sense="maximize",
        shadow_available=True,
        current_success=True,
        shadow_success=False,
        current_objective_value=1.0,
        shadow_objective_value=None,
        objective_abs_diff=None,
        objective_rel_diff=None,
        within_tolerance=False,
        key_flux_diffs={},
        model_summary={},
        shadow_status="failed",
    )

    assert cobrapy_shadow_baseline._result_status([case]) == "completed_with_differences_or_failures"


def test_default_paths_do_not_import_cobrapy_shadow_baseline() -> None:
    default_path_files = [
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "pipeline.py",
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "simulation" / "__init__.py",
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "screens" / "__init__.py",
        REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "reports" / "__init__.py",
        REPO_ROOT / "app" / "services" / "pichia_request_mapping_service.py",
        REPO_ROOT / "app" / "ui" / "streamlit_app.py",
    ]

    for path in default_path_files:
        text = path.read_text(encoding="utf-8")
        tree = ast.parse(text)
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
        assert "pcsec_pichia.analysis.cobrapy_shadow_baseline" not in imported_modules
        assert "pcsec_pichia.analysis.cobrapy_shadow_baseline" not in imported_from_modules
        assert "cobrapy_shadow_baseline" not in text
