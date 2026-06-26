from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _module_ast(relative_path: str) -> ast.Module:
    return ast.parse((REPO_ROOT / relative_path).read_text(encoding="utf-8"))


def _source(relative_path: str) -> str:
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


def test_real_pipeline_solve_tests_are_env_gated() -> None:
    module_ast = _module_ast("python_pichia/tests/test_pipeline_entrypoints.py")
    ungated: list[str] = []

    for node in module_ast.body:
        if not isinstance(node, ast.FunctionDef) or not node.name.startswith("test_"):
            continue
        calls_pipeline_solve = any(
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Name)
            and child.func.id == "run_pichia_secretion_simulation"
            for child in ast.walk(node)
        )
        has_slow_marker = any(
            isinstance(decorator, ast.Name) and decorator.id == "slow_pipeline"
            for decorator in node.decorator_list
        )
        if calls_pipeline_solve and not has_slow_marker:
            ungated.append(node.name)

    source = _source("python_pichia/tests/test_pipeline_entrypoints.py")
    assert 'os.environ.get("PCSEC_RUN_SLOW_PIPELINE_TESTS") != "1"' in source
    assert ungated == []


def test_probe_migration_regression_tests_are_env_gated() -> None:
    source = _source("python_pichia/tests/test_probe_migration.py")

    assert "pytestmark = pytest.mark.skipif" in source
    assert 'os.environ.get("PCSEC_RUN_SLOW_PROBE_TESTS") != "1"' in source
    assert "slow probe migration solve" in source


def test_screen_solve_tests_are_env_gated() -> None:
    module_ast = _module_ast("python_pichia/tests/test_screens_entrypoints.py")
    ungated: list[str] = []
    solve_calls = {"run_knockout_screen", "run_overexpression_screen", "run_reaction_knockout_screen"}

    for node in module_ast.body:
        if not isinstance(node, ast.FunctionDef) or not node.name.startswith("test_"):
            continue
        calls_screen_solve = any(
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Name)
            and child.func.id in solve_calls
            for child in ast.walk(node)
        )
        has_slow_marker = any(
            isinstance(decorator, ast.Name) and decorator.id == "slow_screen"
            for decorator in node.decorator_list
        )
        if calls_screen_solve and not has_slow_marker:
            ungated.append(node.name)

    source = _source("python_pichia/tests/test_screens_entrypoints.py")
    assert 'os.environ.get("PCSEC_RUN_SLOW_SCREEN_TESTS") != "1"' in source
    assert ungated == []


def test_release_validation_docs_keep_slow_gates_explicit() -> None:
    source = _source("docs/pichia_python_release_validation_2026-06-25.md")

    assert 'PCSEC_RUN_SLOW_PIPELINE_TESTS="1"' in source
    assert 'PCSEC_RUN_SLOW_SCREEN_TESTS="1"' in source
    assert 'PCSEC_RUN_SLOW_PROBE_TESTS="1"' in source
    assert "MATLAB harness / baseline 生成" in source
    assert "全模型 KO/OE 批量筛选" in source
