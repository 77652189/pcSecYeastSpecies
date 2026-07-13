from __future__ import annotations

import ast
from pathlib import Path

from app.ui.common import EXPERIMENT_FEEDBACK_PAGE


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_experiment_feedback_page_is_registered_in_navigation_and_entrypoint() -> None:
    common = (REPO_ROOT / "app" / "ui" / "common.py").read_text(encoding="utf-8")
    entrypoint = (REPO_ROOT / "app" / "ui" / "streamlit_app.py").read_text(encoding="utf-8")

    assert EXPERIMENT_FEEDBACK_PAGE in common
    assert "render_experiment_feedback" in entrypoint
    assert "elif page == EXPERIMENT_FEEDBACK_PAGE" in entrypoint


def test_experiment_feedback_view_uses_service_only_and_stable_session_keys() -> None:
    view_path = REPO_ROOT / "app" / "ui" / "views" / "experiment_feedback.py"
    source = view_path.read_text(encoding="utf-8")
    module_ast = ast.parse(source)
    imported_modules: set[str] = set()
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)

    assert "app.services.pichia_experiment_feedback_service" in imported_modules
    assert not any(module.startswith(("pcsec_pichia", "python_pichia")) for module in imported_modules)
    assert "experiment_feedback_selected_run" in source
    assert "experiment_feedback_selected_run_pending" in source
    assert "experiment_feedback_last_import" in source
    assert "experiment_feedback_experiment_upload" in source
    assert "experiment_feedback_prediction_upload" in source
    assert 'type=["csv", "xlsx", "jsonl"]' in source
    assert "st.cache_data" not in source
    assert "st.cache_resource" not in source
    assert "submit_experiment_feedback_import" in source
    assert "validate_experiment_bundle" not in source
    assert "build_prediction_index" not in source
    assert "build_calibration_summary" not in source


def test_experiment_feedback_view_exposes_all_required_states_and_exports() -> None:
    source = (
        REPO_ROOT / "app" / "ui" / "views" / "experiment_feedback.py"
    ).read_text(encoding="utf-8")

    for text in (
        "Validation / Conflicts",
        "Linkage",
        "Calibration",
        "ambiguous",
        "missing_prediction",
        "context_mismatch",
        "不可校准",
        "export_experiment_feedback_issues",
        "export_experiment_feedback_report",
        "下载 prediction-vs-experiment 报告",
        "ranking_assessment",
        "comparable_rank_pair_count",
        "排序证据不足",
        "仅作描述性展示",
    ):
        assert text in source
    assert ".read_bytes()" not in source
