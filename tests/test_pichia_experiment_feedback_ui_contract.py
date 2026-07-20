from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace

from app.ui.common import EXPERIMENT_FEEDBACK_PAGE
from app.ui.views import experiment_feedback as experiment_feedback_view


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
    assert "experiment_feedback_target_metadata" in source
    assert "experiment_feedback_batch_metadata" in source
    assert "experiment_feedback_import_form_state" in source
    assert "_restore_import_form_state()" in source
    assert "on_change=_sync_import_form_field" in source
    assert "experiment_metadata=" in source
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
        "数据校验",
        "预测匹配",
        "历史数据核对",
        "有歧义",
        "无对应预测",
        "条件不匹配",
        "不可核对",
        "export_experiment_feedback_issues",
        "export_experiment_feedback_report",
        "下载预测 vs 实验核对报告",
        "ranking_assessment",
        "comparable_rank_pair_count",
        "样本量不足",
        "仅作描述性展示",
    ):
        assert text in source
    assert ".read_bytes()" not in source


def test_experiment_feedback_import_form_state_survives_widget_cleanup(monkeypatch) -> None:
    session_state: dict[str, object] = {}
    monkeypatch.setattr(
        experiment_feedback_view,
        "st",
        SimpleNamespace(session_state=session_state),
    )

    experiment_feedback_view._restore_import_form_state()
    session_state[experiment_feedback_view.TARGET_METADATA_KEY] = "hLF"
    session_state[experiment_feedback_view.BATCH_METADATA_KEY] = "B01"
    experiment_feedback_view._sync_import_form_field(
        experiment_feedback_view.TARGET_METADATA_KEY
    )
    experiment_feedback_view._sync_import_form_field(
        experiment_feedback_view.BATCH_METADATA_KEY
    )

    del session_state[experiment_feedback_view.TARGET_METADATA_KEY]
    del session_state[experiment_feedback_view.BATCH_METADATA_KEY]
    experiment_feedback_view._restore_import_form_state()

    assert session_state[experiment_feedback_view.TARGET_METADATA_KEY] == "hLF"
    assert session_state[experiment_feedback_view.BATCH_METADATA_KEY] == "B01"
