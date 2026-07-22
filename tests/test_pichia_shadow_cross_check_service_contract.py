from __future__ import annotations

import inspect
import json
import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from app.services import pichia_shadow_cross_check_service as service

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_shadow_cross_check_service_calls_engine_and_returns_manifest_paths(monkeypatch, tmp_path) -> None:
    captured = {}
    monkeypatch.setattr(service, "SHADOW_CROSS_CHECK_RUNS_DIR", tmp_path)

    def fake_run_shadow_lp_cross_check(request, output_dir):
        captured["request"] = request
        captured["output_dir"] = output_dir
        result = _FakeResult()
        return _FakeOutputs(
            manifest_path=output_dir / "cross_check_manifest.json",
            summary_tsv_path=output_dir / "cross_check_summary.tsv",
            report_path=output_dir / "cross_check_report.md",
            diff_path=output_dir / "reference_vs_shadow_diff.json",
            result=result,
        )

    monkeypatch.setattr(service, "run_shadow_lp_cross_check", fake_run_shadow_lp_cross_check)

    response = service.run_pichia_shadow_cross_check(
        target_id="hLF",
        screen_run_id="screen-1",
        saved_result_path=tmp_path / "saved.json",
        output_dir=tmp_path / "cross_check",
    )

    assert captured["request"].target_id == "hLF"
    assert captured["request"].screen_run_id == "screen-1"
    assert captured["request"].saved_result_path == str(tmp_path / "saved.json")
    assert captured["output_dir"] == tmp_path / "cross_check"
    assert response["submitted"] is True
    assert response["status"] == "ok"
    assert response["manifest_path"].endswith("cross_check_manifest.json")
    assert response["within_tolerance"] is True


def test_shadow_cross_check_service_default_output_dir_is_local_runs(monkeypatch) -> None:
    monkeypatch.setattr(service, "SHADOW_CROSS_CHECK_RUNS_DIR", REPO_ROOT / "local_runs" / "shadow_lp_cross_check")
    def fake_run_shadow_lp_cross_check(_request, output_dir):
        return _FakeOutputs(
            manifest_path=output_dir / "cross_check_manifest.json",
            summary_tsv_path=output_dir / "cross_check_summary.tsv",
            report_path=output_dir / "cross_check_report.md",
            diff_path=output_dir / "reference_vs_shadow_diff.json",
            result=_FakeResult(),
        )

    monkeypatch.setattr(service, "run_shadow_lp_cross_check", fake_run_shadow_lp_cross_check)

    response = service.run_pichia_shadow_cross_check(target_id="OPN_ALPHA_FULL_PROJECT")

    assert response["output_dir"].startswith(str(service.SHADOW_CROSS_CHECK_RUNS_DIR))
    assert "OPN_ALPHA_FULL_PROJECT" in response["output_dir"]


def test_shadow_cross_check_service_rejects_output_outside_local_runs(monkeypatch, tmp_path) -> None:
    allowed = tmp_path / "local_runs" / "shadow_lp_cross_check"
    monkeypatch.setattr(service, "SHADOW_CROSS_CHECK_RUNS_DIR", allowed)

    try:
        service.run_pichia_shadow_cross_check(target_id="hLF", output_dir=tmp_path / "Results")
    except ValueError as exc:
        assert "local_runs/shadow_lp_cross_check" in str(exc)
    else:
        raise AssertionError("output outside the shadow cross-check run directory must be rejected")


def test_shadow_cross_check_service_loads_manifest_from_file_or_directory(tmp_path) -> None:
    manifest_dir = tmp_path / "run"
    manifest_dir.mkdir()
    manifest_path = manifest_dir / "cross_check_manifest.json"
    manifest_path.write_text(json.dumps({"result": {"target_id": "hLF"}}), encoding="utf-8")

    assert service.load_pichia_shadow_cross_check_manifest(manifest_dir)["result"]["target_id"] == "hLF"
    assert service.load_pichia_shadow_cross_check_manifest(manifest_path)["result"]["target_id"] == "hLF"


def test_shadow_cross_check_service_stays_facade_only() -> None:
    source = inspect.getsource(service)

    assert "run_shadow_lp_cross_check" in source
    assert "run_shadow_ladder" not in source
    assert "validate_shadow_ladder_against_reference" not in source
    assert "solve_secretion_capacity" not in source


def test_shadow_cross_check_page_is_registered_in_navigation_and_entrypoint() -> None:
    common_text = (REPO_ROOT / "app" / "ui" / "common.py").read_text(encoding="utf-8")
    entrypoint_text = (REPO_ROOT / "app" / "ui" / "streamlit_app.py").read_text(encoding="utf-8")

    assert "SHADOW_CROSS_CHECK_PAGE" in common_text
    assert "Shadow LP一致性验证" in common_text
    assert "render_shadow_cross_check" in entrypoint_text
    assert "elif page == SHADOW_CROSS_CHECK_PAGE" in entrypoint_text


def test_shadow_cross_check_view_uses_service_facade_not_engine_runtime() -> None:
    view_path = REPO_ROOT / "app" / "ui" / "views" / "shadow_cross_check.py"
    module_ast = ast.parse(view_path.read_text(encoding="utf-8"))
    imported_modules: set[str] = set()
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)

    assert "app.services.pichia_shadow_cross_check_service" in imported_modules
    assert not any(module.startswith(("pcsec_pichia", "python_pichia")) for module in imported_modules)


def test_shadow_cross_check_view_initializes_session_and_triggers_service(monkeypatch, tmp_path) -> None:
    import app.ui.views.shadow_cross_check as view

    calls: list[dict[str, Any]] = []
    fake_st = _FakeStreamlit(
        text_values={
            "shadow_lp_cross_check_screen_run_id": "screen-2",
            "shadow_lp_cross_check_saved_result_path": str(tmp_path / "saved.json"),
            "shadow_lp_cross_check_output_dir": "",
        },
        clicked_keys={"shadow_lp_cross_check_run"},
    )

    def fake_run(**kwargs):
        calls.append(kwargs)
        return {
            "status": "ok",
            "target_id": kwargs["target_id"],
            "screen_run_id": kwargs["screen_run_id"],
            "report_path": str(tmp_path / "out" / "cross_check_report.md"),
            "within_tolerance": True,
            "relative_diff": 0.0,
            "warnings": [],
        }

    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(view, "run_pichia_shadow_cross_check", fake_run)

    view.render_shadow_cross_check()

    assert calls == [
        {
            "target_id": "hLF",
            "screen_run_id": "screen-2",
            "saved_result_path": tmp_path / "saved.json",
            "output_dir": None,
        }
    ]
    assert view.SHADOW_CROSS_CHECK_STATE_KEY in fake_st.session_state
    assert view.SHADOW_CROSS_CHECK_MANIFEST_KEY in fake_st.session_state
    assert fake_st.session_state[view.SHADOW_CROSS_CHECK_STATE_KEY]["status"] == "ok"
    assert fake_st.success_calls == ["一致（shadow LP 与参考结果对齐）"]
    assert any("cross_check_report.md" in item for item in fake_st.code_calls)


@dataclass(frozen=True)
class _FakeResult:
    manifest_status: str = "ok"
    target_id: str = "hLF"
    screen_run_id: str = "screen-1"
    within_tolerance: bool = True
    relative_diff: float = 0.0
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class _FakeOutputs:
    manifest_path: Path
    summary_tsv_path: Path
    report_path: Path
    diff_path: Path
    result: _FakeResult


class _FakeStreamlit:
    def __init__(
        self,
        *,
        text_values: dict[str, str] | None = None,
        clicked_keys: set[str] | None = None,
    ) -> None:
        self.session_state: dict[str, Any] = {}
        self.text_values = text_values or {}
        self.clicked_keys = clicked_keys or set()
        self.success_calls: list[str] = []
        self.warning_calls: list[str] = []
        self.error_calls: list[str] = []
        self.code_calls: list[str] = []

    def header(self, _text: str) -> None:
        return None

    def subheader(self, _text: str) -> None:
        return None

    def tabs(self, labels: list[str]) -> list["_FakeStreamlit"]:
        return [self for _label in labels]

    def __enter__(self) -> "_FakeStreamlit":
        return self

    def __exit__(self, *_args) -> None:
        return None

    def selectbox(self, _label: str, options: list[str], *, key: str) -> str:
        return options[0]

    def text_input(self, _label: str, *, value: str = "", key: str) -> str:
        return self.text_values.get(key, value)

    def button(self, _label: str, *, key: str) -> bool:
        return key in self.clicked_keys

    def success(self, text: str) -> None:
        self.success_calls.append(text)

    def warning(self, text: str) -> None:
        self.warning_calls.append(text)

    def error(self, text: str) -> None:
        self.error_calls.append(text)

    def columns(self, count: int) -> list["_FakeStreamlit"]:
        return [self for _ in range(count)]

    def metric(self, _label: str, _value: str) -> None:
        return None

    def code(self, text: str) -> None:
        self.code_calls.append(text)

    def json(self, _value: object) -> None:
        return None

    def markdown(self, _text: str) -> None:
        return None
