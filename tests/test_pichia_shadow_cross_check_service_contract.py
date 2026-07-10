from __future__ import annotations

import inspect
import json
from dataclasses import dataclass
from pathlib import Path

from app.services import pichia_shadow_cross_check_service as service


def test_shadow_cross_check_service_calls_engine_and_returns_manifest_paths(monkeypatch, tmp_path) -> None:
    captured = {}

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
