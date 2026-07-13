from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from pcsec_pichia.analysis.shadow_lp import (
    ShadowLpCrossCheckRequest,
    load_shadow_cross_check_saved_result,
    run_shadow_lp_cross_check,
)


def test_shadow_lp_cross_check_writes_manifest_summary_report_and_diff(tmp_path) -> None:
    saved_path = tmp_path / "saved_result.json"
    saved_path.write_text(
        json.dumps(
            {
                "target_id": "hLF",
                "screen_run_id": "screen-001",
                "reference_capacity": 1.00000001,
                "warnings": ["saved_result_warning"],
            }
        ),
        encoding="utf-8",
    )
    request = ShadowLpCrossCheckRequest(target_id="hLF", saved_result_path=str(saved_path))

    outputs = run_shadow_lp_cross_check(
        request,
        tmp_path / "cross_check",
        ladder_runner=_fake_ladder_runner,
        reference_validator=_aligned_validator,
    )
    result = outputs.result
    manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    report = outputs.report_path.read_text(encoding="utf-8")
    diff = json.loads(outputs.diff_path.read_text(encoding="utf-8"))

    assert result.target_id == "hLF"
    assert result.screen_run_id == "screen-001"
    assert result.reference_capacity == 1.00000001
    assert result.reference_source == "saved_result"
    assert result.saved_reference_capacity == 1.00000001
    assert result.shadow_capacity == 1.00000001
    assert result.within_tolerance is True
    assert result.constraint_layer == "mitochondrial"
    assert result.backend == "scipy-highs"
    assert result.solver_status == "0"
    assert "saved_result_warning" in result.warnings
    assert manifest["result"]["within_tolerance"] is True
    assert diff["validation"]["final_alignment_status"] == "aligned"
    assert "Shadow LP Cross-check Report" in report
    assert "mg/L" not in report
    assert "experimental success rate" in report
    assert outputs.summary_tsv_path.exists()


def test_shadow_lp_cross_check_uses_saved_capacity_for_alignment(tmp_path) -> None:
    saved_path = tmp_path / "saved_result.json"
    saved_path.write_text(
        json.dumps({"target_id": "hLF", "reference_capacity": 0.5}),
        encoding="utf-8",
    )

    outputs = run_shadow_lp_cross_check(
        ShadowLpCrossCheckRequest(target_id="hLF", saved_result_path=str(saved_path)),
        tmp_path / "cross_check",
        ladder_runner=_fake_ladder_runner,
        reference_validator=_aligned_validator,
    )

    assert outputs.result.reference_capacity == 0.5
    assert outputs.result.reference_source == "saved_result"
    assert outputs.result.within_tolerance is False
    assert outputs.result.manifest_status == "review_required"
    assert "saved_reference_capacity_used" in outputs.result.warnings


def test_shadow_lp_cross_check_marks_review_required_without_experimental_claim(tmp_path) -> None:
    request = ShadowLpCrossCheckRequest(target_id="OPN_ALPHA_FULL_PROJECT", relative_tolerance=1e-4)

    outputs = run_shadow_lp_cross_check(
        request,
        tmp_path / "cross_check",
        ladder_runner=_fake_ladder_runner,
        reference_validator=_review_required_validator,
    )
    report = outputs.report_path.read_text(encoding="utf-8")

    assert outputs.result.within_tolerance is False
    assert outputs.result.alignment_status == "review_required"
    assert outputs.result.manifest_status == "review_required"
    assert "shadow_cross_check_review_required" in outputs.result.warnings
    assert "not an experimental infeasibility call" in report


def test_saved_shadow_cross_check_context_loads_minimal_saved_result(tmp_path) -> None:
    saved_path = tmp_path / "saved_result.json"
    saved_path.write_text('{"target_id": "hLF", "run_id": "run-7", "secretion_capacity": 2.5}', encoding="utf-8")

    context = load_shadow_cross_check_saved_result(saved_path)

    assert context.target_id == "hLF"
    assert context.screen_run_id == "run-7"
    assert context.reference_capacity == 2.5
    assert context.source_path == str(saved_path)


def _fake_ladder_runner(target_id: str, **_kwargs: Any) -> _FakeLadder:
    return _FakeLadder(target_id=target_id)


def _aligned_validator(_ladder: _FakeLadder, **_kwargs: Any) -> dict[str, object]:
    return {
        "reference_objective": 1.0,
        "shadow_objective": 1.00000001,
        "objective_abs_diff": 1e-8,
        "objective_rel_diff": 1e-8,
        "constraint_count_diff": 0,
        "final_alignment_status": "aligned",
        "reference_status": "0",
    }


def _review_required_validator(_ladder: _FakeLadder, **_kwargs: Any) -> dict[str, object]:
    return {
        "reference_objective": 1.0,
        "shadow_objective": 1.2,
        "objective_abs_diff": 0.2,
        "objective_rel_diff": 0.2,
        "constraint_count_diff": 2,
        "final_alignment_status": "review_required",
        "reference_status": "0",
    }


@dataclass(frozen=True)
class _FakeLayer:
    layer_id: str = "mitochondrial"
    status: str = "0"
    objective: float = 1.00000001

    def to_dict(self) -> dict[str, object]:
        return {"layer_id": self.layer_id, "status": self.status, "objective": self.objective}


@dataclass(frozen=True)
class _FakeLadder:
    target_id: str
    backend_name: str = "scipy-highs"
    warnings: tuple[str, ...] = ()

    @property
    def final_layer(self) -> _FakeLayer:
        return _FakeLayer()
