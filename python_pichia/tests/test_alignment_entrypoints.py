from __future__ import annotations

from pathlib import Path

import pytest

from pcsec_pichia.alignment import (
    ALIGNMENT_STATUSES,
    AlignmentSummary,
    build_alignment_summary,
    classify_alignment_status,
    hlf_project_710_known_matlab_compatibility_exceptions,
    load_matlab_alignment_artifact,
    opn_known_matlab_compatibility_exceptions,
    summarize_alignment,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_PATH = REPO_ROOT / "local_runs" / "pichia_hlf_opn_probe" / "matlab_stage3_alignment" / "matlab_stage3_alignment_summary.json"
HLF_PROJECT_710_ARTIFACT_PATH = (
    REPO_ROOT
    / "local_runs"
    / "pichia_hlf_opn_probe"
    / "hlf_project_sequence_matlab_harness_2026-06-26"
    / "hlf_project_sequence_matlab_harness_summary.json"
)


def test_opn_alignment_artifact_is_loaded_but_not_marked_aligned() -> None:
    artifact = load_matlab_alignment_artifact(ARTIFACT_PATH)

    summary = build_alignment_summary(
        "OPN_ALPHA_FULL_PROJECT",
        python_result_status="draft",
        artifact_path=ARTIFACT_PATH,
        artifact=artifact,
        rows_python=24434,
        cols_python=29057,
        objective_python=0.0021130308979533004,
    )
    payload = summarize_alignment(summary)

    assert isinstance(summary, AlignmentSummary)
    assert summary.baseline_available is True
    assert summary.matlab_success is True
    assert summary.rows_matlab == 24435
    assert summary.cols_matlab == 29057
    assert summary.rows_diff == -1
    assert summary.cols_diff == 0
    assert summary.objective_relative_diff == pytest.approx(0.0)
    assert summary.constraint_diff_status == "shape_diff"
    assert summary.matlab_alignment_status == "python_draft"
    assert summary.success is False
    assert payload["rows_diff"] == summary.rows_diff
    assert payload["objective_relative_diff"] == pytest.approx(0.0)


def test_opn_alignment_shape_can_match_when_python_ub_row_is_counted() -> None:
    summary = build_alignment_summary(
        "OPN_ALPHA_FULL_PROJECT",
        python_result_status="draft",
        artifact_path=ARTIFACT_PATH,
        rows_python=24435,
        cols_python=29057,
        objective_python=0.0021130308979533004,
    )

    assert summary.rows_diff == 0
    assert summary.cols_diff == 0
    assert summary.constraint_diff_status == "matched"
    assert summary.matlab_alignment_status == "python_draft"
    assert summary.success is False


def test_hlf_alignment_preserves_matlab_failure_diagnostic() -> None:
    summary = build_alignment_summary(
        "hLF",
        python_result_status="draft",
        artifact_path=ARTIFACT_PATH,
    )

    assert summary.baseline_available is True
    assert summary.matlab_success is False
    assert summary.matlab_alignment_status == "matlab_failed"
    assert "MATLAB:sizeDimensionsMustMatch" in summary.diagnostic_message
    assert "calculateMW.m:35" in summary.diagnostic_message
    assert summary.success is False


def test_missing_artifact_returns_baseline_missing_without_unclear_exception(tmp_path: Path) -> None:
    missing_path = tmp_path / "missing_alignment.json"

    artifact = load_matlab_alignment_artifact(missing_path)
    summary = build_alignment_summary(
        "OPN_ALPHA_FULL_PROJECT",
        artifact_path=missing_path,
        artifact=artifact,
    )

    assert artifact["baseline_available"] is False
    assert summary.baseline_available is False
    assert summary.matlab_alignment_status == "baseline_missing"
    assert "missing" in summary.diagnostic_message.lower()


def test_corrected_condition_is_not_treated_as_old_matlab_baseline_equivalent() -> None:
    summary = build_alignment_summary(
        "OPN_ALPHA_FULL_PROJECT",
        python_result_status="corrected_condition",
        artifact_path=ARTIFACT_PATH,
        rows_python=25136,
        cols_python=29057,
        objective_python=0.0021130308979533004,
        constraint_diff_status="matched",
    )

    assert summary.matlab_alignment_status == "pending"
    assert summary.success is False


def test_opn_known_compatibility_exceptions_are_explicit_and_not_fully_aligned() -> None:
    summary = build_alignment_summary(
        "OPN_ALPHA_FULL_PROJECT",
        python_result_status="corrected_condition",
        artifact_path=ARTIFACT_PATH,
        rows_python=24435,
        cols_python=29057,
        objective_python=0.0021130308979533004,
        constraint_diff_status="known_matlab_compatibility_differences",
        compatibility_exceptions=opn_known_matlab_compatibility_exceptions(),
    )
    payload = summarize_alignment(summary)

    assert "aligned_except_known_matlab_compatibility_differences" in ALIGNMENT_STATUSES
    assert summary.python_result_status == "corrected_condition"
    assert summary.matlab_alignment_status == "aligned_except_known_matlab_compatibility_differences"
    assert summary.success is False
    assert payload["is_fully_aligned"] is False
    assert payload["is_aligned_except_known_matlab_compatibility_differences"] is True
    assert len(payload["compatibility_exceptions"]) == 4
    assert {item["count"] for item in payload["compatibility_exceptions"]} == {9, 1418, 18, 2}


def test_opn_corrected_alignment_regression_records_named_known_exceptions() -> None:
    corrected_opn_objective = 0.0021305196599992996
    artifact = load_matlab_alignment_artifact(ARTIFACT_PATH)

    summary = build_alignment_summary(
        "OPN_ALPHA_FULL_PROJECT",
        python_result_status="corrected_condition",
        artifact_path=ARTIFACT_PATH,
        artifact=artifact,
        rows_python=24435,
        cols_python=29057,
        objective_python=corrected_opn_objective,
        constraint_diff_status="known_matlab_compatibility_differences",
        compatibility_exceptions=opn_known_matlab_compatibility_exceptions(),
    )
    payload = summarize_alignment(summary)
    exceptions = {item["id"]: item for item in payload["compatibility_exceptions"]}

    assert summary.baseline_available is True
    assert summary.matlab_success is True
    assert summary.rows_diff == 0
    assert summary.cols_diff == 0
    assert summary.constraint_diff_status == "known_matlab_compatibility_differences"
    assert summary.objective_relative_diff == pytest.approx(0.008276623906890783)
    assert summary.objective_relative_diff < 0.01
    assert summary.matlab_alignment_status == "aligned_except_known_matlab_compatibility_differences"
    assert summary.success is False
    assert payload["is_fully_aligned"] is False
    assert payload["is_aligned_except_known_matlab_compatibility_differences"] is True
    assert exceptions["corrected_medium_exchange_bounds"]["count"] == 9
    assert exceptions["corrected_medium_exchange_bounds"]["category"] == "bound_difference"
    assert exceptions["misfolding_dilution_bounds"]["count"] == 1418
    assert exceptions["misfolding_dilution_bounds"]["category"] == "bound_difference"
    assert exceptions["target_secretory_coupling_missing_in_matlab_artifact"]["count"] == 18
    assert exceptions["target_secretory_coupling_missing_in_matlab_artifact"]["category"] == "row_coefficient_difference"
    assert exceptions["ribosome_optional_row_order_and_proteins_term"]["count"] == 2
    assert exceptions["ribosome_optional_row_order_and_proteins_term"]["category"] == "row_coefficient_difference"


def test_classifier_only_marks_aligned_when_objective_and_constraints_match() -> None:
    assert (
        classify_alignment_status(
            python_result_status="matlab_aligned",
            baseline_available=True,
            matlab_success=True,
            objective_relative_diff=0.005,
            constraint_diff_status="matched",
        )
        == "aligned"
    )
    assert (
        classify_alignment_status(
            python_result_status="matlab_aligned",
            baseline_available=True,
            matlab_success=True,
            objective_relative_diff=0.02,
            constraint_diff_status="matched",
        )
        == "not_aligned"
    )
    assert (
        classify_alignment_status(
            python_result_status="corrected_condition",
            baseline_available=True,
            matlab_success=True,
            objective_relative_diff=0.0,
            constraint_diff_status="known_matlab_compatibility_differences",
            has_compatibility_exceptions=True,
        )
        == "aligned_except_known_matlab_compatibility_differences"
    )
    assert (
        classify_alignment_status(
            python_result_status="corrected_condition",
            baseline_available=True,
            matlab_success=True,
            objective_relative_diff=0.0,
            constraint_diff_status="row_level_diff_missing",
            has_compatibility_exceptions=True,
        )
        == "pending"
    )
    assert (
        classify_alignment_status(
            python_result_status="matlab_aligned",
            baseline_available=True,
            matlab_success=True,
            objective_relative_diff=0.0,
            constraint_diff_status="row_level_diff_missing",
        )
        == "not_aligned"
    )


def test_hlf_matlab_failure_is_not_marked_aligned_when_historical_probe_artifact_exists() -> None:
    artifact = {
        "artifact_path": "local_runs/pichia_hlf_opn_probe/hlf_matlab_harness_probe_2026-06-24/hlf_matlab_harness_probe_summary.json",
        "baseline_available": True,
        "targets": [
            {
                "target_id": "hLF_CLEAN",
                "success": True,
                "production_ratio_fixed": 0.0010297612727174989,
                "lp_stats": {"bound_lines": 29068, "constraint_lines": 24444},
                "fake_info_summary": {
                    "sequence_length": 691,
                    "dsb": 21,
                    "ng": 4,
                    "sequence_ends_with_stop": False,
                },
                "enzymedataTP_proteinMWs": 76164.54,
            },
            {
                "target_id": "hLF",
                "success": False,
                "error_identifier": "MATLAB:sizeDimensionsMustMatch",
                "error_stack": ["calculateMW.m:35", "addTargetProtein.m:118"],
            },
        ],
    }

    original_hlf = build_alignment_summary(
        "hLF",
        python_result_status="corrected_condition",
        artifact=artifact,
    )

    assert original_hlf.target_id == "hLF"
    assert original_hlf.matlab_alignment_status == "matlab_failed"
    assert "calculateMW.m:35" in original_hlf.diagnostic_message
    assert original_hlf.success is False


def test_hlf_project_710_artifact_can_be_marked_aligned_except_known_compatibility_differences() -> None:
    artifact = load_matlab_alignment_artifact(HLF_PROJECT_710_ARTIFACT_PATH)

    summary = build_alignment_summary(
        "hLF_PROJECT_710",
        python_result_status="corrected_condition",
        artifact_path=HLF_PROJECT_710_ARTIFACT_PATH,
        artifact=artifact,
        rows_python=24444,
        cols_python=29068,
        objective_python=0.001112112385054876,
        constraint_diff_status="known_matlab_compatibility_differences",
        compatibility_exceptions=hlf_project_710_known_matlab_compatibility_exceptions(),
    )
    payload = summarize_alignment(summary)
    exceptions = {item["id"]: item for item in payload["compatibility_exceptions"]}

    assert artifact["baseline_available"] is True
    assert len(artifact["targets"]) == 1
    assert summary.target_id == "hLF_PROJECT_710"
    assert summary.baseline_available is True
    assert summary.matlab_success is True
    assert summary.rows_diff == 0
    assert summary.cols_diff == 0
    assert summary.objective_relative_diff == pytest.approx(0.0)
    assert summary.matlab_alignment_status == "aligned_except_known_matlab_compatibility_differences"
    assert summary.success is False
    assert payload["is_fully_aligned"] is False
    assert payload["is_aligned_except_known_matlab_compatibility_differences"] is True
    assert exceptions["corrected_medium_exchange_bounds"]["count"] == 9
    assert exceptions["misfolding_dilution_bounds"]["count"] == 1418
    assert exceptions["ribosome_optional_row_mapping"]["count"] == 2


def test_hlf_project_710_artifact_does_not_override_original_hlf_status() -> None:
    artifact = load_matlab_alignment_artifact(HLF_PROJECT_710_ARTIFACT_PATH)

    original_hlf = build_alignment_summary(
        "hLF",
        python_result_status="corrected_condition",
        artifact=artifact,
    )
    payload = summarize_alignment(original_hlf)

    assert original_hlf.target_id == "hLF"
    assert original_hlf.matlab_success is None
    assert original_hlf.matlab_alignment_status == "pending"
    assert original_hlf.success is False
    assert "Target is missing" in original_hlf.diagnostic_message
    assert payload["is_fully_aligned"] is False
    assert payload["is_aligned_except_known_matlab_compatibility_differences"] is False


# --- ADR-008 的守卫：与旧 MATLAB 实现「对照」的声称边界 -------------------------
# 这些断言是**纯函数**的，不读任何产物文件——ADR-008 要守的是判定规则本身，
# 而基线产物在 local_runs/（gitignored），依赖它就会在干净 clone 上静默变结论。


def test_corrected_condition_never_yields_an_alignment_claim() -> None:
    """ADR-008 决策 1：修正条件下恒 `pending`，与其余参数无关。

    修正后的培养基与旧 MATLAB 基线不是同一个条件，跨条件比对即使数值接近也不
    构成「对齐」。界面上的「待对齐」是这条边界的表现，不是缺陷——谁想把这行
    「修好」，先撞这条守卫，再读 ADR-008。
    """
    for matlab_success in (True, False, None):
        for objective_relative_diff in (0.0, 1e-9, 0.5, None):
            for constraint_diff_status in ("matched", "pending", "known_matlab_compatibility_differences"):
                status = classify_alignment_status(
                    python_result_status="corrected_condition",
                    baseline_available=True,
                    matlab_success=matlab_success,
                    objective_relative_diff=objective_relative_diff,
                    constraint_diff_status=constraint_diff_status,
                )
                # matlab_failed 优先级更高（基线侧自己没跑成功，也谈不上对照）
                expected = "matlab_failed" if matlab_success is False else "pending"
                assert status == expected, (
                    f"corrected_condition 下得到 {status!r}"
                    f"（matlab_success={matlab_success}, diff={objective_relative_diff},"
                    f" constraints={constraint_diff_status}）——违反 ADR-008 决策 1"
                )


def test_alignment_status_vocabulary_is_a_closed_set() -> None:
    """ADR-008 决策 2：七个状态是封闭集合，各自只说一件事。

    断言的是**集合相等**而不是包含：新增一个状态词必须是有意识的决定，
    并且要同步更新 ADR-008 的状态表，否则界面和下游会各自解释。
    """
    assert ALIGNMENT_STATUSES == {
        "pending",
        "baseline_missing",
        "matlab_failed",
        "python_draft",
        "aligned",
        "not_aligned",
        "aligned_except_known_matlab_compatibility_differences",
    }


def test_known_compatibility_status_requires_registered_exceptions() -> None:
    """ADR-008 决策 2 续：未登记的差异不得归入「已知兼容性差异」。

    否则这个状态会退化成「凡是对不上的都算已知差异」。
    """
    common = dict(
        python_result_status="matlab_aligned",
        baseline_available=True,
        matlab_success=True,
        constraint_diff_status="known_matlab_compatibility_differences",
    )
    assert classify_alignment_status(has_compatibility_exceptions=True, **common) == (
        "aligned_except_known_matlab_compatibility_differences"
    )
    assert classify_alignment_status(has_compatibility_exceptions=False, **common) != (
        "aligned_except_known_matlab_compatibility_differences"
    )
