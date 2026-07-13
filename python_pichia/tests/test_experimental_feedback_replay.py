from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from pcsec_pichia.experimental_feedback import (
    CalibrationConfig,
    run_experiment_feedback_replay,
)


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "experimental_feedback"
REPO_ROOT = Path(__file__).resolve().parents[2]


def test_sanitized_hlf_opn_replay_writes_auditable_report_bundle(tmp_path) -> None:
    prediction_path = FIXTURE_ROOT / "round5_replay_predictions.json"
    prediction_run = json.loads(prediction_path.read_text(encoding="utf-8"))

    result = run_experiment_feedback_replay(
        experiment_path=FIXTURE_ROOT / "round5_replay_experiments.jsonl",
        prediction_runs=(prediction_run,),
        prediction_sources=(prediction_path,),
        output_dir=tmp_path / "replay",
        source_classification="sanitized_fixture",
        config=CalibrationConfig(
            increase_threshold_ratio=1.10,
            decrease_threshold_ratio=0.90,
            top_k=(1, 2, 3),
            baseline_hit_rate=0.25,
            minimum_rank_pairs=2,
        ),
    )

    assert result.validation.is_valid is True
    assert result.linkage.matched_count == 4
    assert result.linkage.ambiguous_count == 1
    by_target = {target.target_id: target for target in result.calibration.targets}
    assert (by_target["hLF"].eligible_count, by_target["hLF"].ineligible_count) == (1, 2)
    assert (by_target["OPN"].eligible_count, by_target["OPN"].ineligible_count) == (1, 1)

    summary = json.loads(result.outputs.summary_path.read_text(encoding="utf-8"))
    assert summary["targets"]["hLF"]["ranking_assessment"] == "insufficient_evidence"
    assert summary["targets"]["hLF"]["confirmed_direction_candidates"][0]["gene_id"] == "G-H1"
    assert summary["targets"]["hLF"]["confirmed_direction_candidates"][0]["intervention_type"] == "KO"
    assert summary["targets"]["OPN"]["direction_discordant_count"] == 1
    assert summary["targets"]["OPN"]["direction_discordant_candidates"][0]["gene_id"] == "G-O1"
    assert summary["targets"]["OPN"]["direction_discordant_candidates"][0]["intervention_type"] == "OE"
    assert summary["preserved_status_counts"] == {
        "assay_failed": 1,
        "candidate_measurement_not_evaluable": 1,
        "control_match_missing": 1,
        "negative_observation": 1,
        "prediction_link:ambiguous": 1,
    }

    manifest = json.loads(result.outputs.manifest_path.read_text(encoding="utf-8"))
    assert manifest["source_classification"] == "sanitized_fixture"
    assert manifest["uses_real_experiment_data"] is False
    assert manifest["promotes_curated_data"] is False
    assert manifest["mutates_recommendation_tier"] is False
    assert manifest["mutates_model_constraints"] is False
    assert manifest["files"]["report"] == "prediction_experiment_report.md"

    report = result.outputs.report_path.read_text(encoding="utf-8")
    assert "# Prediction vs Experiment 回放报告" in report
    assert "## hLF" in report
    assert "## OPN" in report
    assert "G-H1" in report
    assert "G-O1" in report
    assert "assay_failed" in report
    assert "control_match_missing" in report
    assert "prediction_link:ambiguous" in report
    assert "未使用真实实验数据" in report
    assert "不会自动修改 recommendation tier 或模型约束" in report


def test_round5_replay_smoke_cli_uses_sanitized_import_fixture(tmp_path) -> None:
    output_dir = tmp_path / "cli-replay"
    completed = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "smoke_experimental_feedback_calibration.py"),
            "--output-dir",
            str(output_dir),
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert (output_dir / "report" / "prediction_experiment_manifest.json").exists()
    assert (output_dir / "report" / "prediction_experiment_summary.json").exists()
    assert (output_dir / "report" / "prediction_experiment_report.md").exists()
    assert "sanitized_fixture" in completed.stdout
