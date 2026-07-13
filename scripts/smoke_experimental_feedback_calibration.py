from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PICHIA_SRC = REPO_ROOT / "python_pichia" / "src"
if str(PYTHON_PICHIA_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_PICHIA_SRC))

from pcsec_pichia.experimental_feedback import (  # noqa: E402
    CalibrationConfig,
    run_experiment_feedback_replay,
)


FIXTURE_ROOT = REPO_ROOT / "python_pichia" / "tests" / "fixtures" / "experimental_feedback"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("local_runs/experiment_feedback/round5_replay"),
    )
    args = parser.parse_args()
    prediction_path = FIXTURE_ROOT / "round5_replay_predictions.json"
    prediction_run = json.loads(prediction_path.read_text(encoding="utf-8"))
    result = run_experiment_feedback_replay(
        experiment_path=FIXTURE_ROOT / "round5_replay_experiments.jsonl",
        prediction_runs=(prediction_run,),
        prediction_sources=(prediction_path,),
        output_dir=args.output_dir,
        source_classification="sanitized_fixture",
        config=CalibrationConfig(
            increase_threshold_ratio=1.10,
            decrease_threshold_ratio=0.90,
            top_k=(1, 2, 3),
            baseline_hit_rate=0.25,
            minimum_rank_pairs=2,
        ),
    )
    by_target = {target.target_id: target for target in result.calibration.targets}
    preserved_reasons = {
        reason
        for record in result.calibration.records
        for reason in record.ineligibility_reasons
    }
    if set(by_target) != {"hLF", "OPN"}:
        return 1
    if result.linkage.ambiguous_count != 1:
        return 1
    if not any(record.hit is False for record in result.calibration.records):
        return 1
    if not any("assay_failed" in record.measurement_statuses for record in result.calibration.records):
        return 1
    if "control_match_missing" not in preserved_reasons:
        return 1
    print(
        json.dumps(
            {
                "source_classification": "sanitized_fixture",
                "manifest_path": str(result.outputs.manifest_path),
                "summary_path": str(result.outputs.summary_path),
                "report_path": str(result.outputs.report_path),
                "hLF": {
                    "eligible": by_target["hLF"].eligible_count,
                    "ineligible": by_target["hLF"].ineligible_count,
                },
                "OPN": {
                    "eligible": by_target["OPN"].eligible_count,
                    "ineligible": by_target["OPN"].ineligible_count,
                },
                "ambiguous_linkage": result.linkage.ambiguous_count,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
