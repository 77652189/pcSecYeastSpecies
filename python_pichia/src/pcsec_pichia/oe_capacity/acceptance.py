from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence


EXPECTED_TARGET_IDS = ("hLF", "OPN_ALPHA_FULL_PROJECT")


@dataclass(frozen=True)
class OECapacityAcceptanceObservation:
    target_id: str
    gene_id: str
    case_kind: str
    elapsed_seconds: float
    screen_status: str
    execution_status: str
    baseline_objective: float | None
    proxy_objective: float | None
    gene_capacity_objective: float | None
    protein_resource_cost_delta: float | None = None
    missing_information: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    output_dir: str = ""

    def validate(self) -> None:
        if not self.target_id.strip() or not self.gene_id.strip():
            raise ValueError("acceptance observation target_id and gene_id are required.")
        if self.case_kind not in {"executable", "boundary"}:
            raise ValueError("case_kind must be executable or boundary.")
        if not math.isfinite(self.elapsed_seconds) or self.elapsed_seconds < 0:
            raise ValueError("elapsed_seconds must be finite and non-negative.")
        for name, value in (
            ("baseline_objective", self.baseline_objective),
            ("proxy_objective", self.proxy_objective),
            ("gene_capacity_objective", self.gene_capacity_objective),
            ("protein_resource_cost_delta", self.protein_resource_cost_delta),
        ):
            if value is not None and not math.isfinite(value):
                raise ValueError(f"{name} must be finite when present.")


@dataclass(frozen=True)
class OECapacityRegressionCheck:
    check_id: str
    passed: bool
    evidence: str

    def validate(self) -> None:
        if not self.check_id.strip() or not self.evidence.strip():
            raise ValueError("regression check id and evidence are required.")


def build_oe_capacity_acceptance_summary(
    observations: Sequence[OECapacityAcceptanceObservation],
    *,
    coverage_by_target: Mapping[str, Mapping[str, object]],
    regression_checks: Sequence[OECapacityRegressionCheck],
    expected_target_ids: Sequence[str] = EXPECTED_TARGET_IDS,
) -> dict[str, object]:
    for observation in observations:
        observation.validate()
    for check in regression_checks:
        check.validate()

    expected = tuple(
        dict.fromkeys(
            str(value).strip() for value in expected_target_ids if str(value).strip()
        )
    )
    target_rows: list[dict[str, object]] = []
    for target_id in expected:
        target_observations = [row for row in observations if row.target_id == target_id]
        coverage = dict(coverage_by_target.get(target_id) or {})
        executable = [row for row in target_observations if row.case_kind == "executable"]
        boundary = [row for row in target_observations if row.case_kind == "boundary"]
        executable_passed = any(
            row.screen_status == "completed"
            and row.execution_status == "gene_level_executable"
            and not row.missing_information
            and row.baseline_objective is not None
            and row.proxy_objective is not None
            and row.gene_capacity_objective is not None
            for row in executable
        )
        boundary_passed = any(
            row.screen_status == "completed" and bool(row.missing_information)
            for row in boundary
        )
        target_rows.append(
            {
                "target_id": target_id,
                "coverage": coverage,
                "coverage_present": bool(coverage),
                "executable_passed": executable_passed,
                "boundary_passed": boundary_passed,
                "observations": [asdict(row) for row in target_observations],
                "max_elapsed_seconds": max(
                    (row.elapsed_seconds for row in target_observations),
                    default=0.0,
                ),
            }
        )

    checks = [asdict(check) for check in regression_checks]
    passed = bool(target_rows) and all(
        row["coverage_present"]
        and row["executable_passed"]
        and row["boundary_passed"]
        for row in target_rows
    ) and bool(checks) and all(check["passed"] for check in checks)
    return {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "phase_2_gene_level_oe",
        "passed": passed,
        "expected_target_ids": list(expected),
        "targets": target_rows,
        "regression_checks": checks,
        "model_relative_only": True,
        "predicts_absolute_yield": False,
        "mutates_recommendation_tier": False,
    }


def write_oe_capacity_acceptance_outputs(
    summary: Mapping[str, object],
    output_dir: str | Path,
) -> tuple[Path, Path]:
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    json_path = root / "phase2_acceptance.json"
    markdown_path = root / "phase2_acceptance.md"
    json_path.write_text(
        json.dumps(dict(summary), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    targets = summary.get("targets") or []
    checks = summary.get("regression_checks") or []
    lines = [
        "# Phase 2 gene-level OE capacity acceptance",
        "",
        f"- Passed: `{bool(summary.get('passed'))}`",
        "- Scope: model-relative baseline / reaction-proxy / gene-capacity comparison",
        "- Absolute yield or experimental success prediction: `false`",
        "",
        "## Target smoke",
        "",
    ]
    for target in targets:
        if not isinstance(target, Mapping):
            continue
        lines.append(
            f"- `{target.get('target_id')}`: executable=`{target.get('executable_passed')}`, "
            f"boundary=`{target.get('boundary_passed')}`, "
            f"coverage=`{target.get('coverage_present')}`, "
            f"max_elapsed_seconds=`{target.get('max_elapsed_seconds')}`"
        )
    lines.extend(["", "## Regression and safety checks", ""])
    for check in checks:
        if not isinstance(check, Mapping):
            continue
        lines.append(
            f"- `{check.get('check_id')}`: `{check.get('passed')}` — {check.get('evidence')}"
        )
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, markdown_path


__all__ = [
    "EXPECTED_TARGET_IDS",
    "OECapacityAcceptanceObservation",
    "OECapacityRegressionCheck",
    "build_oe_capacity_acceptance_summary",
    "write_oe_capacity_acceptance_outputs",
]
