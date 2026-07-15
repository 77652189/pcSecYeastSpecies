from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Mapping

from pcsec_pichia.oe_capacity.schema import (
    OECapacityOutputs,
    OECapacityScreenResult,
    OECapacityScreenRow,
)


def write_oe_capacity_outputs(
    result: OECapacityScreenResult,
    output_dir: str | Path,
    *,
    run_identity: Mapping[str, Any] | None = None,
    capacity_asset: Mapping[str, Any] | None = None,
) -> OECapacityOutputs:
    result.validate()
    resolved = Path(output_dir)
    resolved.mkdir(parents=True, exist_ok=True)
    outputs = OECapacityOutputs(
        output_dir=str(resolved),
        rows_path=str(resolved / "oe_capacity_rows.jsonl"),
        manifest_path=str(resolved / "oe_capacity_manifest.json"),
        report_path=str(resolved / "oe_capacity_report.md"),
    )
    outputs.validate()
    rows_path = Path(outputs.rows_path)
    with rows_path.open("w", encoding="utf-8", newline="\n") as handle:
        for screen_status, rows in (
            ("completed", result.rows),
            ("failed", result.failures),
        ):
            for row in rows:
                payload = _json_ready(asdict(row))
                payload["screen_status"] = str(
                    getattr(row, "screen_status", "") or screen_status
                )
                handle.write(
                    json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n"
                )
    Path(outputs.report_path).write_text(
        _markdown_report(result),
        encoding="utf-8",
    )
    report_path = Path(outputs.report_path)
    row_payloads = [
        {
            **_json_ready(asdict(row)),
            "screen_status": str(getattr(row, "screen_status", "") or status),
        }
        for status, rows in (("completed", result.rows), ("failed", result.failures))
        for row in rows
    ]
    target_ids = sorted({str(row.get("target_id") or "") for row in row_payloads} - {""})
    context_ids = sorted({str(row.get("context_id") or "") for row in row_payloads} - {""})
    gene_ids = sorted({str(row.get("gene_id") or "") for row in row_payloads} - {""})
    required_scenarios = [scenario.value for scenario in result.config.parameter_scenarios]
    complete_scenarios = _scenario_evidence_complete(
        row_payloads,
        required_scenarios,
        feature_enabled=result.config.feature_enabled,
    )
    relative_complete = _relative_scenario_evidence_complete(
        row_payloads,
        required_scenarios,
    )
    identity = {
        "run_id": resolved.name,
        "target_ids": target_ids,
        "context_ids": context_ids,
        "gene_ids": gene_ids,
        "case_kind": "screen",
    }
    if run_identity:
        identity.update(_json_ready(dict(run_identity)))
    asset = {
        "path": "",
        "version": "",
        "sha256": "",
        "reviewed": False,
    }
    if capacity_asset:
        asset.update(_json_ready(dict(capacity_asset)))
    state = (
        "partial_failure"
        if result.rows and result.failures
        else "failed"
        if result.failures
        else "completed"
    )
    execution_status_counts: dict[str, int] = {}
    product_state_counts: dict[str, int] = {}
    for payload in row_payloads:
        key = str(payload.get("execution_status") or "unknown")
        execution_status_counts[key] = execution_status_counts.get(key, 0) + 1
        product_key = str(payload.get("product_state") or "unknown")
        product_state_counts[product_key] = product_state_counts.get(product_key, 0) + 1
    manifest = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_identity": identity,
        "model": {"fingerprint": result.model_fingerprint},
        "capacity_asset": asset,
        "model_fingerprint": result.model_fingerprint,
        "config": _json_ready(asdict(result.config)),
        "completed_count": len(result.rows),
        "failure_count": len(result.failures),
        "status": {
            "state": state,
            "completed_count": len(result.rows),
            "failure_count": len(result.failures),
        },
        "files": {
            "rows": {
                "path": rows_path.name,
                "sha256": _sha256(rows_path),
            },
            "report": {
                "path": report_path.name,
                "sha256": _sha256(report_path),
            },
        },
        "coverage": {
            "total_rows": len(row_payloads),
            "by_execution_status": execution_status_counts,
            "by_product_state": product_state_counts,
        },
        "scenario_completeness": {
            "required": required_scenarios if result.config.feature_enabled else [],
            "complete": complete_scenarios,
        },
        "absolute_scenario_completeness": {
            "required": required_scenarios if result.config.feature_enabled else [],
            "complete": complete_scenarios,
        },
        "relative_scenario_definition_complete": relative_complete,
        "execution_modes": sorted(
            {row.execution_mode.value for row in (*result.rows, *result.failures)}
        ),
        "product_states": sorted(
            {row.product_state.value for row in (*result.rows, *result.failures)}
        ),
        "absolute_capacity_available": any(
            row.absolute_solver_allowed for row in (*result.rows, *result.failures)
        ),
        "model_relative_only": not any(
            row.absolute_solver_allowed for row in (*result.rows, *result.failures)
        ),
        "predicts_absolute_yield": False,
        "mutates_recommendation_tier": False,
        "mutates_model_assets": False,
        "warnings": list(result.warnings),
    }
    manifest_path = Path(outputs.manifest_path)
    temporary_manifest_path = manifest_path.with_name(f".{manifest_path.name}.tmp")
    temporary_manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary_manifest_path.replace(manifest_path)
    return outputs


def _scenario_evidence_complete(
    rows: list[dict[str, Any]],
    required: list[str],
    *,
    feature_enabled: bool,
) -> bool:
    if not feature_enabled or not required:
        return True
    executable = [
        row
        for row in rows
        if row.get("product_state") == "absolute_available"
        and row.get("screen_status") == "completed"
    ]
    if not executable:
        return False
    return all(_row_has_scenarios(row, required) for row in executable)


def _row_has_scenarios(row: Mapping[str, Any], required: list[str]) -> bool:
    raw = row.get("scenario_results") or row.get("gene_capacity_scenario_results")
    if raw is None:
        # Schema v1 rows expose only scenario labels. Keep the manifest honest: the
        # labels describe requested scenarios, not solver evidence.
        return False
    found: set[str] = set()
    for item in raw if isinstance(raw, list) else ():
        if not isinstance(item, Mapping):
            continue
        scenario = item.get("scenario") or item.get("parameter_scenario")
        if isinstance(scenario, Mapping):
            scenario = scenario.get("value")
        baseline = item.get("baseline") or item.get("baseline_snapshot")
        perturbed = item.get("perturbed") or item.get("perturbed_snapshot")
        if (
            str(scenario) in required
            and isinstance(baseline, Mapping)
            and isinstance(perturbed, Mapping)
        ):
            found.add(str(scenario))
    return set(required).issubset(found)


def _relative_scenario_evidence_complete(
    rows: list[dict[str, Any]],
    required: list[str],
) -> bool:
    relative_rows = [
        row for row in rows if row.get("product_state") == "relative_uncalibrated"
    ]
    if not relative_rows:
        return True
    for row in relative_rows:
        factors = {
            str(item[0])
            for item in row.get("relative_capacity_factors") or []
            if isinstance(item, list) and len(item) == 2
        }
        evidence = {
            str(item.get("parameter_scenario") or item.get("scenario") or "")
            for item in row.get("relative_scenario_results") or []
            if isinstance(item, Mapping)
        }
        if not set(required).issubset(factors) or not set(required).issubset(evidence):
            return False
    return True


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _markdown_report(result: OECapacityScreenResult) -> str:
    lines = [
        "# Gene-level OE capacity comparison",
        "",
        "This report separates reaction proxy, relative uncalibrated OE, absolute "
        "capacity availability, and not-executable outcomes. It does not predict mg/L, true expression "
        "fold-change, experimental success probability, or recommendation tier.",
        "",
        f"- Model fingerprint: `{result.model_fingerprint}`",
        f"- Completed rows: {len(result.rows)}",
        f"- Failed or non-executable rows: {len(result.failures)}",
        f"- Gene-capacity feature enabled: {result.config.feature_enabled}",
        f"- Legacy proxy comparison enabled: {result.config.compare_proxy}",
        "",
        "| Status | Gene | Target | Context | Product state | Calibration | Absolute availability | Mode | Execution status | Baseline | Proxy | "
        "Relative objective | Relative delta | Gene capacity | Nominal capacity | Delta vs baseline | Delta vs proxy | Resource cost delta | "
        "Missing information |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for status, rows in (("completed", result.rows), ("failed", result.failures)):
        for row in rows:
            lines.append(_markdown_row(row.screen_status or status, row))
    lines.extend(("", "## Mapping and parameter traceability", ""))
    for status, rows in (("completed", result.rows), ("failed", result.failures)):
        for row in rows:
            identity = f"{row.target_id}/{row.gene_id}/{row.context_id}/{row.dose_id}"
            lines.extend(
                (
                    f"### `{identity}` ({row.screen_status or status})",
                    "",
                    f"- Product mode/state: `{row.product_mode.value}` / `{row.product_state.value}`",
                    f"- Calibration: `{row.calibration_status.value}`",
                    f"- Absolute availability: `{row.absolute_capacity_availability.value}`",
                    f"- Absolute solver allowed: `{row.absolute_solver_allowed}`",
                    f"- Model fingerprint: `{row.model_fingerprint}`",
                    f"- Mapping IDs: {', '.join(row.mapping_ids) or 'none'}",
                    f"- Parameter sources: {', '.join(row.parameter_sources) or 'none'}",
                    "- Parameter confidence: "
                    + (
                        row.parameter_confidence.value
                        if row.parameter_confidence is not None
                        else "not available"
                    ),
                    "- Uncertainty scenarios: "
                    + (
                        ", ".join(
                            scenario.value for scenario in row.uncertainty_scenarios
                        )
                        or "none"
                    ),
                    f"- Warnings: {', '.join(row.warnings) or 'none'}",
                    f"- Limitations: {', '.join(row.limitations) or 'none'}",
                    "",
                )
            )
            if row.scenario_results:
                lines.extend(("#### Scenario solver evidence", ""))
                for scenario_result in row.scenario_results:
                    lines.append(
                        "- "
                        f"`{scenario_result.parameter_scenario.value}`: "
                        f"baseline=`{scenario_result.baseline.solver_status}` "
                        f"(success={scenario_result.baseline.success}); "
                        f"perturbed=`{scenario_result.perturbed.solver_status}` "
                        f"(success={scenario_result.perturbed.success}); "
                        f"failure=`{scenario_result.failure_reason or 'none'}`; "
                        f"message=`{scenario_result.perturbed.message or scenario_result.baseline.message or 'none'}`"
                    )
                lines.append("")
            if row.relative_scenario_results:
                lines.extend(("#### Relative scenario solver evidence", ""))
                for scenario_result in row.relative_scenario_results:
                    lines.append(
                        "- "
                        f"`{scenario_result.parameter_scenario.value}`: "
                        f"baseline=`{scenario_result.baseline.solver_status}` "
                        f"(success={scenario_result.baseline.success}); "
                        f"perturbed=`{scenario_result.perturbed.solver_status}` "
                        f"(success={scenario_result.perturbed.success}); "
                        f"delta=`{_number(scenario_result.objective_delta)}`; "
                        f"failure=`{scenario_result.failure_reason or 'none'}`; "
                        f"message=`{scenario_result.perturbed.message or scenario_result.baseline.message or 'none'}`"
                    )
                lines.append("")
            if row.proxy_attempts:
                lines.extend(("#### Reaction proxy attempts", ""))
                for attempt in row.proxy_attempts:
                    lines.append(
                        "- "
                        f"`{attempt.attempt_id or 'unnamed'}`: "
                        f"status=`{attempt.solver_status}`, success=`{attempt.success}`, "
                        f"objective=`{_number(attempt.secretion_objective)}`, "
                        f"message=`{attempt.message or 'none'}`"
                    )
                lines.append("")
    if result.warnings:
        lines.extend(("", "## Screen warnings", ""))
        lines.extend(f"- {warning}" for warning in result.warnings)
    lines.append("")
    return "\n".join(lines)


def _markdown_row(status: str, row: OECapacityScreenRow) -> str:
    values = (
        status,
        row.gene_id,
        row.target_id,
        row.context_id,
        row.product_state.value,
        row.calibration_status.value,
        row.absolute_capacity_availability.value,
        row.execution_mode.value,
        row.execution_status.value,
        _number(row.baseline_objective),
        _number(row.proxy_objective),
        _number(row.relative_objective),
        _number(row.relative_vs_baseline_delta),
        _number(row.gene_capacity_objective),
        _number(row.nominal_capacity),
        _number(row.gene_capacity_vs_baseline_delta),
        _number(row.gene_capacity_vs_proxy_delta),
        _number(row.protein_resource_cost_delta),
        ", ".join(row.missing_information) or "none",
    )
    return "| " + " | ".join(_escape_markdown(value) for value in values) + " |"


def _number(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.8g}"


def _escape_markdown(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _json_ready(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = ["write_oe_capacity_outputs"]
