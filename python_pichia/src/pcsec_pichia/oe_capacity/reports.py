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
                payload["screen_status"] = screen_status
                handle.write(
                    json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n"
                )
    Path(outputs.report_path).write_text(
        _markdown_report(result),
        encoding="utf-8",
    )
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "model_fingerprint": result.model_fingerprint,
        "config": _json_ready(asdict(result.config)),
        "completed_count": len(result.rows),
        "failure_count": len(result.failures),
        "files": {
            "rows": Path(outputs.rows_path).name,
            "report": Path(outputs.report_path).name,
        },
        "rows_sha256": hashlib.sha256(rows_path.read_bytes()).hexdigest(),
        "execution_modes": sorted(
            {row.execution_mode.value for row in (*result.rows, *result.failures)}
        ),
        "model_relative_only": True,
        "predicts_absolute_yield": False,
        "mutates_recommendation_tier": False,
        "mutates_model_assets": False,
        "warnings": list(result.warnings),
    }
    Path(outputs.manifest_path).write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return outputs


def _markdown_report(result: OECapacityScreenResult) -> str:
    lines = [
        "# Gene-level OE capacity comparison",
        "",
        "This report compares model-relative baseline, legacy reaction proxy, and "
        "gene-enzyme capacity scenarios. It does not predict mg/L, true expression "
        "fold-change, experimental success probability, or recommendation tier.",
        "",
        f"- Model fingerprint: `{result.model_fingerprint}`",
        f"- Completed rows: {len(result.rows)}",
        f"- Failed or non-executable rows: {len(result.failures)}",
        f"- Gene-capacity feature enabled: {result.config.feature_enabled}",
        f"- Legacy proxy comparison enabled: {result.config.compare_proxy}",
        "",
        "| Status | Gene | Target | Mode | Execution status | Baseline | Proxy | "
        "Gene capacity | Delta vs baseline | Delta vs proxy | Resource cost delta | "
        "Missing information |",
        "| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for status, rows in (("completed", result.rows), ("failed", result.failures)):
        for row in rows:
            lines.append(_markdown_row(status, row))
    lines.extend(("", "## Mapping and parameter traceability", ""))
    for status, rows in (("completed", result.rows), ("failed", result.failures)):
        for row in rows:
            identity = f"{row.target_id}/{row.gene_id}/{row.context_id}/{row.dose_id}"
            lines.extend(
                (
                    f"### `{identity}` ({status})",
                    "",
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
                    "",
                )
            )
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
        row.execution_mode.value,
        row.execution_status.value,
        _number(row.baseline_objective),
        _number(row.proxy_objective),
        _number(row.gene_capacity_objective),
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
