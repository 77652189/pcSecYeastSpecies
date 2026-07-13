from __future__ import annotations

import json
import re
from dataclasses import asdict
from enum import Enum
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping, Sequence

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.loading import load_pcsec_pichia_inputs
from pcsec_pichia.oe_capacity import (
    OECapacityScreenConfig,
    OECapacityScreenRequest,
    OEExecutionMode,
    ParameterPolicy,
    ParameterScenario,
    build_current_model_parameter_policy,
    build_gene_enzyme_reaction_catalog,
    build_oe_dose_spec,
    run_gene_level_oe_screen,
    summarize_gene_capacity_catalog,
    write_oe_capacity_outputs,
)
from pcsec_pichia.screens import prepare_screen_inputs
from pcsec_pichia.targets import load_builtin_targets


DEFAULT_OE_CAPACITY_ROOT = Path("local_runs") / "oe_capacity" / "ui_runs"
DEFAULT_TARGET_IDS = ("hLF", "OPN_ALPHA_FULL_PROJECT")


def list_oe_capacity_targets() -> list[dict[str, str]]:
    targets = _target_lookup()
    rows: list[dict[str, str]] = []
    for target_id in DEFAULT_TARGET_IDS:
        target = targets.get(target_id)
        if target is None:
            continue
        rows.append(
            {
                "target_id": target_id,
                "label": "hLF" if target_id == "hLF" else "OPN",
                "protein_id": str(target.protein_id),
                "source": str(target.source),
            }
        )
    return rows


def preview_oe_capacity_candidate(
    *,
    target_id: str,
    gene_id: str,
    growth_rate: float = 0.1,
    carbon_source_id: str = "glucose",
    relative_uncertainty: float = 0.2,
) -> dict[str, Any]:
    normalized_gene = _required_text(gene_id, "gene_id")
    runtime = _prepare_runtime(
        target_id,
        float(growth_rate),
        carbon_source_id,
        float(relative_uncertainty),
    )
    mappings = [
        mapping
        for mapping in runtime.gene_capacity_catalog.mappings
        if mapping.gene_id == normalized_gene
    ]
    parameter_sets = [
        parameter_set
        for parameter_set in runtime.parameter_policy.parameter_sets
        if parameter_set.gene_id == normalized_gene
    ]
    coverage = summarize_gene_capacity_catalog(runtime.gene_capacity_catalog)
    return {
        "target_id": target_id,
        "gene_id": normalized_gene,
        "context_id": _context_id(carbon_source_id, growth_rate),
        "mapping_count": len(mappings),
        "parameter_set_count": len(parameter_sets),
        "mappings": [_json_ready(asdict(mapping)) for mapping in mappings],
        "parameter_sets": [
            _json_ready(asdict(parameter_set)) for parameter_set in parameter_sets
        ],
        "coverage": _json_ready(asdict(coverage)),
        "executable_mapping_count": dict(coverage.by_status).get(
            "gene_level_executable",
            0,
        ),
        "execution_boundary": (
            "Only current-model gene/enzyme/reaction mappings with canonical local "
            "parameters can execute gene-level capacity scenarios."
        ),
    }


def submit_oe_capacity_screen(
    *,
    target_id: str,
    gene_ids: Sequence[str],
    dose_payload: Mapping[str, Any],
    parameter_scenarios: Sequence[str] = ("low", "nominal", "high"),
    execution_mode: str = "comparison",
    feature_enabled: bool = True,
    compare_proxy: bool = True,
    growth_rate: float = 0.1,
    carbon_source_id: str = "glucose",
    relative_uncertainty: float = 0.2,
    run_name: str = "",
    output_root: str | Path = DEFAULT_OE_CAPACITY_ROOT,
) -> dict[str, Any]:
    normalized_genes = _dedupe_gene_ids(gene_ids)
    if not normalized_genes:
        raise ValueError("at least one gene_id is required.")
    runtime = _prepare_runtime(
        target_id,
        float(growth_rate),
        carbon_source_id,
        float(relative_uncertainty),
    )
    dose = build_oe_dose_spec(dict(dose_payload))
    try:
        mode = OEExecutionMode(execution_mode)
    except ValueError as exc:
        raise ValueError(f"unsupported execution_mode: {execution_mode}") from exc
    scenarios = tuple(ParameterScenario(value) for value in parameter_scenarios)
    context_id = _context_id(carbon_source_id, growth_rate)
    requests = tuple(
        OECapacityScreenRequest(
            gene_id=gene_id,
            target_id=target_id,
            context_id=context_id,
            dose=dose,
            execution_mode=mode,
        )
        for gene_id in normalized_genes
    )
    config = OECapacityScreenConfig(
        feature_enabled=feature_enabled,
        compare_proxy=compare_proxy,
        parameter_scenarios=scenarios,
        growth_rate=float(growth_rate),
        solver_options=(("time_limit_seconds", "600"),),
    )
    safe_run_name = _safe_run_name(run_name or f"{target_id}-{context_id}")
    run_dir = Path(output_root) / safe_run_name
    if run_dir.exists():
        raise FileExistsError(f"OE capacity run already exists: {run_dir}")
    result = run_gene_level_oe_screen(runtime, requests, config)
    outputs = write_oe_capacity_outputs(result, run_dir)
    summary = {
        "run_name": safe_run_name,
        "run_dir": str(run_dir),
        "target_id": target_id,
        "context_id": context_id,
        "completed_count": len(result.rows),
        "failure_count": len(result.failures),
        "rows": [_json_ready(asdict(row)) for row in result.rows],
        "failures": [_json_ready(asdict(row)) for row in result.failures],
        "warnings": list(result.warnings),
        "paths": _json_ready(asdict(outputs)),
        "model_relative_only": True,
        "mutates_recommendation_tier": False,
    }
    (run_dir / "ui_run_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def list_oe_capacity_runs(
    output_root: str | Path = DEFAULT_OE_CAPACITY_ROOT,
    *,
    target_id: str = "",
) -> list[dict[str, Any]]:
    root = Path(output_root)
    if not root.exists():
        return []
    rows = [
        load_oe_capacity_run(path)
        for path in root.iterdir()
        if path.is_dir()
    ]
    if target_id:
        rows = [row for row in rows if row.get("target_id") == target_id]
    return sorted(rows, key=lambda row: str(row.get("run_name") or ""), reverse=True)


def load_oe_capacity_run(run_dir: str | Path) -> dict[str, Any]:
    root = Path(run_dir)
    summary_path = root / "ui_run_summary.json"
    if not summary_path.is_file():
        return {
            "run_name": root.name,
            "run_dir": str(root),
            "available": False,
            "rows": [],
            "failures": [],
            "warnings": ["ui_run_summary.json is missing"],
            "paths": {},
        }
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    payload["available"] = True
    return payload


def export_oe_capacity_report(run_dir: str | Path) -> bytes:
    path = Path(run_dir) / "oe_capacity_report.md"
    return path.read_bytes() if path.is_file() else b""


@lru_cache(maxsize=8)
def _prepare_runtime(
    target_id: str,
    growth_rate: float,
    carbon_source_id: str,
    relative_uncertainty: float,
) -> SimpleNamespace:
    targets = _target_lookup()
    try:
        target = targets[target_id]
    except KeyError as exc:
        raise KeyError(f"unknown OE capacity target: {target_id}") from exc
    inputs = load_pcsec_pichia_inputs(
        _repo_root(),
        carbon_source_id=carbon_source_id,
    )
    prepared = prepare_screen_inputs(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_rate,
    )
    if not prepared.get("baseline_success"):
        raise RuntimeError(
            "target baseline preparation failed: "
            + str(prepared.get("baseline_status") or prepared.get("build_status") or "unknown")
        )
    catalog = build_gene_enzyme_reaction_catalog(
        prepared["fixed_model"],
        inputs.metabolic,
        prepared["combined"],
    )
    parameter_policy = build_current_model_parameter_policy(
        catalog,
        prepared["combined"],
        relative_uncertainty=relative_uncertainty,
    )
    return SimpleNamespace(
        target_id=target_id,
        fixed_model=prepared["fixed_model"],
        exchange_reaction_id=prepared["exchange_reaction_id"],
        metabolic=inputs.metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        gene_capacity_catalog=catalog,
        parameter_policy=ParameterPolicy(
            parameter_sets=parameter_policy.parameter_sets,
            scenarios=parameter_policy.scenarios,
            strict_conflicts=parameter_policy.strict_conflicts,
        ),
    )


@lru_cache(maxsize=1)
def _target_lookup() -> dict[str, Any]:
    return {
        target.target_id: target
        for target in load_builtin_targets(_repo_root())
        if target.target_id in DEFAULT_TARGET_IDS
    }


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _context_id(carbon_source_id: str, growth_rate: float) -> str:
    return f"{str(carbon_source_id).strip().lower()}_mu_{float(growth_rate):g}"


def _required_text(value: str, field_name: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"{field_name} must be non-empty.")
    return text


def _dedupe_gene_ids(values: Sequence[str]) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            text
            for value in values
            if (text := str(value or "").strip())
        )
    )


def _safe_run_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", str(value or "").strip()).strip(".-")
    if not cleaned:
        raise ValueError("run_name must contain letters, numbers, dot, underscore, or hyphen.")
    return cleaned


def _json_ready(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "DEFAULT_OE_CAPACITY_ROOT",
    "export_oe_capacity_report",
    "list_oe_capacity_runs",
    "list_oe_capacity_targets",
    "load_oe_capacity_run",
    "preview_oe_capacity_candidate",
    "submit_oe_capacity_screen",
]
