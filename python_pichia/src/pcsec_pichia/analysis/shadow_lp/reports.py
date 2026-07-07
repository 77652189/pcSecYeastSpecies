from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from pcsec_pichia.analysis.shadow_lp.ladder import ShadowLadderResult


NO_ABSOLUTE_YIELD_STATEMENT = (
    "Shadow LP outputs are model-internal relative capacity estimates; "
    "they do not predict absolute fermentation titer or experimental success rate."
)


@dataclass(frozen=True)
class ShadowHardcodeAuditResult:
    production_shadow_path_has_no_reference_objective_literals: bool
    production_shadow_path_has_no_magic_factor_literal: bool
    production_shadow_path_has_no_forbidden_reference_solver_call: bool
    hlf_opn_use_shared_builders: bool
    reference_values_only_used_after_solve_for_comparison: bool
    validation_only_files: tuple[str, ...]
    scanned_files: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def passed(self) -> bool:
        return all(
            (
                self.production_shadow_path_has_no_reference_objective_literals,
                self.production_shadow_path_has_no_magic_factor_literal,
                self.production_shadow_path_has_no_forbidden_reference_solver_call,
                self.hlf_opn_use_shared_builders,
                self.reference_values_only_used_after_solve_for_comparison,
            )
        )

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["passed"] = self.passed
        return payload


def render_shadow_ladder_report_payload(
    ladders: Iterable[ShadowLadderResult],
    audit: ShadowHardcodeAuditResult | None = None,
    comparisons: Iterable[Any] = (),
    validation_matrix: Any | None = None,
) -> dict[str, Any]:
    resolved_ladders = tuple(ladders)
    resolved_audit = audit or run_shadow_hardcode_audit()
    comparison_payloads = tuple(_to_dict(item) for item in comparisons)
    matrix_payload = None if validation_matrix is None else _to_dict(validation_matrix)
    status_mismatches = tuple(
        row
        for row in comparison_payloads
        if row.get("reference_status_category") != row.get("shadow_status_category")
    )
    if matrix_payload is not None:
        status_mismatches = status_mismatches + tuple(
            row
            for row in matrix_payload.get("cases", ())
            if row.get("reference_status") != row.get("shadow_status")
        )
    return {
        "summary": {
            "target_ids": tuple(ladder.target_id for ladder in resolved_ladders),
            "default_large_model_backend": "ScipyHighsBackend",
            "canonical_final_layer": "mitochondrial",
            "production_default_solver_mode": "reference",
            "no_absolute_yield_statement": NO_ABSOLUTE_YIELD_STATEMENT,
        },
        "targets": [ladder.to_dict() for ladder in resolved_ladders],
        "final_objectives": {
            ladder.target_id: {
                "shadow_objective": ladder.final_layer.objective,
                "reference_objective": (ladder.reference_validation or {}).get("reference_objective"),
                "objective_abs_diff": (ladder.reference_validation or {}).get("objective_abs_diff"),
                "objective_rel_diff": (ladder.reference_validation or {}).get("objective_rel_diff"),
                "constraint_count_diff": (ladder.reference_validation or {}).get("constraint_count_diff"),
                "alignment_status": (ladder.reference_validation or {}).get("final_alignment_status"),
            }
            for ladder in resolved_ladders
        },
        "per_layer_objective_changes": {
            ladder.target_id: _objective_changes(ladder)
            for ladder in resolved_ladders
        },
        "constraint_counts": {
            ladder.target_id: {
                layer.layer_id: {
                    "constraint_count": layer.constraint_count,
                    "eq_constraint_count": layer.eq_constraint_count,
                    "ub_constraint_count": layer.ub_constraint_count,
                }
                for layer in ladder.layers
            }
            for ladder in resolved_ladders
        },
        "backend_metadata": {
            ladder.target_id: {
                layer.layer_id: layer.backend_metadata
                for layer in ladder.layers
            }
            for ladder in resolved_ladders
        },
        "skipped_layers": {
            ladder.target_id: {
                layer.layer_id: dict(layer.skipped_layers)
                for layer in ladder.layers
                if layer.skipped_layers
            }
            for ladder in resolved_ladders
        },
        "warnings": {
            ladder.target_id: tuple(
                warning
                for layer in ladder.layers
                for warning in layer.warnings
            )
            for ladder in resolved_ladders
        },
        "compare_mode_summary": comparison_payloads,
        "validation_matrix": matrix_payload,
        "status_mismatches": status_mismatches,
        "no_hardcode_audit": resolved_audit.to_dict(),
    }


def write_shadow_ladder_report(
    ladders: Iterable[ShadowLadderResult],
    output_dir: Path,
    stem: str = "shadow_lp_ladder_report",
    audit: ShadowHardcodeAuditResult | None = None,
    comparisons: Iterable[Any] = (),
    validation_matrix: Any | None = None,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = render_shadow_ladder_report_payload(
        ladders,
        audit=audit,
        comparisons=comparisons,
        validation_matrix=validation_matrix,
    )
    json_path = output_dir / f"{stem}.json"
    markdown_path = output_dir / f"{stem}.md"
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    markdown_path.write_text(render_shadow_ladder_markdown(payload), encoding="utf-8")
    return json_path, markdown_path


def render_shadow_ladder_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Shadow LP Ladder Report",
        "",
        payload["summary"]["no_absolute_yield_statement"],
        "",
        "## Final Objectives",
        "",
        "| target | shadow objective | reference objective | abs diff | rel diff | alignment |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for target_id, row in payload["final_objectives"].items():
        lines.append(
            f"| {target_id} | {_fmt(row['shadow_objective'])} | {_fmt(row['reference_objective'])} | "
            f"{_fmt(row['objective_abs_diff'])} | {_fmt(row['objective_rel_diff'])} | {row['alignment_status']} |"
        )
    lines.extend(
        [
            "",
            "## Constraint Counts",
            "",
            "| target | layer | constraints | equalities | inequalities |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for target_id, layers in payload["constraint_counts"].items():
        for layer_id, row in layers.items():
            lines.append(
                f"| {target_id} | {layer_id} | {row['constraint_count']} | "
                f"{row['eq_constraint_count']} | {row['ub_constraint_count']} |"
            )
    lines.extend(
        [
            "",
            "## Skipped Layers",
            "",
        ]
    )
    for target_id, skipped in payload["skipped_layers"].items():
        for layer_id, reason in skipped.items():
            lines.append(f"- {target_id} `{layer_id}`: {reason.get(layer_id, reason)}")
    lines.extend(
        [
            "",
            "## Backend Metadata",
            "",
            f"- Default large-model backend: {payload['summary']['default_large_model_backend']}",
            f"- Production default solver mode: {payload['summary']['production_default_solver_mode']}",
            f"- No-hardcode audit passed: {payload['no_hardcode_audit']['passed']}",
            "",
            "## Compare Mode Summary",
            "",
        ]
    )
    comparisons = payload.get("compare_mode_summary") or ()
    if comparisons:
        lines.extend(
            [
                "| target | growth rate | rel diff | constraint diff | within tolerance |",
                "|---|---:|---:|---:|---|",
            ]
        )
        for row in comparisons:
            lines.append(
                f"| {row['target_id']} | {_fmt(row['growth_rate'])} | "
                f"{_fmt(row['objective_rel_diff'])} | {_fmt(row['constraint_count_diff'])} | "
                f"{row['within_tolerance']} |"
            )
    else:
        lines.append("- No compare-mode rows provided.")
    lines.extend(
        [
            "",
            "## Validation Matrix",
            "",
        ]
    )
    matrix = payload.get("validation_matrix")
    if matrix:
        lines.extend(
            [
                "| target | growth rate | reference status | shadow status | rel diff | alignment |",
                "|---|---:|---|---|---:|---|",
            ]
        )
        for row in matrix.get("cases", ()):
            lines.append(
                f"| {row['target_id']} | {_fmt(row['growth_rate'])} | {row['reference_status']} | "
                f"{row['shadow_status']} | {_fmt(row['objective_rel_diff'])} | {row['alignment_status']} |"
            )
    else:
        lines.append("- No validation matrix provided.")
    lines.extend(
        [
            "",
            "## Status Mismatches",
            "",
        ]
    )
    mismatches = payload.get("status_mismatches") or ()
    if mismatches:
        for row in mismatches:
            target_id = row.get("target_id", "unknown")
            growth_rate = row.get("growth_rate", "")
            reference_status = row.get("reference_status_category", row.get("reference_status"))
            shadow_status = row.get("shadow_status_category", row.get("shadow_status"))
            lines.append(f"- {target_id} at {_fmt(growth_rate)}: reference={reference_status}, shadow={shadow_status}")
    else:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Warnings",
            "",
        ]
    )
    for target_id, warnings in payload["warnings"].items():
        for warning in warnings:
            lines.append(f"- {target_id}: {warning}")
    return "\n".join(lines) + "\n"


def run_shadow_hardcode_audit(root: Path | None = None) -> ShadowHardcodeAuditResult:
    repo_root = root or Path(__file__).resolve().parents[5]
    shadow_dir = repo_root / "python_pichia" / "src" / "pcsec_pichia" / "analysis" / "shadow_lp"
    validation_only = ("validation.py", "comparison.py")
    audit_only = ("reports.py",)
    production_files = tuple(
        sorted(
            path
            for path in shadow_dir.glob("*.py")
            if path.name not in (*validation_only, *audit_only) and path.name != "__pycache__"
        )
    )
    production_text = "\n".join(path.read_text(encoding="utf-8") for path in production_files)
    warnings: list[str] = []
    objective_literals = (
        "0.003285010027" + "0232106",
        "0.006572021526" + "431409",
    )
    has_no_objective_literals = not any(literal in production_text for literal in objective_literals)
    if not has_no_objective_literals:
        warnings.append("Production shadow path contains reference objective literals.")
    has_no_magic_factor = "target_factor =" not in production_text and "magic factor" not in production_text.lower()
    if not has_no_magic_factor:
        warnings.append("Production shadow path contains a magic factor marker.")
    has_no_forbidden_solver = (
        "solve_pcsec_maximize" not in production_text
        and "solve_secretion_capacity" not in production_text
    )
    if not has_no_forbidden_solver:
        warnings.append("Production shadow solver path calls a forbidden reference solver.")
    validation_text = (shadow_dir / "validation.py").read_text(encoding="utf-8") if (shadow_dir / "validation.py").exists() else ""
    reference_after_solve = "validate_shadow_ladder_against_reference" in validation_text and "solve_secretion_capacity" in validation_text
    return ShadowHardcodeAuditResult(
        production_shadow_path_has_no_reference_objective_literals=has_no_objective_literals,
        production_shadow_path_has_no_magic_factor_literal=has_no_magic_factor,
        production_shadow_path_has_no_forbidden_reference_solver_call=has_no_forbidden_solver,
        hlf_opn_use_shared_builders="build_shadow_constraint_blocks" in production_text,
        reference_values_only_used_after_solve_for_comparison=reference_after_solve,
        validation_only_files=validation_only,
        scanned_files=tuple(str(path.relative_to(repo_root)) for path in production_files),
        warnings=tuple(warnings),
    )


def _objective_changes(ladder: ShadowLadderResult) -> tuple[dict[str, float | None | str], ...]:
    rows: list[dict[str, float | None | str]] = []
    previous: float | None = None
    for layer in ladder.layers:
        delta = None if previous is None or layer.objective is None else layer.objective - previous
        rows.append({"layer_id": layer.layer_id, "objective": layer.objective, "delta_from_previous": delta})
        if layer.objective is not None:
            previous = layer.objective
    return tuple(rows)


def _to_dict(item: Any) -> dict[str, Any]:
    if hasattr(item, "to_dict"):
        return item.to_dict()
    if isinstance(item, Mapping):
        return dict(item)
    raise TypeError(f"Report payload item is not serializable: {type(item)!r}")


def _fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)
