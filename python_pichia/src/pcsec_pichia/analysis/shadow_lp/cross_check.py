from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Mapping

from pcsec_pichia.analysis.shadow_lp.constraint_spec import ShadowConstraintConfig
from pcsec_pichia.analysis.shadow_lp.ladder import ShadowLadderResult, run_shadow_ladder
from pcsec_pichia.analysis.shadow_lp.validation import (
    ReferenceValidationResult,
    validate_shadow_ladder_against_reference,
)


CROSS_CHECK_MANIFEST_FILENAME = "cross_check_manifest.json"
CROSS_CHECK_SUMMARY_TSV_FILENAME = "cross_check_summary.tsv"
CROSS_CHECK_REPORT_FILENAME = "cross_check_report.md"
CROSS_CHECK_DIFF_FILENAME = "reference_vs_shadow_diff.json"
NO_EXPERIMENTAL_CLAIM_STATEMENT = (
    "Shadow LP cross-checks compare model solver paths only; they do not predict absolute secretion titer, "
    "absolute yield, or experimental success rate."
)

ShadowLadderRunner = Callable[..., Any]
ReferenceValidator = Callable[..., Any]


@dataclass(frozen=True)
class ShadowLpCrossCheckRequest:
    target_id: str
    screen_run_id: str = ""
    saved_result_path: str = ""
    relative_tolerance: float = 1e-4

    def validate(self) -> None:
        if not self.target_id.strip():
            raise ValueError("target_id must be non-empty.")


@dataclass(frozen=True)
class SavedShadowCrossCheckContext:
    target_id: str = ""
    screen_run_id: str = ""
    reference_capacity: float | None = None
    source_path: str = ""
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ShadowLpCrossCheckResult:
    target_id: str
    screen_run_id: str
    reference_capacity: float | None
    shadow_capacity: float | None
    absolute_diff: float | None
    relative_diff: float | None
    within_tolerance: bool
    constraint_layer: str
    backend: str
    solver_status: str
    alignment_status: str
    manifest_status: str
    saved_result_path: str = ""
    saved_reference_capacity: float | None = None
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ShadowLpCrossCheckOutputs:
    manifest_path: Path
    summary_tsv_path: Path
    report_path: Path
    diff_path: Path
    result: ShadowLpCrossCheckResult


def load_shadow_cross_check_saved_result(path: Path) -> SavedShadowCrossCheckContext:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("saved result must be a JSON object.")
    target_id = str(payload.get("target_id", ""))
    screen_run_id = str(payload.get("screen_run_id", payload.get("run_id", "")))
    reference_capacity = _optional_float(
        payload.get("reference_capacity", payload.get("reference_objective", payload.get("secretion_capacity")))
    )
    warnings = _tuple_from_json(payload.get("warnings", ()))
    return SavedShadowCrossCheckContext(
        target_id=target_id,
        screen_run_id=screen_run_id,
        reference_capacity=reference_capacity,
        source_path=str(path),
        warnings=warnings,
    )


def run_shadow_lp_cross_check(
    request: ShadowLpCrossCheckRequest,
    output_dir: Path,
    *,
    root: Path | None = None,
    config: ShadowConstraintConfig | None = None,
    ladder_runner: ShadowLadderRunner | None = None,
    reference_validator: ReferenceValidator | None = None,
) -> ShadowLpCrossCheckOutputs:
    request.validate()
    saved_context = (
        load_shadow_cross_check_saved_result(Path(request.saved_result_path))
        if request.saved_result_path
        else SavedShadowCrossCheckContext()
    )
    target_id = saved_context.target_id or request.target_id
    screen_run_id = request.screen_run_id or saved_context.screen_run_id
    context_warnings = list(saved_context.warnings)
    if saved_context.target_id and saved_context.target_id != request.target_id:
        context_warnings.append("saved_result_target_id_overrode_request_target_id")

    resolved_runner = ladder_runner or run_shadow_ladder
    resolved_validator = reference_validator or validate_shadow_ladder_against_reference
    ladder = resolved_runner(target_id, root=root, config=config)
    validation = resolved_validator(
        ladder,
        root=root,
        config=config,
        relative_tolerance=request.relative_tolerance,
    )
    result = _cross_check_result(
        request=request,
        ladder=ladder,
        validation=validation,
        screen_run_id=screen_run_id,
        saved_context=saved_context,
        context_warnings=tuple(context_warnings),
    )
    return write_shadow_lp_cross_check_outputs(result, output_dir, ladder=ladder, validation=validation)


def write_shadow_lp_cross_check_outputs(
    result: ShadowLpCrossCheckResult,
    output_dir: Path,
    *,
    ladder: Any | None = None,
    validation: Any | None = None,
) -> ShadowLpCrossCheckOutputs:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / CROSS_CHECK_MANIFEST_FILENAME
    summary_tsv_path = output_dir / CROSS_CHECK_SUMMARY_TSV_FILENAME
    report_path = output_dir / CROSS_CHECK_REPORT_FILENAME
    diff_path = output_dir / CROSS_CHECK_DIFF_FILENAME
    manifest_payload = {
        "result": result.to_dict(),
        "no_experimental_claim_statement": NO_EXPERIMENTAL_CLAIM_STATEMENT,
    }
    manifest_path.write_text(json.dumps(manifest_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_summary_tsv(result, summary_tsv_path)
    report_path.write_text(render_shadow_lp_cross_check_report(result), encoding="utf-8")
    diff_path.write_text(
        json.dumps(
            {
                "result": result.to_dict(),
                "validation": _to_mapping(validation),
                "final_layer": _final_layer_payload(ladder),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return ShadowLpCrossCheckOutputs(
        manifest_path=manifest_path,
        summary_tsv_path=summary_tsv_path,
        report_path=report_path,
        diff_path=diff_path,
        result=result,
    )


def render_shadow_lp_cross_check_report(result: ShadowLpCrossCheckResult) -> str:
    lines = [
        "# Shadow LP Cross-check Report",
        "",
        NO_EXPERIMENTAL_CLAIM_STATEMENT,
        "",
        "## Summary",
        "",
        f"- target_id: {result.target_id}",
        f"- screen_run_id: {result.screen_run_id}",
        f"- alignment_status: {result.alignment_status}",
        f"- within_tolerance: {result.within_tolerance}",
        f"- backend: {result.backend}",
        f"- solver_status: {result.solver_status}",
        f"- constraint_layer: {result.constraint_layer}",
        "",
        "## Capacity Comparison",
        "",
        "| reference capacity | shadow capacity | absolute diff | relative diff |",
        "|---:|---:|---:|---:|",
        (
            f"| {_fmt(result.reference_capacity)} | {_fmt(result.shadow_capacity)} | "
            f"{_fmt(result.absolute_diff)} | {_fmt(result.relative_diff)} |"
        ),
        "",
        "## Boundary",
        "",
        "- Cross-check failure means model or solver-path consistency needs review.",
        "- Cross-check failure is not an experimental infeasibility call.",
        "- This output does not modify KO/OE recommendation tiers.",
        "",
        "## Warnings",
        "",
    ]
    if result.warnings:
        lines.extend(f"- {warning}" for warning in result.warnings)
    else:
        lines.append("- None.")
    return "\n".join(lines) + "\n"


def _cross_check_result(
    *,
    request: ShadowLpCrossCheckRequest,
    ladder: Any,
    validation: Any,
    screen_run_id: str,
    saved_context: SavedShadowCrossCheckContext,
    context_warnings: tuple[str, ...],
) -> ShadowLpCrossCheckResult:
    final_layer = getattr(ladder, "final_layer")
    validation_payload = _to_mapping(validation)
    relative_diff = _optional_float(validation_payload.get("objective_rel_diff"))
    absolute_diff = _optional_float(validation_payload.get("objective_abs_diff"))
    reference_capacity = _optional_float(validation_payload.get("reference_objective"))
    shadow_capacity = _optional_float(validation_payload.get("shadow_objective", getattr(final_layer, "objective", None)))
    alignment_status = str(validation_payload.get("final_alignment_status", "review_required"))
    within_tolerance = alignment_status == "aligned" or (
        relative_diff is not None and relative_diff <= request.relative_tolerance
    )
    warnings = tuple(
        dict.fromkeys(
            (
                *context_warnings,
                *tuple(getattr(ladder, "warnings", ()) or ()),
                *(() if within_tolerance else ("shadow_cross_check_review_required",)),
            )
        )
    )
    return ShadowLpCrossCheckResult(
        target_id=str(getattr(ladder, "target_id", request.target_id)),
        screen_run_id=screen_run_id,
        reference_capacity=reference_capacity,
        shadow_capacity=shadow_capacity,
        absolute_diff=absolute_diff,
        relative_diff=relative_diff,
        within_tolerance=within_tolerance,
        constraint_layer=str(getattr(final_layer, "layer_id", "")),
        backend=str(getattr(ladder, "backend_name", "")),
        solver_status=str(getattr(final_layer, "status", "")),
        alignment_status="aligned" if within_tolerance else "review_required",
        manifest_status="ok" if within_tolerance else "review_required",
        saved_result_path=saved_context.source_path,
        saved_reference_capacity=saved_context.reference_capacity,
        warnings=warnings,
    )


def _write_summary_tsv(result: ShadowLpCrossCheckResult, path: Path) -> None:
    fieldnames = tuple(ShadowLpCrossCheckResult.__dataclass_fields__)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({key: _tsv_value(value) for key, value in result.to_dict().items()})


def _final_layer_payload(ladder: Any | None) -> Mapping[str, Any]:
    if ladder is None:
        return {}
    final_layer = getattr(ladder, "final_layer", None)
    return _to_mapping(final_layer)


def _to_mapping(value: Any) -> Mapping[str, Any]:
    if value is None:
        return {}
    if hasattr(value, "to_dict"):
        return value.to_dict()
    if hasattr(value, "__dataclass_fields__"):
        return asdict(value)
    if isinstance(value, Mapping):
        return dict(value)
    return {}


def _tuple_from_json(value: object) -> tuple[str, ...]:
    if isinstance(value, tuple):
        return tuple(str(item) for item in value)
    if isinstance(value, list):
        return tuple(str(item) for item in value)
    if value in (None, ""):
        return ()
    raise ValueError("warnings must be a JSON list.")


def _optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _tsv_value(value: object) -> str:
    if isinstance(value, tuple):
        return "|".join(str(item) for item in value)
    if isinstance(value, list):
        return "|".join(str(item) for item in value)
    return "" if value is None else str(value)


def _fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


__all__ = [
    "CROSS_CHECK_DIFF_FILENAME",
    "CROSS_CHECK_MANIFEST_FILENAME",
    "CROSS_CHECK_REPORT_FILENAME",
    "CROSS_CHECK_SUMMARY_TSV_FILENAME",
    "NO_EXPERIMENTAL_CLAIM_STATEMENT",
    "SavedShadowCrossCheckContext",
    "ShadowLpCrossCheckOutputs",
    "ShadowLpCrossCheckRequest",
    "ShadowLpCrossCheckResult",
    "load_shadow_cross_check_saved_result",
    "render_shadow_lp_cross_check_report",
    "run_shadow_lp_cross_check",
    "write_shadow_lp_cross_check_outputs",
]
