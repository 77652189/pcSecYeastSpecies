from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from pcsec_pichia.adapters.cobrapy_shadow import compare_shadow_fba, solve_cobrapy_shadow_fba
from pcsec_pichia.adapters.lp_solver import LinearProgrammingResult, ScipyHiGHSSolver
from pcsec_pichia.adapters.mat_loader import MatStructLoader
from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.core.pichia_model import PichiaModel, default_key_reactions, summary_to_dict


DEFAULT_OBJECTIVE_REACTIONS: tuple[str, ...] = (
    "BIOMASS",
    "Ex_glc_D",
    "Ex_glyc",
    "Ex_meoh",
    "Ex_o2",
)

DEFAULT_OUTPUT_SUBDIR = "cobrapy_shadow_baseline"


@dataclass(frozen=True)
class CobraPyShadowBaselineCase:
    case_id: str
    objective_reaction: str
    sense: str
    shadow_available: bool
    current_success: bool
    shadow_success: bool
    current_objective_value: float | None
    shadow_objective_value: float | None
    objective_abs_diff: float | None
    objective_rel_diff: float | None
    within_tolerance: bool
    key_flux_diffs: dict[str, float | None]
    model_summary: dict[str, object]
    warnings: tuple[str, ...] = ()
    current_status: str = ""
    shadow_status: str = ""
    comparison_status: str = ""

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CobraPyShadowBaselineRun:
    result_status: str
    output_dir: str
    summary_json_path: str | None
    report_markdown_path: str | None
    cases: tuple[CobraPyShadowBaselineCase, ...]
    model_summary: dict[str, object]
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "result_status": self.result_status,
            "output_dir": self.output_dir,
            "summary_json_path": self.summary_json_path,
            "report_markdown_path": self.report_markdown_path,
            "model_summary": dict(self.model_summary),
            "warnings": list(self.warnings),
            "cases": [case.to_dict() for case in self.cases],
        }


def run_cobrapy_shadow_baseline(
    paths: ProjectPaths | None = None,
    objective_reactions: tuple[str, ...] | list[str] = DEFAULT_OBJECTIVE_REACTIONS,
    output_dir: Path | None = None,
    sense: str = "maximize",
    write_artifacts: bool = True,
) -> CobraPyShadowBaselineRun:
    paths = paths or ProjectPaths.discover()
    resolved_output_dir = _resolve_output_dir(paths, output_dir)
    model = MatStructLoader(paths).load_pcsec_pichia_model()
    model_summary = summary_to_dict(model.summary())
    key_reactions = _existing_reactions(model, default_key_reactions())
    cases: list[CobraPyShadowBaselineCase] = []
    warnings = [
        "COBRApy shadow baseline is an opt-in developer parity harness.",
        "This run converts only base GEM stoichiometry, bounds, and objective reactions.",
        "This run does not convert pcSec protein/secretion constraints, KO/OE planning, phenotype evidence, recommendation tiers, mg/L output, or experiment success probabilities.",
    ]

    for objective_reaction in _existing_reactions(model, objective_reactions):
        case_key_reactions = tuple(dict.fromkeys((*key_reactions, objective_reaction)))
        cases.append(
            _run_case(
                model=model,
                model_summary=model_summary,
                objective_reaction=objective_reaction,
                sense=sense,
                key_reactions=case_key_reactions,
            )
        )

    if not cases:
        warnings.append("No configured objective reactions were present in the pcSecPichia model.")

    result_status = _result_status(cases)
    run = CobraPyShadowBaselineRun(
        result_status=result_status,
        output_dir=str(resolved_output_dir),
        summary_json_path=None,
        report_markdown_path=None,
        cases=tuple(cases),
        model_summary=model_summary,
        warnings=tuple(warnings),
    )
    if not write_artifacts:
        return run

    resolved_output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = resolved_output_dir / "cobrapy_shadow_baseline_summary.json"
    report_path = resolved_output_dir / "cobrapy_shadow_baseline_report.md"
    payload = run.to_dict() | {
        "summary_json_path": str(summary_path),
        "report_markdown_path": str(report_path),
    }
    summary_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    report_path.write_text(_render_markdown(payload), encoding="utf-8")
    return CobraPyShadowBaselineRun(
        result_status=result_status,
        output_dir=str(resolved_output_dir),
        summary_json_path=str(summary_path),
        report_markdown_path=str(report_path),
        cases=tuple(cases),
        model_summary=model_summary,
        warnings=tuple(warnings),
    )


def _run_case(
    model: PichiaModel,
    model_summary: dict[str, object],
    objective_reaction: str,
    sense: str,
    key_reactions: tuple[str, ...],
) -> CobraPyShadowBaselineCase:
    case_warnings: list[str] = []
    current_result = _solve_current(model, objective_reaction, sense=sense, key_reactions=key_reactions)
    if current_result is None:
        case_warnings.append("Current SciPy/HiGHS baseline solve raised an exception.")

    shadow_result = solve_cobrapy_shadow_fba(
        model,
        objective_reaction,
        sense=sense,  # type: ignore[arg-type]
        key_reactions=key_reactions,
    )
    if not shadow_result.available:
        case_warnings.append("COBRApy shadow FBA is unavailable; optional dependency is not importable.")
    if shadow_result.available and not shadow_result.success:
        case_warnings.append("COBRApy shadow FBA did not solve successfully.")

    comparison = compare_shadow_fba(current_result, shadow_result, key_reactions=key_reactions) if current_result else None
    return CobraPyShadowBaselineCase(
        case_id=f"{objective_reaction}_{sense}",
        objective_reaction=objective_reaction,
        sense=sense,
        shadow_available=shadow_result.available,
        current_success=bool(current_result.success) if current_result else False,
        shadow_success=shadow_result.success,
        current_objective_value=current_result.objective_value if current_result else None,
        shadow_objective_value=shadow_result.objective_value,
        objective_abs_diff=comparison.objective_abs_diff if comparison else None,
        objective_rel_diff=comparison.objective_rel_diff if comparison else None,
        within_tolerance=bool(comparison.within_tolerance) if comparison else False,
        key_flux_diffs=dict(comparison.key_flux_diffs) if comparison else {},
        model_summary=model_summary,
        warnings=tuple(case_warnings),
        current_status=str(current_result.status) if current_result else "exception",
        shadow_status=shadow_result.status,
        comparison_status=comparison.status if comparison else "not_compared",
    )


def _solve_current(
    model: PichiaModel,
    objective_reaction: str,
    sense: str,
    key_reactions: tuple[str, ...],
) -> LinearProgrammingResult | None:
    try:
        return ScipyHiGHSSolver().solve(
            model,
            objective_reaction,
            sense=sense,  # type: ignore[arg-type]
            key_reactions=key_reactions,
        )
    except Exception:
        return None


def _existing_reactions(model: PichiaModel, reaction_ids: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    reaction_index = model.reaction_index
    return tuple(reaction_id for reaction_id in reaction_ids if reaction_id in reaction_index)


def _resolve_output_dir(paths: ProjectPaths, output_dir: Path | None) -> Path:
    local_runs_dir = paths.local_runs_dir.resolve()
    resolved = (output_dir or (local_runs_dir / DEFAULT_OUTPUT_SUBDIR)).resolve()
    if not _is_relative_to(resolved, local_runs_dir):
        raise ValueError(f"COBRApy shadow baseline artifacts must be written under {local_runs_dir}.")
    return resolved


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _result_status(cases: list[CobraPyShadowBaselineCase]) -> str:
    if not cases:
        return "skipped_no_objective_reactions"
    if all(not case.shadow_available for case in cases):
        return "completed_shadow_unavailable"
    if any(not case.current_success or not case.shadow_success for case in cases):
        return "completed_with_differences_or_failures"
    if any(not case.shadow_available for case in cases):
        return "completed_with_differences_or_failures"
    if all(case.within_tolerance for case in cases):
        return "completed_within_tolerance"
    return "completed_with_differences_or_failures"


def _render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# COBRApy Shadow Baseline",
        "",
        f"Result status: `{payload['result_status']}`",
        "",
        "This is an opt-in developer parity harness for base GEM FBA only.",
        "It does not migrate pcSec protein/secretion constraints, KO/OE planning, recommendation tiers, mg/L output, or experiment success probabilities.",
        "",
        "## Cases",
        "",
        "| case_id | objective | shadow_available | current_success | shadow_success | abs_diff | rel_diff | within_tolerance |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for case in payload["cases"]:
        lines.append(
            "| {case_id} | {objective_reaction} | {shadow_available} | {current_success} | {shadow_success} | {objective_abs_diff} | {objective_rel_diff} | {within_tolerance} |".format(
                **case
            )
        )
    lines.extend(
        [
            "",
            "## Warnings",
            "",
            *[f"- {warning}" for warning in payload["warnings"]],
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run opt-in COBRApy shadow FBA parity checks for the real pcSecPichia model.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory under local_runs/. Defaults to local_runs/cobrapy_shadow_baseline.",
    )
    args = parser.parse_args(argv)
    run = run_cobrapy_shadow_baseline(output_dir=args.output_dir)
    print(json.dumps(run.to_dict(), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "CobraPyShadowBaselineCase",
    "CobraPyShadowBaselineRun",
    "DEFAULT_OBJECTIVE_REACTIONS",
    "DEFAULT_OUTPUT_SUBDIR",
    "run_cobrapy_shadow_baseline",
]
