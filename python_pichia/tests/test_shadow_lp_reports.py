from __future__ import annotations

import json
import shutil
from pathlib import Path

from pcsec_pichia.analysis.shadow_lp import (
    ShadowHardcodeAuditResult,
    ShadowLadderLayerResult,
    ShadowLadderResult,
    render_shadow_ladder_markdown,
    render_shadow_ladder_report_payload,
    write_shadow_ladder_report,
)


def test_shadow_ladder_report_payload_contains_required_sections() -> None:
    ladder = _synthetic_ladder()
    audit = _passing_audit()

    payload = render_shadow_ladder_report_payload((ladder,), audit=audit)

    assert payload["summary"]["default_large_model_backend"] == "ScipyHighsBackend"
    assert payload["summary"]["production_default_solver_mode"] == "reference"
    assert payload["summary"]["canonical_final_layer"] == "mitochondrial"
    assert payload["final_objectives"]["hLF"]["objective_rel_diff"] == 0.0
    assert "constraint_counts" in payload
    assert "backend_metadata" in payload
    assert "ribosome_translation" in payload["skipped_layers"]["hLF"]
    assert "compare_mode_summary" in payload
    assert "validation_matrix" in payload
    assert "status_mismatches" in payload
    assert payload["no_hardcode_audit"]["passed"] is True


def test_shadow_ladder_json_and_markdown_reports_are_written_under_local_runs() -> None:
    output_dir = Path("local_runs") / "shadow_lp_report_tests"
    if output_dir.exists():
        shutil.rmtree(output_dir)

    json_path, markdown_path = write_shadow_ladder_report((_synthetic_ladder(),), output_dir, audit=_passing_audit())

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    markdown = markdown_path.read_text(encoding="utf-8")
    assert payload["final_objectives"]["hLF"]["shadow_objective"] == 0.1
    assert "Shadow LP Ladder Report" in markdown
    assert "constraint" in markdown.lower()
    assert "backend" in markdown.lower()
    assert "mg/L" not in markdown
    assert "absolute fermentation titer" in markdown
    assert "Compare Mode Summary" in markdown
    assert "Validation Matrix" in markdown
    assert "Status Mismatches" in markdown


def test_shadow_ladder_report_payload_accepts_compare_and_validation_matrix_sections() -> None:
    comparison = {
        "target_id": "hLF",
        "growth_rate": 0.1,
        "objective_rel_diff": 0.0,
        "constraint_count_diff": 0,
        "within_tolerance": True,
        "reference_status_category": "optimal",
        "shadow_status_category": "optimal",
    }
    matrix = {
        "cases": (
            {
                "target_id": "hLF",
                "growth_rate": 0.1,
                "reference_status": "optimal",
                "shadow_status": "optimal",
                "objective_rel_diff": 0.0,
                "alignment_status": "aligned",
            },
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "growth_rate": 0.15,
                "reference_status": "optimal",
                "shadow_status": "exception",
                "objective_rel_diff": None,
                "alignment_status": "review_required",
            },
        ),
        "all_required_defaults_aligned": True,
    }

    payload = render_shadow_ladder_report_payload(
        (_synthetic_ladder(),),
        audit=_passing_audit(),
        comparisons=(comparison,),
        validation_matrix=matrix,
    )
    markdown = render_shadow_ladder_markdown(payload)

    assert payload["compare_mode_summary"][0]["within_tolerance"] is True
    assert payload["validation_matrix"]["all_required_defaults_aligned"] is True
    assert len(payload["status_mismatches"]) == 1
    assert "reference=optimal, shadow=exception" in markdown


def _synthetic_ladder() -> ShadowLadderResult:
    layers = (
        ShadowLadderLayerResult(
            target_id="hLF",
            layer_id="target_extension",
            success=True,
            status="implemented_prerequisite",
            message="built",
            objective=None,
            key_fluxes={},
            variable_count=2,
            constraint_count=1,
            eq_constraint_count=1,
            ub_constraint_count=0,
            backend_metadata={"backend": "scipy-highs"},
            timings={"total_seconds": 0.0},
            enabled_layers=("target_extension",),
            skipped_layers={},
        ),
        ShadowLadderLayerResult(
            target_id="hLF",
            layer_id="ribosome_translation",
            success=False,
            status="skipped",
            message="disabled",
            objective=None,
            key_fluxes={},
            variable_count=2,
            constraint_count=1,
            eq_constraint_count=1,
            ub_constraint_count=0,
            backend_metadata={"backend": "scipy-highs"},
            timings={"total_seconds": 0.0},
            enabled_layers=(),
            skipped_layers={"ribosome_translation": "disabled by default"},
            warnings=("disabled by default",),
        ),
        ShadowLadderLayerResult(
            target_id="hLF",
            layer_id="mitochondrial",
            success=True,
            status="0",
            message="optimal",
            objective=0.1,
            key_fluxes={"BIOMASS": 0.1, "hLF exchange": 0.1},
            variable_count=2,
            constraint_count=3,
            eq_constraint_count=2,
            ub_constraint_count=1,
            backend_metadata={"backend": "scipy-highs", "constraint_count": 3},
            timings={"solve_seconds": 0.0},
            enabled_layers=("mitochondrial",),
            skipped_layers={},
        ),
    )
    return ShadowLadderResult(
        target_id="hLF",
        exchange_reaction_id="hLF exchange",
        backend_name="scipy-highs",
        layers=layers,
        reference_validation={
            "reference_objective": 0.1,
            "shadow_objective": 0.1,
            "objective_abs_diff": 0.0,
            "objective_rel_diff": 0.0,
            "constraint_count_diff": 0,
            "final_alignment_status": "aligned",
        },
    )


def _passing_audit() -> ShadowHardcodeAuditResult:
    return ShadowHardcodeAuditResult(
        production_shadow_path_has_no_reference_objective_literals=True,
        production_shadow_path_has_no_magic_factor_literal=True,
        production_shadow_path_has_no_forbidden_reference_solver_call=True,
        hlf_opn_use_shared_builders=True,
        reference_values_only_used_after_solve_for_comparison=True,
        validation_only_files=("validation.py",),
        scanned_files=("python_pichia/src/pcsec_pichia/analysis/shadow_lp/ladder.py",),
        warnings=(),
    )
