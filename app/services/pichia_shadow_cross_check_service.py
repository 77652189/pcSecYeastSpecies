from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.analysis.shadow_lp import (
    CROSS_CHECK_MANIFEST_FILENAME,
    ShadowLpCrossCheckRequest,
    run_shadow_lp_cross_check,
)
from pcsec_pichia.core.paths import ProjectPaths


SHADOW_CROSS_CHECK_RUNS_DIR = (
    ProjectPaths.discover(Path(__file__)).repo_root / "local_runs" / "shadow_lp_cross_check"
)


def run_pichia_shadow_cross_check(
    *,
    target_id: str,
    screen_run_id: str = "",
    saved_result_path: Path | str | None = None,
    output_dir: Path | str | None = None,
    relative_tolerance: float = 1e-4,
) -> dict[str, Any]:
    """Facade for the python_pichia Shadow LP cross-check engine."""

    resolved_output_dir = _resolve_output_dir(output_dir, target_id=target_id)
    request = ShadowLpCrossCheckRequest(
        target_id=target_id,
        screen_run_id=screen_run_id,
        saved_result_path="" if saved_result_path is None else str(saved_result_path),
        relative_tolerance=float(relative_tolerance),
    )
    outputs = run_shadow_lp_cross_check(request, resolved_output_dir)
    result = outputs.result
    return {
        "submitted": True,
        "status": result.manifest_status,
        "target_id": result.target_id,
        "screen_run_id": result.screen_run_id,
        "output_dir": str(resolved_output_dir),
        "manifest_path": str(outputs.manifest_path),
        "summary_tsv_path": str(outputs.summary_tsv_path),
        "report_path": str(outputs.report_path),
        "diff_path": str(outputs.diff_path),
        "within_tolerance": result.within_tolerance,
        "relative_diff": result.relative_diff,
        "warnings": list(result.warnings),
    }


def load_pichia_shadow_cross_check_manifest(path: Path | str) -> dict[str, Any]:
    resolved = Path(path)
    manifest_path = resolved / CROSS_CHECK_MANIFEST_FILENAME if resolved.is_dir() else resolved
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def _resolve_output_dir(output_dir: Path | str | None, *, target_id: str) -> Path:
    base_dir = SHADOW_CROSS_CHECK_RUNS_DIR.resolve()
    if output_dir is not None:
        candidate = Path(output_dir)
        if not candidate.is_absolute():
            candidate = ProjectPaths.discover(Path(__file__)).repo_root / candidate
        resolved = candidate.resolve()
        if resolved != base_dir and base_dir not in resolved.parents:
            raise ValueError("output_dir must stay under local_runs/shadow_lp_cross_check.")
        return resolved
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    safe_target = "".join(character if character.isalnum() or character in {"-", "_"} else "_" for character in target_id)
    return base_dir / f"{stamp}_{safe_target or 'target'}"


__all__ = [
    "SHADOW_CROSS_CHECK_RUNS_DIR",
    "load_pichia_shadow_cross_check_manifest",
    "run_pichia_shadow_cross_check",
]
