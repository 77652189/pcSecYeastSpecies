"""Submits and polls the genome-wide KO/OE tradeoff screen.

Launched as a detached OS subprocess rather than a Streamlit-hosted daemon
thread: this job runs for hours, and the existing background-task pattern
(app/services/pichia_background_tasks.py) uses an in-process thread that
would be killed by a Streamlit restart/auto-reload, silently losing however
much of the run had completed. A detached subprocess survives that.
"""

from __future__ import annotations

import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from pcsec_pichia.core.paths import ProjectPaths

from app.services.genome_wide_screen_registry import RunInfo, list_active_runs

_SCRIPT_RELATIVE_PATH = Path("python_pichia") / "tools" / "run_genome_wide_ko_oe_screen_parallel.py"
DEFAULT_WORKERS = 6

# Detached process creation flags (Windows only; ignored elsewhere).
_DETACHED_FLAGS = 0x00000008 | 0x00000200  # DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP


@dataclass(frozen=True)
class SubmitResult:
    run_name: str
    log_path: Path
    status_path: Path


def check_for_conflicts(paths: ProjectPaths | None = None) -> list[RunInfo]:
    """Return currently active runs; the UI should confirm with the user before submitting if this is non-empty."""
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    return list_active_runs(resolved_paths)


def submit_screen(
    targets: list[str],
    mode: str = "fast",
    workers: int = DEFAULT_WORKERS,
    gene_limit: int | None = None,
    paths: ProjectPaths | None = None,
    scope: str = "gene",
) -> SubmitResult:
    """Launch the genome-wide screen as a detached subprocess. Does not check for conflicts;
    call check_for_conflicts() first and let the caller decide (force/queue/cancel).

    scope="gene" (default) screens all ~1025 model genes (hour-scale); scope="catalog"
    screens the ~30 unique reactions named in the curated SECRETION_GENE_CATALOG literature
    shortlist (minute-scale) - gene_limit is ignored for catalog scope, it has no equivalent.
    """
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    run_name = f"ui_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"
    out_dir = resolved_paths.local_runs_dir / run_name
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run.log"
    status_path = out_dir / "status.json"

    script_path = resolved_paths.repo_root / _SCRIPT_RELATIVE_PATH
    command = [
        sys.executable,
        str(script_path),
        "--targets", ",".join(targets),
        "--mode", mode,
        "--workers", str(workers),
        "--run-name", run_name,
        "--scope", scope,
    ]
    if gene_limit is not None and scope == "gene":
        command.extend(["--limit", str(gene_limit)])

    with log_path.open("w", encoding="utf-8") as log_file:
        popen_kwargs: dict[str, object] = {
            "cwd": str(resolved_paths.repo_root / "python_pichia"),
            "stdout": log_file,
            "stderr": subprocess.STDOUT,
        }
        if sys.platform == "win32":
            popen_kwargs["creationflags"] = _DETACHED_FLAGS
        else:
            popen_kwargs["start_new_session"] = True
        subprocess.Popen(command, **popen_kwargs)  # noqa: S603 - fixed argv, no shell, no user-controlled binary path

    return SubmitResult(run_name=run_name, log_path=log_path, status_path=status_path)


def poll_screen(status_path: Path) -> dict[str, object]:
    """Read the current status.json for a submitted run. Returns {"status": "lost"} if not found."""
    import json

    if not status_path.exists():
        return {"status": "lost"}
    try:
        return json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"status": "lost"}


__all__ = [
    "DEFAULT_WORKERS",
    "SubmitResult",
    "check_for_conflicts",
    "poll_screen",
    "submit_screen",
]
