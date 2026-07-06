"""Tracks active/queued genome-wide KO/OE screen runs so the UI can detect
collisions before launching a new (multi-hour) run.

A run is considered "active" if its status.json exists, its status is
"starting" or "running", and its heartbeat was updated recently. A stale
heartbeat (process likely crashed without writing an error status) is not
treated as active, so a stuck entry cannot permanently block new runs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from pcsec_pichia.core.paths import ProjectPaths

HEARTBEAT_STALE_SECONDS = 30 * 60
REGISTRY_DIRNAME = "genome_wide_ko_oe_screen"


@dataclass(frozen=True)
class RunInfo:
    run_name: str
    status: str
    done: int
    total: int
    targets: tuple[str, ...]
    mode: str
    pid: int | None
    updated_at: str | None
    is_stale: bool
    scope: str = "gene"
    csv_path: str | None = None

    @property
    def progress_label(self) -> str:
        if self.total <= 0:
            return self.status
        return f"{self.done}/{self.total} ({self.status})"


def registry_dir(paths: ProjectPaths) -> Path:
    directory = paths.local_runs_dir / REGISTRY_DIRNAME
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def list_runs(paths: ProjectPaths) -> list[RunInfo]:
    """List every run this service has ever launched, newest first.

    Also recognizes runs produced before status.json existed (e.g. the
    2026-07-02 overnight full-genome run): any local_runs/<name>/ with a
    gene_tradeoff_rows.csv but no status.json gets a status.json backfilled
    from the CSV so it shows up here too, instead of being silently invisible.
    """
    _backfill_missing_status_files(paths)
    runs: list[RunInfo] = []
    for status_path in sorted(
        paths.local_runs_dir.glob("*/status.json"), key=lambda path: path.stat().st_mtime, reverse=True
    ):
        info = _read_run_info(status_path)
        if info is not None:
            runs.append(info)
    return runs


def _backfill_missing_status_files(paths: ProjectPaths) -> None:
    for csv_path in paths.local_runs_dir.glob("*/gene_tradeoff_rows.csv"):
        status_path = csv_path.with_name("status.json")
        if status_path.exists():
            continue
        try:
            _write_backfilled_status(csv_path, status_path)
        except (OSError, ValueError):
            continue  # unreadable/partial CSV from an interrupted legacy run; skip rather than fail the whole list


def _write_backfilled_status(csv_path: Path, status_path: Path) -> None:
    import csv as csv_module

    targets: set[str] = set()
    row_count = 0
    with csv_path.open(encoding="utf-8") as csv_file:
        for row in csv_module.DictReader(csv_file):
            row_count += 1
            target_id = row.get("target_id")
            if target_id:
                targets.add(target_id)
    task_count = row_count // 2  # two rows (KO + OE) per (target, gene) task
    updated_at = datetime.fromtimestamp(csv_path.stat().st_mtime).isoformat()
    payload = {
        "status": "done",
        "done": task_count,
        "total": task_count,
        "targets": sorted(targets),
        "mode": "unknown (backfilled from a run predating status tracking)",
        "pid": None,
        "updated_at": updated_at,
        "backfilled": True,
        "scope": "gene",  # catalog scope did not exist yet when any backfill-eligible run was produced
        "csv_path": str(csv_path),
    }
    status_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def list_active_runs(paths: ProjectPaths) -> list[RunInfo]:
    """List runs that are currently in progress (not stale, not finished)."""
    return [run for run in list_runs(paths) if run.status in {"starting", "running"} and not run.is_stale]


def _read_run_info(status_path: Path) -> RunInfo | None:
    try:
        payload = json.loads(status_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return None
    updated_at = payload.get("updated_at")
    is_stale = _is_stale(updated_at)
    return RunInfo(
        run_name=status_path.parent.name,
        status=str(payload.get("status", "unknown")),
        done=int(payload.get("done") or 0),
        total=int(payload.get("total") or 0),
        targets=tuple(payload.get("targets") or ()),
        mode=str(payload.get("mode") or ""),
        pid=payload.get("pid"),
        updated_at=updated_at,
        is_stale=is_stale,
        scope=str(payload.get("scope") or "gene"),
        csv_path=payload.get("csv_path"),
    )


def _is_stale(updated_at: str | None) -> bool:
    if not updated_at:
        return False
    try:
        age_seconds = (datetime.now() - datetime.fromisoformat(updated_at)).total_seconds()
    except ValueError:
        return False
    return age_seconds > HEARTBEAT_STALE_SECONDS


__all__ = [
    "HEARTBEAT_STALE_SECONDS",
    "RunInfo",
    "list_active_runs",
    "list_runs",
    "registry_dir",
]
