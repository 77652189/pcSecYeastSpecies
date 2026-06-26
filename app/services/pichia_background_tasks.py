from __future__ import annotations

import json
import threading
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any

from pcsec_pichia.core.paths import ProjectPaths

from app.services.pichia_secretion_schema import SecretionRunRequest


def save_last_result(result_dict: dict[str, object], paths: ProjectPaths) -> None:
    """Save the most recent simulation result to a local cache file."""
    cache_path = _last_result_dir(paths) / "result.json"
    cache_path.write_text(
        json.dumps(result_dict, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def load_last_result(paths: ProjectPaths) -> dict[str, object] | None:
    """Load the cached simulation result, or None if no cache exists."""
    cache_path = _last_result_dir(paths) / "result.json"
    if not cache_path.exists():
        return None
    try:
        data = json.loads(cache_path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else None
    except (json.JSONDecodeError, OSError):
        return None


def clear_last_result(paths: ProjectPaths) -> None:
    """Delete the cached simulation result."""
    cache_path = _last_result_dir(paths) / "result.json"
    if cache_path.exists():
        cache_path.unlink()


def submit_background_simulation(
    request: SecretionRunRequest,
    paths: ProjectPaths | None = None,
) -> tuple[str, Path]:
    """Submit a secretion simulation to run in a background daemon thread."""
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    task_id = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    task_dir = _background_dir(resolved_paths) / task_id
    task_dir.mkdir(parents=True, exist_ok=True)
    status_path = task_dir / "status.json"

    _write_status(status_path, "pending", message="任务已提交，等待执行……")

    thread = threading.Thread(
        target=_run_background_task,
        args=(request, resolved_paths, status_path),
        daemon=True,
    )
    with _BACKGROUND_LOCK:
        thread.start()

    return task_id, status_path


def poll_background_simulation(status_path: Path) -> tuple[str, str, dict[str, Any] | None]:
    """Poll the status of a background simulation task.

    Returns (status, message, result_dict_or_None).
    Status is one of: "pending", "running", "done", "error", "lost".
    """
    if not status_path.exists():
        return "lost", "状态文件不存在", None
    try:
        data = json.loads(status_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return "lost", "无法读取状态文件", None
    status: str = data.get("status", "unknown")
    message: str = data.get("message", "")
    result: dict[str, Any] | None = data.get("result")
    return status, message, result


def status_path_for_background_task(task_id: str, paths: ProjectPaths | None = None) -> Path:
    """Return the status path for a submitted background simulation task."""
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    return _background_dir(resolved_paths) / task_id / "status.json"


def _last_result_dir(paths: ProjectPaths) -> Path:
    directory = paths.local_runs_dir / "streamlit_pichia_runs" / ".last_result"
    directory.mkdir(parents=True, exist_ok=True)
    return directory


_BACKGROUND_DIR: Path | None = None
_BACKGROUND_LOCK = threading.Lock()


def _background_dir(paths: ProjectPaths) -> Path:
    global _BACKGROUND_DIR
    if _BACKGROUND_DIR is None:
        _BACKGROUND_DIR = paths.local_runs_dir / "streamlit_pichia_runs" / ".background_tasks"
        _BACKGROUND_DIR.mkdir(parents=True, exist_ok=True)
    return _BACKGROUND_DIR


def _write_status(path: Path, status: str, *, message: str = "", result: dict[str, Any] | None = None) -> None:
    data: dict[str, Any] = {
        "status": status,
        "message": message,
        "updated_at": datetime.now().isoformat(),
    }
    if result is not None:
        data["result"] = result
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _run_background_task(request: SecretionRunRequest, paths: ProjectPaths, status_path: Path) -> None:
    from app.services.pichia_response_summary_service import response_to_summary
    from app.services.pichia_secretion_service import run_pichia_secretion_draft

    try:
        _write_status(status_path, "running", message="仿真正在运行中（约需 1~2 分钟）……")
        response = run_pichia_secretion_draft(request, paths)
        summary = response_to_summary(response)
        _write_status(status_path, "done", message="仿真完成。", result=summary)
    except Exception as exc:
        _write_status(status_path, "error", message=f"{type(exc).__name__}: {exc}")


__all__ = [
    "clear_last_result",
    "load_last_result",
    "poll_background_simulation",
    "save_last_result",
    "status_path_for_background_task",
    "submit_background_simulation",
]
