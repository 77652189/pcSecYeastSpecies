from __future__ import annotations

import csv
import io
import json
from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.homology.crosswalk import load_name_audit_cache, load_rule_transfer_audit_cache


HOMOLOGY_AUDIT_CACHE_DIR = Path("local_runs") / "pichia_homology_cache"
NAME_AUDIT_JSONL = "sce_to_pichia_name_audit.jsonl"
RULE_TRANSFER_JSONL = "sce_to_pichia_rule_transfer_audit.jsonl"
SUMMARY_JSON = "homology_audit_summary.json"
REQUIRED_CACHE_FILES: tuple[str, ...] = (
    NAME_AUDIT_JSONL,
    RULE_TRANSFER_JSONL,
    SUMMARY_JSON,
)


def load_homology_audit_browser_data(
    *,
    query: str = "",
    review_status: str | None = None,
    name_consistency_status: str | None = None,
    rule_transfer_status: str | None = None,
    is_rbh: bool | None = None,
    in_model_gene_index: bool | None = None,
    min_identity: float | None = None,
    paths: ProjectPaths | None = None,
    cache_root: Path | None = None,
) -> dict[str, Any]:
    """Load cached homology audit rows for Streamlit display.

    This service is a facade over cache files. It does not run BLAST or perform
    scientific classification; those decisions live in ``pcsec_pichia.homology``.
    """

    status = homology_audit_cache_status(paths=paths, cache_root=cache_root)
    if not status["cache_available"]:
        return {
            "cache_status": status,
            "summary": {},
            "name_audit_rows": [],
            "rule_transfer_audit_rows": [],
        }

    run_dir = Path(str(status["cache_root"]))
    name_rows = [_row_dict(row) for row in load_name_audit_cache(run_dir / NAME_AUDIT_JSONL)]
    rule_rows = [_row_dict(row) for row in load_rule_transfer_audit_cache(run_dir / RULE_TRANSFER_JSONL)]
    filtered_name_rows = _filter_rows(
        name_rows,
        query=query,
        review_status=review_status,
        name_consistency_status=name_consistency_status,
        is_rbh=is_rbh,
        in_model_gene_index=in_model_gene_index,
        min_identity=min_identity,
    )
    filtered_rule_rows = _filter_rows(
        rule_rows,
        query=query,
        review_status=review_status,
        rule_transfer_status=rule_transfer_status,
        is_rbh=is_rbh,
        in_model_gene_index=in_model_gene_index,
        min_identity=min_identity,
    )
    return {
        "cache_status": status,
        "summary": _read_summary(run_dir / SUMMARY_JSON),
        "name_audit_rows": filtered_name_rows,
        "rule_transfer_audit_rows": filtered_rule_rows,
    }


def homology_audit_cache_status(
    *,
    paths: ProjectPaths | None = None,
    cache_root: Path | None = None,
) -> dict[str, Any]:
    root = _cache_base(paths=paths, cache_root=cache_root)
    run_dir = _latest_valid_run(root)
    if run_dir is None:
        missing = list(REQUIRED_CACHE_FILES)
        return {
            "cache_available": False,
            "cache_root": str(root),
            "generated_at": "",
            "row_count": 0,
            "missing_files": missing,
            "recommended_build_command": _recommended_build_command(root),
        }

    missing_files = _missing_files(run_dir)
    summary = _read_summary(run_dir / SUMMARY_JSON)
    return {
        "cache_available": not missing_files,
        "cache_root": str(run_dir),
        "generated_at": datetime.fromtimestamp(run_dir.stat().st_mtime).isoformat(timespec="seconds"),
        "row_count": int(
            summary.get("name_audit_row_count")
            or summary.get("rule_transfer_row_count")
            or summary.get("homology_row_count")
            or 0
        ),
        "missing_files": missing_files,
        "recommended_build_command": _recommended_build_command(root),
    }


def export_homology_audit_rows(
    rows: list[dict[str, Any]],
    *,
    file_format: str = "tsv",
) -> bytes:
    delimiter = "," if file_format.lower() == "csv" else "\t"
    fieldnames = _export_fieldnames(rows)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: _export_value(row.get(key)) for key in fieldnames})
    return buffer.getvalue().encode("utf-8")


def _cache_base(paths: ProjectPaths | None = None, cache_root: Path | None = None) -> Path:
    if cache_root is not None:
        root = Path(cache_root)
        if all((root / file_name).exists() for file_name in REQUIRED_CACHE_FILES):
            return root
        return root
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    return resolved_paths.repo_root / HOMOLOGY_AUDIT_CACHE_DIR


def _latest_valid_run(root: Path) -> Path | None:
    if all((root / file_name).exists() for file_name in REQUIRED_CACHE_FILES):
        return root
    if not root.exists():
        return None
    candidates = [
        path
        for path in root.iterdir()
        if path.is_dir() and not _missing_files(path)
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _missing_files(run_dir: Path) -> list[str]:
    return [file_name for file_name in REQUIRED_CACHE_FILES if not (run_dir / file_name).exists()]


def _read_summary(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, json.JSONDecodeError, UnicodeDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _filter_rows(
    rows: list[dict[str, Any]],
    *,
    query: str = "",
    review_status: str | None = None,
    name_consistency_status: str | None = None,
    rule_transfer_status: str | None = None,
    is_rbh: bool | None = None,
    in_model_gene_index: bool | None = None,
    min_identity: float | None = None,
) -> list[dict[str, Any]]:
    query_text = str(query or "").strip().lower()
    filtered: list[dict[str, Any]] = []
    for row in rows:
        if query_text and not _row_matches_query(row, query_text):
            continue
        if review_status and row.get("review_status") != review_status and row.get("homology_review_status") != review_status:
            continue
        if name_consistency_status and row.get("name_consistency_status") != name_consistency_status:
            continue
        if rule_transfer_status and row.get("rule_transfer_status") != rule_transfer_status:
            continue
        if is_rbh is not None and bool(row.get("is_rbh")) is not is_rbh:
            continue
        if in_model_gene_index is not None and bool(row.get("in_model_gene_index")) is not in_model_gene_index:
            continue
        if min_identity is not None and not _passes_min_identity(row, min_identity):
            continue
        filtered.append(row)
    return filtered


def _row_matches_query(row: dict[str, Any], query_text: str) -> bool:
    searchable = (
        row.get("internal_common_name"),
        row.get("query_symbol"),
        row.get("internal_sequence_id"),
        row.get("sce_orf"),
        row.get("external_gene_name"),
        row.get("external_locus_tag"),
        row.get("pichia_gene_id"),
        row.get("pichia_model_gene_id"),
        row.get("external_accession"),
        row.get("external_crosscheck_status"),
        row.get("external_crosscheck_sources"),
        row.get("external_crosscheck_warnings"),
    )
    return any(query_text in str(value or "").lower() for value in searchable)


def _passes_min_identity(row: dict[str, Any], min_identity: float) -> bool:
    try:
        return float(row.get("identity_pct")) >= min_identity
    except (TypeError, ValueError):
        return False


def _row_dict(row: Any) -> dict[str, Any]:
    payload = asdict(row) if is_dataclass(row) else dict(row)
    for key, value in list(payload.items()):
        if isinstance(value, tuple):
            payload[key] = list(value)
    return payload


def _export_fieldnames(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return ["internal_common_name", "query_symbol", "review_status"]
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    return keys


def _export_value(value: Any) -> Any:
    if isinstance(value, list):
        return ";".join(str(item) for item in value)
    if value is None:
        return ""
    return value


def _recommended_build_command(root: Path) -> str:
    output_dir = root if root.name != "pichia_homology_cache" else root / "smoke"
    return f"python scripts\\build_pichia_homology_cache.py --catalog-only --output-dir {output_dir}"


__all__ = [
    "HOMOLOGY_AUDIT_CACHE_DIR",
    "NAME_AUDIT_JSONL",
    "RULE_TRANSFER_JSONL",
    "SUMMARY_JSON",
    "export_homology_audit_rows",
    "homology_audit_cache_status",
    "load_homology_audit_browser_data",
]
