from __future__ import annotations

import argparse
import csv
import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python_pichia" / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from app.services.pichia_gene_catalog_service import load_pichia_gene_id_standardization  # noqa: E402
from pcsec_pichia.core.paths import ProjectPaths  # noqa: E402
from pcsec_pichia.services.gene_id_standardization import (  # noqa: E402
    STANDARD_NAME_FIELD_NAMES,
    build_standard_name_lookup,
    enrich_gene_standard_name_fields,
    load_pichia_gene_id_standardization_cache,
    standard_name_fields_for_csv,
)

SCIENTIFIC_NUMERIC_FIELDS = {
    "objective_value",
    "growth_rate",
    "secretion_flux",
    "delta_objective",
    "ranking",
    "secretion_ratio_vs_wildtype",
    "growth_retention_ratio",
    "max_feasible_mu",
    "secretion_at_max_feasible_mu",
    "wildtype_max_feasible_mu",
    "wildtype_secretion_at_max_feasible_mu",
}
REPORT_PATH = REPO_ROOT / "local_runs" / "standard_name_backfill_report.json"
_CSV_FIELD_SIZE_LIMIT = sys.maxsize
while True:
    try:
        csv.field_size_limit(_CSV_FIELD_SIZE_LIMIT)
        break
    except OverflowError:
        _CSV_FIELD_SIZE_LIMIT = int(_CSV_FIELD_SIZE_LIMIT / 10)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.apply == args.dry_run:
        raise SystemExit("Choose exactly one of --dry-run or --apply.")
    report = run_backfill(dry_run=args.dry_run)
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


def run_backfill(*, dry_run: bool) -> dict[str, Any]:
    lookup = build_standard_name_lookup(_load_standard_name_rows(dry_run=dry_run))
    report: dict[str, Any] = {
        "scanned_file_count": 0,
        "updated_file_count": 0,
        "updated_row_count": 0,
        "missing_standard_name_count": 0,
        "not_gene_candidate_count": 0,
        "skipped_files": [],
        "warnings": [],
        "dry_run": dry_run,
        "protected_results_detected": _protected_results_detected(),
        "numeric_invariance_status": "passed",
    }
    local_runs = REPO_ROOT / "local_runs"
    if not local_runs.exists():
        report["warnings"].append("local_runs directory does not exist.")
        return report

    for path in sorted(local_runs.rglob("*")):
        if not path.is_file() or path == REPORT_PATH:
            continue
        if path.suffix.lower() == ".csv":
            changed, rows_updated = _backfill_csv(path, lookup, dry_run=dry_run, report=report)
        elif path.suffix.lower() == ".json":
            changed, rows_updated = _backfill_json(path, lookup, dry_run=dry_run, report=report)
        else:
            continue
        report["scanned_file_count"] += 1
        if changed:
            report["updated_file_count"] += 1
            report["updated_row_count"] += rows_updated

    if not dry_run:
        REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
        REPORT_PATH.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    return report


def _load_standard_name_rows(*, dry_run: bool) -> list[dict[str, Any]]:
    cache_path = (
        REPO_ROOT
        / "local_runs"
        / "streamlit_pichia_runs"
        / "gene_catalog_cache"
        / "pichia_gene_id_standardization.json"
    )
    cached_rows = load_pichia_gene_id_standardization_cache(cache_path) if cache_path.exists() else ()
    if cached_rows:
        return [row.to_dict() for row in cached_rows]
    if dry_run:
        return []
    return load_pichia_gene_id_standardization(paths=ProjectPaths(repo_root=REPO_ROOT))


def _backfill_csv(
    path: Path,
    lookup: dict[str, Any],
    *,
    dry_run: bool,
    report: dict[str, Any],
) -> tuple[bool, int]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        report["skipped_files"].append({"path": str(path), "reason": f"csv_read_failed: {exc}"})
        return False, 0
    if not rows or not _looks_like_candidate_csv(fieldnames):
        return False, 0

    before_numeric = _csv_numeric_snapshot(rows)
    changed = False
    rows_updated = 0
    enriched_rows: list[dict[str, Any]] = []
    for row in rows:
        enriched = enrich_gene_standard_name_fields(dict(row), lookup)
        csv_fields = standard_name_fields_for_csv(enriched)
        existing = {key: row.get(key, "") for key in STANDARD_NAME_FIELD_NAMES}
        for key, value in csv_fields.items():
            row[key] = value
        _count_status(row.get("standard_name_status"), report)
        if existing != {key: row.get(key, "") for key in STANDARD_NAME_FIELD_NAMES}:
            changed = True
            rows_updated += 1
        enriched_rows.append(row)

    after_numeric = _csv_numeric_snapshot(enriched_rows)
    if before_numeric != after_numeric:
        report["numeric_invariance_status"] = "failed"
        report["warnings"].append(f"numeric invariance failed for {path}")
        return False, 0
    if changed and not dry_run:
        output_fields = [*fieldnames, *(field for field in STANDARD_NAME_FIELD_NAMES if field not in fieldnames)]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=output_fields, extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            writer.writerows(enriched_rows)
    return changed, rows_updated


def _backfill_json(
    path: Path,
    lookup: dict[str, Any],
    *,
    dry_run: bool,
    report: dict[str, Any],
) -> tuple[bool, int]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        report["skipped_files"].append({"path": str(path), "reason": f"json_read_failed: {exc}"})
        return False, 0
    before_numeric = _json_numeric_snapshot(payload)
    updated, rows_updated = _enrich_json_value(payload, lookup, report)
    after_numeric = _json_numeric_snapshot(payload)
    if before_numeric != after_numeric:
        report["numeric_invariance_status"] = "failed"
        report["warnings"].append(f"numeric invariance failed for {path}")
        return False, 0
    if updated and not dry_run:
        path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=False), encoding="utf-8")
    return updated, rows_updated


def _enrich_json_value(value: Any, lookup: dict[str, Any], report: dict[str, Any]) -> tuple[bool, int]:
    changed = False
    rows_updated = 0
    if isinstance(value, list):
        for item in value:
            item_changed, item_rows = _enrich_json_value(item, lookup, report)
            changed = changed or item_changed
            rows_updated += item_rows
        return changed, rows_updated
    if not isinstance(value, dict):
        return False, 0

    if _looks_like_candidate_json_row(value):
        before = {key: deepcopy(value.get(key)) for key in STANDARD_NAME_FIELD_NAMES}
        enriched = enrich_gene_standard_name_fields(value, lookup)
        for key in STANDARD_NAME_FIELD_NAMES:
            value[key] = enriched.get(key)
        _count_status(value.get("standard_name_status"), report)
        after = {key: value.get(key) for key in STANDARD_NAME_FIELD_NAMES}
        changed = before != after
        rows_updated += 1 if changed else 0

    for item in value.values():
        item_changed, item_rows = _enrich_json_value(item, lookup, report)
        changed = changed or item_changed
        rows_updated += item_rows
    return changed, rows_updated


def _looks_like_candidate_csv(fieldnames: list[str]) -> bool:
    fields = set(fieldnames)
    return bool(fields & {"gene_id", "canonical_gene_id", "input_gene_id"}) and bool(
        fields & {"intervention_type", "screen_type", "candidate_kind", "target_id"}
    )


def _looks_like_candidate_json_row(row: dict[str, Any]) -> bool:
    return bool({"gene_id", "canonical_gene_id", "input_gene_id"} & set(row)) and bool(
        {"intervention_type", "screen_type", "candidate_kind", "target_id"} & set(row)
    )


def _csv_numeric_snapshot(rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    return [
        {field: str(row.get(field, "")) for field in SCIENTIFIC_NUMERIC_FIELDS if field in row}
        for row in rows
    ]


def _json_numeric_snapshot(value: Any) -> Any:
    values: list[tuple[str, Any]] = []
    _collect_json_numeric_values(value, values)
    return values


def _collect_json_numeric_values(value: Any, values: list[tuple[str, Any]]) -> None:
    if isinstance(value, list):
        for item in value:
            _collect_json_numeric_values(item, values)
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if key in SCIENTIFIC_NUMERIC_FIELDS:
                values.append((key, deepcopy(item)))
            else:
                _collect_json_numeric_values(item, values)


def _count_status(status: Any, report: dict[str, Any]) -> None:
    if status == "missing_standard_name":
        report["missing_standard_name_count"] += 1
    elif status == "not_gene_candidate":
        report["not_gene_candidate_count"] += 1


def _protected_results_detected() -> bool:
    results_dir = REPO_ROOT / "Results"
    if not results_dir.exists():
        return False
    return any(path.suffix.lower() in {".csv", ".json"} for path in results_dir.rglob("*") if path.is_file())


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Backfill display-only Pichia standard gene name fields in local_runs outputs.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dry-run", action="store_true", help="Scan and report without writing files.")
    group.add_argument("--apply", action="store_true", help="Append standard naming fields to supported local_runs files.")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
