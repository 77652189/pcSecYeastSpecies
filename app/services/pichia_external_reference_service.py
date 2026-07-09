from __future__ import annotations

import csv
import io
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.external_refs import (
    DEFAULT_MANIFEST_FILENAME,
    DEFAULT_RECORDS_FILENAME,
    ExternalGeneFunctionEvidence,
    ExternalGprCandidateEvidence,
    ExternalReactionAssociation,
    ExternalReferenceRecord,
    load_external_reference_cache,
    load_external_reference_manifest,
)


EXTERNAL_REFERENCE_CACHE_DIR = Path("local_runs") / "pichia_external_reference_cache"


def load_external_reference_status(cache_root: Path) -> dict[str, Any]:
    """Return local external reference cache status without performing network IO."""

    root = Path(cache_root)
    records_path = _records_path(root)
    manifest_path = _manifest_path(root)
    base = {
        "cache_available": False,
        "cache_root": str(root),
        "records_path": str(records_path),
        "manifest_path": str(manifest_path),
        "record_count": 0,
        "source_counts": {},
        "record_type_counts": {},
        "retrieved_at_range": {"first": "", "last": ""},
        "failed_query_count": 0,
        "missing_files": _missing_files(root),
        "warnings": [],
        "recommended_refresh_command": _recommended_refresh_command(root),
    }
    if not records_path.exists():
        return base
    try:
        records = load_external_reference_cache(records_path)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        return {**base, "warnings": [f"failed to read external reference cache: {type(exc).__name__}: {exc}"]}
    manifest_payload: dict[str, Any] = {}
    if manifest_path.exists():
        try:
            manifest = load_external_reference_manifest(manifest_path)
            manifest_payload = {
                "failed_query_count": manifest.failed_query_count,
                "warnings": list(manifest.warnings),
            }
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            manifest_payload = {"warnings": [f"failed to read external reference manifest: {type(exc).__name__}: {exc}"]}
    retrieved_at_values = sorted(
        str(record.provenance.retrieved_at)
        for record in records
        if str(record.provenance.retrieved_at or "").strip()
    )
    return {
        **base,
        "cache_available": True,
        "record_count": len(records),
        "source_counts": _source_counts(records),
        "record_type_counts": _record_type_counts(records),
        "retrieved_at_range": {
            "first": retrieved_at_values[0] if retrieved_at_values else "",
            "last": retrieved_at_values[-1] if retrieved_at_values else "",
        },
        "failed_query_count": int(manifest_payload.get("failed_query_count", 0)),
        "missing_files": [],
        "warnings": list(manifest_payload.get("warnings", [])),
    }


def load_external_reference_browser_rows(
    cache_root: Path,
    *,
    query: str = "",
    evidence_kind: str | None = None,
    source_database: str | None = None,
    manual_review_only: bool = False,
) -> list[dict[str, Any]]:
    """Flatten local external reference records for service/UI display."""

    records_path = _records_path(Path(cache_root))
    if not records_path.exists():
        return []
    records = load_external_reference_cache(records_path)
    rows = [_browser_row(record) for record in records]
    return _filter_rows(
        rows,
        query=query,
        evidence_kind=evidence_kind,
        source_database=source_database,
        manual_review_only=manual_review_only,
    )


def submit_external_reference_refresh(
    *,
    homology_run_dir: Path,
    sources: tuple[str, ...],
    limit: int | None = None,
) -> dict[str, Any]:
    """Return an explicit cache-build command; do not launch network work from UI."""

    source_text = ",".join(sources)
    output_dir = Path("local_runs") / "pichia_external_reference_cache" / Path(homology_run_dir).name
    command = (
        "python scripts\\build_pichia_external_reference_cache.py "
        f"--sources {source_text} "
        f"--homology-run-dir {homology_run_dir} "
        f"--output-dir {output_dir}"
    )
    if limit is not None:
        command += f" --limit {int(limit)}"
    return {
        "submitted": False,
        "network_performed": False,
        "command": command,
        "message": "Run this command explicitly to refresh the external reference cache.",
    }


def export_external_reference_rows(
    rows: Iterable[Mapping[str, Any]],
    *,
    file_format: str = "tsv",
) -> bytes:
    resolved_rows = [dict(row) for row in rows]
    delimiter = "," if file_format.lower() == "csv" else "\t"
    fieldnames = _export_fieldnames(resolved_rows)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    for row in resolved_rows:
        writer.writerow({key: _export_value(row.get(key)) for key in fieldnames})
    return buffer.getvalue().encode("utf-8")


def _records_path(cache_root: Path) -> Path:
    if _looks_like_records_file(cache_root):
        return cache_root
    return cache_root / DEFAULT_RECORDS_FILENAME


def _manifest_path(cache_root: Path) -> Path:
    if _looks_like_records_file(cache_root):
        return cache_root.with_name(DEFAULT_MANIFEST_FILENAME)
    return cache_root / DEFAULT_MANIFEST_FILENAME


def _looks_like_records_file(cache_root: Path) -> bool:
    return cache_root.name == DEFAULT_RECORDS_FILENAME or cache_root.suffix.lower() == ".jsonl"


def _missing_files(cache_root: Path) -> list[str]:
    missing: list[str] = []
    if not _records_path(cache_root).exists():
        missing.append(DEFAULT_RECORDS_FILENAME)
    return missing


def _browser_row(record: Any) -> dict[str, Any]:
    row = _record_common_fields(record)
    if isinstance(record, ExternalReferenceRecord):
        row.update(
            {
                "gene_id": record.gene_id or "",
                "gene_name": record.gene_name or "",
                "locus_tag": record.locus_tag or "",
                "primary_accession": record.primary_accession,
                "protein_name": record.protein_name or "",
                "name_conflict_status": _name_conflict_status(record),
                "manual_review_reasons": _manual_review_reasons_from_warnings(record.provenance.warnings),
            }
        )
    elif isinstance(record, ExternalGeneFunctionEvidence):
        row.update(
            {
                "gene_id": record.gene_id,
                "protein_name": record.protein_name or "",
                "function_description": record.function_description or "",
                "ec_numbers": list(record.ec_numbers),
                "go_terms": list(record.go_terms),
                "pathways": list(record.pathways),
                "orthology": list(record.orthology),
                "evidence_confidence": record.evidence_scope,
                "manual_review_reasons": _manual_review_reasons_for_gene_function(record),
            }
        )
    elif isinstance(record, ExternalReactionAssociation):
        row.update(
            {
                "source_model_id": record.external_model_id,
                "source_reaction_id": record.external_reaction_id,
                "source_reaction_name": record.external_reaction_name or "",
                "source_gene_rule": record.gene_rule or "",
                "source_gene_ids": list(record.external_gene_ids),
                "gpr_transfer_status": record.association_status,
                "manual_review_reasons": _manual_review_reasons_from_warnings(record.provenance.warnings),
            }
        )
    elif isinstance(record, ExternalGprCandidateEvidence):
        row.update(
            {
                "pichia_gene_id": record.pichia_gene_id or "",
                "query_gene_id": record.query_gene_id or "",
                "source_model_id": record.external_model_id,
                "source_reaction_id": record.external_reaction_id,
                "source_gene_rule": record.external_gene_rule or "",
                "mapped_model_reaction_id": record.mapped_pichia_reaction_id or "",
                "mapped_pichia_gene_ids": list(record.mapped_pichia_gene_ids),
                "gene_mapping_status": record.gene_mapping_status,
                "reaction_mapping_status": record.reaction_mapping_status,
                "gpr_transfer_status": record.gpr_transfer_status,
                "evidence_confidence": record.confidence,
                "supporting_gene_evidence": list(record.supporting_gene_evidence),
                "manual_review_reasons": list(record.blocking_reasons),
            }
        )
    return row


def _record_common_fields(record: Any) -> dict[str, Any]:
    provenance = record.provenance
    return {
        "evidence_kind": str(record.record_type),
        "source_database": provenance.source_database,
        "source_version": provenance.source_version,
        "source_url": provenance.source_url,
        "source_query": provenance.source_query,
        "retrieved_at": provenance.retrieved_at,
        "warnings": list(provenance.warnings),
    }


def _filter_rows(
    rows: list[dict[str, Any]],
    *,
    query: str,
    evidence_kind: str | None,
    source_database: str | None,
    manual_review_only: bool,
) -> list[dict[str, Any]]:
    query_text = str(query or "").strip().lower()
    kind = str(evidence_kind or "").strip().lower()
    source = str(source_database or "").strip().lower()
    filtered: list[dict[str, Any]] = []
    for row in rows:
        if query_text and not _row_matches_query(row, query_text):
            continue
        if kind and str(row.get("evidence_kind") or "").lower() != kind:
            continue
        if source and str(row.get("source_database") or "").lower() != source:
            continue
        if manual_review_only and not row.get("manual_review_reasons"):
            continue
        filtered.append(row)
    return filtered


def _row_matches_query(row: Mapping[str, Any], query_text: str) -> bool:
    return any(query_text in str(value or "").lower() for value in row.values())


def _source_counts(records: Iterable[Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        source = str(record.provenance.source_database or "unknown")
        counts[source] = counts.get(source, 0) + 1
    return dict(sorted(counts.items()))


def _record_type_counts(records: Iterable[Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        record_type = str(record.record_type)
        counts[record_type] = counts.get(record_type, 0) + 1
    return dict(sorted(counts.items()))


def _name_conflict_status(record: ExternalReferenceRecord) -> str:
    warnings = " ".join(str(item).lower() for item in record.provenance.warnings)
    return "external_conflict" if "conflict" in warnings else "not_classified"


def _manual_review_reasons_for_gene_function(record: ExternalGeneFunctionEvidence) -> list[str]:
    if record.evidence_scope == "manual_review_required":
        return ["external gene function annotation requires manual review"]
    return []


def _manual_review_reasons_from_warnings(warnings: Iterable[str]) -> list[str]:
    return [str(item) for item in warnings if str(item).strip()]


def _export_fieldnames(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return ["evidence_kind", "source_database", "manual_review_reasons"]
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    return keys


def _export_value(value: Any) -> Any:
    if is_dataclass(value):
        return _export_value(asdict(value))
    if isinstance(value, Mapping):
        return ";".join(f"{key}={_export_value(item)}" for key, item in value.items())
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(item) for item in value)
    if value is None:
        return ""
    return value


def _recommended_refresh_command(cache_root: Path) -> str:
    output_dir = cache_root if cache_root.name else EXTERNAL_REFERENCE_CACHE_DIR / "manual_refresh"
    return (
        "python scripts\\build_pichia_external_reference_cache.py "
        "--sources uniprot,sgd "
        f"--output-dir {output_dir}"
    )


__all__ = [
    "EXTERNAL_REFERENCE_CACHE_DIR",
    "export_external_reference_rows",
    "load_external_reference_browser_rows",
    "load_external_reference_status",
    "submit_external_reference_refresh",
]
