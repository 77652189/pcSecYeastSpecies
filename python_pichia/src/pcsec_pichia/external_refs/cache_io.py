from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from pcsec_pichia.external_refs.schema import (
    ExternalCacheRecord,
    ExternalReferenceCacheManifest,
    ExternalReferenceSchemaError,
    build_external_reference_manifest,
    canonical_json,
    manifest_from_dict,
    manifest_to_dict,
    record_from_dict,
    record_to_dict,
    sha256_text,
    validate_external_cache_record,
    validate_no_duplicate_cache_keys,
)


DEFAULT_RECORDS_FILENAME = "external_reference_records.jsonl"
DEFAULT_MANIFEST_FILENAME = "external_reference_manifest.json"


def write_external_reference_cache(
    records: Iterable[ExternalCacheRecord],
    output_path: Path,
    *,
    manifest: ExternalReferenceCacheManifest | None = None,
    allow_duplicate_keys: bool = False,
) -> Path:
    """Write external reference records as local JSONL without performing network IO."""

    resolved_records = tuple(records)
    for record in resolved_records:
        validate_external_cache_record(record)
    if not allow_duplicate_keys:
        validate_no_duplicate_cache_keys(resolved_records)
    if manifest is not None:
        manifest.validate()
        if manifest.record_count != len(resolved_records):
            raise ExternalReferenceSchemaError(
                f"Manifest record_count={manifest.record_count} does not match JSONL record count={len(resolved_records)}."
            )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in resolved_records:
            handle.write(json.dumps(record_to_dict(record), ensure_ascii=False, sort_keys=True))
            handle.write("\n")
    return output_path


def load_external_reference_cache(path: Path) -> tuple[ExternalCacheRecord, ...]:
    """Load and validate external reference records from JSONL."""

    records: list[ExternalCacheRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                payload = json.loads(stripped)
            except json.JSONDecodeError as exc:
                raise ExternalReferenceSchemaError(f"Invalid JSONL at {path}:{line_number}.") from exc
            records.append(record_from_dict(payload))
    return tuple(records)


def write_external_reference_manifest(
    manifest: ExternalReferenceCacheManifest,
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(manifest_to_dict(manifest), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def load_external_reference_manifest(path: Path) -> ExternalReferenceCacheManifest:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ExternalReferenceSchemaError(f"Invalid manifest JSON: {path}.") from exc
    if not isinstance(payload, dict):
        raise ExternalReferenceSchemaError("External reference manifest must be a JSON object.")
    return manifest_from_dict(payload)


def write_external_reference_cache_bundle(
    records: Iterable[ExternalCacheRecord],
    output_dir: Path,
    *,
    query_count: int | None = None,
    failed_query_count: int = 0,
    input_cache_fingerprint: str | None = None,
    warnings: tuple[str, ...] = (),
    records_filename: str = DEFAULT_RECORDS_FILENAME,
    manifest_filename: str = DEFAULT_MANIFEST_FILENAME,
) -> ExternalReferenceCacheManifest:
    resolved_records = tuple(records)
    validate_no_duplicate_cache_keys(resolved_records)
    manifest = build_external_reference_manifest(
        resolved_records,
        query_count=query_count,
        failed_query_count=failed_query_count,
        input_cache_fingerprint=input_cache_fingerprint,
        warnings=warnings,
    )
    records_path = output_dir / records_filename
    write_external_reference_cache(resolved_records, records_path, manifest=manifest)
    write_external_reference_manifest(manifest, output_dir / manifest_filename)
    return manifest


def validate_external_reference_cache(
    path: Path,
    *,
    manifest_path: Path | None = None,
    allow_duplicate_keys: bool = False,
) -> ExternalReferenceCacheManifest:
    """Validate JSONL records and return a manifest derived from the current file."""

    records = load_external_reference_cache(path)
    if not allow_duplicate_keys:
        validate_no_duplicate_cache_keys(records)
    generated_manifest = build_external_reference_manifest(
        records,
        input_cache_fingerprint=cache_file_fingerprint(path),
    )
    if manifest_path is not None:
        stored_manifest = load_external_reference_manifest(manifest_path)
        if stored_manifest.record_count != generated_manifest.record_count:
            raise ExternalReferenceSchemaError(
                "Stored manifest record_count does not match validated JSONL record count."
            )
        if dict(stored_manifest.source_counts) != dict(generated_manifest.source_counts):
            raise ExternalReferenceSchemaError(
                "Stored manifest source_counts does not match validated JSONL source counts."
            )
        if dict(stored_manifest.record_type_counts) != dict(generated_manifest.record_type_counts):
            raise ExternalReferenceSchemaError(
                "Stored manifest record_type_counts does not match validated JSONL record type counts."
            )
        if stored_manifest.duplicate_key_count != generated_manifest.duplicate_key_count:
            raise ExternalReferenceSchemaError(
                "Stored manifest duplicate_key_count does not match validated JSONL duplicate count."
            )
    return generated_manifest


def cache_file_fingerprint(path: Path) -> str:
    return sha256_text(path.read_text(encoding="utf-8"))


def records_fingerprint(records: Iterable[ExternalCacheRecord]) -> str:
    payload = [record_to_dict(record) for record in records]
    return sha256_text(canonical_json({"records": payload}))


__all__ = [
    "DEFAULT_MANIFEST_FILENAME",
    "DEFAULT_RECORDS_FILENAME",
    "cache_file_fingerprint",
    "load_external_reference_cache",
    "load_external_reference_manifest",
    "records_fingerprint",
    "validate_external_reference_cache",
    "write_external_reference_cache",
    "write_external_reference_cache_bundle",
    "write_external_reference_manifest",
]
