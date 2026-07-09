from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Iterable

from pcsec_pichia.external_refs.cache_io import (
    DEFAULT_MANIFEST_FILENAME,
    DEFAULT_RECORDS_FILENAME,
    write_external_reference_cache_bundle,
)
from pcsec_pichia.external_refs.clients import (
    ExternalFetchConfig,
    ExternalFetchResult,
    ExternalHttpResponse,
    ExternalReferenceClient,
    default_external_reference_clients,
    fetch_external_references,
)
from pcsec_pichia.external_refs.queries import ExternalReferenceQuery
from pcsec_pichia.external_refs.schema import (
    ExternalReferenceCacheManifest,
    ExternalReferenceRecord,
    canonical_json,
    sha256_text,
)


FAILED_QUERIES_FILENAME = "failed_queries.jsonl"
SUMMARY_FILENAME = "external_reference_summary.md"


def build_external_reference_cache(
    queries: Iterable[ExternalReferenceQuery],
    output_dir: Path,
    *,
    config: ExternalFetchConfig | None = None,
    clients: Iterable[ExternalReferenceClient] | None = None,
    http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
    sleep: Callable[[float], None] | None = None,
) -> ExternalReferenceCacheManifest:
    resolved_queries = tuple(queries)
    resolved_config = config or ExternalFetchConfig()
    resolved_clients = tuple(clients or default_external_reference_clients())
    results = fetch_external_references(
        resolved_queries,
        resolved_clients,
        resolved_config,
        http_get=http_get,
        sleep=_sleep(sleep),
    )
    records = _dedupe_records(record for result in results for record in result.records)
    failed = tuple(result for result in results if result.failed)
    failed_query_count = _failed_query_count(resolved_queries, results)
    warnings = tuple(
        dict.fromkeys(
            warning
            for result in results
            for warning in result.warnings
            if str(warning).strip()
        )
    )
    manifest = write_external_reference_cache_bundle(
        records,
        output_dir,
        query_count=len(resolved_queries),
        failed_query_count=failed_query_count,
        input_cache_fingerprint=_queries_fingerprint(resolved_queries),
        warnings=warnings,
        records_filename=DEFAULT_RECORDS_FILENAME,
        manifest_filename=DEFAULT_MANIFEST_FILENAME,
    )
    write_failed_queries(failed, output_dir / FAILED_QUERIES_FILENAME)
    write_external_reference_summary(
        output_dir / SUMMARY_FILENAME,
        manifest=manifest,
        result_count=len(results),
        failed_result_count=len(failed),
    )
    return manifest


def write_failed_queries(results: Iterable[ExternalFetchResult], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for result in results:
            payload = result.to_failure_dict()
            if payload is not None:
                handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True))
                handle.write("\n")
    return path


def write_external_reference_summary(
    path: Path,
    *,
    manifest: ExternalReferenceCacheManifest,
    result_count: int,
    failed_result_count: int,
) -> Path:
    lines = [
        "# External Reference Cache Summary",
        "",
        f"- schema_version: `{manifest.cache_schema_version}`",
        f"- query_count: {manifest.query_count}",
        f"- fetch_result_count: {result_count}",
        f"- record_count: {manifest.record_count}",
        f"- failed_query_count: {manifest.failed_query_count}",
        f"- failed_result_count: {failed_result_count}",
        "",
        "## Source Counts",
        "",
    ]
    if manifest.source_counts:
        for source, count in sorted(manifest.source_counts.items()):
            lines.append(f"- {source}: {count}")
    else:
        lines.append("- none: 0")
    if manifest.warnings:
        lines.extend(["", "## Warnings", ""])
        for warning in manifest.warnings:
            lines.append(f"- {warning}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _dedupe_records(records: Iterable[ExternalReferenceRecord]) -> tuple[ExternalReferenceRecord, ...]:
    by_key: dict[str, ExternalReferenceRecord] = {}
    for record in records:
        by_key.setdefault(record.cache_key, record)
    return tuple(sorted(by_key.values(), key=lambda record: record.cache_key))


def _failed_query_count(
    queries: tuple[ExternalReferenceQuery, ...],
    results: tuple[ExternalFetchResult, ...],
) -> int:
    successful = {result.query.query_fingerprint for result in results if result.success}
    attempted = {result.query.query_fingerprint for result in results}
    return sum(
        1
        for query in queries
        if query.query_fingerprint in attempted and query.query_fingerprint not in successful
    )


def _queries_fingerprint(queries: tuple[ExternalReferenceQuery, ...]) -> str | None:
    if not queries:
        return None
    return sha256_text(canonical_json({"queries": [query.to_dict() for query in queries]}))


def _sleep(value: Callable[[float], None] | None) -> Callable[[float], None]:
    if value is not None:
        return value
    import time

    return time.sleep


__all__ = [
    "FAILED_QUERIES_FILENAME",
    "SUMMARY_FILENAME",
    "build_external_reference_cache",
    "write_external_reference_summary",
    "write_failed_queries",
]
