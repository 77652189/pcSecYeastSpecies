from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any, Iterable, Mapping

from pcsec_pichia.external_refs.clients import ExternalFetchResult
from pcsec_pichia.external_refs.name_resolution import (
    classify_external_name_consistency,
    select_external_records_for_name_audit_row,
)
from pcsec_pichia.external_refs.schema import ExternalReferenceRecord


def merge_external_fetch_results(
    results: Iterable[ExternalFetchResult],
) -> tuple[ExternalReferenceRecord, ...]:
    """Collect successful external records from fetch results with stable dedupe."""

    by_key: dict[str, ExternalReferenceRecord] = {}
    for result in results:
        for record in result.records:
            by_key.setdefault(record.cache_key, record)
    return tuple(sorted(by_key.values(), key=lambda record: record.cache_key))


def attach_external_references_to_name_audit(
    name_audit_rows: Iterable[Mapping[str, Any] | Any],
    records: Iterable[ExternalReferenceRecord],
) -> tuple[Mapping[str, Any], ...]:
    """Attach external name-resolution fields while preserving original row facts."""

    resolved_records = tuple(records)
    merged: list[Mapping[str, Any]] = []
    for row in name_audit_rows:
        payload = _row_payload(row)
        matched_records = select_external_records_for_name_audit_row(payload, resolved_records)
        candidate = classify_external_name_consistency(
            internal_gene_id=str(payload.get("internal_gene_id") or ""),
            internal_common_name=_optional_text(payload.get("internal_common_name")),
            internal_aliases=_aliases_from_payload(payload),
            external_records=matched_records,
        )
        merged.append(
            {
                **payload,
                "external_name_status": candidate.external_name_status,
                "external_reference_record_count": candidate.matched_record_count,
                "external_reference_sources": candidate.source_databases,
                "external_reference_accessions": candidate.source_accessions,
                "external_reference_gene_names": candidate.matched_gene_names,
                "external_reference_locus_tags": candidate.matched_locus_tags,
                "external_reference_warnings": candidate.warnings,
                "external_manual_review_reasons": candidate.manual_review_reasons,
            }
        )
    return tuple(merged)


def _row_payload(row: Mapping[str, Any] | Any) -> dict[str, Any]:
    if is_dataclass(row):
        return asdict(row)
    if isinstance(row, Mapping):
        return dict(row)
    payload: dict[str, Any] = {}
    for key in dir(row):
        if key.startswith("_"):
            continue
        value = getattr(row, key, None)
        if not callable(value):
            payload[key] = value
    return payload


def _aliases_from_payload(payload: Mapping[str, Any]) -> tuple[str, ...]:
    aliases = []
    for key in ("external_aliases", "internal_aliases", "aliases"):
        aliases.extend(_tuple_values(payload.get(key)))
    for key in ("external_gene_name", "external_locus_tag", "external_accession"):
        value = str(payload.get(key) or "").strip()
        if value:
            aliases.append(value)
    return tuple(dict.fromkeys(aliases))


def _tuple_values(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (tuple, list, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    text = str(value).strip()
    if not text:
        return ()
    return tuple(part.strip() for part in text.split(";") if part.strip())


def _optional_text(value: object) -> str | None:
    text = str(value or "").strip()
    return text or None


__all__ = [
    "attach_external_references_to_name_audit",
    "merge_external_fetch_results",
]
