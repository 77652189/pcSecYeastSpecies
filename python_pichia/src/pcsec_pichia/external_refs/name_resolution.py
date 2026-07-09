from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping

from pcsec_pichia.external_refs.schema import ExternalReferenceRecord


EXTERNAL_MATCH_CONFIRMED = "external_match_confirmed"
EXTERNAL_ALIAS_CONFIRMED = "external_alias_confirmed"
EXTERNAL_LOCUS_CONFIRMED = "external_locus_confirmed"
EXTERNAL_CONFLICT = "external_conflict"
EXTERNAL_REFERENCE_MISSING = "external_reference_missing"


@dataclass(frozen=True)
class NameResolutionCandidate:
    internal_gene_id: str
    internal_common_name: str | None
    external_name_status: str
    matched_record_count: int
    source_databases: tuple[str, ...]
    source_accessions: tuple[str, ...]
    matched_gene_names: tuple[str, ...]
    matched_locus_tags: tuple[str, ...]
    warnings: tuple[str, ...] = ()
    manual_review_reasons: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def classify_external_name_consistency(
    *,
    internal_gene_id: str,
    internal_common_name: str | None,
    internal_aliases: Iterable[str],
    external_records: Iterable[ExternalReferenceRecord],
) -> NameResolutionCandidate:
    """Classify external database name agreement without changing model facts."""

    records = tuple(external_records)
    source_databases = _unique(record.provenance.source_database for record in records)
    source_accessions = _unique(record.primary_accession for record in records)
    gene_names = _unique(record.gene_name or "" for record in records)
    locus_tags = _unique(record.locus_tag or "" for record in records)
    warnings = _unique(
        warning
        for record in records
        for warning in record.provenance.warnings
    )
    if not records:
        return NameResolutionCandidate(
            internal_gene_id=internal_gene_id,
            internal_common_name=internal_common_name,
            external_name_status=EXTERNAL_REFERENCE_MISSING,
            matched_record_count=0,
            source_databases=(),
            source_accessions=(),
            matched_gene_names=(),
            matched_locus_tags=(),
            manual_review_reasons=("no external reference matched the current name-audit row",),
        )

    internal_name_tokens = _name_tokens(internal_common_name)
    alias_tokens = _name_tokens(*internal_aliases)
    locus_tokens = _locus_tokens(internal_gene_id)
    gene_name_tokens = _name_tokens(*(record.gene_name or "" for record in records))
    external_alias_tokens = _name_tokens(*(alias for record in records for alias in record.aliases))
    external_locus_tokens = _locus_tokens(
        *(value for record in records for value in (record.gene_id, record.locus_tag, record.primary_accession))
    )

    if internal_name_tokens and internal_name_tokens & gene_name_tokens:
        status = EXTERNAL_MATCH_CONFIRMED
        review_reasons: tuple[str, ...] = ()
    elif internal_name_tokens and internal_name_tokens & external_alias_tokens:
        status = EXTERNAL_ALIAS_CONFIRMED
        review_reasons = ()
    elif locus_tokens and locus_tokens & external_locus_tokens:
        status = EXTERNAL_LOCUS_CONFIRMED
        review_reasons = ()
    else:
        status = EXTERNAL_CONFLICT
        review_reasons = (
            "external reference names, aliases, or loci do not confirm the current name-audit row",
        )

    if alias_tokens and alias_tokens & gene_name_tokens and status != EXTERNAL_MATCH_CONFIRMED:
        status = EXTERNAL_ALIAS_CONFIRMED
        review_reasons = ()

    return NameResolutionCandidate(
        internal_gene_id=internal_gene_id,
        internal_common_name=internal_common_name,
        external_name_status=status,
        matched_record_count=len(records),
        source_databases=source_databases,
        source_accessions=source_accessions,
        matched_gene_names=gene_names,
        matched_locus_tags=locus_tags,
        warnings=warnings,
        manual_review_reasons=review_reasons,
    )


def select_external_records_for_name_audit_row(
    row: Mapping[str, Any],
    records: Iterable[ExternalReferenceRecord],
) -> tuple[ExternalReferenceRecord, ...]:
    """Return records with accession/name/locus evidence for one name-audit row."""

    row_tokens = _row_match_tokens(row)
    matched: list[ExternalReferenceRecord] = []
    for record in records:
        if row_tokens & _record_match_tokens(record):
            matched.append(record)
    return tuple(sorted(matched, key=lambda record: record.cache_key))


def _row_match_tokens(row: Mapping[str, Any]) -> set[str]:
    values: list[Any] = [
        row.get("internal_gene_id"),
        row.get("internal_common_name"),
        row.get("internal_sequence_id"),
        row.get("external_accession"),
        row.get("external_gene_name"),
        row.get("external_locus_tag"),
    ]
    values.extend(_iter_values(row.get("external_aliases")))
    return _match_tokens(*values)


def _record_match_tokens(record: ExternalReferenceRecord) -> set[str]:
    return _match_tokens(
        record.primary_accession,
        record.gene_id,
        record.gene_name,
        record.locus_tag,
        record.provenance.source_query,
        *record.aliases,
    )


def _match_tokens(*values: object) -> set[str]:
    return {token for value in values for token in (_normalize(value), *_split_tokens(value)) if token}


def _name_tokens(*values: object) -> set[str]:
    return {
        token
        for value in values
        for token in (_normalize(value), *_split_tokens(value))
        if token
    }


def _locus_tokens(*values: object) -> set[str]:
    return {
        token
        for value in values
        for token in (_normalize(value), *_split_tokens(value))
        if token and _looks_like_locus_or_accession(token)
    }


def _split_tokens(value: object) -> tuple[str, ...]:
    text = str(value or "")
    for sep in ("/", ",", ";", "|"):
        text = text.replace(sep, " ")
    return tuple(_normalize(part) for part in text.split() if _normalize(part))


def _normalize(value: object) -> str:
    return str(value or "").strip().upper().replace(" ", "").replace("_", "").replace("-", "")


def _looks_like_locus_or_accession(token: str) -> bool:
    return token.startswith(("PAS", "Y")) or any(char.isdigit() for char in token)


def _iter_values(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    text = str(value).strip()
    if not text:
        return ()
    return tuple(part.strip() for part in text.split(";") if part.strip())


def _unique(values: Iterable[str]) -> tuple[str, ...]:
    result: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in result:
            result.append(text)
    return tuple(result)


__all__ = [
    "EXTERNAL_ALIAS_CONFIRMED",
    "EXTERNAL_CONFLICT",
    "EXTERNAL_LOCUS_CONFIRMED",
    "EXTERNAL_MATCH_CONFIRMED",
    "EXTERNAL_REFERENCE_MISSING",
    "NameResolutionCandidate",
    "classify_external_name_consistency",
    "select_external_records_for_name_audit_row",
]
