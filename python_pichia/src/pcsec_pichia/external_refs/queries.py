from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping

from pcsec_pichia.external_refs.schema import stable_cache_key


ExternalReferenceQueryType = Literal[
    "pichia_gene",
    "sce_homolog",
    "model_gene",
    "external_accession",
]


@dataclass(frozen=True)
class ExternalReferenceQuery:
    query_type: ExternalReferenceQueryType
    query_value: str
    source_context: str
    source_id: str
    source_row_id: str = ""
    preferred_sources: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def query_fingerprint(self) -> str:
        return external_reference_query_fingerprint(self)

    def to_dict(self) -> dict[str, Any]:
        return {
            "query_type": self.query_type,
            "query_value": self.query_value,
            "query_fingerprint": self.query_fingerprint,
            "source_context": self.source_context,
            "source_id": self.source_id,
            "source_row_id": self.source_row_id,
            "preferred_sources": list(self.preferred_sources),
            "warnings": list(self.warnings),
            "metadata": _json_ready(dict(self.metadata)),
        }


def external_reference_query_fingerprint(query: ExternalReferenceQuery) -> str:
    payload = {
        "query_type": query.query_type,
        "query_value": normalize_external_query_name(query.query_value),
        "source_context": query.source_context,
        "source_id": query.source_id,
    }
    return stable_cache_key(json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")))


def normalize_external_query_name(name: str) -> str:
    """Normalize external lookup names without guessing biological identity."""

    return str(name or "").strip().upper()


def dedupe_external_reference_queries(
    queries: Iterable[ExternalReferenceQuery],
) -> tuple[ExternalReferenceQuery, ...]:
    """Stable-deduplicate queries while preserving source/provenance warnings."""

    merged: dict[str, ExternalReferenceQuery] = {}
    for query in queries:
        if not query.query_value.strip():
            continue
        key = query.query_fingerprint
        previous = merged.get(key)
        if previous is None:
            merged[key] = query
            continue
        merged[key] = _merge_queries(previous, query)
    return tuple(sorted(merged.values(), key=_query_sort_key))


def build_external_reference_queries_from_homology_cache(
    homology_cache: Path | Iterable[Any],
    *,
    include_model_genes: bool = True,
    include_sce_queries: bool = True,
    include_pichia_genes: bool = True,
    include_external_accessions: bool = True,
) -> tuple[ExternalReferenceQuery, ...]:
    rows = _load_homology_rows(homology_cache)
    queries: list[ExternalReferenceQuery] = []
    for index, row in enumerate(rows):
        payload = _row_payload(row)
        row_id = _row_id(payload, index, "homology")
        metadata = _metadata_subset(
            payload,
            (
                "internal_common_name",
                "query_symbol",
                "review_status",
                "is_rbh",
                "in_model_gene_index",
            ),
        )
        if include_sce_queries:
            sce_value = payload.get("sce_orf")
            sce_warnings = _row_warnings(payload)
            if not sce_value and payload.get("query_symbol"):
                sce_warnings = (*sce_warnings, "homology cache row has no SCE ORF; falling back to query symbol.")
            queries.extend(
                _query_for_value(
                    sce_value or payload.get("query_symbol"),
                    "sce_homolog",
                    "homology_cache",
                    "homology_cache",
                    row_id,
                    metadata=metadata,
                    warnings=sce_warnings,
                )
            )
        if include_pichia_genes:
            queries.extend(
                _query_for_value(
                    payload.get("pichia_gene_id"),
                    "pichia_gene",
                    "homology_cache",
                    "homology_cache",
                    row_id,
                    metadata=metadata,
                    warnings=_row_warnings(payload),
                )
            )
        if include_model_genes:
            queries.extend(
                _query_for_value(
                    payload.get("pichia_model_gene_id"),
                    "model_gene",
                    "homology_cache",
                    "homology_cache",
                    row_id,
                    metadata=metadata,
                    warnings=_row_warnings(payload),
                )
            )
        if include_external_accessions:
            queries.extend(
                _query_for_value(
                    payload.get("external_accession"),
                    "external_accession",
                    "homology_cache",
                    "homology_cache",
                    row_id,
                    metadata=metadata,
                    warnings=_row_warnings(payload),
                )
            )
    return dedupe_external_reference_queries(queries)


def build_external_reference_queries_from_name_audit(
    name_audit: Path | Iterable[Any],
) -> tuple[ExternalReferenceQuery, ...]:
    rows = _load_name_audit_rows(name_audit)
    queries: list[ExternalReferenceQuery] = []
    for index, row in enumerate(rows):
        payload = _row_payload(row)
        row_id = _row_id(payload, index, "name_audit")
        metadata = _metadata_subset(
            payload,
            (
                "internal_common_name",
                "name_consistency_status",
                "external_crosscheck_status",
                "review_status",
            ),
        )
        queries.extend(
            _query_for_value(
                payload.get("internal_sequence_id"),
                "sce_homolog",
                "name_audit",
                "name_audit",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
        queries.extend(
            _query_for_value(
                payload.get("internal_gene_id"),
                "model_gene",
                "name_audit",
                "name_audit",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
        queries.extend(
            _query_for_value(
                payload.get("external_locus_tag"),
                "pichia_gene",
                "name_audit",
                "name_audit",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
        queries.extend(
            _query_for_value(
                payload.get("external_accession"),
                "external_accession",
                "name_audit",
                "name_audit",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
    return dedupe_external_reference_queries(queries)


def build_external_reference_queries_from_gene_catalog(
    gene_catalog_rows: Iterable[Any],
    *,
    source_context: str = "gene_catalog",
) -> tuple[ExternalReferenceQuery, ...]:
    queries: list[ExternalReferenceQuery] = []
    for index, row in enumerate(gene_catalog_rows):
        payload = _row_payload(row)
        row_id = _row_id(payload, index, source_context)
        metadata = _metadata_subset(
            payload,
            (
                "common_name",
                "category",
                "mapping_status",
                "recommended_use",
                "ko_support_status",
                "oe_support_status",
            ),
        )
        model_gene = payload.get("canonical_gene_id") or payload.get("mapped_model_gene_id") or payload.get("gene_id")
        queries.extend(
            _query_for_value(
                model_gene,
                "model_gene",
                source_context,
                "gene_catalog",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
        pichia_name = payload.get("standard_gene_symbol") or payload.get("display_name") or payload.get("common_name")
        queries.extend(
            _query_for_value(
                pichia_name,
                "pichia_gene",
                source_context,
                "gene_catalog",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
        queries.extend(
            _queries_from_external_ids(
                payload.get("external_ids"),
                source_context,
                "gene_catalog",
                row_id,
                metadata=metadata,
            )
        )
    return dedupe_external_reference_queries(queries)


def build_external_reference_queries_from_ko_oe_candidate_rows(
    candidate_rows: Iterable[Any],
    *,
    source_context: str = "ko_oe_candidate_rows",
) -> tuple[ExternalReferenceQuery, ...]:
    queries: list[ExternalReferenceQuery] = []
    for index, row in enumerate(candidate_rows):
        payload = _row_payload(row)
        row_id = _row_id(payload, index, source_context)
        metadata = _metadata_subset(
            payload,
            (
                "target_id",
                "screen_type",
                "intervention_type",
                "candidate_kind",
                "recommendation_tier",
            ),
        )
        model_gene = (
            payload.get("canonical_gene_id")
            or payload.get("mapped_model_gene_id")
            or payload.get("gene_id")
            or payload.get("input_gene_id")
        )
        queries.extend(
            _query_for_value(
                model_gene,
                "model_gene",
                source_context,
                "ko_oe_candidate_rows",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
        pichia_name = payload.get("common_name") or payload.get("display_name") or payload.get("standard_gene_symbol")
        queries.extend(
            _query_for_value(
                pichia_name,
                "pichia_gene",
                source_context,
                "ko_oe_candidate_rows",
                row_id,
                metadata=metadata,
                warnings=_row_warnings(payload),
            )
        )
        queries.extend(
            _queries_from_external_ids(
                payload.get("external_ids"),
                source_context,
                "ko_oe_candidate_rows",
                row_id,
                metadata=metadata,
            )
        )
    return dedupe_external_reference_queries(queries)


def build_external_reference_queries(
    *,
    homology_cache: Path | Iterable[Any] | None = None,
    name_audit: Path | Iterable[Any] | None = None,
    gene_catalog_rows: Iterable[Any] = (),
    ko_oe_candidate_rows: Iterable[Any] = (),
) -> tuple[ExternalReferenceQuery, ...]:
    queries: list[ExternalReferenceQuery] = []
    if homology_cache is not None:
        queries.extend(build_external_reference_queries_from_homology_cache(homology_cache))
    if name_audit is not None:
        queries.extend(build_external_reference_queries_from_name_audit(name_audit))
    queries.extend(build_external_reference_queries_from_gene_catalog(gene_catalog_rows))
    queries.extend(build_external_reference_queries_from_ko_oe_candidate_rows(ko_oe_candidate_rows))
    return dedupe_external_reference_queries(queries)


def _load_homology_rows(source: Path | Iterable[Any]) -> tuple[Any, ...]:
    if isinstance(source, Path):
        from pcsec_pichia.homology.crosswalk import load_homology_cache

        return load_homology_cache(source)
    return tuple(source)


def _load_name_audit_rows(source: Path | Iterable[Any]) -> tuple[Any, ...]:
    if isinstance(source, Path):
        from pcsec_pichia.homology.crosswalk import load_name_audit_cache

        return load_name_audit_cache(source)
    return tuple(source)


def _query_for_value(
    value: Any,
    query_type: ExternalReferenceQueryType,
    source_context: str,
    source_id: str,
    source_row_id: str,
    *,
    metadata: Mapping[str, Any],
    warnings: tuple[str, ...] = (),
) -> tuple[ExternalReferenceQuery, ...]:
    query_value = str(value or "").strip()
    if not query_value:
        return ()
    preferred_sources = _preferred_sources_for(query_type)
    return (
        ExternalReferenceQuery(
            query_type=query_type,
            query_value=query_value,
            source_context=source_context,
            source_id=source_id,
            source_row_id=source_row_id,
            preferred_sources=preferred_sources,
            warnings=tuple(str(warning) for warning in warnings if str(warning).strip()),
            metadata=dict(metadata),
        ),
    )


def _queries_from_external_ids(
    external_ids: Any,
    source_context: str,
    source_id: str,
    source_row_id: str,
    *,
    metadata: Mapping[str, Any],
) -> tuple[ExternalReferenceQuery, ...]:
    if not isinstance(external_ids, Mapping):
        return ()
    queries: list[ExternalReferenceQuery] = []
    for database, value in sorted(external_ids.items(), key=lambda item: str(item[0])):
        for accession in _values_from_external_id(value):
            queries.extend(
                _query_for_value(
                    accession,
                    "external_accession",
                    source_context,
                    source_id,
                    source_row_id,
                    metadata={**dict(metadata), "external_id_database": str(database)},
                    warnings=(),
                )
            )
    return tuple(queries)


def _values_from_external_id(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    text = str(value).strip()
    if not text:
        return ()
    if ";" in text:
        return tuple(part.strip() for part in text.split(";") if part.strip())
    return (text,)


def _merge_queries(left: ExternalReferenceQuery, right: ExternalReferenceQuery) -> ExternalReferenceQuery:
    return ExternalReferenceQuery(
        query_type=left.query_type,
        query_value=left.query_value,
        source_context=left.source_context,
        source_id=left.source_id,
        source_row_id=_merge_tokens(left.source_row_id, right.source_row_id),
        preferred_sources=tuple(dict.fromkeys((*left.preferred_sources, *right.preferred_sources))),
        warnings=tuple(dict.fromkeys((*left.warnings, *right.warnings))),
        metadata={
            **dict(left.metadata),
            "merged_source_row_ids": _merge_tokens(
                str(left.metadata.get("merged_source_row_ids") or left.source_row_id),
                str(right.metadata.get("merged_source_row_ids") or right.source_row_id),
            ),
        },
    )


def _merge_tokens(*values: str) -> str:
    tokens: list[str] = []
    for value in values:
        for token in str(value or "").split(";"):
            if token.strip() and token.strip() not in tokens:
                tokens.append(token.strip())
    return ";".join(tokens)


def _query_sort_key(query: ExternalReferenceQuery) -> tuple[str, str, str, str]:
    return (
        query.query_type,
        normalize_external_query_name(query.query_value),
        query.source_context,
        query.source_id,
    )


def _row_payload(row: Any) -> dict[str, Any]:
    if is_dataclass(row):
        return asdict(row)
    if isinstance(row, Mapping):
        return dict(row)
    payload: dict[str, Any] = {}
    for key in dir(row):
        if key.startswith("_"):
            continue
        try:
            value = getattr(row, key)
        except AttributeError:
            continue
        if not callable(value):
            payload[key] = value
    return payload


def _row_id(payload: Mapping[str, Any], index: int, prefix: str) -> str:
    for key in (
        "cache_key",
        "row_id",
        "gene_id",
        "canonical_gene_id",
        "pichia_model_gene_id",
        "pichia_gene_id",
        "internal_gene_id",
        "internal_sequence_id",
        "external_accession",
        "input_gene_id",
        "common_name",
    ):
        value = str(payload.get(key) or "").strip()
        if value:
            return f"{prefix}:{value}"
    return f"{prefix}:row:{index}"


def _row_warnings(payload: Mapping[str, Any]) -> tuple[str, ...]:
    warnings: list[str] = [str(warning) for warning in _values_from_external_id(payload.get("warnings"))]
    return tuple(dict.fromkeys(warnings))


def _metadata_subset(payload: Mapping[str, Any], keys: tuple[str, ...]) -> dict[str, Any]:
    return {
        key: _json_ready(payload[key])
        for key in keys
        if key in payload and payload[key] not in (None, "", (), [])
    }


def _preferred_sources_for(query_type: ExternalReferenceQueryType) -> tuple[str, ...]:
    if query_type == "sce_homolog":
        return ("sgd", "uniprot")
    if query_type == "external_accession":
        return ("uniprot", "ncbi", "sgd")
    if query_type == "model_gene":
        return ("uniprot", "ncbi")
    return ("uniprot", "ncbi", "sgd")


def _json_ready(value: Any) -> Any:
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "ExternalReferenceQuery",
    "ExternalReferenceQueryType",
    "build_external_reference_queries",
    "build_external_reference_queries_from_gene_catalog",
    "build_external_reference_queries_from_homology_cache",
    "build_external_reference_queries_from_ko_oe_candidate_rows",
    "build_external_reference_queries_from_name_audit",
    "dedupe_external_reference_queries",
    "external_reference_query_fingerprint",
    "normalize_external_query_name",
]
