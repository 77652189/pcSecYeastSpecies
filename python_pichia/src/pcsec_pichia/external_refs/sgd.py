from __future__ import annotations

import time
import urllib.parse
from dataclasses import dataclass
from typing import Any, Callable, Mapping

from pcsec_pichia.external_refs.clients import (
    ExternalFetchConfig,
    ExternalFetchResult,
    ExternalHttpResponse,
    failed_fetch_result,
    request_json,
)
from pcsec_pichia.external_refs.queries import ExternalReferenceQuery
from pcsec_pichia.external_refs.schema import ExternalReferenceProvenance, ExternalReferenceRecord, utc_now_iso


SGD_LOCUS_URL = "https://www.yeastgenome.org/backend/locus"
SGD_SOURCE_VERSION = "sgd-backend"


@dataclass(frozen=True)
class SgdReferenceClient:
    source_database: str = "sgd"

    def fetch(
        self,
        query: ExternalReferenceQuery,
        config: ExternalFetchConfig,
        *,
        http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> ExternalFetchResult:
        return fetch_sgd_reference(query, config, http_get=http_get, sleep=sleep)


def fetch_sgd_reference(
    query: ExternalReferenceQuery,
    config: ExternalFetchConfig | None = None,
    *,
    http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
    sleep: Callable[[float], None] = time.sleep,
) -> ExternalFetchResult:
    resolved_config = config or ExternalFetchConfig()
    url = build_sgd_url(query)
    payload, response, warnings, attempts = request_json(
        url,
        resolved_config,
        http_get=http_get,
        sleep=sleep,
    )
    if not isinstance(payload, Mapping):
        return failed_fetch_result(
            source_database="sgd",
            query=query,
            source_url=url,
            response=response,
            warnings=warnings or ("SGD returned no JSON object.",),
            attempts=attempts,
            error_type="invalid_response",
        )
    record = _record_from_sgd(payload, query, url, response)
    if not record.primary_accession and not record.gene_name and not record.locus_tag:
        return failed_fetch_result(
            source_database="sgd",
            query=query,
            source_url=url,
            response=response,
            warnings=(*warnings, f"SGD returned no usable locus fields for {query.query_value}."),
            attempts=attempts,
            error_type="no_records",
        )
    return ExternalFetchResult(
        source_database="sgd",
        query=query,
        success=True,
        records=(record,),
        source_url=url,
        retrieved_at=record.provenance.retrieved_at,
        http_status=None if response is None else response.status_code,
        raw_record_sha256=record.provenance.raw_record_sha256,
        warnings=warnings,
        attempts=attempts,
    )


def build_sgd_url(query: ExternalReferenceQuery) -> str:
    return f"{SGD_LOCUS_URL}/{urllib.parse.quote(query.query_value.strip())}"


def _record_from_sgd(
    payload: Mapping[str, Any],
    query: ExternalReferenceQuery,
    source_url: str,
    response: ExternalHttpResponse | None,
) -> ExternalReferenceRecord:
    gene_name = str(payload.get("format_name") or payload.get("display_name") or payload.get("gene_name") or "")
    locus_tag = str(payload.get("systematic_name") or payload.get("locus_tag") or "")
    accession = str(payload.get("sgdid") or payload.get("id") or locus_tag or gene_name)
    aliases = _aliases(payload.get("aliases"))
    retrieved_at = utc_now_iso()
    return ExternalReferenceRecord(
        provenance=ExternalReferenceProvenance(
            source_database="sgd",
            source_version=_header(response.headers if response else {}, "last-modified") or SGD_SOURCE_VERSION,
            source_url=source_url,
            source_query=query.query_value,
            retrieved_at=retrieved_at,
            raw_record_sha256="" if response is None else response.raw_record_sha256,
            warnings=query.warnings,
        ),
        taxon_id="559292",
        organism="Saccharomyces cerevisiae",
        primary_accession=accession,
        gene_id=locus_tag or gene_name or None,
        gene_name=gene_name or None,
        locus_tag=locus_tag or None,
        aliases=_dedupe(aliases, exclude=(gene_name, locus_tag, accession)),
        protein_name=str(payload.get("name_description") or payload.get("description") or "") or None,
    )


def _aliases(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    result: list[str] = []
    for item in value:
        if isinstance(item, Mapping):
            result.append(str(item.get("display_name") or item.get("name") or ""))
        else:
            result.append(str(item))
    return [item.strip() for item in result if item.strip()]


def _dedupe(values: list[str], *, exclude: tuple[str, ...] = ()) -> tuple[str, ...]:
    excluded = {value for value in exclude if value}
    result: list[str] = []
    for value in values:
        text = str(value).strip()
        if text and text not in excluded and text not in result:
            result.append(text)
    return tuple(result)


def _header(headers: Mapping[str, str], key: str) -> str:
    target = key.lower()
    for header_key, value in headers.items():
        if header_key.lower() == target:
            return str(value)
    return ""


__all__ = [
    "SGD_LOCUS_URL",
    "SgdReferenceClient",
    "build_sgd_url",
    "fetch_sgd_reference",
]
