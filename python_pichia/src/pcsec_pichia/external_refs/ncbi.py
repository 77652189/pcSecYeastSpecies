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


NCBI_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_SOURCE_VERSION = "ncbi-entrez-gene"


@dataclass(frozen=True)
class NcbiGeneReferenceClient:
    source_database: str = "ncbi"

    def fetch(
        self,
        query: ExternalReferenceQuery,
        config: ExternalFetchConfig,
        *,
        http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> ExternalFetchResult:
        return fetch_ncbi_gene_reference(query, config, http_get=http_get, sleep=sleep)


def fetch_ncbi_gene_reference(
    query: ExternalReferenceQuery,
    config: ExternalFetchConfig | None = None,
    *,
    http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
    sleep: Callable[[float], None] = time.sleep,
) -> ExternalFetchResult:
    resolved_config = config or ExternalFetchConfig()
    search_url = build_ncbi_esearch_url(query, resolved_config)
    search_payload, search_response, search_warnings, search_attempts = request_json(
        search_url,
        resolved_config,
        http_get=http_get,
        sleep=sleep,
    )
    if not isinstance(search_payload, Mapping):
        return failed_fetch_result(
            source_database="ncbi",
            query=query,
            source_url=search_url,
            response=search_response,
            warnings=search_warnings or ("NCBI esearch returned no JSON object.",),
            attempts=search_attempts,
            error_type="invalid_response",
        )
    ids = _as_list(_as_mapping(search_payload.get("esearchresult")).get("idlist"))
    if not ids:
        return failed_fetch_result(
            source_database="ncbi",
            query=query,
            source_url=search_url,
            response=search_response,
            warnings=(*search_warnings, f"NCBI returned no Gene IDs for {query.query_value}."),
            attempts=search_attempts,
            error_type="no_records",
        )
    gene_uid = str(ids[0])
    if resolved_config.min_interval_seconds > 0:
        sleep(resolved_config.min_interval_seconds)
    summary_url = build_ncbi_esummary_url(gene_uid, resolved_config)
    summary_payload, summary_response, summary_warnings, summary_attempts = request_json(
        summary_url,
        resolved_config,
        http_get=http_get,
        sleep=sleep,
    )
    warnings = (*search_warnings, *summary_warnings)
    if not isinstance(summary_payload, Mapping):
        return failed_fetch_result(
            source_database="ncbi",
            query=query,
            source_url=summary_url,
            response=summary_response,
            warnings=warnings or ("NCBI esummary returned no JSON object.",),
            attempts=search_attempts + summary_attempts,
            error_type="invalid_response",
        )
    record_payload = _as_mapping(_as_mapping(summary_payload.get("result")).get(gene_uid))
    if not record_payload:
        return failed_fetch_result(
            source_database="ncbi",
            query=query,
            source_url=summary_url,
            response=summary_response,
            warnings=(*warnings, f"NCBI summary did not include Gene ID {gene_uid}."),
            attempts=search_attempts + summary_attempts,
            error_type="missing_summary_record",
        )
    record = _record_from_ncbi(record_payload, gene_uid, query, summary_url, summary_response)
    return ExternalFetchResult(
        source_database="ncbi",
        query=query,
        success=True,
        records=(record,),
        source_url=summary_url,
        retrieved_at=record.provenance.retrieved_at,
        http_status=None if summary_response is None else summary_response.status_code,
        raw_record_sha256=record.provenance.raw_record_sha256,
        warnings=warnings,
        attempts=search_attempts + summary_attempts,
    )


def build_ncbi_esearch_url(query: ExternalReferenceQuery, config: ExternalFetchConfig) -> str:
    params = _ncbi_params(
        config,
        {
            "db": "gene",
            "term": query.query_value.strip(),
            "retmode": "json",
        },
    )
    return f"{NCBI_EUTILS_BASE_URL}/esearch.fcgi?{urllib.parse.urlencode(params)}"


def build_ncbi_esummary_url(gene_uid: str, config: ExternalFetchConfig) -> str:
    params = _ncbi_params(
        config,
        {
            "db": "gene",
            "id": gene_uid,
            "retmode": "json",
        },
    )
    return f"{NCBI_EUTILS_BASE_URL}/esummary.fcgi?{urllib.parse.urlencode(params)}"


def _record_from_ncbi(
    payload: Mapping[str, Any],
    gene_uid: str,
    query: ExternalReferenceQuery,
    source_url: str,
    response: ExternalHttpResponse | None,
) -> ExternalReferenceRecord:
    organism = _as_mapping(payload.get("organism"))
    gene_name = str(payload.get("nomenclaturesymbol") or payload.get("name") or payload.get("title") or "")
    locus_tag = str(payload.get("locus_tag") or payload.get("maplocation") or "")
    accession = str(payload.get("accessionversion") or payload.get("caption") or gene_uid)
    aliases = _split_aliases(payload.get("otheraliases"))
    aliases.extend(_split_aliases(payload.get("otherdesignations")))
    retrieved_at = utc_now_iso()
    return ExternalReferenceRecord(
        provenance=ExternalReferenceProvenance(
            source_database="ncbi",
            source_version=NCBI_SOURCE_VERSION,
            source_url=source_url,
            source_query=query.query_value,
            retrieved_at=retrieved_at,
            raw_record_sha256="" if response is None else response.raw_record_sha256,
            warnings=query.warnings,
        ),
        taxon_id=str(organism.get("taxid") or organism.get("tax_id") or ""),
        organism=str(organism.get("scientificname") or organism.get("name") or ""),
        primary_accession=accession,
        gene_id=str(payload.get("uid") or gene_uid),
        gene_name=gene_name or None,
        locus_tag=locus_tag or None,
        aliases=_dedupe(aliases, exclude=(gene_name, locus_tag, accession)),
        protein_name=str(payload.get("description") or "") or None,
    )


def _ncbi_params(config: ExternalFetchConfig, params: Mapping[str, object]) -> dict[str, object]:
    merged: dict[str, object] = {"tool": config.ncbi_tool, **dict(params)}
    email = config.ncbi_email()
    api_key = config.ncbi_api_key()
    if email:
        merged["email"] = email
    if api_key:
        merged["api_key"] = api_key
    return {key: value for key, value in merged.items() if value not in ("", None)}


def _as_mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _as_list(value: object) -> list[Any]:
    return list(value) if isinstance(value, list) else []


def _split_aliases(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).replace(";", ",")
    return [part.strip() for part in text.split(",") if part.strip()]


def _dedupe(values: list[str], *, exclude: tuple[str, ...] = ()) -> tuple[str, ...]:
    excluded = {value for value in exclude if value}
    result: list[str] = []
    for value in values:
        text = str(value).strip()
        if text and text not in excluded and text not in result:
            result.append(text)
    return tuple(result)


__all__ = [
    "NCBI_EUTILS_BASE_URL",
    "NcbiGeneReferenceClient",
    "build_ncbi_esearch_url",
    "build_ncbi_esummary_url",
    "fetch_ncbi_gene_reference",
]
