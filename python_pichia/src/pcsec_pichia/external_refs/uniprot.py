from __future__ import annotations

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


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_SOURCE_VERSION = "uniprot-rest"


@dataclass(frozen=True)
class UniProtReferenceClient:
    source_database: str = "uniprot"

    def fetch(
        self,
        query: ExternalReferenceQuery,
        config: ExternalFetchConfig,
        *,
        http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
        sleep: Callable[[float], None] | None = None,
    ) -> ExternalFetchResult:
        return fetch_uniprot_reference(
            query,
            config,
            http_get=http_get,
            sleep=sleep,
        )


def fetch_uniprot_reference(
    query: ExternalReferenceQuery,
    config: ExternalFetchConfig | None = None,
    *,
    http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
    sleep: Callable[[float], None] | None = None,
) -> ExternalFetchResult:
    resolved_config = config or ExternalFetchConfig()
    url = build_uniprot_url(query)
    payload, response, warnings, attempts = request_json(
        url,
        resolved_config,
        http_get=http_get,
        sleep=_sleep(sleep),
    )
    if not isinstance(payload, Mapping):
        return failed_fetch_result(
            source_database="uniprot",
            query=query,
            source_url=url,
            response=response,
            warnings=warnings or ("UniProt returned no JSON object.",),
            attempts=attempts,
            error_type="invalid_response",
        )
    results = _as_list(payload.get("results"))
    if not results:
        return failed_fetch_result(
            source_database="uniprot",
            query=query,
            source_url=url,
            response=response,
            warnings=(*warnings, f"UniProt returned no records for {query.query_value}."),
            attempts=attempts,
            error_type="no_records",
        )
    record = _record_from_uniprot(results[0], query, url, response)
    return ExternalFetchResult(
        source_database="uniprot",
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


def build_uniprot_url(query: ExternalReferenceQuery) -> str:
    term = query.query_value.strip()
    if query.query_type == "external_accession":
        query_text = f'accession:"{term}"'
    elif query.query_type == "sce_homolog":
        query_text = f'(gene_exact:"{term}" OR xref:SGD-{term})'
    else:
        query_text = f'(gene_exact:"{term}" OR gene:"{term}") AND (organism_id:4922 OR organism_id:644223)'
    params = {
        "query": query_text,
        "format": "json",
        "size": "1",
        "fields": "accession,id,gene_names,organism_name,organism_id,protein_name,reviewed",
    }
    return f"{UNIPROT_SEARCH_URL}?{urllib.parse.urlencode(params)}"


def _record_from_uniprot(
    payload: object,
    query: ExternalReferenceQuery,
    source_url: str,
    response: ExternalHttpResponse | None,
) -> ExternalReferenceRecord:
    record = _as_mapping(payload)
    genes = _as_list(record.get("genes"))
    first_gene = _as_mapping(genes[0]) if genes else {}
    gene_name = _field_value(first_gene.get("geneName")) or str(record.get("uniProtkbId") or "")
    locus_tag = _first_field_value(first_gene, ("orderedLocusNames", "orfNames"))
    aliases: list[str] = []
    for gene in genes:
        gene_payload = _as_mapping(gene)
        aliases.extend(_field_values(gene_payload.get("synonyms")))
        aliases.extend(_field_values(gene_payload.get("orderedLocusNames")))
        aliases.extend(_field_values(gene_payload.get("orfNames")))
    organism = _as_mapping(record.get("organism"))
    protein_name = _protein_name(record.get("proteinDescription"))
    retrieved_at = utc_now_iso()
    return ExternalReferenceRecord(
        provenance=ExternalReferenceProvenance(
            source_database="uniprot",
            source_version=_header(response.headers if response else {}, "x-uniprot-release") or UNIPROT_SOURCE_VERSION,
            source_url=source_url,
            source_query=query.query_value,
            retrieved_at=retrieved_at,
            raw_record_sha256="" if response is None else response.raw_record_sha256,
            warnings=query.warnings,
        ),
        taxon_id=str(organism.get("taxonId") or ""),
        organism=str(organism.get("scientificName") or organism.get("commonName") or ""),
        primary_accession=str(record.get("primaryAccession") or ""),
        gene_id=locus_tag or gene_name or None,
        gene_name=gene_name or None,
        locus_tag=locus_tag or None,
        aliases=_dedupe(aliases, exclude=(gene_name, locus_tag)),
        protein_name=protein_name or None,
        reviewed="reviewed" in str(record.get("entryType") or "").lower(),
    )


def _protein_name(value: object) -> str:
    payload = _as_mapping(value)
    recommended = _as_mapping(payload.get("recommendedName"))
    full_name = _as_mapping(recommended.get("fullName"))
    return str(full_name.get("value") or "")


def _as_mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _as_list(value: object) -> list[Any]:
    return list(value) if isinstance(value, list) else []


def _field_value(value: object) -> str:
    return str(_as_mapping(value).get("value") or "")


def _field_values(value: object) -> list[str]:
    return [_field_value(item) for item in _as_list(value) if _field_value(item)]


def _first_field_value(payload: Mapping[str, Any], keys: tuple[str, ...]) -> str:
    for key in keys:
        values = _field_values(payload.get(key))
        if values:
            return values[0]
    return ""


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


def _sleep(value: Callable[[float], None] | None) -> Callable[[float], None]:
    if value is not None:
        return value
    import time

    return time.sleep


__all__ = [
    "UNIPROT_SEARCH_URL",
    "UniProtReferenceClient",
    "build_uniprot_url",
    "fetch_uniprot_reference",
]
