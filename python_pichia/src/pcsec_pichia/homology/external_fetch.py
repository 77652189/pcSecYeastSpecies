from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable, Mapping

from pcsec_pichia.homology.cache_schema import ExternalNameReference


DEFAULT_USER_AGENT = "pcSecYeastSpecies-homology-cache/0.1"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
NCBI_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SGD_LOCUS_URL = "https://www.yeastgenome.org/backend/locus"


@dataclass(frozen=True)
class ExternalFetchConfig:
    timeout_seconds: float = 10.0
    retry_count: int = 1
    delay_seconds: float = 0.34
    user_agent: str = DEFAULT_USER_AGENT
    ncbi_email: str = ""
    ncbi_api_key: str = ""
    ncbi_tool: str = "pcSecYeastSpecies"
    enabled_sources: tuple[str, ...] = ("uniprot", "ncbi", "sgd")

    @classmethod
    def from_env(cls, **overrides: object) -> "ExternalFetchConfig":
        values: dict[str, object] = {
            "ncbi_email": os.environ.get("NCBI_EMAIL", ""),
            "ncbi_api_key": os.environ.get("NCBI_API_KEY", ""),
            "ncbi_tool": os.environ.get("NCBI_TOOL", "pcSecYeastSpecies"),
            "user_agent": os.environ.get("PCSEC_HOMOLOGY_USER_AGENT", DEFAULT_USER_AGENT),
        }
        values.update({key: value for key, value in overrides.items() if value is not None})
        return cls(**values)


@dataclass(frozen=True)
class HttpResponse:
    status_code: int
    body: str
    headers: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class ExternalFetchResult:
    source_database: str
    query: str
    success: bool
    references: tuple[ExternalNameReference, ...] = ()
    warnings: tuple[str, ...] = ()
    raw_status_code: int | None = None
    error_summary: str = ""


HttpGetter = Callable[[str, ExternalFetchConfig], HttpResponse]
SleepFn = Callable[[float], None]


def fetch_uniprot_name_reference(
    query: str,
    config: ExternalFetchConfig | None = None,
    *,
    http_get: HttpGetter | None = None,
    sleep: SleepFn = time.sleep,
) -> ExternalFetchResult:
    config = config or ExternalFetchConfig()
    if not _source_enabled("uniprot", config):
        return _skipped_result("UniProt", query)
    params = {"query": query, "format": "json", "size": "1"}
    response = _request_json(_url(UNIPROT_SEARCH_URL, params), config, http_get=http_get, sleep=sleep)
    if response.error_summary or response.payload is None:
        return _failed_result("UniProt", query, response)
    results = _list(response.payload.get("results"))
    if not results:
        return ExternalFetchResult(
            source_database="UniProt",
            query=query,
            success=False,
            warnings=(*response.warnings, f"UniProt returned no records for query `{query}`"),
            raw_status_code=response.status_code,
            error_summary="no_records",
        )
    return ExternalFetchResult(
        source_database="UniProt",
        query=query,
        success=True,
        references=(_uniprot_reference(results[0], response.headers),),
        warnings=response.warnings,
        raw_status_code=response.status_code,
    )


def fetch_ncbi_name_reference(
    query: str,
    config: ExternalFetchConfig | None = None,
    *,
    database: str = "gene",
    http_get: HttpGetter | None = None,
    sleep: SleepFn = time.sleep,
) -> ExternalFetchResult:
    config = config or ExternalFetchConfig()
    if not _source_enabled("ncbi", config):
        return _skipped_result("NCBI", query)
    search_params = _ncbi_params(
        config,
        {
            "db": database,
            "term": query,
            "retmode": "json",
        },
    )
    search = _request_json(
        _url(f"{NCBI_EUTILS_BASE_URL}/esearch.fcgi", search_params),
        config,
        http_get=http_get,
        sleep=sleep,
    )
    if search.error_summary or search.payload is None:
        return _failed_result("NCBI", query, search)
    ids = _list(_dict(search.payload.get("esearchresult")).get("idlist"))
    if not ids:
        return ExternalFetchResult(
            source_database="NCBI",
            query=query,
            success=False,
            warnings=(*search.warnings, f"NCBI returned no {database} ids for query `{query}`"),
            raw_status_code=search.status_code,
            error_summary="no_records",
        )

    uid = str(ids[0])
    summary_params = _ncbi_params(
        config,
        {
            "db": database,
            "id": uid,
            "retmode": "json",
        },
    )
    summary = _request_json(
        _url(f"{NCBI_EUTILS_BASE_URL}/esummary.fcgi", summary_params),
        config,
        http_get=http_get,
        sleep=sleep,
    )
    warnings = (*search.warnings, *summary.warnings)
    if summary.error_summary or summary.payload is None:
        return _failed_result("NCBI", query, summary, warnings=warnings)
    record = _dict(_dict(summary.payload.get("result")).get(uid))
    if not record:
        return ExternalFetchResult(
            source_database="NCBI",
            query=query,
            success=False,
            warnings=(*warnings, f"NCBI summary did not include id `{uid}`"),
            raw_status_code=summary.status_code,
            error_summary="missing_summary_record",
        )
    return ExternalFetchResult(
        source_database="NCBI",
        query=query,
        success=True,
        references=(_ncbi_reference(record, uid, database),),
        warnings=warnings,
        raw_status_code=summary.status_code,
    )


def fetch_sgd_name_reference(
    query: str,
    config: ExternalFetchConfig | None = None,
    *,
    http_get: HttpGetter | None = None,
    sleep: SleepFn = time.sleep,
) -> ExternalFetchResult:
    config = config or ExternalFetchConfig()
    if not _source_enabled("sgd", config):
        return _skipped_result("SGD", query)
    url = f"{SGD_LOCUS_URL}/{urllib.parse.quote(query.strip())}"
    response = _request_json(url, config, http_get=http_get, sleep=sleep)
    if response.error_summary or response.payload is None:
        return _failed_result("SGD", query, response)
    return ExternalFetchResult(
        source_database="SGD",
        query=query,
        success=True,
        references=(_sgd_reference(response.payload, response.headers),),
        warnings=response.warnings,
        raw_status_code=response.status_code,
    )


def fetch_external_name_references(
    query: str,
    config: ExternalFetchConfig | None = None,
    *,
    http_get: HttpGetter | None = None,
    sleep: SleepFn = time.sleep,
) -> tuple[ExternalFetchResult, ...]:
    config = config or ExternalFetchConfig()
    results: list[ExternalFetchResult] = []
    if _source_enabled("uniprot", config):
        results.append(fetch_uniprot_name_reference(query, config, http_get=http_get, sleep=sleep))
    if _source_enabled("ncbi", config):
        results.append(fetch_ncbi_name_reference(query, config, http_get=http_get, sleep=sleep))
    if _source_enabled("sgd", config):
        results.append(fetch_sgd_name_reference(query, config, http_get=http_get, sleep=sleep))
    return tuple(results)


@dataclass(frozen=True)
class _JsonResponse:
    payload: dict[str, Any] | None
    status_code: int | None
    headers: Mapping[str, str]
    warnings: tuple[str, ...] = ()
    error_summary: str = ""


def _request_json(
    url: str,
    config: ExternalFetchConfig,
    *,
    http_get: HttpGetter | None,
    sleep: SleepFn,
) -> _JsonResponse:
    getter = http_get or _default_http_get
    attempts = max(0, int(config.retry_count)) + 1
    warnings: list[str] = []
    last_status: int | None = None
    headers: Mapping[str, str] = {}
    for attempt in range(attempts):
        try:
            response = getter(url, config)
            last_status = response.status_code
            headers = response.headers
            if 200 <= response.status_code < 300:
                return _JsonResponse(
                    payload=json.loads(response.body or "{}"),
                    status_code=response.status_code,
                    headers=headers,
                    warnings=tuple(warnings),
                )
            warnings.append(f"HTTP {response.status_code} for {url}")
            if response.status_code not in {408, 429} and response.status_code < 500:
                break
        except (TimeoutError, urllib.error.URLError, OSError, json.JSONDecodeError, ValueError) as exc:
            warnings.append(f"{type(exc).__name__}: {exc}")
        if attempt < attempts - 1 and config.delay_seconds > 0:
            sleep(config.delay_seconds)
    return _JsonResponse(
        payload=None,
        status_code=last_status,
        headers=headers,
        warnings=tuple(warnings),
        error_summary=warnings[-1] if warnings else "request_failed",
    )


def _default_http_get(url: str, config: ExternalFetchConfig) -> HttpResponse:
    request = urllib.request.Request(url, headers={"User-Agent": config.user_agent})
    with urllib.request.urlopen(request, timeout=config.timeout_seconds) as handle:  # noqa: S310 - official HTTPS APIs only
        body = handle.read().decode("utf-8")
        return HttpResponse(
            status_code=int(getattr(handle, "status", 200)),
            body=body,
            headers=dict(handle.headers.items()),
        )


def _uniprot_reference(record: object, headers: Mapping[str, str]) -> ExternalNameReference:
    payload = _dict(record)
    genes = _list(payload.get("genes"))
    first_gene = _dict(genes[0]) if genes else {}
    gene_name = _field_value(first_gene.get("geneName")) or str(payload.get("uniProtkbId") or "")
    locus_tag = _first_field_value(first_gene, ("orderedLocusNames", "orfNames"))
    aliases: list[str] = []
    for gene in genes:
        gene_payload = _dict(gene)
        aliases.extend(_field_values(gene_payload.get("synonyms")))
        aliases.extend(_field_values(gene_payload.get("orderedLocusNames")))
        aliases.extend(_field_values(gene_payload.get("orfNames")))
    return ExternalNameReference(
        source_database="UniProt",
        source_version=_header(headers, "x-uniprot-release"),
        taxon=str(_dict(payload.get("organism")).get("scientificName") or ""),
        accession=str(payload.get("primaryAccession") or ""),
        gene_name=gene_name,
        locus_tag=locus_tag,
        aliases=_dedupe_aliases(aliases, exclude=(gene_name, locus_tag)),
        retrieved_at=_utc_now(),
    )


def _ncbi_reference(record: Mapping[str, Any], uid: str, database: str) -> ExternalNameReference:
    organism = _dict(record.get("organism"))
    gene_name = str(record.get("nomenclaturesymbol") or record.get("name") or record.get("title") or "")
    locus_tag = str(record.get("locus_tag") or record.get("maplocation") or "")
    accession = str(record.get("accessionversion") or record.get("caption") or uid)
    aliases = _split_aliases(record.get("otheraliases"))
    aliases.extend(_split_aliases(record.get("otherdesignations")))
    aliases.extend(_split_aliases(record.get("description")))
    return ExternalNameReference(
        source_database="NCBI",
        source_version=database,
        taxon=str(organism.get("scientificname") or organism.get("name") or ""),
        accession=accession,
        gene_name=gene_name,
        locus_tag=locus_tag,
        aliases=_dedupe_aliases(aliases, exclude=(gene_name, locus_tag, accession)),
        retrieved_at=_utc_now(),
    )


def _sgd_reference(payload: Mapping[str, Any], headers: Mapping[str, str]) -> ExternalNameReference:
    gene_name = str(payload.get("format_name") or payload.get("display_name") or payload.get("gene_name") or "")
    locus_tag = str(payload.get("systematic_name") or payload.get("locus_tag") or "")
    aliases = _list(payload.get("aliases"))
    if aliases and isinstance(aliases[0], Mapping):
        aliases = [str(_dict(item).get("display_name") or _dict(item).get("name") or "") for item in aliases]
    return ExternalNameReference(
        source_database="SGD",
        source_version=_header(headers, "last-modified"),
        taxon=str(payload.get("organism") or "Saccharomyces cerevisiae"),
        accession=str(payload.get("sgdid") or payload.get("id") or ""),
        gene_name=gene_name,
        locus_tag=locus_tag,
        aliases=_dedupe_aliases([str(item) for item in aliases], exclude=(gene_name, locus_tag)),
        retrieved_at=_utc_now(),
    )


def _failed_result(
    source_database: str,
    query: str,
    response: _JsonResponse,
    *,
    warnings: tuple[str, ...] | None = None,
) -> ExternalFetchResult:
    return ExternalFetchResult(
        source_database=source_database,
        query=query,
        success=False,
        warnings=warnings or response.warnings,
        raw_status_code=response.status_code,
        error_summary=response.error_summary or "request_failed",
    )


def _skipped_result(source_database: str, query: str) -> ExternalFetchResult:
    return ExternalFetchResult(
        source_database=source_database,
        query=query,
        success=False,
        warnings=(f"{source_database} source disabled",),
        error_summary="source_disabled",
    )


def _source_enabled(source: str, config: ExternalFetchConfig) -> bool:
    return source.lower() in {item.lower() for item in config.enabled_sources}


def _ncbi_params(config: ExternalFetchConfig, params: dict[str, object]) -> dict[str, object]:
    merged = {"tool": config.ncbi_tool, **params}
    if config.ncbi_email:
        merged["email"] = config.ncbi_email
    if config.ncbi_api_key:
        merged["api_key"] = config.ncbi_api_key
    return merged


def _url(base_url: str, params: Mapping[str, object]) -> str:
    return f"{base_url}?{urllib.parse.urlencode({key: value for key, value in params.items() if value not in {None, ''}})}"


def _dict(value: object) -> dict[str, Any]:
    return dict(value) if isinstance(value, Mapping) else {}


def _list(value: object) -> list[Any]:
    return list(value) if isinstance(value, list) else []


def _field_value(value: object) -> str:
    return str(_dict(value).get("value") or "")


def _field_values(value: object) -> list[str]:
    return [_field_value(item) for item in _list(value) if _field_value(item)]


def _first_field_value(payload: Mapping[str, Any], keys: tuple[str, ...]) -> str:
    for key in keys:
        values = _field_values(payload.get(key))
        if values:
            return values[0]
    return ""


def _split_aliases(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).replace(";", ",")
    return [part.strip() for part in text.split(",") if part.strip()]


def _dedupe_aliases(values: list[str], *, exclude: tuple[str, ...] = ()) -> tuple[str, ...]:
    excluded = {value.strip() for value in exclude if value.strip()}
    aliases: list[str] = []
    for value in values:
        text = str(value).strip()
        if text and text not in excluded and text not in aliases:
            aliases.append(text)
    return tuple(aliases)


def _header(headers: Mapping[str, str], key: str) -> str:
    target = key.lower()
    for header_key, value in headers.items():
        if header_key.lower() == target:
            return str(value)
    return ""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


__all__ = [
    "ExternalFetchConfig",
    "ExternalFetchResult",
    "HttpResponse",
    "fetch_external_name_references",
    "fetch_ncbi_name_reference",
    "fetch_sgd_name_reference",
    "fetch_uniprot_name_reference",
]
