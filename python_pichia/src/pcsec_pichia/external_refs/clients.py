from __future__ import annotations

import os
import json
import time
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Iterable, Mapping, Protocol, runtime_checkable

import httpx

from pcsec_pichia.external_refs.queries import ExternalReferenceQuery
from pcsec_pichia.external_refs.schema import ExternalReferenceRecord, sha256_text, utc_now_iso


DEFAULT_USER_AGENT = "pcSecPichia-external-reference-cache/0.1"
RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}


@dataclass(frozen=True)
class ExternalFetchConfig:
    sources: tuple[str, ...] = ("uniprot", "sgd", "ncbi")
    timeout_seconds: float = 20.0
    retry_attempts: int = 3
    min_interval_seconds: float = 0.2
    user_agent: str = DEFAULT_USER_AGENT
    offline_cache_dir: str | None = None
    ncbi_api_key_env: str = "NCBI_API_KEY"
    ncbi_email_env: str = "NCBI_EMAIL"
    ncbi_tool: str = "pcSecPichia"

    @classmethod
    def from_env(cls, **overrides: Any) -> "ExternalFetchConfig":
        values: dict[str, Any] = {
            "user_agent": os.environ.get("PCSEC_EXTERNAL_REFS_USER_AGENT", DEFAULT_USER_AGENT),
            "ncbi_api_key_env": os.environ.get("PCSEC_NCBI_API_KEY_ENV", "NCBI_API_KEY"),
            "ncbi_email_env": os.environ.get("PCSEC_NCBI_EMAIL_ENV", "NCBI_EMAIL"),
            "ncbi_tool": os.environ.get("NCBI_TOOL", "pcSecPichia"),
        }
        values.update({key: value for key, value in overrides.items() if value is not None})
        return cls(**values)

    def enabled(self, source_database: str) -> bool:
        return source_database.lower() in {source.lower() for source in self.sources}

    def ncbi_api_key(self) -> str:
        return os.environ.get(self.ncbi_api_key_env, "")

    def ncbi_email(self) -> str:
        return os.environ.get(self.ncbi_email_env, "")


@dataclass(frozen=True)
class ExternalHttpResponse:
    status_code: int
    text: str
    url: str
    headers: Mapping[str, str] = field(default_factory=dict)

    @property
    def raw_record_sha256(self) -> str:
        return sha256_text(self.text)

    def json(self) -> Mapping[str, Any] | list[Any]:
        payload = json.loads(self.text or "{}")
        if isinstance(payload, (dict, list)):
            return payload
        raise ValueError("External HTTP response JSON must be an object or array.")


@dataclass(frozen=True)
class ExternalFetchFailure:
    source_database: str
    source_query: str
    source_url: str
    retrieved_at: str
    http_status: int | None
    raw_record_sha256: str | None
    error_type: str
    error_message: str
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalFetchResult:
    source_database: str
    query: ExternalReferenceQuery
    success: bool
    records: tuple[ExternalReferenceRecord, ...] = ()
    failure: ExternalFetchFailure | None = None
    source_url: str = ""
    retrieved_at: str = ""
    http_status: int | None = None
    raw_record_sha256: str | None = None
    warnings: tuple[str, ...] = ()
    attempts: int = 1

    @property
    def failed(self) -> bool:
        return not self.success

    def to_failure_dict(self) -> dict[str, Any] | None:
        if self.failure is None:
            return None
        payload = self.failure.to_dict()
        payload["query"] = self.query.to_dict()
        payload["attempts"] = self.attempts
        return payload


@runtime_checkable
class ExternalReferenceClient(Protocol):
    source_database: str

    def fetch(
        self,
        query: ExternalReferenceQuery,
        config: ExternalFetchConfig,
        *,
        http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> ExternalFetchResult:
        """Fetch external reference records for one query."""


def fetch_external_references(
    queries: Iterable[ExternalReferenceQuery],
    clients: Iterable[ExternalReferenceClient],
    config: ExternalFetchConfig | None = None,
    *,
    http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
    sleep: Callable[[float], None] = time.sleep,
) -> tuple[ExternalFetchResult, ...]:
    resolved_config = config or ExternalFetchConfig()
    enabled_clients = tuple(
        client
        for client in clients
        if resolved_config.enabled(client.source_database)
    )
    results: list[ExternalFetchResult] = []
    dispatched = 0
    for query in queries:
        preferred = set(query.preferred_sources)
        selected = tuple(
            client
            for client in enabled_clients
            if not preferred or client.source_database.lower() in preferred
        )
        for client in selected:
            if dispatched and resolved_config.min_interval_seconds > 0:
                sleep(resolved_config.min_interval_seconds)
            results.append(client.fetch(query, resolved_config, http_get=http_get, sleep=sleep))
            dispatched += 1
    return tuple(results)


def default_external_reference_clients() -> tuple[ExternalReferenceClient, ...]:
    from pcsec_pichia.external_refs.ncbi import NcbiGeneReferenceClient
    from pcsec_pichia.external_refs.sgd import SgdReferenceClient
    from pcsec_pichia.external_refs.uniprot import UniProtReferenceClient

    return (
        UniProtReferenceClient(),
        SgdReferenceClient(),
        NcbiGeneReferenceClient(),
    )


def request_json(
    url: str,
    config: ExternalFetchConfig,
    *,
    http_get: Callable[[str, ExternalFetchConfig], ExternalHttpResponse] | None = None,
    sleep: Callable[[float], None] = time.sleep,
) -> tuple[Mapping[str, Any] | list[Any] | None, ExternalHttpResponse | None, tuple[str, ...], int]:
    getter = http_get or default_http_get
    attempts = max(1, int(config.retry_attempts))
    warnings: list[str] = []
    last_response: ExternalHttpResponse | None = None
    for attempt in range(1, attempts + 1):
        if attempt > 1 and config.min_interval_seconds > 0:
            sleep(config.min_interval_seconds)
        try:
            response = getter(url, config)
            last_response = response
            if 200 <= response.status_code < 300:
                return response.json(), response, tuple(warnings), attempt
            warnings.append(f"HTTP {response.status_code} for {url}")
            if response.status_code not in RETRYABLE_STATUS_CODES:
                break
        except (httpx.HTTPError, TimeoutError, OSError, ValueError) as exc:
            warnings.append(f"{type(exc).__name__}: {exc}")
    return None, last_response, tuple(warnings), attempts


def default_http_get(url: str, config: ExternalFetchConfig) -> ExternalHttpResponse:
    with httpx.Client(timeout=config.timeout_seconds, headers={"User-Agent": config.user_agent}) as client:
        response = client.get(url)
    return ExternalHttpResponse(
        status_code=int(response.status_code),
        text=response.text,
        url=str(response.url),
        headers=dict(response.headers),
    )


def failed_fetch_result(
    *,
    source_database: str,
    query: ExternalReferenceQuery,
    source_url: str,
    response: ExternalHttpResponse | None,
    warnings: tuple[str, ...],
    attempts: int,
    error_type: str = "request_failed",
    error_message: str = "",
) -> ExternalFetchResult:
    retrieved_at = utc_now_iso()
    failure = ExternalFetchFailure(
        source_database=source_database,
        source_query=query.query_value,
        source_url=source_url,
        retrieved_at=retrieved_at,
        http_status=None if response is None else response.status_code,
        raw_record_sha256=None if response is None else response.raw_record_sha256,
        error_type=error_type,
        error_message=error_message or (warnings[-1] if warnings else error_type),
        warnings=warnings,
    )
    return ExternalFetchResult(
        source_database=source_database,
        query=query,
        success=False,
        failure=failure,
        source_url=source_url,
        retrieved_at=retrieved_at,
        http_status=failure.http_status,
        raw_record_sha256=failure.raw_record_sha256,
        warnings=warnings,
        attempts=attempts,
    )


__all__ = [
    "DEFAULT_USER_AGENT",
    "ExternalFetchConfig",
    "ExternalFetchFailure",
    "ExternalFetchResult",
    "ExternalHttpResponse",
    "ExternalReferenceClient",
    "default_http_get",
    "default_external_reference_clients",
    "failed_fetch_result",
    "fetch_external_references",
    "request_json",
]
