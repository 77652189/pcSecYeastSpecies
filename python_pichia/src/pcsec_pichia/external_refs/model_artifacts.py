from __future__ import annotations

import hashlib
import json
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable

from pcsec_pichia.external_refs.model_inventory import ExternalModelInventoryRecord
from pcsec_pichia.external_refs.schema import utc_now_iso


ARTIFACT_MANIFEST_FILENAME = "external_model_artifact_manifest.json"
ARTIFACT_FAILURES_FILENAME = "external_model_artifact_failures.jsonl"

ArtifactFetcher = Callable[[str, float], bytes]


@dataclass(frozen=True)
class ExternalModelArtifactRequest:
    model_id: str
    artifact_url: str
    artifact_type: str
    filename: str
    expected_sha256: str = ""
    source_page_url: str = ""
    requires_manual_access: bool = False
    notes: str = ""
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        for field_name in ("model_id", "artifact_type", "filename"):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"{field_name} must be non-empty.")
        if not self.requires_manual_access and not self.artifact_url:
            raise ValueError("artifact_url is required unless requires_manual_access is true.")
        if self.expected_sha256 and (
            len(self.expected_sha256) != 64
            or any(char not in "0123456789abcdef" for char in self.expected_sha256.lower())
        ):
            raise ValueError("expected_sha256 must be a 64-character hex digest.")


@dataclass(frozen=True)
class ExternalModelArtifactResult:
    model_id: str
    artifact_url: str
    artifact_type: str
    filename: str
    download_status: str
    local_path: str
    checksum_sha256: str
    bytes_written: int
    source_page_url: str = ""
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class ExternalModelArtifactManifest:
    generated_at: str
    request_count: int
    downloaded_count: int
    failed_count: int
    manual_required_count: int
    checksum_mismatch_count: int
    results: tuple[ExternalModelArtifactResult, ...]


@dataclass(frozen=True)
class ExternalModelArtifactOutputs:
    manifest_path: Path
    failures_path: Path
    manifest: ExternalModelArtifactManifest


def build_artifact_requests_from_inventory(
    records: Iterable[ExternalModelInventoryRecord],
) -> tuple[ExternalModelArtifactRequest, ...]:
    """Convert inventory rows into explicit artifact requests.

    Round B does not infer direct downloadable files from publication or
    repository landing pages. Those rows are kept as manual-required records
    until a direct artifact URL is curated.
    """

    requests: list[ExternalModelArtifactRequest] = []
    for record in records:
        artifact_type = _preferred_artifact_type(record)
        requests.append(
            ExternalModelArtifactRequest(
                model_id=record.model_id,
                artifact_url="",
                artifact_type=artifact_type,
                filename=_default_filename(record.model_id, artifact_type),
                source_page_url=record.source_url or record.publication_url,
                requires_manual_access=True,
                notes=record.notes,
                warnings=tuple(dict.fromkeys((*record.warnings, "direct_artifact_url_required"))),
            )
        )
    return tuple(requests)


def cache_external_model_artifacts(
    requests: Iterable[ExternalModelArtifactRequest],
    output_dir: Path,
    *,
    fetcher: ArtifactFetcher | None = None,
    timeout_seconds: float = 30.0,
) -> ExternalModelArtifactOutputs:
    resolved_requests = tuple(requests)
    for request in resolved_requests:
        request.validate()
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved_fetcher = fetcher or _fetch_url
    results = tuple(
        _cache_one_artifact(
            request,
            output_dir=output_dir,
            fetcher=resolved_fetcher,
            timeout_seconds=timeout_seconds,
        )
        for request in resolved_requests
    )
    manifest = _manifest_from_results(results)
    manifest_path = output_dir / ARTIFACT_MANIFEST_FILENAME
    failures_path = output_dir / ARTIFACT_FAILURES_FILENAME
    manifest_path.write_text(
        json.dumps(_json_ready(asdict(manifest)), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_failures(results, failures_path)
    return ExternalModelArtifactOutputs(
        manifest_path=manifest_path,
        failures_path=failures_path,
        manifest=manifest,
    )


def _cache_one_artifact(
    request: ExternalModelArtifactRequest,
    *,
    output_dir: Path,
    fetcher: ArtifactFetcher,
    timeout_seconds: float,
) -> ExternalModelArtifactResult:
    if request.requires_manual_access:
        return ExternalModelArtifactResult(
            model_id=request.model_id,
            artifact_url=request.artifact_url,
            artifact_type=request.artifact_type,
            filename=request.filename,
            download_status="manual_download_required",
            local_path="",
            checksum_sha256="",
            bytes_written=0,
            source_page_url=request.source_page_url,
            warnings=tuple(dict.fromkeys((*request.warnings, "manual_download_required"))),
        )
    try:
        payload = fetcher(request.artifact_url, timeout_seconds)
    except Exception as exc:  # pragma: no cover - exact network errors are environment-dependent.
        return ExternalModelArtifactResult(
            model_id=request.model_id,
            artifact_url=request.artifact_url,
            artifact_type=request.artifact_type,
            filename=request.filename,
            download_status="download_failed",
            local_path="",
            checksum_sha256="",
            bytes_written=0,
            source_page_url=request.source_page_url,
            warnings=tuple(dict.fromkeys((*request.warnings, f"{type(exc).__name__}: {exc}"))),
        )
    checksum = hashlib.sha256(payload).hexdigest()
    if request.expected_sha256 and checksum.lower() != request.expected_sha256.lower():
        return ExternalModelArtifactResult(
            model_id=request.model_id,
            artifact_url=request.artifact_url,
            artifact_type=request.artifact_type,
            filename=request.filename,
            download_status="checksum_mismatch",
            local_path="",
            checksum_sha256=checksum,
            bytes_written=0,
            source_page_url=request.source_page_url,
            warnings=tuple(dict.fromkeys((*request.warnings, "checksum_mismatch"))),
        )
    artifact_path = output_dir / "artifacts" / _safe_path_part(request.model_id) / _safe_path_part(request.filename)
    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    artifact_path.write_bytes(payload)
    return ExternalModelArtifactResult(
        model_id=request.model_id,
        artifact_url=request.artifact_url,
        artifact_type=request.artifact_type,
        filename=request.filename,
        download_status="downloaded",
        local_path=str(artifact_path),
        checksum_sha256=checksum,
        bytes_written=len(payload),
        source_page_url=request.source_page_url,
        warnings=request.warnings,
    )


def _manifest_from_results(
    results: tuple[ExternalModelArtifactResult, ...],
) -> ExternalModelArtifactManifest:
    downloaded = sum(1 for result in results if result.download_status == "downloaded")
    manual = sum(1 for result in results if result.download_status == "manual_download_required")
    checksum_mismatch = sum(1 for result in results if result.download_status == "checksum_mismatch")
    failed = sum(1 for result in results if result.download_status in {"download_failed", "checksum_mismatch"})
    return ExternalModelArtifactManifest(
        generated_at=utc_now_iso(),
        request_count=len(results),
        downloaded_count=downloaded,
        failed_count=failed,
        manual_required_count=manual,
        checksum_mismatch_count=checksum_mismatch,
        results=results,
    )


def _write_failures(results: tuple[ExternalModelArtifactResult, ...], path: Path) -> None:
    failures = tuple(result for result in results if result.download_status != "downloaded")
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for result in failures:
            handle.write(json.dumps(_json_ready(asdict(result)), ensure_ascii=False, sort_keys=True))
            handle.write("\n")


def _fetch_url(url: str, timeout_seconds: float) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "pcSecPichia external model artifact cache"})
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:  # noqa: S310 - explicit user-triggered cache builder.
        return response.read()


def _preferred_artifact_type(record: ExternalModelInventoryRecord) -> str:
    for candidate in ("SBML", "model archive", "MATLAB", "repository"):
        if candidate in record.available_artifact_types:
            return candidate
    return record.available_artifact_types[0] if record.available_artifact_types else "manual"


def _default_filename(model_id: str, artifact_type: str) -> str:
    extension = {
        "SBML": ".xml",
        "MATLAB": ".mat",
        "model archive": ".zip",
        "repository": ".zip",
    }.get(artifact_type, ".artifact")
    return f"{_safe_path_part(model_id)}{extension}"


def _safe_path_part(value: str) -> str:
    cleaned = [character if character.isalnum() or character in {"-", "_", "."} else "_" for character in value]
    return "".join(cleaned).strip("._") or "artifact"


def _json_ready(value: object) -> object:
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "ARTIFACT_FAILURES_FILENAME",
    "ARTIFACT_MANIFEST_FILENAME",
    "ArtifactFetcher",
    "ExternalModelArtifactManifest",
    "ExternalModelArtifactOutputs",
    "ExternalModelArtifactRequest",
    "ExternalModelArtifactResult",
    "build_artifact_requests_from_inventory",
    "cache_external_model_artifacts",
]
