from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import asdict, dataclass, replace
from enum import Enum
from pathlib import Path
from typing import Any, Mapping

from pcsec_pichia.external_refs.clients import (
    ExternalFetchConfig,
    ExternalHttpResponse,
    default_http_get,
)
from pcsec_pichia.external_refs.queries import ExternalReferenceQuery
from pcsec_pichia.external_refs.uniprot import build_uniprot_url, fetch_uniprot_reference
from pcsec_pichia.errors import OECapacityValidationError


ExternalCapacitySourceValidationError = OECapacityValidationError


class ExternalCapacitySourceType(str, Enum):
    QUANTITATIVE_PROTEOMICS = "quantitative_proteomics"
    EXTERNAL_ENZYME_MODEL = "external_enzyme_model"
    KINETICS_DATABASE = "kinetics_database"
    LITERATURE = "literature"
    IDENTITY_REFERENCE = "identity_reference"


class RetrievalMode(str, Enum):
    ONLINE = "online"
    MANUAL_IMPORT = "manual_import"
    OFFLINE_REPLAY = "offline_replay"


@dataclass(frozen=True)
class ExternalCapacitySource:
    source_id: str
    source_type: ExternalCapacitySourceType
    source_version: str
    source_url: str
    retrieved_at: str
    query: str
    raw_sha256: str
    license_id: str
    retrieval_mode: RetrievalMode
    cache_path: str
    license_url: str = ""
    adapter_id: str = "pcsec_pichia.external_capacity"
    adapter_version: str = "1"
    terms_reviewed: bool = False
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        for name in (
            "source_id",
            "source_version",
            "retrieved_at",
            "query",
            "raw_sha256",
            "license_id",
            "cache_path",
            "adapter_id",
            "adapter_version",
        ):
            _require_text(getattr(self, name), name)
        _require_sha256(self.raw_sha256, "raw_sha256")


def cache_uniprot_identity_source(
    gene_id: str,
    output_dir: str | Path,
    *,
    http_get: Any = None,
    retrieval_mode: RetrievalMode = RetrievalMode.ONLINE,
) -> ExternalCapacitySource:
    query = ExternalReferenceQuery(
        query_type="pichia_gene",
        query_value=gene_id,
        source_context="oe_capacity_round_6a",
        source_id="g6pdh2_identity",
        preferred_sources=("uniprot",),
    )
    root = Path(output_dir)
    raw_dir = root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_path = raw_dir / f"uniprot-{_safe_name(gene_id)}.json"
    metadata_path = raw_dir / f"uniprot-{_safe_name(gene_id)}.source.json"
    config = ExternalFetchConfig(sources=("uniprot",))
    cached_source: ExternalCapacitySource | None = None
    if retrieval_mode is RetrievalMode.OFFLINE_REPLAY:
        if not raw_path.is_file():
            raise _validation_error(f"offline UniProt raw cache is missing: {raw_path}")
        if not metadata_path.is_file():
            raise _validation_error(
                f"offline UniProt source metadata is missing: {metadata_path}"
            )
        cached_source = _source_from_dict(_load_json_object(metadata_path))
        cached_source.validate()
        if _sha256_file(raw_path) != cached_source.raw_sha256:
            raise _validation_error("offline UniProt raw cache sha256 mismatch.")
        response = ExternalHttpResponse(
            status_code=200,
            text=raw_path.read_text(encoding="utf-8"),
            url=build_uniprot_url(query),
            headers={},
        )
    else:
        getter = http_get or default_http_get
        response = getter(build_uniprot_url(query), config)
        if not 200 <= response.status_code < 300:
            raise _validation_error(
                f"UniProt identity fetch failed with HTTP {response.status_code}."
            )
        _atomic_write_text(raw_path, response.text)
    result = fetch_uniprot_reference(
        query,
        config,
        http_get=lambda url, cfg: response,
        sleep=lambda _: None,
    )
    if not result.success or not result.records:
        raise _validation_error("UniProt identity response did not yield a record.")
    record = result.records[0]
    if gene_id not in {record.gene_id, record.locus_tag, *record.aliases}:
        raise _validation_error(
            "UniProt identity response does not confirm the requested gene."
        )
    source = ExternalCapacitySource(
        source_id=f"uniprot:{record.primary_accession}",
        source_type=ExternalCapacitySourceType.IDENTITY_REFERENCE,
        source_version=record.provenance.source_version,
        source_url=result.source_url,
        retrieved_at=result.retrieved_at,
        query=json.dumps(query.to_dict(), ensure_ascii=False, sort_keys=True),
        raw_sha256=_sha256_file(raw_path),
        license_id="UniProt terms of use",
        license_url="https://www.uniprot.org/help/license",
        retrieval_mode=retrieval_mode,
        cache_path=str(raw_path),
        terms_reviewed=True,
        warnings=("identity_only_not_capacity_evidence",),
    )
    if cached_source is not None:
        source = replace(
            cached_source,
            retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
            cache_path=str(raw_path),
            raw_sha256=_sha256_file(raw_path),
        )
    else:
        _atomic_write_json(metadata_path, _json_ready(asdict(source)))
    source.validate()
    return source


def _source_from_dict(item: Mapping[str, Any]) -> ExternalCapacitySource:
    return ExternalCapacitySource(
        source_id=str(item.get("source_id") or ""),
        source_type=ExternalCapacitySourceType(str(item.get("source_type") or "")),
        source_version=str(item.get("source_version") or ""),
        source_url=str(item.get("source_url") or ""),
        retrieved_at=str(item.get("retrieved_at") or ""),
        query=str(item.get("query") or ""),
        raw_sha256=str(item.get("raw_sha256") or ""),
        license_id=str(item.get("license_id") or ""),
        retrieval_mode=RetrievalMode(str(item.get("retrieval_mode") or "")),
        cache_path=str(item.get("cache_path") or ""),
        license_url=str(item.get("license_url") or ""),
        adapter_id=str(item.get("adapter_id") or ""),
        adapter_version=str(item.get("adapter_version") or ""),
        terms_reviewed=bool(item.get("terms_reviewed")),
        warnings=tuple(str(value) for value in item.get("warnings") or ()),
    )


def _load_json_object(path: Path) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise _validation_error(f"failed to load JSON object: {path}") from exc
    if not isinstance(payload, Mapping):
        raise _validation_error(f"{path} must be an object.")
    return payload


def _atomic_write_json(path: Path, payload: object) -> None:
    _atomic_write_text(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_name, path)
    except Exception:
        try:
            os.unlink(temp_name)
        except OSError:
            pass
        raise


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_ready(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


def _safe_name(value: str) -> str:
    safe = "".join(
        character
        if character.isalnum() or character in {"-", "_", "."}
        else "_"
        for character in value
    )
    return safe.strip("._") or "artifact"


def _validation_error(message: str) -> Exception:
    return ExternalCapacitySourceValidationError(message)


def _require_text(value: object, field_name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise _validation_error(f"{field_name} must be non-empty text.")


def _require_sha256(value: str, field_name: str) -> None:
    if len(value) != 64 or any(
        character not in "0123456789abcdef" for character in value.lower()
    ):
        raise _validation_error(
            f"{field_name} must be a 64-character hex digest."
        )


__all__ = [
    "ExternalCapacitySource",
    "ExternalCapacitySourceValidationError",
    "ExternalCapacitySourceType",
    "RetrievalMode",
    "cache_uniprot_identity_source",
]
