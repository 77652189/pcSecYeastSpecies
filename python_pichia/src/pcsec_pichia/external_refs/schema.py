from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import UTC, datetime
from typing import Any, Literal, Mapping, TypeAlias


CACHE_SCHEMA_VERSION = "external_refs.v1"

ExternalRecordType = Literal[
    "external_reference",
    "gene_function",
    "reaction_association",
    "gpr_candidate",
]


class ExternalReferenceSchemaError(ValueError):
    """Raised when an external reference cache record violates the local schema."""


@dataclass(frozen=True)
class ExternalReferenceProvenance:
    source_database: str
    source_version: str
    source_url: str
    source_query: str
    retrieved_at: str
    raw_record_sha256: str
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        missing = [
            field_name
            for field_name in (
                "source_database",
                "source_version",
                "source_url",
                "source_query",
                "retrieved_at",
                "raw_record_sha256",
            )
            if not str(getattr(self, field_name, "")).strip()
        ]
        if missing:
            raise ExternalReferenceSchemaError(f"Missing provenance field(s): {', '.join(missing)}")
        if len(self.raw_record_sha256) != 64 or any(char not in "0123456789abcdef" for char in self.raw_record_sha256.lower()):
            raise ExternalReferenceSchemaError("raw_record_sha256 must be a 64-character hex SHA-256 digest.")
        _parse_timestamp(self.retrieved_at, "retrieved_at")


@dataclass(frozen=True)
class ExternalReferenceRecord:
    provenance: ExternalReferenceProvenance
    taxon_id: str
    organism: str
    primary_accession: str
    gene_id: str | None = None
    gene_name: str | None = None
    locus_tag: str | None = None
    aliases: tuple[str, ...] = ()
    protein_name: str | None = None
    reviewed: bool | None = None
    sequence_accession: str | None = None
    protein_sequence_sha256: str | None = None
    record_type: ExternalRecordType = "external_reference"

    @property
    def cache_key(self) -> str:
        return stable_cache_key(
            self.record_type,
            self.provenance.source_database,
            self.taxon_id,
            self.primary_accession,
        )

    def validate(self) -> None:
        _expect_record_type(self.record_type, "external_reference")
        self.provenance.validate()
        _require_nonempty("taxon_id", self.taxon_id)
        _require_nonempty("organism", self.organism)
        _require_nonempty("primary_accession", self.primary_accession)


@dataclass(frozen=True)
class ExternalGeneFunctionEvidence:
    provenance: ExternalReferenceProvenance
    gene_id: str
    protein_name: str | None = None
    function_description: str | None = None
    ec_numbers: tuple[str, ...] = ()
    go_terms: tuple[str, ...] = ()
    pathways: tuple[str, ...] = ()
    orthology: tuple[str, ...] = ()
    reviewed: bool | None = None
    evidence_scope: str = "annotation_only"
    record_type: ExternalRecordType = "gene_function"

    @property
    def cache_key(self) -> str:
        return stable_cache_key(
            self.record_type,
            self.provenance.source_database,
            self.gene_id,
            self.provenance.source_query,
        )

    def validate(self) -> None:
        _expect_record_type(self.record_type, "gene_function")
        self.provenance.validate()
        _require_nonempty("gene_id", self.gene_id)
        if not any((self.protein_name, self.function_description, self.ec_numbers, self.go_terms, self.pathways, self.orthology)):
            raise ExternalReferenceSchemaError("ExternalGeneFunctionEvidence needs at least one annotation field.")


@dataclass(frozen=True)
class ExternalReactionAssociation:
    provenance: ExternalReferenceProvenance
    external_model_id: str
    external_reaction_id: str
    external_reaction_name: str | None = None
    external_gene_ids: tuple[str, ...] = ()
    gene_rule: str | None = None
    ec_numbers: tuple[str, ...] = ()
    mapped_pichia_reaction_id: str | None = None
    mapped_pichia_gene_ids: tuple[str, ...] = ()
    association_status: str = "external_gpr_candidate"
    record_type: ExternalRecordType = "reaction_association"

    @property
    def cache_key(self) -> str:
        return stable_cache_key(
            self.record_type,
            self.provenance.source_database,
            self.external_model_id,
            self.external_reaction_id,
            self.gene_rule or "",
        )

    def validate(self) -> None:
        _expect_record_type(self.record_type, "reaction_association")
        self.provenance.validate()
        _require_nonempty("external_model_id", self.external_model_id)
        _require_nonempty("external_reaction_id", self.external_reaction_id)
        if not self.external_gene_ids and not self.gene_rule and not self.ec_numbers:
            raise ExternalReferenceSchemaError("ExternalReactionAssociation needs gene IDs, a gene rule, or EC numbers.")


@dataclass(frozen=True)
class ExternalGprCandidateEvidence:
    provenance: ExternalReferenceProvenance
    external_model_id: str
    external_reaction_id: str
    external_gene_rule: str
    candidate_status: str = "external_gpr_candidate"
    mapped_pichia_reaction_id: str | None = None
    mapped_pichia_gene_ids: tuple[str, ...] = ()
    mapping_warnings: tuple[str, ...] = ()
    record_type: ExternalRecordType = "gpr_candidate"

    @property
    def cache_key(self) -> str:
        return stable_cache_key(
            self.record_type,
            self.provenance.source_database,
            self.external_model_id,
            self.external_reaction_id,
            self.external_gene_rule,
        )

    def validate(self) -> None:
        _expect_record_type(self.record_type, "gpr_candidate")
        self.provenance.validate()
        _require_nonempty("external_model_id", self.external_model_id)
        _require_nonempty("external_reaction_id", self.external_reaction_id)
        _require_nonempty("external_gene_rule", self.external_gene_rule)
        if self.candidate_status == "model_gpr_executable" and (
            not self.mapped_pichia_reaction_id or not self.mapped_pichia_gene_ids
        ):
            raise ExternalReferenceSchemaError(
                "model_gpr_executable requires mapped Pichia reaction and gene IDs."
            )


ExternalCacheRecord: TypeAlias = (
    ExternalReferenceRecord
    | ExternalGeneFunctionEvidence
    | ExternalReactionAssociation
    | ExternalGprCandidateEvidence
)


@dataclass(frozen=True)
class ExternalReferenceCacheManifest:
    generated_at: str
    cache_schema_version: str
    query_count: int
    record_count: int
    failed_query_count: int
    source_counts: Mapping[str, int]
    record_type_counts: Mapping[str, int] = field(default_factory=dict)
    duplicate_key_count: int = 0
    input_cache_fingerprint: str | None = None
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _parse_timestamp(self.generated_at, "generated_at")
        if self.cache_schema_version != CACHE_SCHEMA_VERSION:
            raise ExternalReferenceSchemaError(
                f"Unsupported cache schema version: {self.cache_schema_version!r}."
            )
        for field_name in ("query_count", "record_count", "failed_query_count", "duplicate_key_count"):
            if int(getattr(self, field_name)) < 0:
                raise ExternalReferenceSchemaError(f"{field_name} must be non-negative.")
        _validate_count_mapping("source_counts", self.source_counts)
        _validate_count_mapping("record_type_counts", self.record_type_counts)
        if sum(int(value) for value in self.source_counts.values()) != int(self.record_count):
            raise ExternalReferenceSchemaError("source_counts total must match record_count.")
        if self.record_type_counts and sum(int(value) for value in self.record_type_counts.values()) != int(self.record_count):
            raise ExternalReferenceSchemaError("record_type_counts total must match record_count.")


RECORD_CLASS_BY_TYPE: Mapping[ExternalRecordType, type[ExternalCacheRecord]] = {
    "external_reference": ExternalReferenceRecord,
    "gene_function": ExternalGeneFunctionEvidence,
    "reaction_association": ExternalReactionAssociation,
    "gpr_candidate": ExternalGprCandidateEvidence,
}


def utc_now_iso() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def stable_cache_key(*parts: object) -> str:
    normalized = "\x1f".join(str(part or "").strip().lower() for part in parts)
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def record_to_dict(record: ExternalCacheRecord) -> dict[str, Any]:
    validate_external_cache_record(record)
    return _json_ready(asdict(record))


def record_from_dict(payload: Mapping[str, Any]) -> ExternalCacheRecord:
    record_type = str(payload.get("record_type", "external_reference"))
    if record_type not in RECORD_CLASS_BY_TYPE:
        raise ExternalReferenceSchemaError(f"Unsupported external record_type: {record_type!r}.")
    normalized = dict(payload)
    provenance_payload = normalized.get("provenance")
    if not isinstance(provenance_payload, Mapping):
        raise ExternalReferenceSchemaError("External cache record requires a provenance mapping.")
    normalized["provenance"] = ExternalReferenceProvenance(**_tuple_normalized_dict(provenance_payload))
    cls = RECORD_CLASS_BY_TYPE[record_type]  # type: ignore[index]
    record = cls(**_tuple_normalized_dict(normalized))
    validate_external_cache_record(record)
    return record


def validate_external_cache_record(record: ExternalCacheRecord) -> None:
    if not is_dataclass(record) or not hasattr(record, "validate"):
        raise ExternalReferenceSchemaError(f"Unsupported external cache record object: {type(record)!r}.")
    record.validate()


def validate_no_duplicate_cache_keys(records: tuple[ExternalCacheRecord, ...]) -> None:
    seen: dict[str, ExternalCacheRecord] = {}
    duplicates: list[str] = []
    for record in records:
        key = record.cache_key
        if key in seen:
            duplicates.append(key)
        seen[key] = record
    if duplicates:
        raise ExternalReferenceSchemaError(f"Duplicate external reference cache key(s): {', '.join(duplicates)}")


def build_external_reference_manifest(
    records: tuple[ExternalCacheRecord, ...],
    *,
    query_count: int | None = None,
    failed_query_count: int = 0,
    input_cache_fingerprint: str | None = None,
    warnings: tuple[str, ...] = (),
) -> ExternalReferenceCacheManifest:
    for record in records:
        validate_external_cache_record(record)
    source_counts: dict[str, int] = {}
    record_type_counts: dict[str, int] = {}
    seen: set[str] = set()
    duplicate_count = 0
    for record in records:
        source = record.provenance.source_database
        source_counts[source] = source_counts.get(source, 0) + 1
        record_type_counts[record.record_type] = record_type_counts.get(record.record_type, 0) + 1
        if record.cache_key in seen:
            duplicate_count += 1
        seen.add(record.cache_key)
    manifest = ExternalReferenceCacheManifest(
        generated_at=utc_now_iso(),
        cache_schema_version=CACHE_SCHEMA_VERSION,
        query_count=len(records) if query_count is None else int(query_count),
        record_count=len(records),
        failed_query_count=int(failed_query_count),
        source_counts=source_counts,
        record_type_counts=record_type_counts,
        duplicate_key_count=duplicate_count,
        input_cache_fingerprint=input_cache_fingerprint,
        warnings=warnings,
    )
    manifest.validate()
    return manifest


def manifest_to_dict(manifest: ExternalReferenceCacheManifest) -> dict[str, Any]:
    manifest.validate()
    return _json_ready(asdict(manifest))


def manifest_from_dict(payload: Mapping[str, Any]) -> ExternalReferenceCacheManifest:
    manifest = ExternalReferenceCacheManifest(**_tuple_normalized_dict(payload))
    manifest.validate()
    return manifest


def _require_nonempty(field_name: str, value: str | None) -> None:
    if not str(value or "").strip():
        raise ExternalReferenceSchemaError(f"{field_name} must be non-empty.")


def _expect_record_type(actual: str, expected: str) -> None:
    if actual != expected:
        raise ExternalReferenceSchemaError(f"record_type must be {expected!r}, got {actual!r}.")


def _parse_timestamp(value: str, field_name: str) -> None:
    try:
        datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError as exc:
        raise ExternalReferenceSchemaError(f"{field_name} must be an ISO-8601 timestamp.") from exc


def _validate_count_mapping(field_name: str, counts: Mapping[str, int]) -> None:
    for key, value in counts.items():
        if not str(key).strip():
            raise ExternalReferenceSchemaError(f"{field_name} keys must be non-empty.")
        if int(value) < 0:
            raise ExternalReferenceSchemaError(f"{field_name} values must be non-negative.")


def _tuple_normalized_dict(payload: Mapping[str, Any]) -> dict[str, Any]:
    normalized: dict[str, Any] = {}
    tuple_fields = {
        "aliases",
        "warnings",
        "ec_numbers",
        "go_terms",
        "pathways",
        "orthology",
        "external_gene_ids",
        "mapped_pichia_gene_ids",
        "mapping_warnings",
    }
    for key, value in payload.items():
        if key in tuple_fields and isinstance(value, list):
            normalized[key] = tuple(value)
        else:
            normalized[key] = value
    return normalized


def _json_ready(value: Any) -> Any:
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


def canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
