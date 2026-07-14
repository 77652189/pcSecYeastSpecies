from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Mapping

from pcsec_pichia.external_refs.schema import utc_now_iso
from pcsec_pichia.oe_capacity.external_candidate_schema import (
    CANDIDATE_MANIFEST_FILENAME,
    CANDIDATE_RECORDS_FILENAME,
    RAW_MEASUREMENTS_FILENAME,
    CapacityApplicabilityScope,
    CapacityCandidateStatus,
    CapacityConfidence,
    CapacityConversionStep,
    CapacityModelBinding,
    CapacityParameterKind,
    ExternalCapacityCandidate,
    ExternalCapacityCandidateBundle,
    ExternalCapacitySource,
    ExternalCapacitySourceType,
    HostCondition,
    RawCapacityMeasurement,
    RetrievalMode,
)
from pcsec_pichia.oe_capacity.schema import OECapacityValidationError


@dataclass(frozen=True)
class ExternalCapacityCandidateOutputs:
    records_path: Path
    measurements_path: Path
    manifest_path: Path
    bundle_sha256: str


@dataclass(frozen=True)
class ExternalCapacityCandidateSnapshot:
    bundle: ExternalCapacityCandidateBundle
    manifest_sha256: str


def write_external_capacity_candidate_cache(
    bundle: ExternalCapacityCandidateBundle,
    output_dir: str | Path,
) -> ExternalCapacityCandidateOutputs:
    bundle.validate()
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    records_path = root / CANDIDATE_RECORDS_FILENAME
    measurements_path = root / RAW_MEASUREMENTS_FILENAME
    _atomic_write_jsonl(records_path, (_json_ready(asdict(item)) for item in bundle.candidates))
    _atomic_write_jsonl(measurements_path, (_json_ready(asdict(item)) for item in bundle.measurements))
    records_sha = _sha256_file(records_path)
    measurements_sha = _sha256_file(measurements_path)
    manifest = {
        "schema_version": bundle.schema_version,
        "generated_at": bundle.generated_at,
        "model_fingerprints": list(bundle.model_fingerprints),
        "sources": [_json_ready(asdict(item)) for item in bundle.sources],
        "records_file": records_path.name,
        "records_sha256": records_sha,
        "records_count": len(bundle.candidates),
        "measurements_file": measurements_path.name,
        "measurements_sha256": measurements_sha,
        "measurements_count": len(bundle.measurements),
    }
    manifest_path = root / CANDIDATE_MANIFEST_FILENAME
    _atomic_write_json(manifest_path, manifest)
    return ExternalCapacityCandidateOutputs(
        records_path=records_path,
        measurements_path=measurements_path,
        manifest_path=manifest_path,
        bundle_sha256=_sha256_file(manifest_path),
    )


def load_external_capacity_candidate_bundle(
    source: str | Path,
    *,
    expected_manifest_sha256: str | None = None,
) -> ExternalCapacityCandidateBundle:
    return load_external_capacity_candidate_snapshot(
        source,
        expected_manifest_sha256=expected_manifest_sha256,
    ).bundle


def load_external_capacity_candidate_snapshot(
    source: str | Path,
    *,
    expected_manifest_sha256: str | None = None,
) -> ExternalCapacityCandidateSnapshot:
    manifest_path = Path(source)
    if manifest_path.is_dir():
        manifest_path = manifest_path / CANDIDATE_MANIFEST_FILENAME
    payload, manifest_sha256 = _load_json_object_snapshot_with_sha256(
        manifest_path,
        expected_sha256=expected_manifest_sha256,
        label="manifest",
    )
    records_path = _resolve_manifest_artifact(
        manifest_path,
        payload.get("records_file"),
        "candidate records",
    )
    measurements_path = _resolve_manifest_artifact(
        manifest_path,
        payload.get("measurements_file"),
        "capacity measurements",
    )
    sources = tuple(_source_from_dict(item) for item in _as_object_list(payload.get("sources"), "sources"))
    measurements = tuple(
        _measurement_from_dict(item)
        for item in _read_verified_jsonl_snapshot(
            measurements_path,
            str(payload.get("measurements_sha256") or ""),
            "capacity measurements",
        )
    )
    candidates = tuple(
        _candidate_from_dict(item)
        for item in _read_verified_jsonl_snapshot(
            records_path,
            str(payload.get("records_sha256") or ""),
            "candidate records",
        )
    )
    if len(measurements) != int(payload.get("measurements_count", -1)):
        raise OECapacityValidationError("capacity measurement count mismatch.")
    if len(candidates) != int(payload.get("records_count", -1)):
        raise OECapacityValidationError("capacity candidate count mismatch.")
    bundle = ExternalCapacityCandidateBundle(
        schema_version=int(payload.get("schema_version", 0)),
        generated_at=str(payload.get("generated_at") or ""),
        model_fingerprints=tuple(
            str(value) for value in payload.get("model_fingerprints") or ()
        ),
        sources=sources,
        measurements=measurements,
        candidates=candidates,
    )
    bundle.validate()
    return ExternalCapacityCandidateSnapshot(
        bundle=bundle,
        manifest_sha256=manifest_sha256,
    )


def import_capacity_measurements(
    source_path: str | Path,
    *,
    source_id: str,
    source_type: ExternalCapacitySourceType,
    source_version: str,
    source_url: str,
    license_id: str,
    query: str,
    output_dir: str | Path,
    expected_sha256: str = "",
    license_url: str = "",
    terms_reviewed: bool = False,
) -> tuple[ExternalCapacitySource, tuple[RawCapacityMeasurement, ...]]:
    path = Path(source_path)
    if not path.is_file():
        raise OECapacityValidationError(f"capacity import file does not exist: {path}")
    checksum = _sha256_file(path)
    if expected_sha256:
        _require_sha256(expected_sha256, "expected_sha256")
        if checksum.lower() != expected_sha256.lower():
            raise OECapacityValidationError("capacity import checksum mismatch.")
    root = Path(output_dir)
    raw_dir = root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    copied = raw_dir / _safe_name(path.name)
    if copied.exists():
        raise FileExistsError(f"raw capacity import already exists: {copied}")
    shutil.copyfile(path, copied)
    rows = _load_measurement_rows(copied)
    measurements = tuple(_measurement_from_import_row(row, source_id=source_id) for row in rows)
    for measurement in measurements:
        measurement.validate()
    warnings = tuple(
        item
        for item, missing in (
            ("source_version_review_required", not source_version.strip()),
            ("license_review_required", not license_id.strip()),
            ("expected_sha256_not_predeclared", not expected_sha256),
        )
        if missing
    )
    source = ExternalCapacitySource(
        source_id=source_id,
        source_type=source_type,
        source_version=source_version or "unversioned",
        source_url=source_url,
        retrieved_at=utc_now_iso(),
        query=query,
        raw_sha256=checksum,
        license_id=license_id or "unreviewed",
        license_url=license_url,
        retrieval_mode=RetrievalMode.MANUAL_IMPORT,
        cache_path=str(copied),
        terms_reviewed=terms_reviewed,
        warnings=warnings,
    )
    source.validate()
    return source, measurements

def _source_artifact_matches(source: ExternalCapacitySource) -> bool:
    path = Path(source.cache_path)
    return path.is_file() and _sha256_file(path) == source.raw_sha256

def _measurement_from_import_row(row: Mapping[str, Any], *, source_id: str) -> RawCapacityMeasurement:
    condition_payload = row.get("condition") if isinstance(row.get("condition"), Mapping) else row
    return RawCapacityMeasurement(
        measurement_id=str(row.get("measurement_id") or ""),
        source_id=source_id,
        parameter_kind=CapacityParameterKind(str(row.get("parameter_kind") or "")),
        nominal_value=_float(row.get("nominal_value"), "nominal_value"),
        lower_bound=_float(row.get("lower_bound", row.get("nominal_value")), "lower_bound"),
        upper_bound=_float(row.get("upper_bound", row.get("nominal_value")), "upper_bound"),
        unit=str(row.get("unit") or ""),
        condition=_condition_from_dict(condition_payload),
        external_gene_id=str(row.get("external_gene_id") or ""),
        external_protein_id=str(row.get("external_protein_id") or ""),
        external_enzyme_id=str(row.get("external_enzyme_id") or ""),
        biomass_basis=str(row.get("biomass_basis") or ""),
        notes=str(row.get("notes") or ""),
    )


def _load_measurement_rows(path: Path) -> list[Mapping[str, Any]]:
    suffix = path.suffix.lower()
    if suffix == ".jsonl":
        return _read_jsonl(path)
    if suffix == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, Mapping):
            payload = payload.get("measurements")
        return _as_object_list(payload, "measurements")
    if suffix in {".csv", ".tsv"}:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t" if suffix == ".tsv" else ","))
    raise OECapacityValidationError("capacity import supports JSON, JSONL, CSV, or TSV.")


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


def _condition_from_dict(item: Mapping[str, Any]) -> HostCondition:
    return HostCondition(
        species=str(item.get("species") or ""),
        strain=str(item.get("strain") or ""),
        medium=str(item.get("medium") or ""),
        carbon_source=str(item.get("carbon_source") or ""),
        culture_mode=str(item.get("culture_mode") or ""),
        growth_rate_per_h=_float(item.get("growth_rate_per_h"), "growth_rate_per_h"),
        temperature_c=_optional_float(item.get("temperature_c"), "temperature_c"),
        ph=_optional_float(item.get("ph"), "ph"),
        oxygen_condition=str(item.get("oxygen_condition") or ""),
        biomass_basis=str(item.get("biomass_basis") or "gDW"),
    )


def _measurement_from_dict(item: Mapping[str, Any]) -> RawCapacityMeasurement:
    return RawCapacityMeasurement(
        measurement_id=str(item.get("measurement_id") or ""),
        source_id=str(item.get("source_id") or ""),
        parameter_kind=CapacityParameterKind(str(item.get("parameter_kind") or "")),
        nominal_value=_float(item.get("nominal_value"), "nominal_value"),
        lower_bound=_float(item.get("lower_bound"), "lower_bound"),
        upper_bound=_float(item.get("upper_bound"), "upper_bound"),
        unit=str(item.get("unit") or ""),
        condition=_condition_from_dict(_as_object(item.get("condition"), "condition")),
        external_gene_id=str(item.get("external_gene_id") or ""),
        external_protein_id=str(item.get("external_protein_id") or ""),
        external_enzyme_id=str(item.get("external_enzyme_id") or ""),
        biomass_basis=str(item.get("biomass_basis") or ""),
        notes=str(item.get("notes") or ""),
    )


def _binding_from_dict(item: Mapping[str, Any]) -> CapacityModelBinding:
    return CapacityModelBinding(
        target_id=str(item.get("target_id") or ""),
        context_id=str(item.get("context_id") or ""),
        mapping_id=str(item.get("mapping_id") or ""),
        model_fingerprint=str(item.get("model_fingerprint") or ""),
        gene_id=str(item.get("gene_id") or ""),
        enzyme_id=str(item.get("enzyme_id") or ""),
        reaction_id=str(item.get("reaction_id") or ""),
        formation_or_dilution_reaction_id=str(item.get("formation_or_dilution_reaction_id") or ""),
        mapping_evidence=tuple(str(value) for value in item.get("mapping_evidence") or ()),
        external_gene_id=str(item.get("external_gene_id") or ""),
        external_protein_id=str(item.get("external_protein_id") or ""),
        external_enzyme_id=str(item.get("external_enzyme_id") or ""),
    )


def _step_from_dict(item: Mapping[str, Any]) -> CapacityConversionStep:
    return CapacityConversionStep(
        step_id=str(item.get("step_id") or ""),
        input_value=_float(item.get("input_value"), "input_value"),
        input_unit=str(item.get("input_unit") or ""),
        output_value=_float(item.get("output_value"), "output_value"),
        output_unit=str(item.get("output_unit") or ""),
        formula=str(item.get("formula") or ""),
        factor=_float(item.get("factor"), "factor"),
        source_ref=str(item.get("source_ref") or ""),
        missing_metadata=tuple(str(value) for value in item.get("missing_metadata") or ()),
    )


def _candidate_from_dict(item: Mapping[str, Any]) -> ExternalCapacityCandidate:
    return ExternalCapacityCandidate(
        candidate_id=str(item.get("candidate_id") or ""),
        applicability_scope=CapacityApplicabilityScope(str(item.get("applicability_scope") or "")),
        source_ids=tuple(str(value) for value in item.get("source_ids") or ()),
        measurement_ids=tuple(str(value) for value in item.get("measurement_ids") or ()),
        model_bindings=tuple(
            _binding_from_dict(value)
            for value in _as_object_list(item.get("model_bindings"), "model_bindings")
        ),
        condition=_condition_from_dict(_as_object(item.get("condition"), "condition")),
        nominal_capacity=_optional_float(item.get("nominal_capacity"), "nominal_capacity"),
        lower_capacity=_optional_float(item.get("lower_capacity"), "lower_capacity"),
        upper_capacity=_optional_float(item.get("upper_capacity"), "upper_capacity"),
        unit=str(item.get("unit") or ""),
        confidence=CapacityConfidence(str(item.get("confidence") or "")),
        status=CapacityCandidateStatus(str(item.get("status") or "")),
        conversion_steps=tuple(_step_from_dict(value) for value in _as_object_list(item.get("conversion_steps"), "conversion_steps")),
        target_id=str(item.get("target_id") or ""),
        conflicts=tuple(str(value) for value in item.get("conflicts") or ()),
        missing_information=tuple(str(value) for value in item.get("missing_information") or ()),
        rejection_reasons=tuple(str(value) for value in item.get("rejection_reasons") or ()),
        warnings=tuple(str(value) for value in item.get("warnings") or ()),
    )


def _read_jsonl(path: Path) -> list[Mapping[str, Any]]:
    try:
        snapshot = path.read_bytes()
    except OSError as exc:
        raise OECapacityValidationError(f"failed to read JSONL: {path}") from exc
    return _read_jsonl_snapshot(snapshot, path)


def _read_jsonl_snapshot(snapshot: bytes, path: Path) -> list[Mapping[str, Any]]:
    rows: list[Mapping[str, Any]] = []
    try:
        text = snapshot.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise OECapacityValidationError(f"invalid UTF-8 JSONL: {path}") from exc
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            raise OECapacityValidationError(f"invalid JSONL at {path}:{line_number}.") from exc
        rows.append(_as_object(payload, f"JSONL row {line_number}"))
    return rows


def _load_json_object(path: Path) -> Mapping[str, Any]:
    return _load_json_object_snapshot(path)


def _load_json_object_snapshot(
    path: Path,
    *,
    expected_sha256: str | None = None,
    label: str = "JSON object",
) -> Mapping[str, Any]:
    payload, _ = _load_json_object_snapshot_with_sha256(
        path,
        expected_sha256=expected_sha256,
        label=label,
    )
    return payload


def _load_json_object_snapshot_with_sha256(
    path: Path,
    *,
    expected_sha256: str | None = None,
    label: str = "JSON object",
) -> tuple[Mapping[str, Any], str]:
    try:
        snapshot = path.read_bytes()
    except OSError as exc:
        raise OECapacityValidationError(f"failed to load JSON object: {path}") from exc
    snapshot_sha256 = _sha256_bytes(snapshot)
    if expected_sha256 is not None:
        _require_sha256(expected_sha256, f"{label} sha256")
        if snapshot_sha256 != expected_sha256.lower():
            raise OECapacityValidationError(f"{label} sha256 mismatch.")
    try:
        payload = _as_object(json.loads(snapshot.decode("utf-8")), str(path))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise OECapacityValidationError(f"failed to load JSON object: {path}") from exc
    return payload, snapshot_sha256


def _resolve_manifest_artifact(manifest_path: Path, value: object, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise OECapacityValidationError(f"{label} file must be a plain basename.")
    relative_path = Path(value)
    if relative_path.is_absolute() or relative_path.name != value or len(relative_path.parts) != 1:
        raise OECapacityValidationError(f"{label} file must be a plain basename.")
    manifest_dir = manifest_path.parent.resolve()
    resolved_path = (manifest_dir / relative_path).resolve()
    if resolved_path.parent != manifest_dir:
        raise OECapacityValidationError(f"{label} file must be a plain basename.")
    return resolved_path


def _read_verified_jsonl_snapshot(
    path: Path,
    expected_sha256: str,
    label: str,
) -> list[Mapping[str, Any]]:
    _require_sha256(expected_sha256, f"{label} sha256")
    try:
        snapshot = path.read_bytes()
    except OSError as exc:
        raise OECapacityValidationError(f"{label} file is missing: {path}") from exc
    if _sha256_bytes(snapshot) != expected_sha256.lower():
        raise OECapacityValidationError(f"{label} sha256 mismatch.")
    return _read_jsonl_snapshot(snapshot, path)


def _atomic_write_json(path: Path, payload: object) -> None:
    _atomic_write_text(path, json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n")


def _atomic_write_jsonl(path: Path, rows: Iterable[object]) -> None:
    text = "".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in rows)
    _atomic_write_text(path, text)


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
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


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_json(payload: object) -> str:
    return hashlib.sha256(json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


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
    safe = "".join(character if character.isalnum() or character in {"-", "_", "."} else "_" for character in value)
    return safe.strip("._") or "artifact"


def _as_object(value: object, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise OECapacityValidationError(f"{label} must be an object.")
    return value


def _as_object_list(value: object, label: str) -> list[Mapping[str, Any]]:
    if not isinstance(value, list):
        raise OECapacityValidationError(f"{label} must be an array.")
    return [_as_object(item, label) for item in value]


def _require_sha256(value: str, field_name: str) -> None:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value.lower()):
        raise OECapacityValidationError(f"{field_name} must be a 64-character hex digest.")


def _float(value: object, field_name: str) -> float:
    if isinstance(value, bool):
        raise OECapacityValidationError(f"{field_name} must be numeric.")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise OECapacityValidationError(f"{field_name} must be numeric.") from exc
    if not math.isfinite(parsed):
        raise OECapacityValidationError(f"{field_name} must be finite.")
    return parsed


def _optional_float(value: object, field_name: str) -> float | None:
    if value in (None, ""):
        return None
    return _float(value, field_name)
