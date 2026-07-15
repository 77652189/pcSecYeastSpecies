from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import tempfile
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

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


@dataclass(frozen=True)
class PrideMaxQuantEvidence:
    source: ExternalCapacitySource
    measurement: RawCapacityMeasurement
    raw_values: tuple[float, ...]
    sample_ids: tuple[str, ...]
    mapping_evidence: tuple[str, ...]
    source_bundle_sha256: str
    artifact_sha256s: Mapping[str, str]
    quantitative_boundary: Mapping[str, Any]


@dataclass(frozen=True)
class EcPichiaG6PDH2Evidence:
    source: ExternalCapacitySource
    gene_id: str
    enzyme_id: str
    reaction_id: str
    ec_number: str
    molecular_weight_g_per_mol: float
    kcat_per_s: float
    kcat_source_label: str
    reaction_protein_coefficient: float
    reported_concentration: float
    usage_lower_bound: float
    protein_pool_lower_bound: float
    reported_concentration_unit: str
    gecko_expected_concentration_unit: str
    mapping_evidence: tuple[str, ...]
    conflicts: tuple[str, ...]
    missing_information: tuple[str, ...]


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


def parse_pride_maxquant_g6pdh2_evidence(
    source: ExternalCapacitySource,
) -> PrideMaxQuantEvidence:
    """Parse a cached PRIDE MaxQuant bundle without treating iBAQ as absolute."""

    source.validate()
    if source.source_type is not ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS:
        raise OECapacityValidationError(
            "PRIDE MaxQuant parsing requires a quantitative_proteomics source."
        )
    if source.adapter_id != "pcsec_pichia.pride.maxquant":
        raise OECapacityValidationError("unsupported PRIDE quantitative adapter_id.")
    bundle_path = Path(source.cache_path)
    payload, bundle_sha256 = _load_json_object_snapshot_with_sha256(
        bundle_path,
        expected_sha256=source.raw_sha256,
        label="PRIDE source bundle",
    )
    target = _as_object(payload.get("target"), "PRIDE target")
    condition_payload = _as_object(payload.get("condition"), "PRIDE condition")
    artifacts = _as_object(payload.get("artifacts"), "PRIDE artifacts")
    protein_groups_path = _resolve_source_bundle_artifact(
        bundle_path,
        _as_object(artifacts.get("protein_groups"), "PRIDE protein_groups"),
    )
    sample_map_path = _resolve_source_bundle_artifact(
        bundle_path,
        _as_object(artifacts.get("sample_map"), "PRIDE sample_map"),
    )
    database_fasta_path = _resolve_source_bundle_artifact(
        bundle_path,
        _as_object(artifacts.get("database_fasta"), "PRIDE database_fasta"),
    )
    query_fasta_path = _resolve_source_bundle_artifact(
        bundle_path,
        _as_object(artifacts.get("query_fasta"), "PRIDE query_fasta"),
    )
    sample_ids = tuple(str(value) for value in target.get("sample_ids") or ())
    if not sample_ids:
        raise OECapacityValidationError("PRIDE target requires sample_ids.")
    metric_name = str(target.get("metric_name") or "")
    external_protein_id = str(target.get("external_protein_id") or "")
    query_protein_id = str(target.get("query_protein_id") or "")
    external_gene_id = str(target.get("external_gene_id") or "")
    for value, name in (
        (metric_name, "metric_name"),
        (external_protein_id, "external_protein_id"),
        (query_protein_id, "query_protein_id"),
    ):
        if not value:
            raise OECapacityValidationError(f"PRIDE target requires {name}.")
    _validate_pride_sample_map(sample_map_path, sample_ids)
    row = _find_maxquant_protein_row(protein_groups_path, external_protein_id)
    raw_values = _maxquant_metric_values(row, metric_name, sample_ids)
    query_header, query_sequence = _fasta_record_by_accession(
        query_fasta_path,
        query_protein_id,
    )
    external_header, external_sequence = _fasta_record_by_accession(
        database_fasta_path,
        external_protein_id,
    )
    if query_sequence != external_sequence:
        raise OECapacityValidationError(
            "PRIDE external protein is not an exact sequence match for the query protein."
        )
    condition = _condition_from_dict(condition_payload)
    nominal = float(statistics.median(raw_values))
    measurement = RawCapacityMeasurement(
        measurement_id=(
            f"{source.source_id}:{external_protein_id}:"
            f"{condition.oxygen_condition}:{metric_name}"
        ),
        source_id=source.source_id,
        parameter_kind=CapacityParameterKind.ABUNDANCE,
        nominal_value=nominal,
        lower_bound=min(raw_values),
        upper_bound=max(raw_values),
        unit=f"{metric_name.strip().lower()}_intensity",
        condition=condition,
        external_gene_id=external_gene_id,
        external_protein_id=external_protein_id,
        biomass_basis=condition.biomass_basis,
        notes=json.dumps(
            {
                "aggregation": "median_positive_replicates",
                "raw_values": list(raw_values),
                "sample_ids": list(sample_ids),
                "query_fasta_header": query_header,
                "external_fasta_header": external_header,
                "absolute_abundance": False,
            },
            ensure_ascii=False,
            sort_keys=True,
        ),
    )
    measurement.validate()
    quantitative_boundary = _as_object(
        payload.get("quantitative_boundary"),
        "PRIDE quantitative_boundary",
    )
    artifact_sha256s = {
        role: str(_as_object(record, f"PRIDE {role}").get("sha256") or "")
        for role, record in artifacts.items()
    }
    for role, sha256 in artifact_sha256s.items():
        _require_sha256(sha256, f"PRIDE {role} sha256")
    return PrideMaxQuantEvidence(
        source=source,
        measurement=measurement,
        raw_values=raw_values,
        sample_ids=sample_ids,
        mapping_evidence=(
            f"pride:{payload.get('project_accession')}",
            f"uniprot:{query_protein_id}",
            f"pride_protein:{external_protein_id}",
            f"exact_sequence_identity:{len(query_sequence)}/{len(query_sequence)}",
        ),
        source_bundle_sha256=bundle_sha256,
        artifact_sha256s=artifact_sha256s,
        quantitative_boundary=quantitative_boundary,
    )


def parse_ecpichia_g6pdh2_source_assessment(
    source: ExternalCapacitySource,
) -> EcPichiaG6PDH2Evidence:
    """Parse ecPichia raw YAML values without creating a capacity measurement."""

    source.validate()
    if source.source_type is not ExternalCapacitySourceType.EXTERNAL_ENZYME_MODEL:
        raise OECapacityValidationError(
            "ecPichia parsing requires an external_enzyme_model source."
        )
    if source.adapter_id != "pcsec_pichia.ecpichia.supplement_yaml_assessment":
        raise OECapacityValidationError("unsupported ecPichia supplement adapter_id.")
    path = Path(source.cache_path)
    if not path.is_file() or _sha256_file(path) != source.raw_sha256:
        raise OECapacityValidationError("ecPichia supplement sha256 mismatch.")
    text = path.read_text(encoding="utf-8")
    reaction = _yaml_omap_entry(text, "id", "G6PDH2", indent=4)
    kinetic = _yaml_omap_entry(text, "id", "G6PDH2", indent=2)
    protein = _yaml_omap_entry(text, "genes", "PAS_chr2-1_0308", indent=2)
    usage = _yaml_omap_entry(text, "id", "usage_prot_C4R099", indent=4)
    pool = _yaml_omap_entry(text, "id", "prot_pool_exchange", indent=4)

    gene_id = _yaml_quoted_value(reaction, "gene_reaction_rule")
    ec_number = _yaml_quoted_value(reaction, "eccodes")
    enzyme_id = _yaml_quoted_value(protein, "enzymes")
    protein_coefficient = _yaml_number(reaction, "prot_C4R099")
    molecular_weight = _yaml_number(protein, "mw")
    concentration = _yaml_number(protein, "concs")
    kcat = _yaml_number(kinetic, "kcat")
    kcat_source = _yaml_quoted_value(kinetic, "source")
    usage_lower_bound = _yaml_number(usage, "lower_bound")
    pool_lower_bound = _yaml_number(pool, "lower_bound")
    if gene_id != "PAS_chr2-1_0308" or enzyme_id != "C4R099":
        raise OECapacityValidationError("ecPichia G6PDH2 gene/enzyme binding mismatch.")
    if ec_number != "1.1.1.49":
        raise OECapacityValidationError("ecPichia G6PDH2 EC number mismatch.")
    for value, label in (
        (molecular_weight, "molecular weight"),
        (concentration, "reported concentration"),
        (kcat, "kcat"),
    ):
        if not math.isfinite(value) or value <= 0:
            raise OECapacityValidationError(f"ecPichia {label} must be finite and positive.")
    if protein_coefficient >= 0 or usage_lower_bound >= 0 or pool_lower_bound >= 0:
        raise OECapacityValidationError(
            "ecPichia enzyme coefficients and pool bounds require negative GECKO signs."
        )
    expected_coefficient = -(molecular_weight / (kcat * 3600.0))
    if not math.isclose(protein_coefficient, expected_coefficient, rel_tol=1e-12):
        raise OECapacityValidationError(
            "ecPichia G6PDH2 protein coefficient is inconsistent with MW/kcat."
        )
    if not math.isclose(-usage_lower_bound, concentration, rel_tol=1e-12):
        raise OECapacityValidationError(
            "ecPichia usage bound is inconsistent with reported concentration."
        )
    return EcPichiaG6PDH2Evidence(
        source=source,
        gene_id=gene_id,
        enzyme_id=enzyme_id,
        reaction_id="G6PDH2",
        ec_number=ec_number,
        molecular_weight_g_per_mol=molecular_weight,
        kcat_per_s=kcat,
        kcat_source_label=kcat_source,
        reaction_protein_coefficient=protein_coefficient,
        reported_concentration=concentration,
        usage_lower_bound=usage_lower_bound,
        protein_pool_lower_bound=pool_lower_bound,
        reported_concentration_unit="unit_not_declared_in_yaml",
        gecko_expected_concentration_unit="mg_per_gDCW_by_GECKO_convention",
        mapping_evidence=(
            "gene:PAS_chr2-1_0308",
            "enzyme:C4R099",
            "reaction:G6PDH2",
            "ec:1.1.1.49",
        ),
        conflicts=(
            "source_unit_missing_and_supplement_header_requires_review",
            "supplement_table_yaml_binding_requires_reconciliation",
            "source_condition_not_verified_for_formal_glucose_mu_0.1",
        ),
        missing_information=(
            "supplement_reuse_license_missing",
            "lfq_to_absolute_conversion_missing",
            "biomass_normalization_and_sample_selection_missing",
            "brenda_record_metadata_missing",
            "condition_metadata_missing",
            "formation_flux_conversion_missing",
            "current_model_capacity_handle_binding_review",
        ),
    )


def _yaml_omap_entry(text: str, key: str, value: str, *, indent: int) -> str:
    prefix = " " * indent
    field_indent = " " * (indent + 2)
    pattern = re.compile(
        rf"(?ms)^{re.escape(prefix)}- !!omap\r?\n"
        rf"{re.escape(field_indent)}- {re.escape(key)}: \"{re.escape(value)}\"\r?\n"
        rf"(?P<body>.*?)(?=^{re.escape(prefix)}- !!omap|\Z)"
    )
    matches = tuple(pattern.finditer(text))
    if len(matches) != 1:
        raise OECapacityValidationError(
            f"ecPichia YAML requires exactly one {key}={value} entry at indent {indent}."
        )
    return matches[0].group(0)


def _yaml_quoted_value(block: str, key: str) -> str:
    matches = re.findall(rf"(?m)^\s+- {re.escape(key)}: \"([^\"]+)\"\s*$", block)
    if len(matches) != 1:
        raise OECapacityValidationError(
            f"ecPichia YAML requires exactly one quoted {key} value."
        )
    return matches[0]


def _yaml_number(block: str, key: str) -> float:
    matches = re.findall(
        rf"(?m)^\s+- {re.escape(key)}: (-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$",
        block,
    )
    if len(matches) != 1:
        raise OECapacityValidationError(
            f"ecPichia YAML requires exactly one numeric {key} value."
        )
    return _float(matches[0], f"ecPichia {key}")

def _source_artifact_matches(source: ExternalCapacitySource) -> bool:
    path = Path(source.cache_path)
    return path.is_file() and _sha256_file(path) == source.raw_sha256


def _resolve_source_bundle_artifact(
    bundle_path: Path,
    record: Mapping[str, Any],
) -> Path:
    filename = str(record.get("filename") or "")
    if not filename or Path(filename).name != filename:
        raise OECapacityValidationError(
            "PRIDE source artifact filename must be a basename."
        )
    root = bundle_path.parent.resolve()
    path = (root / filename).resolve()
    if path.parent != root or not path.is_file():
        raise OECapacityValidationError(
            "PRIDE source artifact is missing or outside the source bundle."
        )
    expected_sha256 = str(record.get("sha256") or "")
    _require_sha256(expected_sha256, "PRIDE source artifact sha256")
    if _sha256_file(path) != expected_sha256.lower():
        raise OECapacityValidationError("PRIDE source artifact sha256 mismatch.")
    return path


def _validate_pride_sample_map(path: Path, sample_ids: Sequence[str]) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = tuple(csv.DictReader(handle, delimiter="\t"))
    available = {str(row.get("protein groups.txt") or "") for row in rows}
    missing = tuple(sample_id for sample_id in sample_ids if sample_id not in available)
    if missing:
        raise OECapacityValidationError(
            "PRIDE sample map is missing selected sample IDs: " + ", ".join(missing)
        )


def _find_maxquant_protein_row(
    path: Path,
    external_protein_id: str,
) -> Mapping[str, Any]:
    matches: list[Mapping[str, Any]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            identifiers = {
                identifier.strip()
                for field_name in ("Protein IDs", "Majority protein IDs")
                for identifier in str(row.get(field_name) or "").split(";")
                if identifier.strip()
            }
            if external_protein_id in identifiers:
                matches.append(row)
    if len(matches) != 1:
        raise OECapacityValidationError(
            "PRIDE MaxQuant table requires exactly one target protein row."
        )
    return matches[0]


def _maxquant_metric_values(
    row: Mapping[str, Any],
    metric_name: str,
    sample_ids: Sequence[str],
) -> tuple[float, ...]:
    values: list[float] = []
    for sample_id in sample_ids:
        column = f"{metric_name} {sample_id}"
        raw = row.get(column)
        if raw in (None, ""):
            raise OECapacityValidationError(
                f"PRIDE MaxQuant table is missing quantitative column {column}."
            )
        value = _float(raw, column)
        if value <= 0:
            raise OECapacityValidationError(
                f"PRIDE quantitative value must be positive for selected sample {sample_id}."
            )
        values.append(value)
    return tuple(values)


def _fasta_record_by_accession(path: Path, accession: str) -> tuple[str, str]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence)))
            header = line[1:].strip()
            sequence = []
        elif header:
            sequence.append(line.strip())
    if header:
        records.append((header, "".join(sequence)))
    matches = tuple(
        (item_header, item_sequence)
        for item_header, item_sequence in records
        if accession in item_header.split()[0].split("|")
    )
    if len(matches) != 1 or not matches[0][1]:
        raise OECapacityValidationError(
            f"FASTA requires exactly one non-empty record for {accession}."
        )
    return matches[0]

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
