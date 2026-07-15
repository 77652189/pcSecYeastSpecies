from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass, replace
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Sequence

from pcsec_pichia.external_refs.clients import (
    ExternalFetchConfig,
    ExternalHttpResponse,
    default_http_get,
)
from pcsec_pichia.external_refs.queries import ExternalReferenceQuery
from pcsec_pichia.external_refs.schema import utc_now_iso
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
class PrideMaxQuantSourceProfile:
    project_accession: str
    protein_groups_filename: str
    sample_map_filename: str
    database_fasta_filename: str
    query_fasta_url: str
    query_protein_id: str
    external_gene_id: str
    external_protein_id: str
    sample_ids: tuple[str, ...]
    metric_name: str
    species: str
    strain: str
    medium: str
    carbon_source: str
    culture_mode: str
    growth_rate_per_h: float
    temperature_c: float
    ph: float
    oxygen_condition: str
    biomass_basis: str
    condition_source: str


@dataclass(frozen=True)
class EcPichiaSupplementSourceProfile:
    source_id: str
    source_url: str
    source_version: str
    doi: str
    artifact_filename: str
    artifact_sha256: str
    upstream_archive_sha256: str
    license_id: str
    license_url: str
    gene_id: str
    enzyme_id: str
    reaction_id: str


PXD055501_G6PDH2_PROFILE = PrideMaxQuantSourceProfile(
    project_accession="PXD055501",
    protein_groups_filename="lfq_proteinGroups_PRIDE.txt",
    sample_map_filename="sample_names.txt",
    database_fasta_filename="database.fasta",
    query_fasta_url="https://rest.uniprot.org/uniprotkb/C4R099.fasta",
    query_protein_id="C4R099",
    external_gene_id="ZWF1",
    external_protein_id="F2QTE5",
    sample_ids=("WT_T1_ft2", "WT_T1_ft3", "WT_T1_ft4"),
    metric_name="iBAQ",
    species="Komagataella phaffii",
    strain="CBS7435 wild-type",
    medium="glucose-limited defined chemostat medium",
    carbon_source="glucose",
    culture_mode="chemostat",
    growth_rate_per_h=0.075,
    temperature_c=25.0,
    ph=5.5,
    oxygen_condition="normoxic T0",
    biomass_basis="relative_iBAQ_intensity_not_biomass_normalized",
    condition_source="PXD055501; doi:10.1111/1751-7915.70106",
)


ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE = EcPichiaSupplementSourceProfile(
    source_id="ecpichia:supplementary-8:g6pdh2",
    source_url=(
        "https://ars.els-cdn.com/content/image/"
        "1-s2.0-S1369703X25003146-mmc1.zip"
    ),
    source_version="doi:10.1016/j.bej.2025.109940;PII:S1369703X25003146",
    doi="10.1016/j.bej.2025.109940",
    artifact_filename="Supplementary 8.yml",
    artifact_sha256="317ab62f77c95feb2758f9ad7ed5efe18ff8430c747fbb880c03bb4d6b943d34",
    upstream_archive_sha256="bea45233dc4feb81295315c4e73ca2ca4c886f648822dda27347be8892a3620c",
    license_id="reuse_terms_not_established",
    license_url="https://www.elsevier.com/tdm/tdmrep-policy.json",
    gene_id="PAS_chr2-1_0308",
    enzyme_id="C4R099",
    reaction_id="G6PDH2",
)


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


def cache_pride_maxquant_source(
    profile: PrideMaxQuantSourceProfile,
    output_dir: str | Path,
    *,
    http_get: Any = None,
    retrieval_mode: RetrievalMode = RetrievalMode.ONLINE,
) -> ExternalCapacitySource:
    """Cache one reviewed PRIDE MaxQuant source without interpreting its values."""

    _validate_pride_profile(profile)
    root = Path(output_dir)
    raw_dir = root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    safe_accession = _safe_name(profile.project_accession)
    bundle_path = raw_dir / f"pride-{safe_accession}.bundle.json"
    metadata_path = root / f"pride-{safe_accession}.source.json"
    if retrieval_mode is RetrievalMode.OFFLINE_REPLAY:
        if not bundle_path.is_file() or not metadata_path.is_file():
            raise _validation_error(
                f"offline PRIDE cache is incomplete for {profile.project_accession}."
            )
        cached = _source_from_dict(_load_json_object(metadata_path))
        cached.validate()
        if _sha256_file(bundle_path) != cached.raw_sha256:
            raise _validation_error("offline PRIDE source bundle sha256 mismatch.")
        _validate_pride_bundle(bundle_path, profile)
        replayed = replace(
            cached,
            retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
            cache_path=str(bundle_path),
            raw_sha256=_sha256_file(bundle_path),
        )
        _validate_pride_source_coherence(replayed, bundle_path, profile)
        return replayed

    config = ExternalFetchConfig(sources=("pride",), timeout_seconds=60.0)
    getter = http_get or default_http_get
    project_url = _pride_project_url(profile.project_accession)
    files_url = _pride_files_url(profile.project_accession)
    project_response = _fetch_text(project_url, getter, config, "PRIDE project")
    files_response = _fetch_text(files_url, getter, config, "PRIDE files")
    project_payload = _as_mapping_json(project_response.text, "PRIDE project")
    files_payload = _as_sequence_json(files_response.text, "PRIDE files")
    if str(project_payload.get("accession") or "") != profile.project_accession:
        raise _validation_error("PRIDE project accession mismatch.")

    project_path = raw_dir / f"{safe_accession}.project.json"
    files_path = raw_dir / f"{safe_accession}.files.json"
    _atomic_write_text(project_path, project_response.text)
    _atomic_write_text(files_path, files_response.text)
    selected: dict[str, dict[str, object]] = {}
    file_roles = {
        "protein_groups": profile.protein_groups_filename,
        "sample_map": profile.sample_map_filename,
        "database_fasta": profile.database_fasta_filename,
    }
    for role, filename in file_roles.items():
        entry = _pride_file_entry(files_payload, filename)
        download_url = _pride_https_download_url(entry)
        response = _fetch_text(download_url, getter, config, f"PRIDE {role}")
        destination = raw_dir / _safe_name(filename)
        _atomic_write_text(destination, response.text)
        expected_size = int(entry.get("fileSizeBytes") or 0)
        if expected_size and destination.stat().st_size != expected_size:
            raise _validation_error(f"PRIDE {role} file size mismatch.")
        selected[role] = {
            "filename": destination.name,
            "download_url": download_url,
            "sha256": _sha256_file(destination),
            "file_size_bytes": destination.stat().st_size,
            "pride_file_accession": str(entry.get("accession") or ""),
        }
    query_response = _fetch_text(
        profile.query_fasta_url,
        getter,
        config,
        "UniProt query FASTA",
    )
    query_path = raw_dir / f"{_safe_name(profile.query_protein_id)}.fasta"
    _atomic_write_text(query_path, query_response.text)
    selected["query_fasta"] = {
        "filename": query_path.name,
        "download_url": profile.query_fasta_url,
        "sha256": _sha256_file(query_path),
        "file_size_bytes": query_path.stat().st_size,
    }

    raw_license = str(project_payload.get("license") or "")
    license_id, license_url, terms_reviewed = _normalize_pride_license(raw_license)
    publication_date = str(project_payload.get("publicationDate") or "").split("T", 1)[0]
    warnings = [
        "relative_ibaq_not_absolute_abundance",
        "source_growth_rate_0.075_not_formal_growth_rate_0.1",
    ]
    if not publication_date:
        warnings.append("source_version_review_required")
    if not terms_reviewed:
        warnings.append("license_review_required")
    condition = {
        "species": profile.species,
        "strain": profile.strain,
        "medium": profile.medium,
        "carbon_source": profile.carbon_source,
        "culture_mode": profile.culture_mode,
        "growth_rate_per_h": profile.growth_rate_per_h,
        "temperature_c": profile.temperature_c,
        "ph": profile.ph,
        "oxygen_condition": profile.oxygen_condition,
        "biomass_basis": profile.biomass_basis,
        "source_ref": profile.condition_source,
    }
    bundle = {
        "schema_version": 1,
        "adapter_id": "pcsec_pichia.pride.maxquant",
        "adapter_version": "1",
        "project_accession": profile.project_accession,
        "project_title": str(project_payload.get("title") or ""),
        "project_publication_date": publication_date,
        "project_license_raw": raw_license,
        "project_metadata": {
            "filename": project_path.name,
            "sha256": _sha256_file(project_path),
            "source_url": project_url,
        },
        "files_metadata": {
            "filename": files_path.name,
            "sha256": _sha256_file(files_path),
            "source_url": files_url,
        },
        "artifacts": selected,
        "target": {
            "query_protein_id": profile.query_protein_id,
            "external_gene_id": profile.external_gene_id,
            "external_protein_id": profile.external_protein_id,
            "metric_name": profile.metric_name,
            "sample_ids": list(profile.sample_ids),
        },
        "condition": condition,
        "quantitative_boundary": {
            "raw_value_available": True,
            "absolute_abundance_available": False,
            "model_flux_conversion_available": False,
            "missing": [
                "absolute abundance calibration",
                "biomass-normalized enzyme amount",
                "paired condition-matched kcat",
                "formal glucose_mu_0.1 condition match",
            ],
        },
    }
    _atomic_write_json(bundle_path, bundle)
    source = ExternalCapacitySource(
        source_id=f"pride:{profile.project_accession}",
        source_type=ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS,
        source_version=(
            f"{profile.project_accession}:{publication_date}"
            if publication_date
            else "unversioned"
        ),
        source_url=project_url,
        retrieved_at=utc_now_iso(),
        query=json.dumps(
            {
                "project_accession": profile.project_accession,
                "query_protein_id": profile.query_protein_id,
                "external_gene_id": profile.external_gene_id,
                "external_protein_id": profile.external_protein_id,
                "metric_name": profile.metric_name,
                "sample_ids": list(profile.sample_ids),
            },
            ensure_ascii=False,
            sort_keys=True,
        ),
        raw_sha256=_sha256_file(bundle_path),
        license_id=license_id,
        license_url=license_url,
        retrieval_mode=RetrievalMode.ONLINE,
        cache_path=str(bundle_path),
        adapter_id="pcsec_pichia.pride.maxquant",
        adapter_version="1",
        terms_reviewed=terms_reviewed,
        warnings=tuple(warnings),
    )
    source.validate()
    _atomic_write_json(metadata_path, _json_ready(asdict(source)))
    _validate_pride_bundle(bundle_path, profile)
    _validate_pride_source_coherence(source, bundle_path, profile)
    return source


def cache_ecpichia_supplement_source(
    profile: EcPichiaSupplementSourceProfile,
    output_dir: str | Path,
    *,
    source_file: str | Path | None = None,
    retrieval_mode: RetrievalMode = RetrievalMode.MANUAL_IMPORT,
) -> ExternalCapacitySource:
    """Cache an ecPichia YAML artifact for assessment-only offline replay."""

    _validate_ecpichia_profile(profile)
    root = Path(output_dir)
    raw_dir = root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    cached_path = raw_dir / "ecpichia-supplementary-8.yml"
    metadata_path = root / "ecpichia-supplementary-8.source.json"
    if retrieval_mode is RetrievalMode.OFFLINE_REPLAY:
        if not cached_path.is_file() or not metadata_path.is_file():
            raise _validation_error("offline ecPichia supplement cache is incomplete.")
        cached = _source_from_dict(_load_json_object(metadata_path))
        cached.validate()
        if _sha256_file(cached_path) != profile.artifact_sha256:
            raise _validation_error("offline ecPichia supplement sha256 mismatch.")
        source = replace(
            cached,
            retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
            cache_path=str(cached_path),
            raw_sha256=_sha256_file(cached_path),
        )
        _validate_ecpichia_source_coherence(source, profile)
        return source
    if retrieval_mode is not RetrievalMode.MANUAL_IMPORT:
        raise _validation_error(
            "ecPichia supplement requires manual_import or offline_replay."
        )
    if source_file is None:
        raise _validation_error("ecPichia manual import requires source_file.")
    source_path = Path(source_file)
    if not source_path.is_file():
        raise _validation_error(f"ecPichia source file is missing: {source_path}")
    if _sha256_file(source_path) != profile.artifact_sha256:
        raise _validation_error("ecPichia supplement sha256 does not match reviewed profile.")
    if source_path.resolve() != cached_path.resolve():
        _atomic_copy_file(source_path, cached_path)
    source = ExternalCapacitySource(
        source_id=profile.source_id,
        source_type=ExternalCapacitySourceType.EXTERNAL_ENZYME_MODEL,
        source_version=profile.source_version,
        source_url=profile.source_url,
        retrieved_at=utc_now_iso(),
        query=json.dumps(
            {
                "gene_id": profile.gene_id,
                "enzyme_id": profile.enzyme_id,
                "reaction_id": profile.reaction_id,
            },
            sort_keys=True,
        ),
        raw_sha256=_sha256_file(cached_path),
        license_id=profile.license_id,
        license_url=profile.license_url,
        retrieval_mode=RetrievalMode.MANUAL_IMPORT,
        cache_path=str(cached_path),
        adapter_id="pcsec_pichia.ecpichia.supplement_yaml_assessment",
        adapter_version="1",
        terms_reviewed=False,
        warnings=(
            "supplement_reuse_license_missing",
            "assessment_only_not_capacity_measurement",
            "formation_flux_conversion_missing",
        ),
    )
    source.validate()
    _validate_ecpichia_source_coherence(source, profile)
    _atomic_write_json(metadata_path, _json_ready(asdict(source)))
    return source


def _validate_pride_profile(profile: PrideMaxQuantSourceProfile) -> None:
    for field_name in (
        "project_accession",
        "protein_groups_filename",
        "sample_map_filename",
        "database_fasta_filename",
        "query_fasta_url",
        "query_protein_id",
        "external_gene_id",
        "external_protein_id",
        "metric_name",
        "species",
        "strain",
        "medium",
        "carbon_source",
        "culture_mode",
        "biomass_basis",
        "condition_source",
    ):
        _require_text(getattr(profile, field_name), field_name)
    if not profile.sample_ids:
        raise _validation_error("PRIDE profile requires sample_ids.")
    if profile.growth_rate_per_h <= 0 or profile.temperature_c <= 0 or profile.ph <= 0:
        raise _validation_error("PRIDE profile condition values must be positive.")


def _validate_ecpichia_profile(profile: EcPichiaSupplementSourceProfile) -> None:
    for field_name in (
        "source_id",
        "source_url",
        "source_version",
        "doi",
        "artifact_filename",
        "license_id",
        "license_url",
        "gene_id",
        "enzyme_id",
        "reaction_id",
    ):
        _require_text(getattr(profile, field_name), field_name)
    _require_sha256(profile.artifact_sha256, "ecPichia artifact_sha256")
    _require_sha256(
        profile.upstream_archive_sha256,
        "ecPichia upstream_archive_sha256",
    )


def _pride_project_url(accession: str) -> str:
    return f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}"


def _pride_files_url(accession: str) -> str:
    return f"{_pride_project_url(accession)}/files"


def _fetch_text(
    url: str,
    getter: Any,
    config: ExternalFetchConfig,
    label: str,
) -> ExternalHttpResponse:
    response = getter(url, config)
    if not 200 <= response.status_code < 300:
        raise _validation_error(
            f"{label} fetch failed with HTTP {response.status_code}."
        )
    if not response.text:
        raise _validation_error(f"{label} response is empty.")
    return response


def _as_mapping_json(text: str, label: str) -> Mapping[str, Any]:
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise _validation_error(f"{label} response is not valid JSON.") from exc
    if not isinstance(payload, Mapping):
        raise _validation_error(f"{label} response must be an object.")
    return payload


def _as_sequence_json(text: str, label: str) -> Sequence[Mapping[str, Any]]:
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise _validation_error(f"{label} response is not valid JSON.") from exc
    if not isinstance(payload, list) or any(
        not isinstance(item, Mapping) for item in payload
    ):
        raise _validation_error(f"{label} response must be an array of objects.")
    return payload


def _pride_file_entry(
    files: Sequence[Mapping[str, Any]],
    filename: str,
) -> Mapping[str, Any]:
    matches = tuple(item for item in files if item.get("fileName") == filename)
    if len(matches) != 1:
        raise _validation_error(
            f"PRIDE project must contain exactly one file named {filename}."
        )
    return matches[0]


def _pride_https_download_url(entry: Mapping[str, Any]) -> str:
    locations = entry.get("publicFileLocations")
    if not isinstance(locations, list):
        raise _validation_error("PRIDE file has no public locations.")
    for location in locations:
        if not isinstance(location, Mapping):
            continue
        value = str(location.get("value") or "")
        if value.startswith("ftp://ftp.pride.ebi.ac.uk/"):
            return "https://ftp.pride.ebi.ac.uk/" + value.split(
                "ftp://ftp.pride.ebi.ac.uk/", 1
            )[1]
    raise _validation_error("PRIDE file has no supported HTTPS download URL.")


def _normalize_pride_license(raw_license: str) -> tuple[str, str, bool]:
    if "CC0" in raw_license.upper() or "PUBLIC DOMAIN" in raw_license.upper():
        return (
            "CC0-1.0",
            "https://creativecommons.org/publicdomain/zero/1.0/",
            True,
        )
    return "unreviewed", "", False


def _validate_pride_bundle(
    bundle_path: Path,
    profile: PrideMaxQuantSourceProfile,
) -> None:
    bundle = _load_json_object(bundle_path)
    if bundle.get("schema_version") != 1:
        raise _validation_error("PRIDE source bundle requires schema_version=1.")
    if bundle.get("project_accession") != profile.project_accession:
        raise _validation_error("PRIDE source bundle project mismatch.")
    root = bundle_path.parent.resolve()
    records: list[Mapping[str, Any]] = []
    for key in ("project_metadata", "files_metadata"):
        value = bundle.get(key)
        if not isinstance(value, Mapping):
            raise _validation_error(f"PRIDE source bundle is missing {key}.")
        records.append(value)
    artifacts = bundle.get("artifacts")
    if not isinstance(artifacts, Mapping):
        raise _validation_error("PRIDE source bundle is missing artifacts.")
    for role in ("protein_groups", "sample_map", "database_fasta", "query_fasta"):
        value = artifacts.get(role)
        if not isinstance(value, Mapping):
            raise _validation_error(f"PRIDE source bundle is missing {role}.")
        records.append(value)
    for record in records:
        filename = str(record.get("filename") or "")
        if not filename or Path(filename).name != filename:
            raise _validation_error("PRIDE source artifact filename must be a basename.")
        path = (root / filename).resolve()
        if path.parent != root or not path.is_file():
            raise _validation_error("PRIDE source artifact is missing or outside cache.")
        expected_sha256 = str(record.get("sha256") or "")
        _require_sha256(expected_sha256, "PRIDE source artifact sha256")
        if _sha256_file(path) != expected_sha256.lower():
            raise _validation_error("PRIDE source artifact sha256 mismatch.")


def _validate_pride_source_coherence(
    source: ExternalCapacitySource,
    bundle_path: Path,
    profile: PrideMaxQuantSourceProfile,
) -> None:
    bundle = _load_json_object(bundle_path)
    raw_license = str(bundle.get("project_license_raw") or "")
    license_id, license_url, terms_reviewed = _normalize_pride_license(raw_license)
    publication_date = str(bundle.get("project_publication_date") or "")
    expected_version = (
        f"{profile.project_accession}:{publication_date}"
        if publication_date
        else "unversioned"
    )
    expected = {
        "source_id": f"pride:{profile.project_accession}",
        "source_type": ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS,
        "source_version": expected_version,
        "source_url": _pride_project_url(profile.project_accession),
        "license_id": license_id,
        "license_url": license_url,
        "adapter_id": "pcsec_pichia.pride.maxquant",
        "adapter_version": "1",
        "terms_reviewed": terms_reviewed,
    }
    for field_name, expected_value in expected.items():
        if getattr(source, field_name) != expected_value:
            raise _validation_error(
                f"PRIDE source metadata does not match bundle field {field_name}."
            )

    target = bundle.get("target")
    if not isinstance(target, Mapping):
        raise _validation_error("PRIDE source bundle is missing target metadata.")
    expected_target = {
        "query_protein_id": profile.query_protein_id,
        "external_gene_id": profile.external_gene_id,
        "external_protein_id": profile.external_protein_id,
        "metric_name": profile.metric_name,
        "sample_ids": list(profile.sample_ids),
    }
    if dict(target) != expected_target:
        raise _validation_error("PRIDE source bundle target metadata mismatch.")

    condition = bundle.get("condition")
    if not isinstance(condition, Mapping):
        raise _validation_error("PRIDE source bundle is missing condition metadata.")
    expected_condition = {
        "species": profile.species,
        "strain": profile.strain,
        "medium": profile.medium,
        "carbon_source": profile.carbon_source,
        "culture_mode": profile.culture_mode,
        "growth_rate_per_h": profile.growth_rate_per_h,
        "temperature_c": profile.temperature_c,
        "ph": profile.ph,
        "oxygen_condition": profile.oxygen_condition,
        "biomass_basis": profile.biomass_basis,
        "source_ref": profile.condition_source,
    }
    if dict(condition) != expected_condition:
        raise _validation_error("PRIDE source bundle condition metadata mismatch.")


def _validate_ecpichia_source_coherence(
    source: ExternalCapacitySource,
    profile: EcPichiaSupplementSourceProfile,
) -> None:
    expected = {
        "source_id": profile.source_id,
        "source_type": ExternalCapacitySourceType.EXTERNAL_ENZYME_MODEL,
        "source_version": profile.source_version,
        "source_url": profile.source_url,
        "raw_sha256": profile.artifact_sha256,
        "license_id": profile.license_id,
        "license_url": profile.license_url,
        "adapter_id": "pcsec_pichia.ecpichia.supplement_yaml_assessment",
        "adapter_version": "1",
        "terms_reviewed": False,
        "query": json.dumps(
            {
                "gene_id": profile.gene_id,
                "enzyme_id": profile.enzyme_id,
                "reaction_id": profile.reaction_id,
            },
            sort_keys=True,
        ),
    }
    for field_name, expected_value in expected.items():
        if getattr(source, field_name) != expected_value:
            raise _validation_error(
                f"ecPichia source metadata does not match reviewed {field_name}."
            )


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


def _atomic_copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=str(destination.parent),
    )
    os.close(fd)
    try:
        shutil.copyfile(source, temp_name)
        os.replace(temp_name, destination)
    except Exception:
        try:
            os.unlink(temp_name)
        except OSError:
            pass
        raise


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
    "ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE",
    "EcPichiaSupplementSourceProfile",
    "ExternalCapacitySource",
    "ExternalCapacitySourceValidationError",
    "ExternalCapacitySourceType",
    "PXD055501_G6PDH2_PROFILE",
    "PrideMaxQuantSourceProfile",
    "RetrievalMode",
    "cache_ecpichia_supplement_source",
    "cache_pride_maxquant_source",
    "cache_uniprot_identity_source",
]
