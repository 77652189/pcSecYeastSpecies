from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Mapping, Sequence

from pcsec_pichia.external_refs.capacity_sources import (
    ECPICHIA_G6PDH2_TABLE_PROFILE,
    ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE,
    ExternalCapacitySourceType,
    PXD055501_G6PDH2_PROFILE,
    RetrievalMode,
    cache_ecpichia_supplement_table_source,
    cache_ecpichia_supplement_source,
    cache_pride_maxquant_source,
    cache_uniprot_identity_source,
)
from pcsec_pichia.loading import load_pcsec_pichia_inputs
from pcsec_pichia.oe_capacity.external_candidate_evaluation import (
    build_capacity_candidate,
    build_capacity_model_binding,
    evaluate_ecpichia_g6pdh2_provenance,
)
from pcsec_pichia.oe_capacity.external_candidate_io import (
    CANDIDATE_MANIFEST_FILENAME,
    import_capacity_measurements,
    load_external_capacity_candidate_bundle,
    load_external_capacity_candidate_snapshot,
    parse_ecpichia_g6pdh2_table_evidence,
    parse_ecpichia_g6pdh2_source_assessment,
    parse_pride_maxquant_g6pdh2_evidence,
    write_external_capacity_candidate_cache,
)
from pcsec_pichia.oe_capacity.external_candidate_promotion import (
    PromotionDecision,
    build_capacity_promotion_manifest,
    promote_capacity_candidates,
)
from pcsec_pichia.oe_capacity.external_candidate_schema import (
    CapacityApplicabilityScope,
    CapacityCandidateStatus,
    CapacityConfidence,
    CapacityParameterKind,
    ExternalCapacityCandidateBundle,
    HostCondition,
)
from pcsec_pichia.oe_capacity.mapping import build_gene_enzyme_reaction_catalog
from pcsec_pichia.oe_capacity.parameters import (
    build_current_model_parameter_policy,
    load_capacity_anchor_catalog,
)
from pcsec_pichia.oe_capacity.schema import CapacityAnchorCatalog, ParameterPolicy
from pcsec_pichia.screens import prepare_screen_inputs
from pcsec_pichia.targets import load_builtin_targets


DEFAULT_TARGET_IDS = ("hLF", "OPN_ALPHA_FULL_PROJECT")
DEFAULT_CAPACITY_ASSET_PATH = Path("Enzymedata") / "oe_capacity_baseline_capacity.json"
G6PDH2_GENE_ID = "PAS_chr2-1_0308"
G6PDH2_ENZYME_ID = "G6PDH2_no_1_fwd_complex"
G6PDH2_FORMATION_HANDLE = "G6PDH2_no_1_fwd_complex_formation"


@dataclass(frozen=True)
class ExternalCapacityAuditRequest:
    repo_root: Path
    output_dir: Path
    offline_replay: bool = False
    identity_cache_dir: Path | None = None
    quantitative_cache_dir: Path | None = None
    ecpichia_cache_dir: Path | None = None
    ecpichia_table_cache_dir: Path | None = None
    pride_pxd055501: bool = False
    ecpichia_assessment: bool = False
    ecpichia_provenance_closure: bool = False
    ecpichia_supplement_file: Path | None = None
    ecpichia_table_file: Path | None = None
    measurement_file: Path | None = None
    source_id: str = ""
    source_type: ExternalCapacitySourceType = (
        ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS
    )
    source_version: str = ""
    source_url: str = ""
    license_id: str = ""
    license_url: str = ""
    query: str = f"{G6PDH2_GENE_ID} glucose mu=0.1"
    expected_sha256: str = ""
    terms_reviewed: bool = False
    target_ids: tuple[str, ...] = DEFAULT_TARGET_IDS
    carbon_source_id: str = "glucose"
    growth_rate: float = 0.1
    relative_uncertainty: float = 0.2


@dataclass(frozen=True)
class ExternalCapacityAuditOutputs:
    audit_json_path: Path
    audit_markdown_path: Path
    candidate_manifest_path: Path
    candidate_count: int
    promotion_ready_count: int
    completion_outcome: str = "in_progress"
    provenance_gap_json_path: Path | None = None
    provenance_gap_markdown_path: Path | None = None

    def summary(self) -> dict[str, object]:
        return {
            "audit": str(self.audit_json_path),
            "candidate_count": self.candidate_count,
            "promotion_ready_count": self.promotion_ready_count,
            "formal_promotion_performed": False,
            "completion_outcome": self.completion_outcome,
            "provenance_gap": (
                str(self.provenance_gap_json_path)
                if self.provenance_gap_json_path is not None
                else None
            ),
        }


def capacity_asset_version(
    repo_root: str | Path,
    asset_path: str | Path = DEFAULT_CAPACITY_ASSET_PATH,
) -> str:
    path = _resolve_repo_path(repo_root, asset_path)
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "missing"


def prepare_external_candidate_runtime(
    repo_root: str | Path,
    *,
    target_id: str,
    growth_rate: float = 0.1,
    carbon_source_id: str = "glucose",
    relative_uncertainty: float = 0.2,
    capacity_asset_path: str | Path = DEFAULT_CAPACITY_ASSET_PATH,
    expected_asset_sha256: str | None = None,
) -> SimpleNamespace:
    root = Path(repo_root).resolve()
    targets = {
        target.target_id: target
        for target in load_builtin_targets(root)
        if target.target_id in DEFAULT_TARGET_IDS
    }
    try:
        target = targets[target_id]
    except KeyError as exc:
        raise KeyError(f"unknown OE capacity target: {target_id}") from exc
    inputs = load_pcsec_pichia_inputs(root, carbon_source_id=carbon_source_id)
    prepared = prepare_screen_inputs(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_rate,
    )
    if not prepared.get("baseline_success"):
        raise RuntimeError(
            "target baseline preparation failed: "
            + str(
                prepared.get("baseline_status")
                or prepared.get("build_status")
                or "unknown"
            )
        )
    catalog = build_gene_enzyme_reaction_catalog(
        prepared["fixed_model"],
        inputs.metabolic,
        prepared["combined"],
    )
    resolved_asset_path = _resolve_repo_path(root, capacity_asset_path)
    anchor_catalog, asset_metadata = load_capacity_asset_snapshot(resolved_asset_path)
    loaded_sha = str(asset_metadata.get("sha256") or "missing")
    if expected_asset_sha256 is not None and loaded_sha != expected_asset_sha256:
        raise RuntimeError(
            "OE capacity asset changed during runtime preparation; retry the run."
        )
    context_id = _context_id(carbon_source_id, growth_rate)
    parameter_policy = build_current_model_parameter_policy(
        catalog,
        prepared["combined"],
        capacity_anchors=anchor_catalog,
        target_id=target_id,
        context_id=context_id,
        relative_uncertainty=relative_uncertainty,
    )
    return SimpleNamespace(
        target_id=target_id,
        fixed_model=prepared["fixed_model"],
        exchange_reaction_id=prepared["exchange_reaction_id"],
        metabolic=inputs.metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        gene_capacity_catalog=catalog,
        parameter_policy=ParameterPolicy(
            parameter_sets=parameter_policy.parameter_sets,
            scenarios=parameter_policy.scenarios,
            strict_conflicts=parameter_policy.strict_conflicts,
        ),
        capacity_asset_version=loaded_sha,
        capacity_asset_metadata=asset_metadata,
        capacity_anchor_catalog=anchor_catalog,
    )


def load_capacity_asset_snapshot(
    path: str | Path,
) -> tuple[CapacityAnchorCatalog, dict[str, Any]]:
    asset_path = Path(path)
    if not asset_path.is_file():
        catalog = CapacityAnchorCatalog(
            model_fingerprint="missing-capacity-asset",
            anchors=(),
            source_ref=str(asset_path),
        )
        return catalog, {
            "path": str(asset_path),
            "version": "missing",
            "sha256": "",
            "reviewed": False,
        }
    raw = asset_path.read_bytes()
    payload = json.loads(raw.decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("OE capacity asset root must be a JSON object.")
    catalog = load_capacity_anchor_catalog(payload)
    return catalog, {
        "path": str(asset_path),
        "version": str(payload.get("asset_version") or ""),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "reviewed": bool(catalog.anchors),
    }


def load_external_candidate_review(
    repo_root: str | Path,
    candidate_root: str | Path,
    *,
    target_id: str,
    growth_rate: float = 0.1,
    carbon_source_id: str = "glucose",
    capacity_asset_path: str | Path = DEFAULT_CAPACITY_ASSET_PATH,
) -> dict[str, Any]:
    root = Path(candidate_root)
    manifest_path = root if root.is_file() else root / CANDIDATE_MANIFEST_FILENAME
    asset_path = _resolve_repo_path(repo_root, capacity_asset_path)
    context_id = _context_id(carbon_source_id, growth_rate)
    asset_sha = capacity_asset_version(repo_root, capacity_asset_path)
    if not manifest_path.is_file():
        return {
            "available": False,
            "candidate_root": str(root),
            "target_id": target_id,
            "context_id": context_id,
            "candidates": [],
            "message": "尚无 Round 6A 外部容量候选 cache；正式容量资产保持不变。",
            "formal_asset_sha256": asset_sha,
        }
    try:
        candidate_snapshot = load_external_capacity_candidate_snapshot(manifest_path)
        bundle = candidate_snapshot.bundle
    except Exception as exc:
        return {
            "available": False,
            "candidate_root": str(root),
            "target_id": target_id,
            "context_id": context_id,
            "candidates": [],
            "message": f"候选 cache 校验失败：{type(exc).__name__}: {exc}",
            "formal_asset_sha256": asset_sha,
        }
    candidates: list[dict[str, Any]] = []
    for candidate in bundle.candidates:
        condition_matches = (
            candidate.condition.carbon_source.strip().lower()
            == str(carbon_source_id).strip().lower()
            and abs(candidate.condition.growth_rate_per_h - float(growth_rate))
            <= 1e-12
        )
        binding_matches = any(
            binding.target_id == target_id and binding.context_id == context_id
            for binding in candidate.model_bindings
        )
        if candidate.applicability_scope is CapacityApplicabilityScope.TARGET_SPECIFIC:
            applicable = candidate.target_id == target_id and condition_matches
        else:
            applicable = condition_matches
        if not (applicable and binding_matches):
            continue
        payload = _json_ready(asdict(candidate))
        payload["promotion_eligible"] = (
            candidate.status is CapacityCandidateStatus.REVIEW_READY
            and candidate.applicability_scope
            is not CapacityApplicabilityScope.HOMOLOG_TRANSFERRED
        )
        candidates.append(payload)
    return {
        "available": True,
        "candidate_root": str(root),
        "candidate_manifest_path": str(manifest_path),
        "candidate_manifest_sha256": candidate_snapshot.manifest_sha256,
        "formal_asset_path": str(asset_path),
        "formal_asset_sha256": asset_sha,
        "model_fingerprints": list(bundle.model_fingerprints),
        "target_id": target_id,
        "context_id": context_id,
        "candidates": candidates,
        "source_count": len(bundle.sources),
        "measurement_count": len(bundle.measurements),
        "sources": [_json_ready(asdict(item)) for item in bundle.sources],
        "message": (
            "外部候选仅供审核，不是 reviewed anchor；只有显式批准并通过 hash/model "
            "校验后才会写入正式资产。"
        ),
    }


def preview_external_candidate_promotion(
    repo_root: str | Path,
    *,
    candidate_root: str | Path,
    candidate_ids: Sequence[str],
    target_id: str,
    expected_candidate_manifest_sha256: str,
    expected_asset_sha256: str,
    capacity_asset_path: str | Path = DEFAULT_CAPACITY_ASSET_PATH,
) -> dict[str, Any]:
    root = Path(candidate_root)
    manifest_path = root if root.is_file() else root / CANDIDATE_MANIFEST_FILENAME
    if capacity_asset_version(repo_root, capacity_asset_path) != expected_asset_sha256:
        raise ValueError("formal capacity asset changed after review; reload before promotion.")
    bundle = load_external_capacity_candidate_bundle(
        manifest_path,
        expected_manifest_sha256=expected_candidate_manifest_sha256,
    )
    promotion = build_capacity_promotion_manifest(
        bundle,
        candidate_ids,
        _resolve_repo_path(repo_root, capacity_asset_path),
        expected_asset_sha256=expected_asset_sha256,
    )
    selected = {
        item.candidate_id: item
        for item in bundle.candidates
        if item.candidate_id in promotion.candidate_ids
    }
    return {
        "decision": promotion.decision.value,
        "candidate_ids": list(promotion.candidate_ids),
        "candidate_manifest_sha256": expected_candidate_manifest_sha256,
        "formal_asset_sha256": expected_asset_sha256,
        "target_id": target_id,
        "model_fingerprints": list(bundle.model_fingerprints),
        "eligible": all(
            item.status is CapacityCandidateStatus.REVIEW_READY
            and item.applicability_scope
            is not CapacityApplicabilityScope.HOMOLOG_TRANSFERRED
            for item in selected.values()
        ),
        "additions": [_json_ready(asdict(item)) for item in selected.values()],
        "warning": "预览不会写入 Enzymedata，也不会启动正式 acceptance。",
    }


def promote_external_candidate_selection(
    repo_root: str | Path,
    *,
    candidate_root: str | Path,
    candidate_ids: Sequence[str],
    reviewer: str,
    expected_candidate_manifest_sha256: str,
    expected_asset_sha256: str,
    explicit_approval: bool,
    capacity_asset_path: str | Path = DEFAULT_CAPACITY_ASSET_PATH,
    runtime_resolver: Callable[[str, float, str, float, str], SimpleNamespace]
    | None = None,
) -> dict[str, Any]:
    if explicit_approval is not True:
        raise ValueError("explicit_approval=True is required for formal capacity promotion.")
    reviewer_text = str(reviewer or "").strip()
    if not reviewer_text:
        raise ValueError("reviewer must be non-empty.")
    root = Path(candidate_root)
    manifest_path = root if root.is_file() else root / CANDIDATE_MANIFEST_FILENAME
    if capacity_asset_version(repo_root, capacity_asset_path) != expected_asset_sha256:
        raise ValueError("formal capacity asset changed after review; promotion refused.")
    bundle = load_external_capacity_candidate_bundle(
        manifest_path,
        expected_manifest_sha256=expected_candidate_manifest_sha256,
    )
    selected_ids = set(candidate_ids)
    selected_candidates = tuple(
        candidate
        for candidate in bundle.candidates
        if candidate.candidate_id in selected_ids
    )
    catalogs: dict[str, Any] = {}
    for candidate in selected_candidates:
        for binding in candidate.model_bindings:
            if runtime_resolver is None:
                runtime = prepare_external_candidate_runtime(
                    repo_root,
                    target_id=binding.target_id,
                    growth_rate=candidate.condition.growth_rate_per_h,
                    carbon_source_id=candidate.condition.carbon_source,
                    relative_uncertainty=0.2,
                    capacity_asset_path=capacity_asset_path,
                    expected_asset_sha256=expected_asset_sha256,
                )
            else:
                runtime = runtime_resolver(
                    binding.target_id,
                    float(candidate.condition.growth_rate_per_h),
                    candidate.condition.carbon_source,
                    0.2,
                    expected_asset_sha256,
                )
            catalogs[binding.target_id] = runtime.gene_capacity_catalog
    promotion = build_capacity_promotion_manifest(
        bundle,
        candidate_ids,
        _resolve_repo_path(repo_root, capacity_asset_path),
        decision=PromotionDecision.APPROVED,
        reviewer=reviewer_text,
        reviewed_at=datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        reason="Explicit Round 6A reviewer approval via public audit API.",
        expected_asset_sha256=expected_asset_sha256,
    )
    promoted = promote_capacity_candidates(bundle, promotion, catalogs=catalogs)
    return {
        "decision": promoted.decision.value,
        "candidate_ids": list(promoted.candidate_ids),
        "asset_path": promoted.asset_path,
        "promoted_asset_sha256": promoted.promoted_asset_sha256,
        "reviewer": promoted.reviewer,
        "reviewed_at": promoted.reviewed_at,
        "acceptance_started": False,
    }


def run_external_capacity_candidate_audit(
    request: ExternalCapacityAuditRequest,
) -> ExternalCapacityAuditOutputs:
    output_dir = Path(request.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if request.pride_pxd055501 and request.measurement_file is not None:
        raise ValueError(
            "pride_pxd055501 and measurement_file are mutually exclusive source paths."
        )
    identity_source = cache_uniprot_identity_source(
        G6PDH2_GENE_ID,
        request.identity_cache_dir or (output_dir / "identity_source"),
        retrieval_mode=(
            RetrievalMode.OFFLINE_REPLAY
            if request.offline_replay
            else RetrievalMode.ONLINE
        ),
    )
    capacity_sources = []
    measurements = ()
    source_assessments: list[dict[str, Any]] = []
    pride_evidence = None
    ecpichia_evidence = None
    ecpichia_table_evidence = None
    provenance_closure = None
    runtime_condition = _default_host_condition(request)
    if request.pride_pxd055501:
        pride_source = cache_pride_maxquant_source(
            PXD055501_G6PDH2_PROFILE,
            request.quantitative_cache_dir
            or (output_dir / "quantitative_source"),
            retrieval_mode=(
                RetrievalMode.OFFLINE_REPLAY
                if request.offline_replay
                else RetrievalMode.ONLINE
            ),
        )
        pride_evidence = parse_pride_maxquant_g6pdh2_evidence(pride_source)
        capacity_sources.append(pride_source)
        measurements = (pride_evidence.measurement,)
        runtime_condition = pride_evidence.measurement.condition
        source_assessments.append(
            {
                "assessment_id": "pride-pxd055501-relative-ibaq",
                "source_kind": "relative_abundance",
                "source_id": pride_source.source_id,
                "project_accession": PXD055501_G6PDH2_PROFILE.project_accession,
                "raw_metric": PXD055501_G6PDH2_PROFILE.metric_name,
                "raw_values": list(pride_evidence.raw_values),
                "raw_unit": pride_evidence.measurement.unit,
                "sample_ids": list(pride_evidence.sample_ids),
                "source_bundle_sha256": pride_evidence.source_bundle_sha256,
                "artifact_sha256s": dict(pride_evidence.artifact_sha256s),
                "source_version": pride_source.source_version,
                "source_url": pride_source.source_url,
                "license_id": pride_source.license_id,
                "license_url": pride_source.license_url,
                "mapping_evidence": list(pride_evidence.mapping_evidence),
                "condition": asdict(pride_evidence.measurement.condition),
                "condition_source_ref": PXD055501_G6PDH2_PROFILE.condition_source,
                "quantitative_boundary": dict(pride_evidence.quantitative_boundary),
                "formal_context_id": _context_id(
                    request.carbon_source_id, request.growth_rate
                ),
                "source_context_id": _context_id(
                    runtime_condition.carbon_source,
                    runtime_condition.growth_rate_per_h,
                ),
                "formal_context_match": (
                    runtime_condition.carbon_source.strip().lower()
                    == request.carbon_source_id.strip().lower()
                    and abs(
                        runtime_condition.growth_rate_per_h - request.growth_rate
                    )
                    <= 1e-12
                ),
                "promotion_ready": False,
                "assessment": "quantitative_relative_only",
                "quantitative_value_available": True,
                "target_value_available": True,
                "capacity_value_available": False,
                "model_flux_conversion_available": False,
                "missing_information": [
                    "absolute_abundance_calibration",
                    "biomass_normalization",
                    "paired_condition_matched_kcat",
                    "formal_glucose_mu_0.1_condition_match",
                ],
            }
        )
    if request.ecpichia_assessment or request.ecpichia_provenance_closure:
        ecpichia_source = cache_ecpichia_supplement_source(
            ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE,
            request.ecpichia_cache_dir or (output_dir / "ecpichia_source"),
            source_file=request.ecpichia_supplement_file,
            retrieval_mode=(
                RetrievalMode.OFFLINE_REPLAY
                if request.offline_replay
                else RetrievalMode.MANUAL_IMPORT
            ),
        )
        ecpichia_evidence = parse_ecpichia_g6pdh2_source_assessment(
            ecpichia_source
        )
        capacity_sources.append(ecpichia_source)
        source_assessments.append(
            {
                "assessment_id": "ecpichia-supplementary-8-g6pdh2",
                "source_kind": "external_enzyme_model_assessment",
                "source_id": ecpichia_source.source_id,
                "source_version": ecpichia_source.source_version,
                "source_url": ecpichia_source.source_url,
                "doi": ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE.doi,
                "license_id": ecpichia_source.license_id,
                "license_url": ecpichia_source.license_url,
                "terms_reviewed": ecpichia_source.terms_reviewed,
                "raw_artifact_sha256": ecpichia_source.raw_sha256,
                "upstream_archive_sha256": (
                    ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE.upstream_archive_sha256
                ),
                "raw_fields": [
                    "gene_id",
                    "enzyme_id",
                    "reaction_id",
                    "ec_number",
                    "molecular_weight_g_per_mol",
                    "kcat_per_s",
                    "kcat_source_label",
                    "reaction_protein_coefficient",
                    "reported_concentration",
                    "usage_lower_bound",
                    "protein_pool_lower_bound",
                ],
                "raw_values": {
                    "gene_id": ecpichia_evidence.gene_id,
                    "enzyme_id": ecpichia_evidence.enzyme_id,
                    "reaction_id": ecpichia_evidence.reaction_id,
                    "ec_number": ecpichia_evidence.ec_number,
                    "molecular_weight_g_per_mol": (
                        ecpichia_evidence.molecular_weight_g_per_mol
                    ),
                    "kcat_per_s": ecpichia_evidence.kcat_per_s,
                    "kcat_source_label": ecpichia_evidence.kcat_source_label,
                    "reaction_protein_coefficient": (
                        ecpichia_evidence.reaction_protein_coefficient
                    ),
                    "reported_concentration": (
                        ecpichia_evidence.reported_concentration
                    ),
                    "usage_lower_bound": ecpichia_evidence.usage_lower_bound,
                    "protein_pool_lower_bound": (
                        ecpichia_evidence.protein_pool_lower_bound
                    ),
                },
                "raw_units": {
                    "molecular_weight": "g_per_mol",
                    "kcat": "per_s",
                    "reported_concentration": (
                        ecpichia_evidence.reported_concentration_unit
                    ),
                    "gecko_expected_concentration": (
                        ecpichia_evidence.gecko_expected_concentration_unit
                    ),
                },
                "condition": {
                    "species": "Komagataella phaffii",
                    "carbon_source": None,
                    "culture_mode": None,
                    "growth_rate_per_h": None,
                    "temperature_c": None,
                    "ph": None,
                    "oxygen_condition": None,
                    "source_ref": ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE.doi,
                    "evidence_status": "bibliographic_not_embedded_in_yaml",
                },
                "mapping_evidence": list(ecpichia_evidence.mapping_evidence),
                "conflicts": list(ecpichia_evidence.conflicts),
                "quantitative_value_available": True,
                "target_value_available": True,
                "capacity_value_available": False,
                "model_flux_conversion_available": False,
                "promotion_ready": False,
                "assessment": "review_required_source_conflicts",
                "missing_information": list(
                    ecpichia_evidence.missing_information
                ),
            }
        )
    if request.ecpichia_provenance_closure:
        ecpichia_table_source = cache_ecpichia_supplement_table_source(
            ECPICHIA_G6PDH2_TABLE_PROFILE,
            request.ecpichia_table_cache_dir or (output_dir / "ecpichia_table_source"),
            source_file=request.ecpichia_table_file,
            retrieval_mode=(
                RetrievalMode.OFFLINE_REPLAY
                if request.offline_replay
                else RetrievalMode.MANUAL_IMPORT
            ),
        )
        ecpichia_table_evidence = parse_ecpichia_g6pdh2_table_evidence(
            ecpichia_table_source
        )
        capacity_sources.append(ecpichia_table_source)
    asset_sha = capacity_asset_version(request.repo_root)
    bindings = []
    fingerprints = []
    catalogs = {}
    for target_id in request.target_ids:
        runtime = prepare_external_candidate_runtime(
            request.repo_root,
            target_id=target_id,
            growth_rate=runtime_condition.growth_rate_per_h,
            carbon_source_id=runtime_condition.carbon_source,
            relative_uncertainty=request.relative_uncertainty,
            expected_asset_sha256=asset_sha,
        )
        binding = build_capacity_model_binding(
            runtime.gene_capacity_catalog,
            target_id=target_id,
            context_id=_context_id(
                runtime_condition.carbon_source,
                runtime_condition.growth_rate_per_h,
            ),
            gene_id=G6PDH2_GENE_ID,
            external_gene_id=(
                pride_evidence.measurement.external_gene_id
                if pride_evidence is not None
                else G6PDH2_GENE_ID
            ),
            external_protein_id=(
                pride_evidence.measurement.external_protein_id
                if pride_evidence is not None
                else "C4R099"
            ),
            mapping_evidence=(
                "uniprot:C4R099",
                "ncbi_gene:8198996",
                "kegg:ppa:PAS_chr2-1_0308",
                "refseq:XP_002491203.1",
                "ec:1.1.1.49",
                *(
                    pride_evidence.mapping_evidence
                    if pride_evidence is not None
                    else ()
                ),
            ),
        )
        bindings.append(binding)
        fingerprints.append(binding.model_fingerprint)
        catalogs[target_id] = runtime.gene_capacity_catalog
    if request.ecpichia_provenance_closure:
        formal_context_id = _context_id(request.carbon_source_id, request.growth_rate)
        closure_mapping_evidence = (
            "ecpichia_yaml_gene:PAS_chr2-1_0308",
            "ecpichia_yaml_enzyme:C4R099",
            "ecpichia_yaml_reaction:G6PDH2",
            "current_model_formation:G6PDH2_no_1_fwd_complex_formation",
        )
        if all(binding.context_id == formal_context_id for binding in bindings):
            formal_bindings = [
                replace(
                    binding,
                    external_gene_id=G6PDH2_GENE_ID,
                    external_protein_id="C4R099",
                    external_enzyme_id="C4R099",
                    mapping_evidence=tuple(
                        dict.fromkeys(
                            (*binding.mapping_evidence, *closure_mapping_evidence)
                        )
                    ),
                )
                for binding in bindings
            ]
        else:
            formal_bindings = []
            for target_id in request.target_ids:
                formal_runtime = prepare_external_candidate_runtime(
                    request.repo_root,
                    target_id=target_id,
                    growth_rate=request.growth_rate,
                    carbon_source_id=request.carbon_source_id,
                    relative_uncertainty=request.relative_uncertainty,
                    expected_asset_sha256=asset_sha,
                )
                formal_bindings.append(
                    build_capacity_model_binding(
                        formal_runtime.gene_capacity_catalog,
                        target_id=target_id,
                        context_id=formal_context_id,
                        gene_id=G6PDH2_GENE_ID,
                        external_gene_id=G6PDH2_GENE_ID,
                        external_protein_id="C4R099",
                        external_enzyme_id="C4R099",
                        mapping_evidence=closure_mapping_evidence,
                    )
                )
        provenance_closure = evaluate_ecpichia_g6pdh2_provenance(
            ecpichia_evidence,
            ecpichia_table_evidence,
            model_bindings=tuple(formal_bindings),
            formal_context_id=formal_context_id,
        )
        source_assessments.append(
            {
                "assessment_id": "a0c-ecpichia-g6pdh2-provenance-closure",
                "assessment": provenance_closure.completion_outcome,
                "source_id": ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE.source_id,
                "source_version": ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE.source_version,
                "source_url": ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE.source_url,
                "license_id": ECPICHIA_G6PDH2_TABLE_PROFILE.license_id,
                "license_url": ECPICHIA_G6PDH2_TABLE_PROFILE.license_url,
                "raw_values": {
                    "yaml": {
                        "gene_id": ecpichia_evidence.gene_id,
                        "enzyme_id": ecpichia_evidence.enzyme_id,
                        "molecular_weight_g_per_mol": ecpichia_evidence.molecular_weight_g_per_mol,
                        "kcat_per_s": ecpichia_evidence.kcat_per_s,
                        "reported_concentration": ecpichia_evidence.reported_concentration,
                    },
                    "published_table": {
                        "gene_id": ecpichia_table_evidence.gene_id,
                        "enzyme_id": ecpichia_table_evidence.enzyme_id,
                        "molecular_weight_g_per_mol": ecpichia_table_evidence.molecular_weight_g_per_mol,
                        "kcat_per_s": ecpichia_table_evidence.kcat_per_s,
                        "reported_concentration": ecpichia_table_evidence.reported_concentration_text,
                    },
                },
                "raw_units": {
                    "yaml_concentration": ecpichia_evidence.reported_concentration_unit,
                    "published_table_concentration": ecpichia_table_evidence.reported_concentration_unit,
                },
                "condition": {
                    "formal_context_id": formal_context_id,
                    "source_condition": "not_embedded_in_published_artifacts",
                },
                "mapping_evidence": [asdict(item) for item in formal_bindings],
                "coefficient_trace": dict(provenance_closure.coefficient_trace),
                "unit_trace": list(provenance_closure.conditional_unit_trace),
                "conflicts": list(provenance_closure.source_conflicts),
                "missing_information": list(provenance_closure.missing_information),
                "nominal_capacity": None,
                "promotion_preview_available": False,
                "promotion_ready": False,
                "capacity_value_available": False,
                "model_flux_conversion_available": False,
            }
        )
    candidates = ()
    if pride_evidence is not None:
        candidate = build_capacity_candidate(
            candidate_id="g6pdh2-pride-pxd055501-t0-ibaq",
            applicability_scope=CapacityApplicabilityScope.HOST_CONDITION,
            model_bindings=tuple(bindings),
            catalogs=catalogs,
            sources=tuple(capacity_sources),
            condition=pride_evidence.measurement.condition,
            abundance=pride_evidence.measurement,
            confidence=CapacityConfidence.LOW,
        )
        formal_context_matches = source_assessments[0]["formal_context_match"] is True
        if not formal_context_matches:
            candidate = replace(
                candidate,
                status=CapacityCandidateStatus.REVIEW_REQUIRED,
                missing_information=tuple(
                    dict.fromkeys(
                        (
                            *candidate.missing_information,
                            "formal_glucose_mu_0.1_condition_match",
                        )
                    )
                ),
                warnings=tuple(
                    dict.fromkeys(
                        (
                            *candidate.warnings,
                            "source_context_glucose_mu_0.075_not_formal_glucose_mu_0.1",
                        )
                    )
                ),
            )
            candidate.validate()
        candidates = (candidate,)
    elif request.measurement_file is not None:
        if not request.source_id.strip():
            raise ValueError("source_id is required with measurement_file.")
        capacity_source, measurements = import_capacity_measurements(
            request.measurement_file,
            source_id=request.source_id,
            source_type=request.source_type,
            source_version=request.source_version,
            source_url=request.source_url,
            license_id=request.license_id,
            license_url=request.license_url,
            query=request.query,
            output_dir=output_dir / "capacity_source",
            expected_sha256=request.expected_sha256,
            terms_reviewed=request.terms_reviewed,
        )
        capacity_sources.append(capacity_source)
        by_kind = {
            kind: tuple(
                item for item in measurements if item.parameter_kind is kind
            )
            for kind in CapacityParameterKind
        }
        duplicates = {
            kind.value: len(rows)
            for kind, rows in by_kind.items()
            if len(rows) > 1
        }
        if duplicates:
            raise ValueError(
                "Round 6A v1 accepts at most one measurement per parameter kind: "
                + json.dumps(duplicates, sort_keys=True)
            )
        candidate = build_capacity_candidate(
            candidate_id=f"g6pdh2-{request.source_id}",
            applicability_scope=CapacityApplicabilityScope.HOST_CONDITION,
            model_bindings=tuple(bindings),
            catalogs=catalogs,
            sources=tuple(capacity_sources),
            condition=(
                measurements[0].condition
                if measurements
                else _default_host_condition(request)
            ),
            abundance=(
                by_kind[CapacityParameterKind.ABUNDANCE][0]
                if by_kind[CapacityParameterKind.ABUNDANCE]
                else None
            ),
            kcat=(
                by_kind[CapacityParameterKind.KCAT][0]
                if by_kind[CapacityParameterKind.KCAT]
                else None
            ),
            direct_capacity=(
                by_kind[CapacityParameterKind.BASELINE_CAPACITY][0]
                if by_kind[CapacityParameterKind.BASELINE_CAPACITY]
                else None
            ),
            confidence=CapacityConfidence.MEDIUM,
        )
        candidates = (candidate,)
    bundle = ExternalCapacityCandidateBundle(
        model_fingerprints=tuple(fingerprints),
        sources=(identity_source, *capacity_sources),
        measurements=measurements,
        candidates=candidates,
    )
    candidate_outputs = write_external_capacity_candidate_cache(
        bundle, output_dir / "candidates"
    )
    promotion_ready_count = sum(
        1 for item in candidates if item.status is CapacityCandidateStatus.REVIEW_READY
    )
    audit = _build_audit_payload(
        request=request,
        identity_source=identity_source,
        capacity_sources=capacity_sources,
        measurements=measurements,
        candidates=candidates,
        bindings=bindings,
        source_assessments=source_assessments,
        candidate_manifest=candidate_outputs.manifest_path,
        candidate_manifest_sha256=candidate_outputs.bundle_sha256,
    )
    audit_json_path = output_dir / "g6pdh2_capacity_candidate_audit.json"
    audit_json_path.write_text(
        json.dumps(
            audit,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            default=_json_default,
        )
        + "\n",
        encoding="utf-8",
    )
    audit_markdown_path = output_dir / "g6pdh2_capacity_candidate_audit.md"
    audit_markdown_path.write_text(
        _render_audit_markdown(candidates, source_assessments), encoding="utf-8"
    )
    gap_json_path = None
    gap_markdown_path = None
    completion_outcome = "in_progress"
    if provenance_closure is not None:
        completion_outcome = provenance_closure.completion_outcome
        gap_json_path = output_dir / "g6pdh2_ecpichia_provenance_gap.json"
        gap_json_path.write_text(
            json.dumps(provenance_closure.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        gap_markdown_path = output_dir / "g6pdh2_ecpichia_provenance_gap.md"
        gap_markdown_path.write_text(
            _render_ecpichia_provenance_gap(provenance_closure.to_dict()),
            encoding="utf-8",
        )
    return ExternalCapacityAuditOutputs(
        audit_json_path=audit_json_path,
        audit_markdown_path=audit_markdown_path,
        candidate_manifest_path=candidate_outputs.manifest_path,
        candidate_count=len(candidates),
        promotion_ready_count=promotion_ready_count,
        completion_outcome=completion_outcome,
        provenance_gap_json_path=gap_json_path,
        provenance_gap_markdown_path=gap_markdown_path,
    )


def _build_audit_payload(
    *,
    request: ExternalCapacityAuditRequest,
    identity_source: Any,
    capacity_sources: Sequence[Any],
    measurements: Sequence[Any],
    candidates: Sequence[Any],
    bindings: Sequence[Any],
    source_assessments: Sequence[Mapping[str, Any]],
    candidate_manifest: Path,
    candidate_manifest_sha256: str,
) -> dict[str, Any]:
    pride_assessment = next(
        (
            item
            for item in source_assessments
            if item.get("assessment") == "quantitative_relative_only"
        ),
        None,
    )
    ecpichia_assessment = next(
        (
            item
            for item in source_assessments
            if item.get("assessment_id") == "ecpichia-supplementary-8-g6pdh2"
        ),
        None,
    )
    assessment_missing = tuple(
        dict.fromkeys(
            str(value)
            for item in source_assessments
            for value in item.get("missing_information") or ()
        )
    )
    candidate_missing = tuple(
        dict.fromkeys(
            value for item in candidates for value in item.missing_information
        )
    )
    assessment_conflicts = tuple(
        dict.fromkeys(
            str(value)
            for item in source_assessments
            for value in item.get("conflicts") or ()
        )
    )
    has_review_ready_candidate = any(
        item.status is CapacityCandidateStatus.REVIEW_READY for item in candidates
    )
    provenance_assessment = next(
        (
            item
            for item in source_assessments
            if item.get("assessment_id")
            == "a0c-ecpichia-g6pdh2-provenance-closure"
        ),
        None,
    )
    return {
        "schema_version": 1,
        "round": "round_6a_external_capacity_candidates",
        "gene_id": G6PDH2_GENE_ID,
        "enzyme_id": G6PDH2_ENZYME_ID,
        "formation_handle": G6PDH2_FORMATION_HANDLE,
        "targets": list(request.target_ids),
        "context_id": _context_id(request.carbon_source_id, request.growth_rate),
        "completion_outcome": (
            provenance_assessment.get("assessment")
            if provenance_assessment is not None
            else "in_progress"
        ),
        "source_assessments": [dict(item) for item in source_assessments],
        "formal_model_bindings": (
            list(provenance_assessment.get("mapping_evidence") or ())
            if provenance_assessment is not None
            else []
        ),
        "sources_checked": [
            {
                **asdict(identity_source),
                "capacity_value_available": False,
                "role": "identity_crosswalk_only",
            },
            {
                "source": "same-host quantitative Pichia proteomics",
                "status": (
                    "parsed_relative_quantitative_evidence"
                    if pride_assessment is not None
                    else "manual_import_required"
                ),
                "quantitative_value_available": pride_assessment is not None,
                "target_value_available": pride_assessment is not None,
                "capacity_value_available": False,
            },
            {
                "source": "iPichia/ecPichia calibrated enzyme capacity",
                "status": (
                    "parsed_raw_values_review_required"
                    if ecpichia_assessment is not None
                    else "formal_file_import_required"
                ),
                "quantitative_value_available": ecpichia_assessment is not None,
                "target_value_available": ecpichia_assessment is not None,
                "capacity_value_available": False,
            },
            {
                "source": "Pichia literature plus BRENDA/SABIO-RK kinetics",
                "status": "abundance_and_condition_pair_missing",
                "capacity_value_available": False,
            },
            *(
                {
                    **asdict(source),
                    "capacity_value_available": any(
                        item.source_id == source.source_id for item in measurements
                    ),
                    "role": "external_source_or_measurement",
                }
                for source in capacity_sources
            ),
        ],
        "model_bindings": [asdict(item) for item in bindings],
        "candidate_count": len(candidates),
        "promotion_ready_count": sum(
            1
            for item in candidates
            if item.status is CapacityCandidateStatus.REVIEW_READY
        ),
        "candidates": [asdict(item) for item in candidates],
        "rejection_reasons": (
            [
                "UniProt confirms identity but contains no baseline capacity value.",
                "No reviewed same-host abundance or direct capacity measurement was supplied.",
                "No traceable abundance x kcat conversion chain is available for glucose mu=0.1.",
                *assessment_conflicts,
            ]
            if not candidates
            else list(
                dict.fromkeys(
                    (
                        *(
                            reason
                            for item in candidates
                            for reason in (*item.rejection_reasons, *item.conflicts)
                        ),
                        *assessment_conflicts,
                    )
                )
            )
        ),
        "missing_information": list(
            dict.fromkeys(
                (
                    *candidate_missing,
                    *assessment_missing,
                    *(
                        ()
                        if has_review_ready_candidate
                        else (
                            "reviewed_absolute_condition_matched_g6pdh2_abundance",
                            "direct_or_reviewed_komagataella_g6pdh2_kcat",
                            "reviewed_model_flux_conversion_chain",
                        )
                    ),
                )
            )
        ),
        "forbidden_fallbacks_used": [],
        "formal_promotion_performed": False,
        "candidate_manifest": str(candidate_manifest),
        "candidate_manifest_sha256": candidate_manifest_sha256,
    }


def _render_audit_markdown(
    candidates: Sequence[Any],
    source_assessments: Sequence[Mapping[str, Any]],
) -> str:
    ready_count = sum(
        1 for item in candidates if item.status is CapacityCandidateStatus.REVIEW_READY
    )
    missing_line = (
        "- Missing: none at candidate-contract level; explicit reviewer approval is still required.\n"
        if candidates and candidates[0].status is CapacityCandidateStatus.REVIEW_READY
        else "- Missing: reviewed condition-matched abundance/direct capacity and conversion metadata.\n"
    )
    assessment_lines = ""
    for assessment in source_assessments:
        source_name = (
            assessment.get("project_accession")
            or assessment.get("assessment_id")
            or assessment.get("source_id")
        )
        assessment_lines += (
            f"\n## Source assessment: {source_name}\n\n"
            f"- Status: {assessment.get('assessment')}.\n"
            f"- Version/license: {assessment.get('source_version')}; "
            f"{assessment.get('license_id')} ({assessment.get('license_url')}).\n"
            f"- Source URL: {assessment.get('source_url')}.\n"
            f"- Raw values: {assessment.get('raw_values')}.\n"
            f"- Raw units: {assessment.get('raw_units') or assessment.get('raw_unit')}.\n"
            f"- Condition: {assessment.get('condition')}.\n"
            f"- Mapping evidence: {assessment.get('mapping_evidence')}.\n"
            f"- Conflicts: {assessment.get('conflicts')}.\n"
            f"- Missing: {assessment.get('missing_information')}.\n"
            "- Capacity closure/promotion: unavailable / false.\n"
        )
        if assessment.get("raw_metric"):
            assessment_lines += (
                "- Unit chain: iBAQ intensity -> retained raw quantitative evidence; "
                "model_flux unavailable.\n"
                "- Boundary: relative quantitative evidence only.\n"
            )
    return (
        "# G6PDH2 external capacity candidate audit\n\n"
        f"- Candidate count: {len(candidates)}.\n"
        f"- Promotion-ready count: {ready_count}.\n"
        "- UniProt: C4R099 confirms PAS_chr2-1_0308 identity only.\n"
        "- hLF and OPN retain separate current-model fingerprints and bindings.\n"
        + assessment_lines
        + missing_line
        + "- Promotion: not performed.\n"
        + "- Forbidden fallbacks: 1000 upper bound, optimal flux, fixed 1.0, fixture were not used.\n"
    )


def _render_ecpichia_provenance_gap(payload: Mapping[str, Any]) -> str:
    trace = payload["coefficient_trace"]
    bindings = payload["model_bindings"]
    return (
        "# A0c ecPichia G6PDH2 provenance gap\n\n"
        f"- Completion outcome: `{payload['completion_outcome']}`.\n"
        f"- GECKO coefficient reproduced: `{trace['matches']}`; "
        f"expected `{trace['expected_coefficient_mg_h_per_mmol']}`.\n"
        f"- Formal current-model bindings: `{len(bindings)}` at `glucose_mu_0.1`.\n"
        "- Nominal capacity: unavailable.\n"
        "- Promotion preview: unavailable.\n\n"
        "## Frozen source artifacts\n\n"
        + "".join(
            f"- `{item['source_id']}`: `{item['raw_sha256']}`; "
            f"license `{item['license_id']}`; retrieved `{item['retrieved_at']}`.\n"
            for item in payload["source_artifacts"]
        )
        + "\n"
        "## Conflicts\n\n"
        + "".join(f"- `{item}`\n" for item in payload["source_conflicts"])
        + "\n## Missing information\n\n"
        + "".join(f"- `{item}`\n" for item in payload["missing_information"])
        + "\n## Boundary\n\n"
        "The conditional catalytic-flux calculation is retained only as an "
        "audit trace. It is not a reviewed abundance, current-model formation "
        "capacity, or promotion candidate.\n"
    )


def _default_host_condition(request: ExternalCapacityAuditRequest) -> HostCondition:
    return HostCondition(
        species="Komagataella phaffii",
        strain="GS115",
        medium="defined minimal",
        carbon_source=request.carbon_source_id,
        culture_mode="chemostat",
        growth_rate_per_h=request.growth_rate,
        biomass_basis="gDW",
    )


def _context_id(carbon_source_id: str, growth_rate: float) -> str:
    return f"{str(carbon_source_id).strip().lower()}_mu_{float(growth_rate):g}"


def _resolve_repo_path(repo_root: str | Path, path: str | Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else Path(repo_root).resolve() / candidate


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


def _json_default(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    return str(value)


__all__ = [
    "DEFAULT_CAPACITY_ASSET_PATH",
    "DEFAULT_TARGET_IDS",
    "ExternalCapacityAuditOutputs",
    "ExternalCapacityAuditRequest",
    "ExternalCapacitySourceType",
    "G6PDH2_GENE_ID",
    "capacity_asset_version",
    "load_capacity_asset_snapshot",
    "load_external_candidate_review",
    "prepare_external_candidate_runtime",
    "preview_external_candidate_promotion",
    "promote_external_candidate_selection",
    "run_external_capacity_candidate_audit",
]
