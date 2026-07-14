from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Sequence

from pcsec_pichia.external_refs.capacity_sources import (
    ExternalCapacitySourceType,
    RetrievalMode,
    cache_uniprot_identity_source,
)
from pcsec_pichia.loading import load_pcsec_pichia_inputs
from pcsec_pichia.oe_capacity.external_candidate_evaluation import (
    build_capacity_candidate,
    build_capacity_model_binding,
)
from pcsec_pichia.oe_capacity.external_candidate_io import (
    CANDIDATE_MANIFEST_FILENAME,
    import_capacity_measurements,
    load_external_capacity_candidate_bundle,
    load_external_capacity_candidate_snapshot,
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

    def summary(self) -> dict[str, object]:
        return {
            "audit": str(self.audit_json_path),
            "candidate_count": self.candidate_count,
            "promotion_ready_count": self.promotion_ready_count,
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
    identity_source = cache_uniprot_identity_source(
        G6PDH2_GENE_ID,
        request.identity_cache_dir or (output_dir / "identity_source"),
        retrieval_mode=(
            RetrievalMode.OFFLINE_REPLAY
            if request.offline_replay
            else RetrievalMode.ONLINE
        ),
    )
    asset_sha = capacity_asset_version(request.repo_root)
    bindings = []
    fingerprints = []
    catalogs = {}
    for target_id in request.target_ids:
        runtime = prepare_external_candidate_runtime(
            request.repo_root,
            target_id=target_id,
            growth_rate=request.growth_rate,
            carbon_source_id=request.carbon_source_id,
            relative_uncertainty=request.relative_uncertainty,
            expected_asset_sha256=asset_sha,
        )
        binding = build_capacity_model_binding(
            runtime.gene_capacity_catalog,
            target_id=target_id,
            context_id=_context_id(request.carbon_source_id, request.growth_rate),
            gene_id=G6PDH2_GENE_ID,
            external_gene_id=G6PDH2_GENE_ID,
            external_protein_id="C4R099",
            mapping_evidence=(
                "uniprot:C4R099",
                "ncbi_gene:8198996",
                "kegg:ppa:PAS_chr2-1_0308",
                "refseq:XP_002491203.1",
                "ec:1.1.1.49",
            ),
        )
        bindings.append(binding)
        fingerprints.append(binding.model_fingerprint)
        catalogs[target_id] = runtime.gene_capacity_catalog
    capacity_sources = []
    measurements = ()
    candidates = ()
    if request.measurement_file is not None:
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
        _render_audit_markdown(candidates), encoding="utf-8"
    )
    return ExternalCapacityAuditOutputs(
        audit_json_path=audit_json_path,
        audit_markdown_path=audit_markdown_path,
        candidate_manifest_path=candidate_outputs.manifest_path,
        candidate_count=len(candidates),
        promotion_ready_count=promotion_ready_count,
    )


def _build_audit_payload(
    *,
    request: ExternalCapacityAuditRequest,
    identity_source: Any,
    capacity_sources: Sequence[Any],
    measurements: Sequence[Any],
    candidates: Sequence[Any],
    bindings: Sequence[Any],
    candidate_manifest: Path,
    candidate_manifest_sha256: str,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "round": "round_6a_external_capacity_candidates",
        "gene_id": G6PDH2_GENE_ID,
        "enzyme_id": G6PDH2_ENZYME_ID,
        "formation_handle": G6PDH2_FORMATION_HANDLE,
        "targets": list(request.target_ids),
        "context_id": _context_id(request.carbon_source_id, request.growth_rate),
        "sources_checked": [
            {
                **asdict(identity_source),
                "capacity_value_available": False,
                "role": "identity_crosswalk_only",
            },
            {
                "source": "same-host quantitative Pichia proteomics",
                "status": "manual_import_required",
                "capacity_value_available": False,
            },
            {
                "source": "iPichia/ecPichia calibrated enzyme capacity",
                "status": "license_and_artifact_review_required",
                "capacity_value_available": False,
            },
            {
                "source": "Pichia literature plus BRENDA/SABIO-RK kinetics",
                "status": "abundance_and_condition_pair_missing",
                "capacity_value_available": False,
            },
            *(
                [
                    {
                        **asdict(capacity_sources[0]),
                        "capacity_value_available": bool(measurements),
                        "role": "capacity_measurement_source",
                    }
                ]
                if capacity_sources
                else []
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
            ]
            if not candidates
            else [
                reason
                for item in candidates
                for reason in (*item.rejection_reasons, *item.conflicts)
            ]
        ),
        "missing_information": (
            [
                "reviewed quantitative abundance or direct baseline capacity",
                "source version/hash/license for the capacity-valued artifact",
                "condition-matched biomass basis and conversion metadata",
            ]
            if not candidates
            else list(
                dict.fromkeys(
                    value
                    for item in candidates
                    for value in item.missing_information
                )
            )
        ),
        "forbidden_fallbacks_used": [],
        "formal_promotion_performed": False,
        "candidate_manifest": str(candidate_manifest),
        "candidate_manifest_sha256": candidate_manifest_sha256,
    }


def _render_audit_markdown(candidates: Sequence[Any]) -> str:
    ready_count = sum(
        1 for item in candidates if item.status is CapacityCandidateStatus.REVIEW_READY
    )
    missing_line = (
        "- Missing: none at candidate-contract level; explicit reviewer approval is still required.\n"
        if candidates and candidates[0].status is CapacityCandidateStatus.REVIEW_READY
        else "- Missing: reviewed condition-matched abundance/direct capacity and conversion metadata.\n"
    )
    return (
        "# G6PDH2 external capacity candidate audit\n\n"
        f"- Candidate count: {len(candidates)}.\n"
        f"- Promotion-ready count: {ready_count}.\n"
        "- UniProt: C4R099 confirms PAS_chr2-1_0308 identity only.\n"
        "- hLF and OPN retain separate current-model fingerprints and bindings.\n"
        + missing_line
        + "- Promotion: not performed.\n"
        + "- Forbidden fallbacks: 1000 upper bound, optimal flux, fixed 1.0, fixture were not used.\n"
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
