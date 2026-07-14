from __future__ import annotations

import hashlib
import json
import os
import tempfile
import time
from contextlib import contextmanager
from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from pcsec_pichia.oe_capacity.external_candidate_evaluation import (
    _validate_binding_against_catalog,
)
from pcsec_pichia.oe_capacity.external_candidate_io import (
    _as_object_list,
    _atomic_write_json,
    _json_ready,
    _sha256_file,
    _sha256_json,
    _source_artifact_matches,
)
from pcsec_pichia.oe_capacity.external_candidate_schema import (
    PROMOTION_MANIFEST_FILENAME,
    CapacityApplicabilityScope,
    CapacityCandidateStatus,
    CapacityPromotionManifest,
    ExternalCapacityCandidateBundle,
    PromotionDecision,
)
from pcsec_pichia.oe_capacity.schema import (
    GeneCapacityCatalog,
    OECapacityValidationError,
)


def build_capacity_promotion_manifest(
    bundle: ExternalCapacityCandidateBundle,
    candidate_ids: Sequence[str],
    asset_path: str | Path,
    *,
    decision: PromotionDecision = PromotionDecision.PENDING,
    reviewer: str = "",
    reviewed_at: str = "",
    reason: str = "",
    expected_asset_sha256: str | None = None,
) -> CapacityPromotionManifest:
    bundle.validate()
    selected = tuple(dict.fromkeys(str(item).strip() for item in candidate_ids if str(item).strip()))
    known = {item.candidate_id: item for item in bundle.candidates}
    if not selected or any(item not in known for item in selected):
        raise OECapacityValidationError("promotion references unknown or empty candidate_ids.")
    if decision is PromotionDecision.APPROVED:
        for candidate_id in selected:
            candidate = known[candidate_id]
            if candidate.status is not CapacityCandidateStatus.REVIEW_READY:
                raise OECapacityValidationError("only review_ready candidates can be approved for promotion.")
            if candidate.applicability_scope is CapacityApplicabilityScope.HOMOLOG_TRANSFERRED:
                raise OECapacityValidationError("homolog_transferred candidate cannot be promoted alone.")
    asset = Path(asset_path)
    current_asset_sha256 = _sha256_file(asset) if asset.is_file() else "missing"
    if (
        expected_asset_sha256 is not None
        and current_asset_sha256 != expected_asset_sha256
    ):
        raise OECapacityValidationError(
            "formal capacity asset changed after promotion review."
        )
    manifest = CapacityPromotionManifest(
        decision=decision,
        candidate_ids=selected,
        model_fingerprints=bundle.model_fingerprints,
        candidate_bundle_sha256=_sha256_json(_json_ready(asdict(bundle))),
        asset_path=str(asset_path),
        expected_asset_sha256=(expected_asset_sha256 or current_asset_sha256),
        reviewer=reviewer,
        reviewed_at=reviewed_at,
        reason=reason,
    )
    manifest.validate()
    return manifest


def promote_capacity_candidates(
    bundle: ExternalCapacityCandidateBundle,
    manifest: CapacityPromotionManifest,
    *,
    catalogs: Mapping[str, GeneCapacityCatalog],
) -> CapacityPromotionManifest:
    bundle.validate()
    manifest.validate()
    with _exclusive_asset_promotion_lock(Path(manifest.asset_path)):
        return _promote_capacity_candidates_locked(
            bundle,
            manifest,
            catalogs=catalogs,
        )


def _promote_capacity_candidates_locked(
    bundle: ExternalCapacityCandidateBundle,
    manifest: CapacityPromotionManifest,
    *,
    catalogs: Mapping[str, GeneCapacityCatalog],
) -> CapacityPromotionManifest:
    bundle.validate()
    manifest.validate()
    if manifest.decision is not PromotionDecision.APPROVED:
        raise OECapacityValidationError("promotion requires an approved manifest.")
    if manifest.model_fingerprints != bundle.model_fingerprints:
        raise OECapacityValidationError("promotion model_fingerprints mismatch.")
    if manifest.candidate_bundle_sha256 != _sha256_json(_json_ready(asdict(bundle))):
        raise OECapacityValidationError("promotion candidate bundle hash mismatch.")
    asset_path = Path(manifest.asset_path)
    current_bytes = asset_path.read_bytes() if asset_path.is_file() else None
    current_hash = (
        hashlib.sha256(current_bytes).hexdigest()
        if current_bytes is not None
        else "missing"
    )
    if current_hash != manifest.expected_asset_sha256:
        raise OECapacityValidationError("formal capacity asset changed after promotion review.")
    existing_payload: Mapping[str, Any] = {}
    if current_bytes is not None:
        decoded = json.loads(current_bytes.decode("utf-8"))
        if not isinstance(decoded, Mapping):
            raise OECapacityValidationError(
                "formal capacity asset root must be a JSON object."
            )
        existing_payload = decoded
    existing_anchors = list(
        _as_object_list(existing_payload.get("anchors", []), "anchors")
    )
    anchors: list[dict[str, Any]] = [dict(item) for item in existing_anchors]
    identities = {
        _anchor_identity(item): item for item in existing_anchors
    }
    source_lookup = {item.source_id: item for item in bundle.sources}
    selected = {
        item.candidate_id: item
        for item in bundle.candidates
        if item.candidate_id in manifest.candidate_ids
    }
    for candidate_id in manifest.candidate_ids:
        candidate = selected[candidate_id]
        for binding in candidate.model_bindings:
            _validate_binding_against_catalog(binding, catalogs.get(binding.target_id))
        for source_id in candidate.source_ids:
            source = source_lookup[source_id]
            if not _source_artifact_matches(source):
                raise OECapacityValidationError(
                    f"capacity source artifact hash verification failed: {source_id}."
                )
        if candidate.nominal_capacity is None:
            raise OECapacityValidationError("promoted candidate requires canonical capacity.")
        for binding in candidate.model_bindings:
            anchor = {
                "anchor_id": f"{candidate.candidate_id}:{binding.target_id}",
                "target_id": binding.target_id,
                "context_id": binding.context_id,
                "gene_id": binding.gene_id,
                "enzyme_id": binding.enzyme_id,
                "formation_or_dilution_reaction_id": binding.formation_or_dilution_reaction_id,
                "model_fingerprint": binding.model_fingerprint,
                "baseline_capacity": candidate.nominal_capacity,
                "lower_capacity": candidate.lower_capacity,
                "upper_capacity": candidate.upper_capacity,
                "unit": "model_flux",
                "source_ref": ",".join(candidate.source_ids),
                "source_version": manifest.candidate_bundle_sha256,
                "reviewed_by": manifest.reviewer,
                "reviewed_at": manifest.reviewed_at,
                "applicability_scope": candidate.applicability_scope.value,
                "host_condition": _json_ready(asdict(candidate.condition)),
                "mapping_id": binding.mapping_id,
                "mapping_evidence": list(binding.mapping_evidence),
                "source_metadata": [
                    _json_ready(asdict(source_lookup[source_id]))
                    for source_id in candidate.source_ids
                ],
                "conversion_steps": [
                    _json_ready(asdict(step)) for step in candidate.conversion_steps
                ],
            }
            identity = _anchor_identity(anchor)
            if identity in identities and dict(identities[identity]) != anchor:
                raise OECapacityValidationError(
                    f"formal capacity asset already contains a conflicting anchor: {identity}."
                )
            if identity not in identities:
                anchors.append(anchor)
                identities[identity] = anchor
    payload = {
        "schema_version": 1,
        "asset_version": (
            f"round6a-{manifest.reviewed_at}-{manifest.candidate_bundle_sha256[:12]}"
        ),
        "model_fingerprint": "multi-target:" + _sha256_json(list(bundle.model_fingerprints))[:16],
        "source_ref": "Round 6A explicitly approved external capacity promotion",
        "anchors": anchors,
    }
    latest_hash = _sha256_file(asset_path) if asset_path.is_file() else "missing"
    if latest_hash != manifest.expected_asset_sha256:
        raise OECapacityValidationError(
            "formal capacity asset changed during promotion execution."
        )
    _atomic_write_json(asset_path, payload)
    promoted_hash = _sha256_file(asset_path)
    promoted = CapacityPromotionManifest(
        **{**asdict(manifest), "promoted_asset_sha256": promoted_hash}
    )
    promotion_path = asset_path.parent / PROMOTION_MANIFEST_FILENAME
    _atomic_write_json(promotion_path, _json_ready(asdict(promoted)))
    return promoted


@contextmanager
def _exclusive_asset_promotion_lock(
    asset_path: Path,
    *,
    timeout_seconds: float = 30.0,
):
    lock_root = Path(tempfile.gettempdir()) / "pcsec_pichia_oe_capacity_locks"
    lock_root.mkdir(parents=True, exist_ok=True)
    lock_key = hashlib.sha256(
        str(asset_path.resolve()).casefold().encode("utf-8")
    ).hexdigest()
    lock_path = lock_root / f"{lock_key}.lock"
    with lock_path.open("a+b") as handle:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"\0")
            handle.flush()
        deadline = time.monotonic() + timeout_seconds
        while True:
            try:
                _lock_file_handle(handle)
                break
            except OSError as exc:
                if time.monotonic() >= deadline:
                    raise OECapacityValidationError(
                        "timed out waiting for the formal capacity asset promotion lock."
                    ) from exc
                time.sleep(0.05)
        try:
            yield
        finally:
            _unlock_file_handle(handle)


def _lock_file_handle(handle: Any) -> None:
    handle.seek(0)
    if os.name == "nt":
        import msvcrt

        msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        return
    import fcntl

    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)


def _unlock_file_handle(handle: Any) -> None:
    handle.seek(0)
    if os.name == "nt":
        import msvcrt

        msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        return
    import fcntl

    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _anchor_identity(item: Mapping[str, Any]) -> tuple[str, str, str, str, str, str]:
    return (
        str(item.get("model_fingerprint") or ""),
        str(item.get("target_id") or ""),
        str(item.get("context_id") or ""),
        str(item.get("gene_id") or ""),
        str(item.get("enzyme_id") or ""),
        str(item.get("formation_or_dilution_reaction_id") or ""),
    )
