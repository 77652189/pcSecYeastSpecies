from __future__ import annotations

import hashlib
import json
import threading
import time
from dataclasses import replace
from pathlib import Path

import pytest

from pcsec_pichia.external_refs.clients import ExternalHttpResponse
from pcsec_pichia.oe_capacity import (
    CapacityApplicabilityScope,
    CapacityCandidateStatus,
    CapacityConfidence,
    ConfidenceLevel,
    EvidenceSourceType,
    GeneCapacityCatalog,
    GeneEnzymeReactionMapping,
    GPRRole,
    CapacityParameterKind,
    ExternalCapacityCandidateBundle,
    ExternalCapacitySource,
    ExternalCapacitySourceType,
    HostCondition,
    OECapacityValidationError,
    OEExecutionStatus,
    PromotionDecision,
    RawCapacityMeasurement,
    RetrievalMode,
    build_capacity_candidate,
    build_capacity_model_binding,
    build_capacity_promotion_manifest,
    cache_uniprot_identity_source,
    import_capacity_measurements,
    load_external_capacity_candidate_bundle,
    promote_capacity_candidates,
    write_external_capacity_candidate_cache,
)
from pcsec_pichia.oe_capacity.acceptance import _validate_capacity_asset
from pcsec_pichia.oe_capacity import external_candidate_audit
from pcsec_pichia.oe_capacity import external_candidate_promotion
from pcsec_pichia.oe_capacity.external_candidate_io import (
    load_external_capacity_candidate_snapshot,
)
from pcsec_pichia.oe_capacity.parameters import load_capacity_anchor_catalog


HLF_FINGERPRINT = "41056e7a1098c43535f380cf700ae007f7e1027f14c44241397f4db451d0c3d5"
OPN_FINGERPRINT = "03df424f9b57535bb6fa8ef8e5cca392ae728f0b3db5fd34cc29cdcb0b1aea2f"


def _condition() -> HostCondition:
    return HostCondition(
        species="Komagataella phaffii",
        strain="GS115",
        medium="defined minimal",
        carbon_source="glucose",
        culture_mode="chemostat",
        growth_rate_per_h=0.1,
        biomass_basis="gDW",
    )


def _source(source_id: str = "proteomics:study-1") -> ExternalCapacitySource:
    artifact_path = Path(__file__).resolve()
    artifact_sha256 = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
    return ExternalCapacitySource(
        source_id=source_id,
        source_type=ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS,
        source_version="study-1-table-s2",
        source_url="https://example.test/study-1",
        retrieved_at="2026-07-14T00:00:00Z",
        query="PAS_chr2-1_0308 glucose mu=0.1",
        raw_sha256=artifact_sha256,
        license_id="CC-BY-4.0",
        license_url="https://creativecommons.org/licenses/by/4.0/",
        retrieval_mode=RetrievalMode.MANUAL_IMPORT,
        cache_path=str(artifact_path),
        terms_reviewed=True,
    )


def _measurement(
    measurement_id: str = "capacity-1",
    *,
    source_id: str = "proteomics:study-1",
    parameter_kind: CapacityParameterKind = CapacityParameterKind.BASELINE_CAPACITY,
    value: float = 0.25,
    unit: str = "model_flux",
) -> RawCapacityMeasurement:
    return RawCapacityMeasurement(
        measurement_id=measurement_id,
        source_id=source_id,
        parameter_kind=parameter_kind,
        nominal_value=value,
        lower_bound=value * 0.8,
        upper_bound=value * 1.2,
        unit=unit,
        condition=_condition(),
        external_gene_id="PAS_chr2-1_0308",
        biomass_basis="gDW",
    )


def _catalog(target_id: str) -> GeneCapacityCatalog:
    fingerprint = HLF_FINGERPRINT if target_id == "hLF" else OPN_FINGERPRINT
    mapping_id = "oe-map-hlf" if target_id == "hLF" else "oe-map-opn"
    return GeneCapacityCatalog(
        model_fingerprint=fingerprint,
        mappings=(
            GeneEnzymeReactionMapping(
                mapping_id=mapping_id,
                model_fingerprint=fingerprint,
                gene_id="PAS_chr2-1_0308",
                enzyme_id="G6PDH2_no_1_fwd_complex",
                reaction_id="G6PDH2_no_1_fwd",
                gpr_rule="PAS_chr2-1_0308",
                gpr_role=GPRRole.SINGLE_GENE,
                enzyme_variable_id="G6PDH2_no_1_fwd_complex_formation",
                formation_or_dilution_reaction_id="G6PDH2_no_1_fwd_complex_formation",
                mapping_source=EvidenceSourceType.CURRENT_MODEL,
                mapping_confidence=ConfidenceLevel.HIGH,
                execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
            ),
        ),
    )


def _catalogs():
    return {
        "hLF": _catalog("hLF"),
        "OPN_ALPHA_FULL_PROJECT": _catalog("OPN_ALPHA_FULL_PROJECT"),
    }


def _bindings(*, target_specific: str = ""):
    targets = (target_specific,) if target_specific else ("hLF", "OPN_ALPHA_FULL_PROJECT")
    catalogs = _catalogs()
    return tuple(
        build_capacity_model_binding(
            catalogs[target_id],
            target_id=target_id,
            context_id="glucose_mu_0.1",
            gene_id="PAS_chr2-1_0308",
            external_gene_id="PAS_chr2-1_0308",
            external_protein_id="C4R099",
            mapping_evidence=("uniprot:C4R099", "ncbi:8198996", "kegg:ppa:PAS_chr2-1_0308"),
        )
        for target_id in targets
    )


def _candidate_kwargs(scope: CapacityApplicabilityScope, source, measurement):
    target_id = "hLF" if scope is CapacityApplicabilityScope.TARGET_SPECIFIC else ""
    return {
        "applicability_scope": scope,
        "target_id": target_id,
        "model_bindings": _bindings(target_specific=target_id),
        "catalogs": _catalogs(),
        "sources": (source,),
        "condition": _condition(),
        "direct_capacity": measurement,
    }


def _invalid_binding():
    binding = _bindings(target_specific="hLF")[0]
    return type(binding)(
        target_id=binding.target_id,
        context_id=binding.context_id,
        mapping_id=binding.mapping_id,
        model_fingerprint=binding.model_fingerprint,
        gene_id="PAS_chr2-1_0308",
        enzyme_id="G6PDH2_no_1_fwd_complex",
        reaction_id="not-current-reaction",
        formation_or_dilution_reaction_id="G6PDH2_no_1_fwd_complex_formation",
        mapping_evidence=binding.mapping_evidence,
        external_gene_id="PAS_chr2-1_0308",
        external_protein_id="C4R099",
    )


def _bundle(*, scope: CapacityApplicabilityScope = CapacityApplicabilityScope.HOST_CONDITION):
    source = _source()
    measurement = _measurement()
    candidate = build_capacity_candidate(
        candidate_id="g6pdh2-host-condition",
        **_candidate_kwargs(scope, source, measurement),
        confidence=CapacityConfidence.HIGH,
    )
    return ExternalCapacityCandidateBundle(
        model_fingerprints=(HLF_FINGERPRINT, OPN_FINGERPRINT),
        sources=(source,),
        measurements=(measurement,),
        candidates=(candidate,),
        generated_at="2026-07-14T00:00:00Z",
    )


def test_host_condition_candidate_reuses_one_provenance_across_targets() -> None:
    bundle = _bundle()
    candidate = bundle.candidates[0]

    assert candidate.target_id == ""
    assert candidate.status is CapacityCandidateStatus.REVIEW_READY
    assert candidate.source_ids == ("proteomics:study-1",)

    manifest = build_capacity_promotion_manifest(
        bundle,
        (candidate.candidate_id,),
        "capacity.json",
        decision=PromotionDecision.APPROVED,
        reviewer="reviewer@example.test",
        reviewed_at="2026-07-14T01:00:00Z",
    )
    assert manifest.decision is PromotionDecision.APPROVED


def test_target_specific_requires_target_and_other_scopes_reject_target() -> None:
    source = _source()
    measurement = _measurement()
    with pytest.raises(OECapacityValidationError, match="requires target_id"):
        build_capacity_candidate(
            candidate_id="bad-target",
            applicability_scope=CapacityApplicabilityScope.TARGET_SPECIFIC,
            model_bindings=_bindings(target_specific="hLF"),
            catalogs=_catalogs(),
            sources=(source,),
            condition=_condition(),
            direct_capacity=measurement,
        )

    with pytest.raises(OECapacityValidationError, match="only target_specific"):
        build_capacity_candidate(
            candidate_id="bad-host",
            applicability_scope=CapacityApplicabilityScope.HOST_CONDITION,
            target_id="hLF",
            model_bindings=_bindings(),
            catalogs=_catalogs(),
            sources=(source,),
            condition=_condition(),
            direct_capacity=measurement,
        )


def test_homolog_transfer_stays_low_confidence_and_cannot_be_promoted() -> None:
    source = _source()
    measurement = _measurement()
    candidate = build_capacity_candidate(
        candidate_id="homolog-only",
        applicability_scope=CapacityApplicabilityScope.HOMOLOG_TRANSFERRED,
        model_bindings=_bindings(),
        catalogs=_catalogs(),
        sources=(source,),
        condition=_condition(),
        direct_capacity=measurement,
        confidence=CapacityConfidence.HIGH,
    )
    assert candidate.confidence is CapacityConfidence.LOW
    assert candidate.status is CapacityCandidateStatus.REVIEW_REQUIRED
    assert "independent_pichia_capacity_evidence" in candidate.missing_information


def test_missing_conversion_never_invents_relative_capacity() -> None:
    source = _source()
    measurement = _measurement(unit="copies/cell")
    candidate = build_capacity_candidate(
        candidate_id="unsupported-unit",
        applicability_scope=CapacityApplicabilityScope.HOST_CONDITION,
        model_bindings=_bindings(),
        catalogs=_catalogs(),
        sources=(source,),
        condition=_condition(),
        direct_capacity=measurement,
    )
    assert candidate.nominal_capacity is None
    assert candidate.lower_capacity is None
    assert candidate.upper_capacity is None
    assert candidate.status is CapacityCandidateStatus.REVIEW_REQUIRED
    assert candidate.missing_information == ("supported_capacity_unit_conversion",)


def test_source_review_status_and_model_binding_are_derived_from_artifacts() -> None:
    measurement = _measurement()
    unreviewed = ExternalCapacitySource(
        **{**_source().__dict__, "terms_reviewed": False}
    )
    candidate = build_capacity_candidate(
        candidate_id="unreviewed-source",
        **_candidate_kwargs(
            CapacityApplicabilityScope.HOST_CONDITION,
            unreviewed,
            measurement,
        ),
    )
    assert candidate.status is CapacityCandidateStatus.REVIEW_REQUIRED
    assert "source_license_version_or_terms_requires_review" in candidate.warnings

    identity_source = ExternalCapacitySource(
        **{
            **_source("identity:C4R099").__dict__,
            "source_type": ExternalCapacitySourceType.IDENTITY_REFERENCE,
        }
    )
    identity_measurement = _measurement(source_id=identity_source.source_id)
    identity_candidate = build_capacity_candidate(
        candidate_id="identity-is-not-capacity",
        **_candidate_kwargs(
            CapacityApplicabilityScope.HOST_CONDITION,
            identity_source,
            identity_measurement,
        ),
    )
    assert identity_candidate.status is CapacityCandidateStatus.REVIEW_REQUIRED
    assert "capacity_valued_source" in identity_candidate.missing_information

    with pytest.raises(OECapacityValidationError, match="identity does not match"):
        build_capacity_candidate(
            candidate_id="forged-binding",
            applicability_scope=CapacityApplicabilityScope.TARGET_SPECIFIC,
            target_id="hLF",
            model_bindings=(_invalid_binding(),),
            catalogs=_catalogs(),
            sources=(_source(),),
            condition=_condition(),
            direct_capacity=measurement,
        )


def test_abundance_times_kcat_conversion_is_traceable() -> None:
    source = _source()
    abundance = _measurement(
        "abundance-1",
        parameter_kind=CapacityParameterKind.ABUNDANCE,
        value=500.0,
        unit="nmol_enzyme/gDW",
    )
    kcat = _measurement(
        "kcat-1",
        parameter_kind=CapacityParameterKind.KCAT,
        value=10.0,
        unit="1/s",
    )
    candidate = build_capacity_candidate(
        candidate_id="derived-capacity",
        applicability_scope=CapacityApplicabilityScope.HOST_CONDITION,
        model_bindings=_bindings(),
        catalogs=_catalogs(),
        sources=(source,),
        condition=_condition(),
        abundance=abundance,
        kcat=kcat,
    )
    assert candidate.nominal_capacity == pytest.approx(18.0)
    assert [step.step_id for step in candidate.conversion_steps] == [
        "abundance-to-mmol-per-gdw",
        "kcat-to-per-hour",
        "abundance-times-kcat",
    ]


def test_candidate_cache_roundtrip_and_hash_tamper_detection(tmp_path: Path) -> None:
    outputs = write_external_capacity_candidate_cache(_bundle(), tmp_path)
    loaded = load_external_capacity_candidate_bundle(outputs.manifest_path)
    assert loaded == _bundle()

    outputs.records_path.write_text(outputs.records_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(OECapacityValidationError, match="sha256 mismatch"):
        load_external_capacity_candidate_bundle(outputs.manifest_path)


@pytest.mark.parametrize("field_name", ["records_file", "measurements_file"])
@pytest.mark.parametrize("path_kind", ["parent", "absolute", "subdirectory"])
def test_candidate_bundle_manifest_rejects_non_basename_artifact_paths(
    tmp_path: Path,
    field_name: str,
    path_kind: str,
) -> None:
    outputs = write_external_capacity_candidate_cache(_bundle(), tmp_path / "bundle")
    manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    source_path = outputs.records_path if field_name == "records_file" else outputs.measurements_path
    if path_kind == "parent":
        artifact_path = tmp_path / f"outside-{source_path.name}"
        manifest[field_name] = f"../{artifact_path.name}"
    elif path_kind == "absolute":
        artifact_path = tmp_path / f"absolute-{source_path.name}"
        manifest[field_name] = str(artifact_path.resolve())
    else:
        artifact_path = outputs.manifest_path.parent / "nested" / source_path.name
        manifest[field_name] = f"nested/{source_path.name}"
    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    artifact_path.write_bytes(source_path.read_bytes())
    outputs.manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(OECapacityValidationError, match="plain basename"):
        load_external_capacity_candidate_bundle(outputs.manifest_path)


def test_candidate_bundle_load_reads_each_artifact_once_and_can_pin_manifest_hash(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    outputs = write_external_capacity_candidate_cache(_bundle(), tmp_path)
    original_open = Path.open
    artifact_paths = {outputs.records_path.resolve(), outputs.measurements_path.resolve()}
    open_counts = {path: 0 for path in artifact_paths}

    def counted_open(path: Path, *args: object, **kwargs: object):
        resolved = path.resolve()
        if resolved in open_counts:
            open_counts[resolved] += 1
            if open_counts[resolved] > 1:
                raise AssertionError(f"artifact was reopened after snapshot: {resolved}")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", counted_open)
    loaded = load_external_capacity_candidate_bundle(
        outputs.manifest_path,
        expected_manifest_sha256=outputs.bundle_sha256,
    )

    assert loaded == _bundle()
    assert set(open_counts.values()) == {1}


def test_candidate_bundle_load_rejects_wrong_expected_manifest_hash(tmp_path: Path) -> None:
    outputs = write_external_capacity_candidate_cache(_bundle(), tmp_path)

    with pytest.raises(OECapacityValidationError, match="manifest sha256 mismatch"):
        load_external_capacity_candidate_bundle(
            outputs.manifest_path,
            expected_manifest_sha256="0" * 64,
        )


def test_manual_import_copies_and_hashes_source_without_promotion(tmp_path: Path) -> None:
    source_file = tmp_path / "measurements.csv"
    source_file.write_text(
        "measurement_id,parameter_kind,nominal_value,lower_bound,upper_bound,unit,external_gene_id,species,strain,medium,carbon_source,culture_mode,growth_rate_per_h,biomass_basis\n"
        "capacity-1,baseline_capacity,0.25,0.2,0.3,model_flux,PAS_chr2-1_0308,Komagataella phaffii,GS115,defined minimal,glucose,chemostat,0.1,gDW\n",
        encoding="utf-8",
    )
    expected = hashlib.sha256(source_file.read_bytes()).hexdigest()
    source, measurements = import_capacity_measurements(
        source_file,
        source_id="manual:study-1",
        source_type=ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS,
        source_version="table-s2",
        source_url="https://example.test/study-1",
        license_id="CC-BY-4.0",
        query="G6PDH2",
        output_dir=tmp_path / "cache",
        expected_sha256=expected,
        terms_reviewed=True,
    )
    assert source.raw_sha256 == expected
    assert source.retrieval_mode is RetrievalMode.MANUAL_IMPORT
    assert measurements[0].condition.growth_rate_per_h == pytest.approx(0.1)
    assert Path(source.cache_path).read_bytes() == source_file.read_bytes()
    ready = build_capacity_candidate(
        candidate_id="manual-ready",
        applicability_scope=CapacityApplicabilityScope.HOST_CONDITION,
        model_bindings=_bindings(),
        catalogs=_catalogs(),
        sources=(source,),
        condition=measurements[0].condition,
        direct_capacity=measurements[0],
    )
    assert ready.status is CapacityCandidateStatus.REVIEW_READY

    unpinned_source, unpinned_measurements = import_capacity_measurements(
        source_file,
        source_id="manual:unpinned",
        source_type=ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS,
        source_version="table-s2",
        source_url="https://example.test/study-1",
        license_id="CC-BY-4.0",
        query="G6PDH2",
        output_dir=tmp_path / "cache-unpinned",
        terms_reviewed=True,
    )
    unpinned = build_capacity_candidate(
        candidate_id="manual-unpinned",
        applicability_scope=CapacityApplicabilityScope.HOST_CONDITION,
        model_bindings=_bindings(),
        catalogs=_catalogs(),
        sources=(unpinned_source,),
        condition=unpinned_measurements[0].condition,
        direct_capacity=unpinned_measurements[0],
    )
    assert unpinned.status is CapacityCandidateStatus.REVIEW_REQUIRED
    assert "expected_sha256_not_predeclared" in unpinned_source.warnings


def test_uniprot_online_cache_and_offline_replay_preserve_identity(tmp_path: Path) -> None:
    raw = json.dumps(
        {
            "results": [
                {
                    "primaryAccession": "C4R099",
                    "uniProtkbId": "A0A1E4RTV1_PICPA",
                    "entryType": "UniProtKB unreviewed (TrEMBL)",
                    "genes": [{"geneName": {"value": "G6PDH2"}, "orderedLocusNames": [{"value": "PAS_chr2-1_0308"}]}],
                    "organism": {"taxonId": 644223, "scientificName": "Komagataella phaffii"},
                    "proteinDescription": {"recommendedName": {"fullName": {"value": "Glucose-6-phosphate 1-dehydrogenase"}}},
                }
            ]
        }
    )

    def fake_get(url, config):
        return ExternalHttpResponse(
            status_code=200,
            text=raw,
            url=url,
            headers={"X-UniProt-Release": "2026_02", "Authorization": "must-not-persist"},
        )

    online = cache_uniprot_identity_source("PAS_chr2-1_0308", tmp_path, http_get=fake_get)
    offline = cache_uniprot_identity_source(
        "PAS_chr2-1_0308",
        tmp_path,
        http_get=lambda *_: pytest.fail("offline replay must not call the network"),
        retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
    )
    assert online.source_id == offline.source_id == "uniprot:C4R099"
    assert online.source_version == "2026_02"
    assert online.raw_sha256 == offline.raw_sha256
    assert "Authorization" not in Path(online.cache_path).read_text(encoding="utf-8")
    assert "identity_only_not_capacity_evidence" in online.warnings

    raw_path = Path(online.cache_path)
    raw_path.write_text(raw_path.read_text(encoding="utf-8") + " ", encoding="utf-8")
    with pytest.raises(OECapacityValidationError, match="sha256 mismatch"):
        cache_uniprot_identity_source(
            "PAS_chr2-1_0308",
            tmp_path,
            http_get=lambda *_: pytest.fail("tampered replay must not call the network"),
            retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
        )


def test_promotion_defaults_to_pending_and_approved_write_is_atomic(tmp_path: Path) -> None:
    bundle = _bundle()
    asset_path = tmp_path / "oe_capacity_baseline_capacity.json"
    pending = build_capacity_promotion_manifest(
        bundle,
        (bundle.candidates[0].candidate_id,),
        asset_path,
    )
    assert pending.decision is PromotionDecision.PENDING
    assert not asset_path.exists()
    with pytest.raises(OECapacityValidationError, match="approved"):
        promote_capacity_candidates(bundle, pending, catalogs=_catalogs())

    previous_anchor = {
        "anchor_id": "keep-existing",
        "target_id": "hLF",
        "context_id": "glucose_mu_0.1",
        "gene_id": "OTHER",
        "enzyme_id": "OTHER_complex",
        "formation_or_dilution_reaction_id": "OTHER_formation",
        "model_fingerprint": HLF_FINGERPRINT,
        "baseline_capacity": 0.5,
        "unit": "model_flux",
        "source_ref": "existing-reviewed-source",
        "source_version": "v1",
        "reviewed_by": "previous-reviewer",
        "reviewed_at": "2026-07-13T00:00:00Z",
    }
    asset_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "asset_version": "existing-v1",
                "model_fingerprint": "multi-target:existing",
                "anchors": [previous_anchor],
            }
        ),
        encoding="utf-8",
    )

    approved = build_capacity_promotion_manifest(
        bundle,
        (bundle.candidates[0].candidate_id,),
        asset_path,
        decision=PromotionDecision.APPROVED,
        reviewer="reviewer@example.test",
        reviewed_at="2026-07-14T01:00:00Z",
    )
    promoted = promote_capacity_candidates(bundle, approved, catalogs=_catalogs())
    payload = json.loads(asset_path.read_text(encoding="utf-8"))
    assert {item["target_id"] for item in payload["anchors"]} == {"hLF", "OPN_ALPHA_FULL_PROJECT"}
    promoted_rows = [item for item in payload["anchors"] if item["gene_id"] == "PAS_chr2-1_0308"]
    assert {item["model_fingerprint"] for item in promoted_rows} == {HLF_FINGERPRINT, OPN_FINGERPRINT}
    assert len({item["source_ref"] for item in promoted_rows}) == 1
    assert all(item["applicability_scope"] == "host_condition" for item in promoted_rows)
    assert all(item["host_condition"]["strain"] == "GS115" for item in promoted_rows)
    assert previous_anchor in payload["anchors"]
    assert payload["asset_version"].startswith("round6a-")
    assert load_capacity_anchor_catalog(asset_path).anchors
    assert _validate_capacity_asset(asset_path)["valid"] is True
    assert promoted.promoted_asset_sha256 == hashlib.sha256(asset_path.read_bytes()).hexdigest()


def test_promotion_manifest_pins_the_reviewed_asset_hash(tmp_path: Path) -> None:
    bundle = _bundle()
    asset_path = tmp_path / "oe_capacity_baseline_capacity.json"
    asset_path.write_text(
        json.dumps({"schema_version": 1, "anchors": []}),
        encoding="utf-8",
    )
    reviewed_sha256 = hashlib.sha256(asset_path.read_bytes()).hexdigest()

    with pytest.raises(OECapacityValidationError, match="changed after promotion review"):
        build_capacity_promotion_manifest(
            bundle,
            (bundle.candidates[0].candidate_id,),
            asset_path,
            expected_asset_sha256="0" * 64,
        )

    approved = build_capacity_promotion_manifest(
        bundle,
        (bundle.candidates[0].candidate_id,),
        asset_path,
        decision=PromotionDecision.APPROVED,
        reviewer="reviewer@example.test",
        reviewed_at="2026-07-14T01:00:00Z",
        expected_asset_sha256=reviewed_sha256,
    )
    assert approved.expected_asset_sha256 == reviewed_sha256

    asset_path.write_text(
        json.dumps({"schema_version": 1, "anchors": [], "changed": True}),
        encoding="utf-8",
    )
    with pytest.raises(OECapacityValidationError, match="changed after promotion review"):
        promote_capacity_candidates(bundle, approved, catalogs=_catalogs())


def test_candidate_review_reports_the_hash_from_the_displayed_snapshot(
    tmp_path: Path,
    monkeypatch,
) -> None:
    outputs = write_external_capacity_candidate_cache(_bundle(), tmp_path / "cache")
    original_loader = load_external_capacity_candidate_snapshot

    def replace_manifest_after_snapshot(path):
        snapshot = original_loader(path)
        Path(path).write_text("{}\n", encoding="utf-8")
        return snapshot

    monkeypatch.setattr(
        external_candidate_audit,
        "load_external_capacity_candidate_snapshot",
        replace_manifest_after_snapshot,
    )
    review = external_candidate_audit.load_external_candidate_review(
        tmp_path,
        outputs.manifest_path,
        target_id="hLF",
        capacity_asset_path=tmp_path / "missing-capacity-asset.json",
    )

    assert review["available"] is True
    assert review["candidate_manifest_sha256"] == outputs.bundle_sha256
    assert review["candidates"]


def test_concurrent_promotions_cannot_silently_overwrite_the_asset(
    tmp_path: Path,
    monkeypatch,
) -> None:
    bundle = _bundle()
    asset_path = tmp_path / "oe_capacity_baseline_capacity.json"
    asset_path.write_text(
        json.dumps({"schema_version": 1, "anchors": []}),
        encoding="utf-8",
    )
    reviewed_sha256 = hashlib.sha256(asset_path.read_bytes()).hexdigest()
    manifests = tuple(
        build_capacity_promotion_manifest(
            bundle,
            (bundle.candidates[0].candidate_id,),
            asset_path,
            decision=PromotionDecision.APPROVED,
            reviewer=f"reviewer-{index}@example.test",
            reviewed_at=f"2026-07-14T01:00:0{index}Z",
            expected_asset_sha256=reviewed_sha256,
        )
        for index in (1, 2)
    )
    original_write = external_candidate_promotion._atomic_write_json
    first_write_entered = threading.Event()
    release_first_write = threading.Event()
    second_finished = threading.Event()
    delayed_once = False
    guard = threading.Lock()
    successes: list[object] = []
    failures: list[Exception] = []

    def delayed_write(path, payload):
        nonlocal delayed_once
        should_delay = False
        if Path(path) == asset_path:
            with guard:
                if not delayed_once:
                    delayed_once = True
                    should_delay = True
        if should_delay:
            first_write_entered.set()
            assert release_first_write.wait(timeout=5)
        original_write(path, payload)

    def run_promotion(manifest, *, mark_finished: bool = False):
        try:
            successes.append(
                promote_capacity_candidates(bundle, manifest, catalogs=_catalogs())
            )
        except Exception as exc:
            failures.append(exc)
        finally:
            if mark_finished:
                second_finished.set()

    monkeypatch.setattr(external_candidate_promotion, "_atomic_write_json", delayed_write)
    first = threading.Thread(target=run_promotion, args=(manifests[0],))
    second = threading.Thread(
        target=run_promotion,
        args=(manifests[1],),
        kwargs={"mark_finished": True},
    )
    first.start()
    assert first_write_entered.wait(timeout=5)
    second.start()
    time.sleep(0.15)
    assert second_finished.is_set() is False
    release_first_write.set()
    first.join(timeout=5)
    second.join(timeout=5)

    assert first.is_alive() is False
    assert second.is_alive() is False
    assert len(successes) == 1
    assert len(failures) == 1
    assert isinstance(failures[0], OECapacityValidationError)
    assert "changed after promotion review" in str(failures[0])


def test_promotion_revalidates_source_artifact_and_current_model_binding(
    tmp_path: Path,
) -> None:
    bundle = _bundle()
    candidate = bundle.candidates[0]
    bad_source = replace(bundle.sources[0], raw_sha256="0" * 64)
    bad_source_bundle = replace(bundle, sources=(bad_source,))
    bad_source_manifest = build_capacity_promotion_manifest(
        bad_source_bundle,
        (candidate.candidate_id,),
        tmp_path / "bad-source.json",
        decision=PromotionDecision.APPROVED,
        reviewer="reviewer",
        reviewed_at="2026-07-14T01:00:00Z",
    )
    with pytest.raises(OECapacityValidationError, match="source artifact hash"):
        promote_capacity_candidates(
            bad_source_bundle,
            bad_source_manifest,
            catalogs=_catalogs(),
        )

    forged_candidate = replace(
        candidate,
        model_bindings=(_invalid_binding(), candidate.model_bindings[1]),
    )
    forged_bundle = replace(bundle, candidates=(forged_candidate,))
    forged_manifest = build_capacity_promotion_manifest(
        forged_bundle,
        (candidate.candidate_id,),
        tmp_path / "forged-binding.json",
        decision=PromotionDecision.APPROVED,
        reviewer="reviewer",
        reviewed_at="2026-07-14T01:00:00Z",
    )
    with pytest.raises(OECapacityValidationError, match="identity does not match"):
        promote_capacity_candidates(
            forged_bundle,
            forged_manifest,
            catalogs=_catalogs(),
        )
