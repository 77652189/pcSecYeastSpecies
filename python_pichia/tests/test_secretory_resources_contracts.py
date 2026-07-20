from __future__ import annotations

import pytest

from pcsec_pichia.secretory_resources import (
    CalibrationMode,
    EvidenceClass,
    ExecutionStatus,
    ResourceApplicability,
    ResourceCategory,
    ResourceSource,
    SecretoryResource,
    SecretoryResourceCatalog,
    SecretoryResourceValidationError,
)


def _source() -> ResourceSource:
    return ResourceSource(
        source_ref="pcsec_pichia.core.target_protein_plan.build_target_protein_plan",
        version="1",
        evidence_class=EvidenceClass.CURRENT_MODEL_HANDLE,
    )


def _applicability(target_id: str = "hLF") -> ResourceApplicability:
    return ResourceApplicability(host="Komagataella phaffii", target_id=target_id, condition="any")


def _executable_resource(**overrides: object) -> SecretoryResource:
    fields = dict(
        resource_id="hLF_er_translocation",
        category=ResourceCategory.ER_TRANSLOCATION,
        canonical_unit="model_flux",
        model_handles=("hLF_Post_translation_PSTA_sec_SEC61SEC63C_complex",),
        source=_source(),
        applicability=_applicability(),
        status=ExecutionStatus.EXECUTABLE,
        calibration_mode=CalibrationMode.RELATIVE_UNCALIBRATED,
    )
    fields.update(overrides)
    return SecretoryResource(**fields)  # type: ignore[arg-type]


def test_valid_executable_resource_passes() -> None:
    _executable_resource().validate()


def test_executable_requires_model_handles() -> None:
    with pytest.raises(SecretoryResourceValidationError, match="model_handle"):
        _executable_resource(model_handles=()).validate()


def test_executable_requires_relative_calibration_mode() -> None:
    with pytest.raises(SecretoryResourceValidationError, match="relative_uncalibrated"):
        _executable_resource(calibration_mode=None).validate()
    with pytest.raises(SecretoryResourceValidationError, match="relative_uncalibrated"):
        _executable_resource(calibration_mode=CalibrationMode.ABSOLUTE_CALIBRATED).validate()


@pytest.mark.parametrize("field_name", ("nominal", "lower", "upper"))
def test_round0_forbids_any_numeric_value_regardless_of_status(field_name: str) -> None:
    with pytest.raises(SecretoryResourceValidationError, match="does not compute"):
        _executable_resource(**{field_name: 1.0}).validate()


def test_evidence_only_must_not_carry_calibration_mode() -> None:
    resource = _executable_resource(
        status=ExecutionStatus.EVIDENCE_ONLY,
        calibration_mode=CalibrationMode.RELATIVE_UNCALIBRATED,
    )
    with pytest.raises(SecretoryResourceValidationError, match="calibration_mode"):
        resource.validate()


def test_evidence_only_allows_handles_without_calibration_mode() -> None:
    resource = _executable_resource(status=ExecutionStatus.EVIDENCE_ONLY, calibration_mode=None)
    resource.validate()


@pytest.mark.parametrize("status", (ExecutionStatus.UNAVAILABLE, ExecutionStatus.NOT_APPLICABLE))
def test_unavailable_and_not_applicable_forbid_model_handles(status: ExecutionStatus) -> None:
    resource = _executable_resource(status=status, calibration_mode=None)
    with pytest.raises(SecretoryResourceValidationError, match="model_handles"):
        resource.validate()


@pytest.mark.parametrize("status", (ExecutionStatus.UNAVAILABLE, ExecutionStatus.NOT_APPLICABLE))
def test_unavailable_and_not_applicable_pass_without_handles(status: ExecutionStatus) -> None:
    resource = _executable_resource(status=status, calibration_mode=None, model_handles=())
    resource.validate()


def test_conflict_requires_at_least_one_conflict_entry() -> None:
    resource = _executable_resource(
        status=ExecutionStatus.CONFLICT, calibration_mode=None, model_handles=(), conflicts=()
    )
    with pytest.raises(SecretoryResourceValidationError, match="conflict"):
        resource.validate()


def test_conflict_forbids_calibration_mode() -> None:
    resource = _executable_resource(
        status=ExecutionStatus.CONFLICT,
        calibration_mode=CalibrationMode.RELATIVE_UNCALIBRATED,
        conflicts=("two sources disagree on gene identity",),
    )
    with pytest.raises(SecretoryResourceValidationError, match="calibration_mode"):
        resource.validate()


def test_valid_conflict_resource_passes() -> None:
    resource = _executable_resource(
        status=ExecutionStatus.CONFLICT,
        calibration_mode=None,
        model_handles=(),
        conflicts=("two sources disagree on gene identity",),
    )
    resource.validate()


def test_empty_canonical_unit_is_rejected() -> None:
    with pytest.raises(SecretoryResourceValidationError, match="canonical_unit"):
        _executable_resource(canonical_unit="").validate()


def test_catalog_feature_off_must_be_empty() -> None:
    catalog = SecretoryResourceCatalog(
        target_id="hLF", host="Komagataella phaffii", feature_enabled=False, resources=(_executable_resource(),)
    )
    with pytest.raises(SecretoryResourceValidationError, match="feature_enabled=False"):
        catalog.validate()


def test_catalog_rejects_applicability_target_id_mismatch() -> None:
    catalog = SecretoryResourceCatalog(
        target_id="OPN_ALPHA_FULL_PROJECT",
        host="Komagataella phaffii",
        feature_enabled=True,
        resources=(_executable_resource(applicability=_applicability("hLF")),),
    )
    with pytest.raises(SecretoryResourceValidationError, match="target_id"):
        catalog.validate()


def test_catalog_rejects_duplicate_resource_ids() -> None:
    catalog = SecretoryResourceCatalog(
        target_id="hLF",
        host="Komagataella phaffii",
        feature_enabled=True,
        resources=(_executable_resource(), _executable_resource()),
    )
    with pytest.raises(SecretoryResourceValidationError, match="duplicate"):
        catalog.validate()


def test_catalog_by_category_filters() -> None:
    translocation = _executable_resource()
    disulfide = _executable_resource(
        resource_id="hLF_disulfide_bond_formation", category=ResourceCategory.DISULFIDE_BOND_FORMATION
    )
    catalog = SecretoryResourceCatalog(
        target_id="hLF", host="Komagataella phaffii", feature_enabled=True, resources=(translocation, disulfide)
    )
    assert catalog.by_category(ResourceCategory.DISULFIDE_BOND_FORMATION) == (disulfide,)
