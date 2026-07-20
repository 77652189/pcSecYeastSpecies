from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest

from pcsec_pichia.secretory_resources import (
    CalibrationMode,
    EvidenceClass,
    ExecutionStatus,
    IN_SCOPE_CATEGORIES,
    ResourceCategory,
    SecretoryResourceCatalog,
    build_secretory_resource_catalog,
    plan_secretory_resource_constraints,
    summarize_secretory_resource_coverage,
    validate_secretory_resource_catalog,
)
from pcsec_pichia.secretory_resources import validation as secretory_resources_validation
from pcsec_pichia.targets import TargetSpec, load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def _builtin_targets() -> dict[str, TargetSpec]:
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}


def _target(target_id: str) -> TargetSpec:
    return _builtin_targets()[target_id]


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_catalog_covers_exactly_the_in_scope_categories(target_id: str) -> None:
    catalog = build_secretory_resource_catalog(_target(target_id))

    assert catalog.target_id == target_id
    assert catalog.feature_enabled is True
    assert {resource.category for resource in catalog.resources} == set(IN_SCOPE_CATEGORIES)
    assert len(catalog.resources) == len(IN_SCOPE_CATEGORIES)
    for resource in catalog.resources:
        assert resource.applicability.target_id == target_id


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_catalog_passes_validation_with_no_issues(target_id: str) -> None:
    catalog = build_secretory_resource_catalog(_target(target_id))
    result = validate_secretory_resource_catalog(catalog)
    assert result.is_valid, result.issues


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_translocation_and_trafficking_track_through_er(target_id: str) -> None:
    target = _target(target_id)
    catalog = build_secretory_resource_catalog(target)
    translocation = catalog.by_category(ResourceCategory.ER_TRANSLOCATION)[0]
    trafficking = catalog.by_category(ResourceCategory.VESICLE_TRAFFICKING)[0]

    if target.through_er:
        for resource in (translocation, trafficking):
            assert resource.status is ExecutionStatus.EXECUTABLE
            assert resource.calibration_mode is CalibrationMode.RELATIVE_UNCALIBRATED
            assert resource.model_handles
    else:
        for resource in (translocation, trafficking):
            assert resource.status is ExecutionStatus.NOT_APPLICABLE
            assert resource.model_handles == ()


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_disulfide_status_tracks_target_structural_profile(target_id: str) -> None:
    target = _target(target_id)
    resource = build_secretory_resource_catalog(target).by_category(ResourceCategory.DISULFIDE_BOND_FORMATION)[0]

    if target.disulfide_sites > 0:
        assert resource.status is ExecutionStatus.EXECUTABLE
        assert resource.model_handles
        assert all("dsb" in handle.lower() for handle in resource.model_handles)
    else:
        assert resource.status is ExecutionStatus.NOT_APPLICABLE
        assert resource.model_handles == ()


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_glycosylation_status_tracks_target_structural_profile(target_id: str) -> None:
    target = _target(target_id)
    resource = build_secretory_resource_catalog(target).by_category(ResourceCategory.GLYCOSYLATION)[0]

    if target.n_glycosylation_sites > 0 or target.o_glycosylation_sites > 0:
        assert resource.status is ExecutionStatus.EXECUTABLE
        assert resource.model_handles
    else:
        assert resource.status is ExecutionStatus.NOT_APPLICABLE
        assert resource.model_handles == ()


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_folding_chaperone_is_executable_with_confirmed_kcat_backing(target_id: str) -> None:
    # Every one of this category's handles is confirmed (not assumed) to
    # match a real, non-placeholder kcat entry in the base secretory enzyme
    # dataset (see catalog._folding_chaperone_resource); it is executable on
    # that evidence, same bar as translocation/disulfide/glycosylation/
    # vesicle_trafficking. evidence_class stays classifier_inferred because
    # the handle *identification* method (name match) is still weaker than a
    # dedicated stage lookup -- a separate concern from whether real
    # capacity data backs the identified handles.
    target = _target(target_id)
    resource = build_secretory_resource_catalog(target).by_category(ResourceCategory.FOLDING_CHAPERONE)[0]

    if target.through_er:
        assert resource.status is ExecutionStatus.EXECUTABLE
        assert resource.calibration_mode is CalibrationMode.RELATIVE_UNCALIBRATED
        assert resource.source.evidence_class is EvidenceClass.CLASSIFIER_INFERRED
        assert any("kcat" in limitation for limitation in resource.limitations)
    else:
        assert resource.status is ExecutionStatus.NOT_APPLICABLE


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_folding_chaperone_kcat_backing_is_verified_complete(target_id: str) -> None:
    # Direct regression guard for the specific evidence the executable
    # promotion above relies on: every handle actually matches a real complex
    # entry, not just "most of them" -- if the base enzyme dataset asset ever
    # changes and this drops below 100%, the promotion needs re-justifying.
    resource = build_secretory_resource_catalog(_target(target_id)).by_category(ResourceCategory.FOLDING_CHAPERONE)[0]
    note = next(limitation for limitation in resource.limitations if "kcat" in limitation)
    handle_count = len(resource.model_handles)
    assert f"{handle_count}/{handle_count} handles are confirmed" in note
    assert "All handles matched." in note


def test_target_specific_cost_kcat_backing_is_verified_absent() -> None:
    # Direct regression guard for the "real gap" claim: translation/
    # degradation handles must show zero matches against the secretory
    # enzyme dataset, confirming the missing piece is target-protein kdeg,
    # not a missed enzyme kcat this catalog failed to wire up.
    resource = build_secretory_resource_catalog(_target("hLF")).by_category(ResourceCategory.TARGET_SPECIFIC_COST)[0]
    note = next(limitation for limitation in resource.limitations if "Checked against" in limitation)
    assert "none of this resource's handles matched" in note


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_er_quality_control_kcat_backing_is_genuinely_partial(target_id: str) -> None:
    # er_quality_control is the one executable category where coverage is
    # NOT expected to be 100%: the "B" branch ERAD reactions emitted when a
    # PTM type is absent (see target_protein_plan.py's _misfolding_reactions)
    # have no distinct catalyzing complex. This guards against silently
    # regressing to "0 matched" (which would mean the join broke) or
    # silently becoming "100% matched" without updating the comment that
    # explains why it isn't.
    resource = build_secretory_resource_catalog(_target(target_id)).by_category(ResourceCategory.ER_QUALITY_CONTROL)[0]
    note = next(limitation for limitation in resource.limitations if "handles are confirmed backed" in limitation)
    matched = int(note.split("/", maxsplit=1)[0])
    assert 0 < matched < len(resource.model_handles)


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_er_quality_control_is_executable_for_every_target(target_id: str) -> None:
    # target_protein_plan.py emits at least one misfolding-stage reaction
    # unconditionally, whether the target goes through the ER or not (see
    # _misfolding_reactions/_non_secretory_folding_reactions), so this
    # category has no not_applicable branch, unlike translocation/disulfide/
    # glycosylation/vesicle_trafficking.
    target = _target(target_id)
    resource = build_secretory_resource_catalog(target).by_category(ResourceCategory.ER_QUALITY_CONTROL)[0]

    assert resource.status is ExecutionStatus.EXECUTABLE
    assert resource.calibration_mode is CalibrationMode.RELATIVE_UNCALIBRATED
    assert resource.model_handles
    assert any("misfold" in handle.lower() or "erad" in handle.lower() for handle in resource.model_handles)


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_target_specific_cost_is_executable_for_every_target(target_id: str) -> None:
    # r_{protein_id}_peptide_translation and r_{protein_id}_subunit_degradation
    # are both unconditional in target_protein_plan.py, so this category is
    # always executable too.
    target = _target(target_id)
    resource = build_secretory_resource_catalog(target).by_category(ResourceCategory.TARGET_SPECIFIC_COST)[0]

    assert resource.status is ExecutionStatus.EXECUTABLE
    assert resource.calibration_mode is CalibrationMode.RELATIVE_UNCALIBRATED
    assert resource.model_handles
    assert any("translation" in handle.lower() for handle in resource.model_handles)
    assert any("degradation" in handle.lower() for handle in resource.model_handles)


def test_target_specific_cost_never_shares_handles_across_targets() -> None:
    # Direct regression guard for the Round 0 contract's explicit
    # requirement (EXECUTION_PLAN.md/handoff.md): hLF's and OPN's
    # target-specific cost must never share a handle, i.e. one target's
    # evidence must never be copied into another target's resource.
    hlf_resource = build_secretory_resource_catalog(_target("hLF")).by_category(ResourceCategory.TARGET_SPECIFIC_COST)[0]
    opn_resource = build_secretory_resource_catalog(_target("OPN_ALPHA_FULL_PROJECT")).by_category(
        ResourceCategory.TARGET_SPECIFIC_COST
    )[0]

    assert set(hlf_resource.model_handles).isdisjoint(set(opn_resource.model_handles))
    assert hlf_resource.resource_id != opn_resource.resource_id


@pytest.mark.parametrize("target_id", ("hLF", "OPN_ALPHA_FULL_PROJECT"))
def test_er_quality_control_and_target_specific_cost_handles_do_not_overlap_other_categories(target_id: str) -> None:
    # Defensive uniqueness check: the two new categories claim the
    # "misfolding" and "translation"/"degradation" plan stages exclusively;
    # they must not silently re-claim a handle another category already
    # owns for this target.
    catalog = build_secretory_resource_catalog(_target(target_id))
    new_handles = set(catalog.by_category(ResourceCategory.ER_QUALITY_CONTROL)[0].model_handles) | set(
        catalog.by_category(ResourceCategory.TARGET_SPECIFIC_COST)[0].model_handles
    )
    other_handles: set[str] = set()
    for category in IN_SCOPE_CATEGORIES:
        if category in (ResourceCategory.ER_QUALITY_CONTROL, ResourceCategory.TARGET_SPECIFIC_COST):
            continue
        for resource in catalog.by_category(category):
            other_handles.update(resource.model_handles)

    assert new_handles.isdisjoint(other_handles)


def test_feature_off_returns_empty_but_valid_catalog() -> None:
    catalog = build_secretory_resource_catalog(_target("hLF"), feature_enabled=False)

    assert catalog.feature_enabled is False
    assert catalog.resources == ()
    catalog.validate()  # baseline/feature-off contract: empty catalog is valid


def test_coverage_summary_counts_match_catalog_size() -> None:
    catalog = build_secretory_resource_catalog(_target("hLF"))
    summary = summarize_secretory_resource_coverage(catalog)

    assert summary.target_id == "hLF"
    assert summary.total == len(IN_SCOPE_CATEGORIES)
    assert (
        summary.executable + summary.evidence_only + summary.unavailable + summary.not_applicable + summary.conflict
        == summary.total
    )
    for resource in catalog.resources:
        if resource.status in (ExecutionStatus.UNAVAILABLE, ExecutionStatus.CONFLICT):
            assert any(resource.resource_id in explanation for explanation in summary.unavailable_or_conflict_explanations)


def test_coverage_summary_on_feature_off_is_all_zero() -> None:
    catalog = build_secretory_resource_catalog(_target("hLF"), feature_enabled=False)
    summary = summarize_secretory_resource_coverage(catalog)
    assert summary.total == 0
    assert summary.unavailable_or_conflict_explanations == ()


def test_plan_never_calls_a_solver_and_marks_executable_ready() -> None:
    catalog = build_secretory_resource_catalog(_target("OPN_ALPHA_FULL_PROJECT"))
    plan = plan_secretory_resource_constraints(catalog)

    assert plan.backend == "none"
    assert len(plan.entries) == len(catalog.resources)
    for entry, resource in zip(plan.entries, catalog.resources):
        assert entry.resource_id == resource.resource_id
        if resource.status is ExecutionStatus.EXECUTABLE:
            assert entry.action == "relative_comparison_ready"
        else:
            assert entry.action == "not_planned"
            assert entry.reason


def test_validation_flags_forbidden_placeholder_handle_as_an_issue() -> None:
    catalog = build_secretory_resource_catalog(_target("hLF"))
    tampered_resources = tuple(
        resource if resource.category is not ResourceCategory.ER_TRANSLOCATION else _with_forbidden_handle(resource)
        for resource in catalog.resources
    )
    tampered = SecretoryResourceCatalog(
        target_id=catalog.target_id,
        host=catalog.host,
        feature_enabled=catalog.feature_enabled,
        resources=tampered_resources,
    )

    result = validate_secretory_resource_catalog(tampered)

    assert not result.is_valid
    assert any(issue.code == "forbidden_handle_literal" for issue in result.issues)


def test_validation_flags_deferred_category_as_out_of_round0_scope(monkeypatch: pytest.MonkeyPatch) -> None:
    # All seven real categories are in scope now, so there is no real
    # ResourceCategory member left to use as a "deferred" example. Simulate
    # the scenario this gate actually protects against instead: a future
    # category added to the enum without a matching IN_SCOPE_CATEGORIES
    # authorization update.
    catalog = build_secretory_resource_catalog(_target("hLF"))
    monkeypatch.setattr(
        secretory_resources_validation,
        "IN_SCOPE_CATEGORIES",
        tuple(category for category in IN_SCOPE_CATEGORIES if category is not ResourceCategory.ER_QUALITY_CONTROL),
    )

    result = validate_secretory_resource_catalog(catalog)

    assert not result.is_valid
    assert any(issue.code == "category_out_of_round0_scope" for issue in result.issues)


def _with_forbidden_handle(resource):
    from dataclasses import replace

    return replace(resource, model_handles=("1000",))
