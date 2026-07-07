from __future__ import annotations

import pytest

from pcsec_pichia.analysis.shadow_lp import run_shadow_validation_matrix


@pytest.fixture(scope="module")
def validation_matrix() -> object:
    return run_shadow_validation_matrix(
        target_ids=("hLF", "OPN_ALPHA_FULL_PROJECT"),
        growth_rates=(0.05, 0.10, 0.15),
    )


def test_shadow_validation_matrix_records_all_target_growth_cases(validation_matrix: object) -> None:
    assert validation_matrix.target_ids == ("hLF", "OPN_ALPHA_FULL_PROJECT")
    assert validation_matrix.growth_rates == (0.05, 0.10, 0.15)
    assert len(validation_matrix.cases) == 6
    assert {
        (case.target_id, case.growth_rate)
        for case in validation_matrix.cases
    } == {
        ("hLF", 0.05),
        ("hLF", 0.10),
        ("hLF", 0.15),
        ("OPN_ALPHA_FULL_PROJECT", 0.05),
        ("OPN_ALPHA_FULL_PROJECT", 0.10),
        ("OPN_ALPHA_FULL_PROJECT", 0.15),
    }


def test_shadow_validation_matrix_requires_default_growth_alignment(validation_matrix: object) -> None:
    defaults = {
        case.target_id: case
        for case in validation_matrix.cases
        if case.growth_rate == 0.10
    }

    assert defaults["hLF"].alignment_status == "aligned"
    assert defaults["hLF"].objective_rel_diff <= 1e-4
    assert defaults["hLF"].constraint_count_diff == 0
    assert defaults["OPN_ALPHA_FULL_PROJECT"].alignment_status == "aligned"
    assert defaults["OPN_ALPHA_FULL_PROJECT"].objective_rel_diff <= 1e-4
    assert defaults["OPN_ALPHA_FULL_PROJECT"].constraint_count_diff == 0
    assert validation_matrix.all_required_defaults_aligned is True


def test_shadow_validation_matrix_non_aligned_cases_are_structured(validation_matrix: object) -> None:
    allowed_statuses = {
        "optimal",
        "infeasible",
        "unbounded",
        "timeout_iteration_limit",
        "unavailable_backend",
        "numerical_failure",
        "exception",
    }
    for case in validation_matrix.cases:
        assert case.reference_status in allowed_statuses
        assert case.shadow_status in allowed_statuses
        assert case.alignment_status in {"aligned", "review_required"}
        assert isinstance(case.key_flux_diffs, dict)
        if not case.within_tolerance:
            assert case.alignment_status == "review_required"
