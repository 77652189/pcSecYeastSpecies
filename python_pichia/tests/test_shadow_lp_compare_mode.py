from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest

from pcsec_pichia.analysis.shadow_lp import compare_secretion_capacity, run_shadow_hardcode_audit
from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.probe import TargetSpec
from pcsec_pichia.targets import load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def _inputs() -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT)


@lru_cache(maxsize=1)
def _builtin_targets() -> dict[str, TargetSpec]:
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}


@pytest.fixture(scope="module")
def comparisons() -> dict[str, object]:
    inputs = _inputs()
    targets = _builtin_targets()
    return {
        target_id: compare_secretion_capacity(
            inputs.prepared_model,
            targets[target_id],
            inputs.amino_acids,
            inputs.metabolic,
            inputs.secretory,
            inputs.combined,
            growth_rate=0.10,
        )
        for target_id in ("hLF", "OPN_ALPHA_FULL_PROJECT")
    }


def test_compare_mode_aligns_hlf_and_opn_within_tolerance(comparisons: dict[str, object]) -> None:
    for comparison in comparisons.values():
        assert comparison.reference_result.solver_mode == "reference"
        assert comparison.shadow_result.solver_mode == "shadow"
        assert comparison.reference_status_category == "optimal"
        assert comparison.shadow_status_category == "optimal"
        assert comparison.status_match is True
        assert comparison.within_tolerance is True
        assert comparison.objective_rel_diff <= 1e-4
        assert comparison.constraint_count_diff == 0


def test_compare_mode_returns_structured_key_flux_diffs(comparisons: dict[str, object]) -> None:
    for comparison in comparisons.values():
        assert comparison.key_flux_diffs
        assert comparison.shadow_result.exchange_reaction_id in comparison.key_flux_diffs
        row = comparison.key_flux_diffs[comparison.shadow_result.exchange_reaction_id]
        assert set(row) == {"reference", "shadow", "abs_diff"}
        assert row["abs_diff"] is not None


def test_reference_solver_is_only_allowed_in_compare_or_validation_boundary() -> None:
    audit = run_shadow_hardcode_audit()

    assert audit.passed
    assert set(audit.validation_only_files) == {"validation.py", "comparison.py"}
    assert audit.production_shadow_path_has_no_forbidden_reference_solver_call
