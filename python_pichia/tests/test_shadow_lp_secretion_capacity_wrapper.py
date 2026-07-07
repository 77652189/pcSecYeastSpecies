from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest

from pcsec_pichia.analysis.shadow_lp import (
    SHADOW_CONSTRAINT_COUNT_KEYS,
    run_shadow_hardcode_audit,
    solve_shadow_secretion_capacity,
)
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
def wrapper_results() -> dict[str, object]:
    inputs = _inputs()
    targets = _builtin_targets()
    return {
        target_id: solve_shadow_secretion_capacity(
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


def test_shadow_secretion_capacity_wrapper_preserves_result_contract(wrapper_results: dict[str, object]) -> None:
    for target_id, result in wrapper_results.items():
        assert result.success is True
        assert result.target_id == target_id
        assert result.objective_value is not None
        assert result.secretion_flux == pytest.approx(result.objective_value)
        assert result.exchange_reaction_id
        assert result.key_fluxes
        assert result.solver_mode == "shadow"
        assert result.result_status == "shadow_lp_capacity"
        assert result.shadow_metadata is not None
        assert result.shadow_metadata["reference_solver_used"] is False
        assert result.shadow_metadata["canonical_final_layer"] == "mitochondrial"
        assert result.objective_value == pytest.approx(result.shadow_metadata["final_layer"]["objective"])


def test_shadow_secretion_capacity_constraint_counts_are_complete(wrapper_results: dict[str, object]) -> None:
    for result in wrapper_results.values():
        assert tuple(result.constraint_counts) == SHADOW_CONSTRAINT_COUNT_KEYS
        assert result.constraint_counts["eq_total"] > 0
        assert result.constraint_counts["ub_total"] > 0
        assert result.constraint_counts["ribosome_translation"] == 0
        assert result.constraint_counts["misfolding"] == 0
        assert result.shadow_metadata is not None
        assert "ribosome_translation" in result.shadow_metadata["skipped_layers"]
        assert "misfolding" in result.shadow_metadata["skipped_layers"]


def test_shadow_secretion_capacity_path_has_no_forbidden_reference_solver_call() -> None:
    audit = run_shadow_hardcode_audit()

    assert audit.passed
    assert audit.production_shadow_path_has_no_forbidden_reference_solver_call
    assert "comparison.py" in audit.validation_only_files
    assert not any(path.endswith("comparison.py") for path in audit.scanned_files)
