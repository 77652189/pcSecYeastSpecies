from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest

from pcsec_pichia.constraints import PcSecConstraintResult, build_pcsec_constraints, summarize_pcsec_constraints
from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.targets import TargetSpec, load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def _inputs() -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT)


@lru_cache(maxsize=1)
def _builtin_targets() -> dict[str, TargetSpec]:
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}


@pytest.mark.parametrize(
    ("target_id", "expected_columns", "expected_counts"),
    (
        (
            "OPN_ALPHA_FULL_PROJECT",
            29057,
            {
                "stoichiometric": 20221,
                "metabolic_coupling": 2732,
                "secretory_coupling": 58,
                "protein_mass": 2,
                "proteasome": 1,
                "ribosome_assembly": 1,
                "ribosome_translation": 0,
                "misfolding": 0,
                "mitochondrial": 1,
                "eq_total": 23015,
                "ub_total": 1,
            },
        ),
        (
            "hLF",
            29068,
            {
                "stoichiometric": 20230,
                "metabolic_coupling": 2732,
                "secretory_coupling": 58,
                "protein_mass": 2,
                "proteasome": 1,
                "ribosome_assembly": 1,
                "ribosome_translation": 0,
                "misfolding": 0,
                "mitochondrial": 1,
                "eq_total": 23024,
                "ub_total": 1,
            },
        ),
    ),
)
def test_builtin_targets_build_default_pcsec_constraint_summary(
    target_id: str,
    expected_columns: int,
    expected_counts: dict[str, int],
) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    summary = summarize_pcsec_constraints(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
    )

    assert summary["target_id"] == target_id
    assert summary["supported"] is True
    assert summary["build_status"] == "stoichiometric_target_model_built"
    assert summary["constraint_counts"] == expected_counts
    assert summary["eq_shape"] == (expected_counts["eq_total"], expected_columns)
    assert summary["ub_shape"] == (expected_counts["ub_total"], expected_columns)
    assert summary["matlab_alignment_status"] == "pending"


@pytest.mark.parametrize(
    ("target_id", "expected_eq_total", "expected_columns"),
    (
        ("OPN_ALPHA_FULL_PROJECT", 24434, 29057),
        ("hLF", 24443, 29068),
    ),
)
def test_builtin_targets_build_optional_pcsec_constraints(
    target_id: str,
    expected_eq_total: int,
    expected_columns: int,
) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    result = build_pcsec_constraints(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        write_ribosome_translation_constraint=True,
        write_misfolding_constraints=True,
    )

    assert isinstance(result, PcSecConstraintResult)
    assert result.supported is True
    assert result.constraint_counts["ribosome_translation"] == 1
    assert result.constraint_counts["misfolding"] == 1418
    assert result.constraint_counts["eq_total"] == expected_eq_total
    assert result.eq_shape == (expected_eq_total, expected_columns)
    assert result.ub_shape == (1, expected_columns)
    assert result.A_eq is not None
    assert result.b_eq is not None
    assert result.A_eq.shape[0] == result.b_eq.shape[0]
    assert result.A_ub is not None
    assert result.b_ub is not None
    assert result.A_ub.shape[0] == result.b_ub.shape[0]
    assert result.matlab_alignment_status == "pending"
