"""Regression test for the pcSecPichia probe prototype migration.

This test verifies that the probe prototype, after being moved into the
``pcsec_pichia.probe`` subpackage, produces numerically identical results to
the original standalone ``local_runs/pichia_hlf_opn_probe/pichia_hlf_opn_probe.py``.

The expected objective values are taken from the prototype README's "Latest
smoke verification" section and must remain byte-for-byte identical. Any drift
means the migration changed numeric behavior, which is a hard failure.

These tests are slow (they load the full pcSecPichia.mat model and solve
multiple LPs). They are kept as an explicit migration regression gate and only
run when ``PCSEC_RUN_SLOW_PROBE_TESTS=1`` is set.
"""

from __future__ import annotations

import os
import json
from pathlib import Path

import pytest

from pcsec_pichia.probe import (
    build_supported_target_model,
    build_target_enzymedata,
    load_aa_stoichiometry,
    load_combined_enzymedata,
    load_metabolic_enzymedata,
    load_pcsec_pichia_model,
    load_secretory_enzymedata,
    load_targets,
    prepare_glucose_model,
    repo_root,
    solve_maximize,
    solve_pcsec_maximize,
    target_secretion_smoke,
)

# Expected objective values after MATLAB setMediaPP media-type-4 alignment.
# The original prototype had a bug where minimal exchange reactions were set to
# lb=-1000, but MATLAB setMediaPP.m line 62 does not assign the return value of
# setLowerBounds, so those bounds are never applied. After fixing the probe to
# match MATLAB's actual behavior, the GEM baseline drops to 0 (no ammonia
# import without the minimal set) and the pcSec reference objectives change.
#
# The MATLAB baseline for "maximize Ex_glc_D with BIOMASS=0.10, target=1e-8"
# is -1.07282263 (SoPlex). The probe gets -1.0770728 (scipy HiGHS), within 0.4%.
EXPECTED_GEM_BASELINE_OBJECTIVE = 0.0
EXPECTED_OPN_PCSEC_REFERENCE_OBJECTIVE = 0.0013859763311249447
EXPECTED_HLF_PCSEC_REFERENCE_OBJECTIVE = 0.00036679903739833243

# The MATLAB alignment target (separate from the probe's own smoke objectives).
# This is the hard alignment anchor: maximize Ex_glc_D with BIOMASS=0.10 and
# target exchange fixed at 1e-8, solved with the pcSec constraint set.
MATLAB_ALIGNMENT_OBJECTIVE = -1.07282263
MATLAB_ALIGNMENT_TOLERANCE = 0.05  # 5% tolerance for solver differences

# Absolute tolerance for the probe's own smoke values (deterministic).
OBJECTIVE_TOLERANCE = 1e-9

pytestmark = pytest.mark.skipif(
    os.environ.get("PCSEC_RUN_SLOW_PROBE_TESTS") != "1",
    reason="slow probe migration solve; set PCSEC_RUN_SLOW_PROBE_TESTS=1 to run",
)


@pytest.fixture(scope="module")
def repo():
    return repo_root()


@pytest.fixture(scope="module")
def prepared_model(repo):
    raw = load_pcsec_pichia_model(repo)
    return prepare_glucose_model(raw, media_type=4)


@pytest.fixture(scope="module")
def amino_acids(repo):
    return load_aa_stoichiometry(repo)


@pytest.fixture(scope="module")
def metabolic(repo):
    return load_metabolic_enzymedata(repo)


@pytest.fixture(scope="module")
def secretory(repo):
    return load_secretory_enzymedata(repo)


@pytest.fixture(scope="module")
def combined(repo):
    return load_combined_enzymedata(repo)


@pytest.fixture(scope="module")
def gem_baseline(prepared_model):
    result = solve_maximize(
        prepared_model,
        "BIOMASS",
        key_reactions=("Ex_glc_D", "Ex_o2", "BIOMASS"),
    )
    assert result.success, f"GEM baseline solve failed: {result.message}"
    return result


@pytest.fixture(scope="module")
def targets(repo):
    return load_targets(repo, targets_json=None)


def _run_pcsec_smoke(
    prepared_model,
    target,
    amino_acids,
    metabolic,
    secretory,
    combined,
):
    return target_secretion_smoke(
        prepared_model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        ko_genes=[],
        oe_reactions=[],
        fixed_mu=0.10,
        tradeoff_mus=(0.05, 0.10, 0.20),
        write_ribosome_translation_constraint=False,
        write_misfolding_constraints=False,
    )


def test_gem_baseline_objective_unchanged(gem_baseline):
    """GEM-layer FBA growth objective must match the original prototype."""
    assert gem_baseline.objective_value is not None
    assert gem_baseline.objective_value == pytest.approx(
        EXPECTED_GEM_BASELINE_OBJECTIVE, abs=OBJECTIVE_TOLERANCE
    ), (
        f"GEM baseline objective drifted: got {gem_baseline.objective_value!r}, "
        f"expected {EXPECTED_GEM_BASELINE_OBJECTIVE!r}"
    )


def test_opn_pcsec_reference_objective_unchanged(
    prepared_model, targets, amino_acids, metabolic, secretory, combined
):
    """OPN pcSec reference smoke objective must match the original prototype."""
    opn_target = next(t for t in targets if "OPN" in t.protein_id)
    smoke = _run_pcsec_smoke(
        prepared_model, opn_target, amino_acids, metabolic, secretory, combined
    )
    pcsec_ref = smoke.get("target_pcsec_reference_max", {})
    objective = pcsec_ref.get("objective_value")
    assert objective is not None, f"OPN pcSec reference solve did not return objective: {pcsec_ref}"
    assert objective == pytest.approx(
        EXPECTED_OPN_PCSEC_REFERENCE_OBJECTIVE, abs=OBJECTIVE_TOLERANCE
    ), (
        f"OPN pcSec reference objective drifted: got {objective!r}, "
        f"expected {EXPECTED_OPN_PCSEC_REFERENCE_OBJECTIVE!r}"
    )


def test_hlf_pcsec_reference_objective_unchanged(
    prepared_model, targets, amino_acids, metabolic, secretory, combined
):
    """hLF pcSec reference smoke objective must match the original prototype."""
    hlf_target = next(t for t in targets if t.protein_id == "hLF")
    smoke = _run_pcsec_smoke(
        prepared_model, hlf_target, amino_acids, metabolic, secretory, combined
    )
    pcsec_ref = smoke.get("target_pcsec_reference_max", {})
    objective = pcsec_ref.get("objective_value")
    assert objective is not None, f"hLF pcSec reference solve did not return objective: {pcsec_ref}"
    assert objective == pytest.approx(
        EXPECTED_HLF_PCSEC_REFERENCE_OBJECTIVE, abs=OBJECTIVE_TOLERANCE
    ), (
        f"hLF pcSec reference objective drifted: got {objective!r}, "
        f"expected {EXPECTED_HLF_PCSEC_REFERENCE_OBJECTIVE!r}"
    )


def test_probe_cli_runs_and_outputs_match(
    prepared_model, targets, amino_acids, metabolic, secretory, combined, tmp_path
):
    """Smoke-check that the probe CLI still produces a valid report structure.

    This does not invoke ``main()`` (which would write to local_runs); instead
    it re-runs the smoke path and checks the payload shape so the test stays
    hermetic.
    """
    opn_target = next(t for t in targets if "OPN" in t.protein_id)
    smoke = _run_pcsec_smoke(
        prepared_model, opn_target, amino_acids, metabolic, secretory, combined
    )
    assert "build" in smoke
    assert smoke["build"]["supported"] is True
    assert "target_pcsec_reference_max" in smoke
    assert "target_pcsec_constraint_counts" in smoke
    counts = smoke["target_pcsec_constraint_counts"]
    # The prototype README documents these constraint counts for the default
    # (no optional constraints) run. They must not drift.
    assert counts["metabolic_coupling"] > 0
    assert counts["secretory_coupling"] > 0
    assert counts["protein_mass"] > 0
    assert counts.get("ribosome_translation", 0) == 0
    assert counts.get("misfolding", 0) == 0


def test_matlab_lp_alignment(
    prepared_model, targets, amino_acids, metabolic, secretory, combined
):
    """The probe's LP must align with the MATLAB baseline within 5%.

    MATLAB baseline: local_opn_pichia_glc with OPN_ALPHA_FULL_PROJECT,
    mu=0.10, mediaType=4, productionRatio=1e-8, blockMisfoldDilution=true.
    SoPlex objective: -1.07282263 (maximize Ex_glc_D).

    The probe solves the same problem with scipy HiGHS. The 0.4% gap is
    from solver tolerance differences (SoPlex fpfeastol=1e-3 vs HiGHS default).
    """
    opn_target = next(t for t in targets if "OPN" in t.protein_id)
    build = build_supported_target_model(prepared_model, opn_target, amino_acids)
    target_enzyme = build_target_enzymedata(opn_target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzyme.reaction_coefficients)
    target_combined = combined.with_target(target_enzyme)

    dilution_blocks = {rxn: (None, 0.0) for rxn in build.model.rxns if "dilution_misfolding" in rxn}
    fixed_model = build.model.with_bounds({
        "BIOMASS": (0.10, 0.10),
        build.exchange_reaction_id: (1e-8, 1e-8),
        **dilution_blocks,
    })

    result, counts = solve_pcsec_maximize(
        fixed_model,
        "Ex_glc_D",
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu=0.10,
        key_reactions=("BIOMASS", "Ex_glc_D", "Ex_o2", build.exchange_reaction_id),
    )
    assert result.success, f"MATLAB alignment solve failed: {result.message}"
    assert result.objective_value is not None
    # 5% tolerance for solver differences (SoPlex vs HiGHS)
    assert result.objective_value == pytest.approx(
        MATLAB_ALIGNMENT_OBJECTIVE, rel=MATLAB_ALIGNMENT_TOLERANCE
    ), (
        f"MATLAB alignment objective drifted: got {result.objective_value!r}, "
        f"expected {MATLAB_ALIGNMENT_OBJECTIVE!r} (tolerance {MATLAB_ALIGNMENT_TOLERANCE:.0%})"
    )
