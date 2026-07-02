from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from pcsec_pichia.alignment import AlignmentSummary
from pcsec_pichia.constraints import build_pcsec_constraint_matrices
from pcsec_pichia.engines.base import PichiaSimulationRequest, PichiaSimulationRunResult
from pcsec_pichia.loading import (
    CobraModel,
    PcSecPichiaInputs,
    load_pcsec_pichia_inputs,
    medium_condition_id_for,
    prepare_carbon_source_model,
    prepare_glucose_model,
    set_media_pp,
)
from pcsec_pichia.reports import build_candidate_table, write_outputs
from pcsec_pichia.screens import run_pcsec_ko_screen, run_pcsec_oe_screen
from pcsec_pichia.secretion_plan import target_reaction_plan
from pcsec_pichia.simulation import build_supported_target_model, solve_pcsec_maximize
from pcsec_pichia.targets import TargetSpec, load_targets


REPO_ROOT = Path(__file__).resolve().parents[2]
PROBE_DIR = REPO_ROOT / "local_runs" / "pichia_hlf_opn_probe"


def _require_local_probe_artifact(path: Path) -> None:
    if not path.exists():
        pytest.skip(f"local probe artifact is not available in this checkout: {path}")
    assert path.stat().st_size > 0, f"empty Stage 3 artifact: {path}"


def _tiny_exchange_model() -> CobraModel:
    rxns = [
        "Ex_nh4",
        "Ex_o2",
        "Ex_pi",
        "Ex_so4",
        "Ex_fe2",
        "Ex_h",
        "Ex_h2o",
        "Ex_na1",
        "Ex_k",
        "Ex_co2",
        "Ex_btn",
        "Ex_arg_L",
        "Ex_glc_D",
        "BIOMASS",
        "Ex_glyc",
        "BIOMASS_glyc",
        "Ex_meoh",
        "BIOMASS_meoh",
    ]
    return CobraModel(
        source_file="tiny",
        rxns=rxns,
        mets=["dummy"],
        genes=[],
        lb=np.full(len(rxns), -999.0),
        ub=np.full(len(rxns), 999.0),
        b=np.zeros(1),
        s_matrix=sparse.csc_matrix((1, len(rxns))),
        rules=[],
        gr_rules=[],
    )


def test_formal_modules_expose_scope_freeze_entrypoints() -> None:
    assert callable(load_targets)
    assert callable(target_reaction_plan)
    assert callable(build_pcsec_constraint_matrices)
    assert callable(build_supported_target_model)
    assert callable(solve_pcsec_maximize)
    assert callable(run_pcsec_ko_screen)
    assert callable(run_pcsec_oe_screen)
    assert callable(build_candidate_table)
    assert callable(write_outputs)
    assert TargetSpec is not None
    assert AlignmentSummary(target_id="OPN", success=True, baseline_source="matlab").success is True
    assert PcSecPichiaInputs is not None
    assert callable(load_pcsec_pichia_inputs)




def test_pichia_public_request_and_result_defaults_are_explicit() -> None:
    request = PichiaSimulationRequest(target_id="OPN", candidate_id="OPN_ALPHA_FULL_PROJECT")
    assert request.compatibility_mode == "corrected"
    assert request.glycosylation_mode == "native"
    assert request.enable_ribosome_translation_constraint is False
    assert request.enable_misfolding_constraint is False

    result = PichiaSimulationRunResult(
        success=True,
        target_id="OPN",
        candidate_id="OPN_ALPHA_FULL_PROJECT",
        mu=0.10,
        production_ratio=1e-8,
        media_type=4,
        message="ok",
    )
    assert result.result_status == "draft"
    assert result.constraint_counts == {}
    assert result.candidate_table_path is None
    assert result.tradeoff_path is None
    assert result.alignment_summary == {}


def test_medium_compatibility_modes_remain_explicit() -> None:
    model = _tiny_exchange_model()
    compat = set_media_pp(model, media_type=4, compatibility_mode="matlab_compat")
    corrected = set_media_pp(model, media_type=4, compatibility_mode="corrected")

    minimal = ("Ex_nh4", "Ex_o2", "Ex_pi", "Ex_so4", "Ex_fe2", "Ex_h", "Ex_h2o", "Ex_na1", "Ex_k", "Ex_co2")
    for reaction_id in minimal:
        assert compat.lb[compat.reaction_index[reaction_id]] == 0.0
        assert corrected.lb[corrected.reaction_index[reaction_id]] == -1000.0

    assert compat.lb[compat.reaction_index["Ex_btn"]] == -2.0
    assert compat.lb[compat.reaction_index["Ex_arg_L"]] == -0.08

    glucose = prepare_glucose_model(model, media_type=4, compatibility_mode="matlab_compat")
    assert glucose.lb[glucose.reaction_index["Ex_glc_D"]] == -1000.0
    assert glucose.lb[glucose.reaction_index["Ex_o2"]] == -1000.0

    default_glucose = prepare_glucose_model(model, media_type=4)
    for reaction_id in minimal:
        assert default_glucose.lb[default_glucose.reaction_index[reaction_id]] == -1000.0


def test_carbon_source_model_preparation_switches_exchange_bounds() -> None:
    model = _tiny_exchange_model()

    expected = {
        "glucose": {
            "Ex_glc_D": (-1000.0, 1000.0),
            "Ex_glyc": (0.0, 1000.0),
            "BIOMASS_glyc": (0.0, 0.0),
            "Ex_meoh": (0.0, 1000.0),
            "BIOMASS_meoh": (0.0, 0.0),
        },
        "glycerol": {
            "Ex_glc_D": (0.0, 1000.0),
            "Ex_glyc": (-1000.0, 1000.0),
            "BIOMASS_glyc": (0.0, 1000.0),
            "Ex_meoh": (0.0, 1000.0),
            "BIOMASS_meoh": (0.0, 0.0),
        },
        "methanol": {
            "Ex_glc_D": (0.0, 1000.0),
            "Ex_glyc": (0.0, 1000.0),
            "BIOMASS_glyc": (0.0, 0.0),
            "Ex_meoh": (-1000.0, 1000.0),
            "BIOMASS_meoh": (0.0, 1000.0),
        },
        "glucose_glycerol": {
            "Ex_glc_D": (-1000.0, 1000.0),
            "Ex_glyc": (-1000.0, 1000.0),
            "BIOMASS_glyc": (0.0, 1000.0),
            "Ex_meoh": (0.0, 1000.0),
            "BIOMASS_meoh": (0.0, 0.0),
        },
        "glycerol_methanol": {
            "Ex_glc_D": (0.0, 1000.0),
            "Ex_glyc": (-1000.0, 1000.0),
            "BIOMASS_glyc": (0.0, 1000.0),
            "Ex_meoh": (-1000.0, 1000.0),
            "BIOMASS_meoh": (0.0, 1000.0),
        },
    }

    for carbon_source_id, bounds in expected.items():
        prepared = prepare_carbon_source_model(model, media_type=4, carbon_source_id=carbon_source_id)
        for reaction_id, (lower_bound, upper_bound) in bounds.items():
            index = prepared.reaction_index[reaction_id]
            assert prepared.lb[index] == lower_bound
            assert prepared.ub[index] == upper_bound


def test_medium_condition_id_tracks_carbon_source_and_supplement_layer() -> None:
    assert medium_condition_id_for(4, "corrected", "methanol") == "methanol_ynb_core_aa_corrected"
    assert medium_condition_id_for(2, "corrected", "glycerol") == "glycerol_ynb_minimal_corrected"
    assert medium_condition_id_for(4, "corrected", "glucose_glycerol") == "glucose_glycerol_ynb_core_aa_corrected"
    assert medium_condition_id_for(5, "corrected", "glycerol_methanol") == "glycerol_methanol_ynb_all_aa_corrected"


def test_stage3_optional_constraint_artifacts_remain_valid() -> None:
    for prefix, expected_rows in (
        ("stage3_builtin_optional", 24),
        ("stage3_custom_input_optional", 14),
    ):
        summary_path = PROBE_DIR / f"{prefix}_summary.json"
        report_path = PROBE_DIR / f"{prefix}_REPORT.md"
        candidates_path = PROBE_DIR / f"{prefix}_candidates.csv"
        tradeoff_path = PROBE_DIR / f"{prefix}_tradeoff.csv"
        for path in (summary_path, report_path, candidates_path, tradeoff_path):
            _require_local_probe_artifact(path)

        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        for target in summary["targets"]:
            smoke = target["target_secretion_smoke"]
            reference = smoke["target_pcsec_reference_max"]
            counts = smoke["target_pcsec_constraint_counts"]
            assert smoke["target_exchange_max"]["success"] is True
            assert reference["success"] is True
            assert counts["ribosome_translation"] == 1
            assert counts["misfolding"] == 1418
            assert smoke["candidate_table"]
            assert smoke["pcsec_growth_tradeoff"]

        with candidates_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        assert len(rows) == expected_rows
        assert "complex_subunit_ids" in (reader.fieldnames or [])
        assert "complex_subunit_stoichiometry" in (reader.fieldnames or [])
        assert any((row.get("complex_subunit_ids") or "").strip() for row in rows)

        report = report_path.read_text(encoding="utf-8")
        assert "MATLAB" in report
        assert "alignment" in report.lower()


def test_matlab_alignment_summary_preserves_hlf_diagnostic() -> None:
    summary_path = PROBE_DIR / "matlab_stage3_alignment" / "matlab_stage3_alignment_summary.json"
    _require_local_probe_artifact(summary_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    by_target = {item["target_id"]: item for item in summary}

    assert by_target["OPN_ALPHA_FULL_PROJECT"]["success"] is True
    assert by_target["OPN_ALPHA_FULL_PROJECT"]["added_reactions"] == 31
    assert by_target["OPN_ALPHA_FULL_PROJECT"]["added_metabolites"] == 26

    hlf = by_target["hLF"]
    assert hlf["success"] is False
    assert hlf["error_identifier"] == "MATLAB:sizeDimensionsMustMatch"
    assert any("calculateMW.m:35" in frame for frame in hlf["error_stack"])
