from __future__ import annotations

import sys
import types

import pytest

import app.services.per_strain_oe_candidate_run as run_mod
from app.services.per_strain_oe_candidate_run import run_next_oe_candidate_analysis


def _seed_fake_engine(monkeypatch, engine_fn) -> None:
    """Install a fake `pcsec_pichia.next_oe_candidates` so the service's lazy import resolves to it.

    Reverted automatically by monkeypatch; also stubs the analysis-reload guard to a no-op so the
    test never touches the real engine.
    """
    monkeypatch.setattr(run_mod, "_ensure_pcsec_pichia_analysis_api", lambda: None)
    monkeypatch.setitem(sys.modules, "pcsec_pichia", sys.modules.get("pcsec_pichia") or types.ModuleType("pcsec_pichia"))
    fake_module = types.ModuleType("pcsec_pichia.next_oe_candidates")
    fake_module.analyze_next_oe_candidates = engine_fn
    monkeypatch.setitem(sys.modules, "pcsec_pichia.next_oe_candidates", fake_module)


def _engine_payload() -> dict[str, object]:
    return {
        "target_id": "hLF",
        "modified_solve_success": True,
        "modified_objective_value": 1.5,
        "lp_attribution_status": "draft_cost_slope_analysis",
        "oe_actionable_bottlenecks": [
            {"reaction_id": "sec_A_complex_formation", "secretory_process": "disulfide_folding", "abs_marginal": 5.0},
        ],
        "floor_constraints_not_oe_addressable": [{"reaction_id": "floor_x"}],
        "dose_response": {
            "baseline_objective": 1.2,
            "shapes_by_reaction": {
                "sec_A_complex_formation": {
                    "reaction_id": "sec_A_complex_formation",
                    "shape": "saturating",
                    "max_relative_gain": 0.2,
                    "half_gain_factor": 2.0,
                }
            },
        },
        "applied_modifications": {"ko_reaction_ids": [], "oe_reaction_ids": ["sec_A_complex_formation"], "oe_factor": 2.0},
        "modification_warnings": ["OE reaction not found in model, skipped: ghost"],
        "carbon_source_id": "glucose",
        "medium_condition_id": "cond_glucose",
        "mu": 0.10,
        "top_n": 6,
    }


def test_service_merges_engine_output_with_c1_readout(monkeypatch) -> None:
    _seed_fake_engine(monkeypatch, lambda **kwargs: _engine_payload())

    readout = run_next_oe_candidate_analysis(target_id="hLF", oe_reaction_ids=("sec_A_complex_formation",))

    # C1 ranked candidates assembled from the engine bottlenecks + dose-response.
    assert [c["reaction"] for c in readout["candidates"]] == ["sec_A_complex_formation"]
    assert readout["dose_response_available"] is True
    assert readout["candidates"][0]["shape"] == "saturating"
    # engine metadata carried through
    assert readout["modified_solve_success"] is True
    assert readout["modified_objective_value"] == 1.5
    assert readout["applied_modifications"]["oe_reaction_ids"] == ["sec_A_complex_formation"]
    assert readout["modification_warnings"] == ["OE reaction not found in model, skipped: ghost"]
    assert readout["floor_constraints_not_oe_addressable"] == [{"reaction_id": "floor_x"}]
    assert readout["error"] is None
    assert readout["caveats"]  # honest caveats always present


def test_service_forwards_call_arguments_to_engine(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_engine(**kwargs):
        captured.update(kwargs)
        return _engine_payload()

    _seed_fake_engine(monkeypatch, fake_engine)
    run_next_oe_candidate_analysis(
        target_id="hLF",
        ko_reaction_ids=("R_ko",),
        oe_reaction_ids=("sec_A_complex_formation",),
        mu=0.2,
        carbon_source_id="glycerol",
        enable_misfolding_constraint=True,
        top_n=4,
    )
    assert captured["target_id"] == "hLF"
    assert captured["ko_reaction_ids"] == ("R_ko",)
    assert captured["oe_reaction_ids"] == ("sec_A_complex_formation",)
    assert captured["mu"] == 0.2
    assert captured["carbon_source_id"] == "glycerol"
    assert captured["enable_misfolding_constraint"] is True
    assert captured["top_n"] == 4


def test_service_returns_error_readout_when_engine_raises(monkeypatch) -> None:
    def boom(**kwargs):
        raise RuntimeError("solver exploded")

    _seed_fake_engine(monkeypatch, boom)
    readout = run_next_oe_candidate_analysis(target_id="hLF")

    assert readout["candidates"] == []
    assert readout["dose_response_available"] is False
    assert readout["modified_solve_success"] is False
    assert readout["error"].startswith("RuntimeError")
    assert "solver exploded" in readout["error"]
