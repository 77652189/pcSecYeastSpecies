from __future__ import annotations

import types

import pcsec_pichia.next_oe_candidates as mod
from pcsec_pichia.next_oe_candidates import analyze_next_oe_candidates


def _bottleneck(reaction_id: str, marginal: float) -> dict[str, object]:
    return {
        "bound_type": "upper",
        "reaction_id": reaction_id,
        "secretory_process": "disulfide_folding",
        "abs_marginal": marginal,
        "marginal": marginal,
        "oe_actionable": True,
    }


def _install(monkeypatch, *, bottlenecks, sim_success=True, sweep_success=True):
    captured: dict[str, object] = {}

    inputs = types.SimpleNamespace(
        prepared_model=types.SimpleNamespace(
            reaction_index={"sec_A_complex_formation": 0, "sec_B_complex_formation": 1, "R_ko": 2},
            rxns=("sec_A_complex_formation", "sec_B_complex_formation", "R_ko"),
        ),
        amino_acids=object(),
        metabolic=object(),
        secretory=object(),
        combined=object(),
        carbon_source_id="glucose",
        medium_condition_id="cond_glucose",
    )
    monkeypatch.setattr(mod, "load_pcsec_pichia_inputs", lambda *a, **k: inputs)
    monkeypatch.setattr(mod, "_resolve_target", lambda request, root: types.SimpleNamespace(target_id=request.target_id))
    monkeypatch.setattr(mod, "build_secretion_plan", lambda target: object())
    monkeypatch.setattr(mod, "build_pcsec_constraints", lambda *a, **k: types.SimpleNamespace(constraint_counts={}))

    def fake_solve(*a, **k):
        captured["solve_strain_modifications"] = k.get("strain_modifications")
        return types.SimpleNamespace(success=sim_success, objective_value=1.23)

    monkeypatch.setattr(mod, "solve_secretion_capacity", fake_solve)
    monkeypatch.setattr(mod, "analyze_target_protein_lp_attribution", lambda *a, **k: object())
    monkeypatch.setattr(
        mod,
        "summarize_protein_lp_attribution",
        lambda lp: {
            "result_status": "draft_cost_slope_analysis",
            "oe_actionable_bottlenecks": bottlenecks,
            "floor_constraints_not_oe_addressable": [{"reaction_id": "floor_x"}],
        },
    )

    def fake_sweep(*a, **k):
        captured["sweep_reactions"] = list(k.get("reactions") or [])
        captured["sweep_strain_modifications"] = k.get("strain_modifications")
        reaction_points = {r: ((1.0, 1.23), (2.0, 1.5)) for r in (k.get("reactions") or [])}
        return types.SimpleNamespace(
            success=sweep_success,
            reaction_points=reaction_points,
            baseline_objective=1.23,
            tested_factors=(1.5, 2.0),
            warnings=(),
        )

    monkeypatch.setattr(mod, "run_oe_dose_response_sweep", fake_sweep)
    monkeypatch.setattr(
        mod,
        "classify_oe_dose_response_sweep",
        lambda reaction_points, baseline: [types.SimpleNamespace(reaction_id=r) for r in reaction_points],
    )
    monkeypatch.setattr(
        mod,
        "summarize_oe_dose_response_shape",
        lambda shape: {"reaction_id": shape.reaction_id, "shape": "saturating", "max_relative_gain": 0.2, "half_gain_factor": 2.0},
    )
    return captured


def test_two_pass_passes_modifications_to_both_solve_and_sweep_and_assembles_shapes(monkeypatch) -> None:
    captured = _install(
        monkeypatch,
        bottlenecks=[_bottleneck("sec_A_complex_formation", 5.0), _bottleneck("sec_B_complex_formation", 3.0)],
    )
    result = analyze_next_oe_candidates(
        target_id="hLF",
        oe_reaction_ids=("sec_A_complex_formation",),
        ko_reaction_ids=("R_ko",),
        top_n=2,
    )

    # Same stacked modifications drive the modified solve AND the dose-response baseline.
    solve_mods = captured["solve_strain_modifications"]
    assert solve_mods is not None
    assert solve_mods.oe_reaction_ids == ("sec_A_complex_formation",)
    assert solve_mods.ko_reaction_ids == ("R_ko",)
    assert captured["sweep_strain_modifications"] is solve_mods
    # The top-N bottleneck reactions (by shadow price) are what gets swept.
    assert captured["sweep_reactions"] == ["sec_A_complex_formation", "sec_B_complex_formation"]
    # dose_response.shapes_by_reaction is keyed by reaction id for C1 to consume.
    dose = result["dose_response"]
    assert set(dose["shapes_by_reaction"]) == {"sec_A_complex_formation", "sec_B_complex_formation"}
    assert result["modified_solve_success"] is True
    assert result["applied_modifications"]["oe_factor"] == 2.0


def test_empty_modifications_solve_with_none_strain_modifications(monkeypatch) -> None:
    captured = _install(monkeypatch, bottlenecks=[_bottleneck("sec_A_complex_formation", 5.0)])
    analyze_next_oe_candidates(target_id="hLF", top_n=1)
    assert captured["solve_strain_modifications"] is None  # wildtype path, byte-identical


def test_unknown_modification_reactions_are_warned(monkeypatch) -> None:
    _install(monkeypatch, bottlenecks=[])
    result = analyze_next_oe_candidates(
        target_id="hLF", oe_reaction_ids=("not_a_reaction",), ko_reaction_ids=("also_missing",)
    )
    warnings = result["modification_warnings"]
    assert any("not_a_reaction" in w for w in warnings)
    assert any("also_missing" in w for w in warnings)


def test_infeasible_modified_solve_yields_no_dose_response(monkeypatch) -> None:
    captured = _install(
        monkeypatch,
        bottlenecks=[_bottleneck("sec_A_complex_formation", 5.0)],
        sim_success=False,
    )
    result = analyze_next_oe_candidates(target_id="hLF", oe_reaction_ids=("sec_A_complex_formation",), top_n=2)
    assert result["modified_solve_success"] is False
    assert result["dose_response"] is None
    assert "sweep_reactions" not in captured  # sweep never called when solve failed


def test_no_bottlenecks_skips_dose_response(monkeypatch) -> None:
    captured = _install(monkeypatch, bottlenecks=[])
    result = analyze_next_oe_candidates(target_id="hLF", top_n=6)
    assert result["dose_response"] is None
    assert "sweep_reactions" not in captured
    assert result["oe_actionable_bottlenecks"] == []
