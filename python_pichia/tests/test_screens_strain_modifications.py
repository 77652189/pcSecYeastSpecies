from __future__ import annotations

from types import SimpleNamespace

import pcsec_pichia.screens as screens
from pcsec_pichia.screens import ScreenResult, run_knockout_screen
from pcsec_pichia.strain_modifications import StrainModifications


def _prepared(success: bool) -> dict[str, object]:
    return {
        "baseline_success": success,
        "baseline": SimpleNamespace(success=success, objective_value=(1.0 if success else None)),
        "constraint_counts": {},
        "fixed_model": "WT_MODEL",
        "secretory": "WT_SEC",
        "combined": "WT_COMB",
        "exchange_reaction_id": "r_target_exchange",
    }


def test_apply_modifications_and_rebaseline_overrides_model_and_resolves_new_baseline(monkeypatch) -> None:
    seen: dict[str, object] = {}

    def fake_apply(model, sec, comb, mods):
        seen["apply_args"] = (model, sec, comb, mods)
        return "MOD_MODEL", "MOD_SEC", "MOD_COMB", ({"reaction_id": "x", "kind": "OE"},), ()

    def fake_solve(model, exchange, **kwargs):
        seen["solve_model"] = model
        seen["solve_sec"] = kwargs.get("secretory")
        seen["solve_comb"] = kwargs.get("combined")
        return SimpleNamespace(success=True, objective_value=0.7), {"eq_total": 1}

    monkeypatch.setattr(screens, "apply_strain_modifications", fake_apply)
    monkeypatch.setattr(screens, "solve_pcsec_maximize", fake_solve)

    prepared = _prepared(True)
    mods = StrainModifications(oe_reaction_ids=("x",))
    out = screens._apply_modifications_and_rebaseline(
        prepared, mods, metabolic=object(), growth_rate=0.1,
        write_ribosome_translation_constraint=False, write_misfolding_constraints=False,
    )

    # 改造后的 model/enzyme data 覆盖进 prepared，baseline 重解到改造后的解。
    assert out["fixed_model"] == "MOD_MODEL"
    assert out["secretory"] == "MOD_SEC"
    assert out["combined"] == "MOD_COMB"
    assert out["baseline"].objective_value == 0.7
    assert out["baseline_success"] is True
    # 基线重解用的是改造后的 model + enzyme data，不是野生型。
    assert seen["solve_model"] == "MOD_MODEL"
    assert seen["solve_sec"] == "MOD_SEC" and seen["solve_comb"] == "MOD_COMB"
    # 返回副本，不就地改原 prepared。
    assert prepared["fixed_model"] == "WT_MODEL"


def test_infeasible_rebaseline_marks_baseline_success_false(monkeypatch) -> None:
    monkeypatch.setattr(
        screens, "apply_strain_modifications",
        lambda model, sec, comb, mods: ("MOD_MODEL", "MOD_SEC", "MOD_COMB", (), ()),
    )
    monkeypatch.setattr(
        screens, "solve_pcsec_maximize",
        lambda model, exchange, **kwargs: (SimpleNamespace(success=False, objective_value=None), {}),
    )
    out = screens._apply_modifications_and_rebaseline(
        _prepared(True), StrainModifications(ko_reaction_ids=("R_ko",)), metabolic=object(),
        growth_rate=0.1, write_ribosome_translation_constraint=False, write_misfolding_constraints=False,
    )
    assert out["baseline_success"] is False


def test_run_knockout_screen_without_modifications_never_rebaselines(monkeypatch) -> None:
    calls: list[int] = []
    monkeypatch.setattr(screens, "_prepare_screen_inputs", lambda *a, **k: _prepared(False))
    monkeypatch.setattr(
        screens, "_apply_modifications_and_rebaseline",
        lambda *a, **k: (calls.append(1), _prepared(True))[1],
    )
    result = run_knockout_screen(
        model="m", target=SimpleNamespace(target_id="hLF"), amino_acids="aa",
        metabolic="met", secretory="S", combined="C", genes=["G1"],
    )
    assert isinstance(result, ScreenResult)
    assert calls == []  # 无改造 → 走野生型路径，绝不重解改造后基线


def test_run_knockout_screen_with_modifications_reroutes_through_rebaseline(monkeypatch) -> None:
    calls: list[object] = []

    def spy_rebaseline(prepared, modifications, *a, **k):
        calls.append(modifications)
        return _prepared(False)  # 改造后不可行 → 优雅早返回，不碰后半段求解

    monkeypatch.setattr(screens, "_prepare_screen_inputs", lambda *a, **k: _prepared(True))
    monkeypatch.setattr(screens, "_apply_modifications_and_rebaseline", spy_rebaseline)

    mods = StrainModifications(oe_reaction_ids=("sec_Pdi1p_complex_formation",))
    result = run_knockout_screen(
        model="m", target=SimpleNamespace(target_id="hLF"), amino_acids="aa",
        metabolic="met", secretory="S", combined="C", genes=["G1"],
        strain_modifications=mods,
    )
    assert isinstance(result, ScreenResult)
    assert len(calls) == 1 and calls[0] is mods  # 改造经由 rebaseline 应用到 KO 基线
