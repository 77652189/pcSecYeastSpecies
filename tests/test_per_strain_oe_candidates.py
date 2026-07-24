from __future__ import annotations

from app.services.per_strain_oe_candidates import build_next_oe_candidates_readout


def _bottleneck(reaction_id: str, layer: str, abs_marginal: float) -> dict[str, object]:
    return {"reaction_id": reaction_id, "secretory_process": layer, "abs_marginal": abs_marginal, "oe_actionable": True}


def test_empty_bottlenecks_gives_empty_readout() -> None:
    readout = build_next_oe_candidates_readout([])
    assert readout["candidates"] == []
    assert readout["dose_response_available"] is False


def test_without_dose_response_ranks_by_shadow_price() -> None:
    bottlenecks = [
        _bottleneck("R_A", "ER 折叠 / DSB", 100.0),
        _bottleneck("R_C", "Golgi", 30.0),
        _bottleneck("R_B", "Golgi", 50.0),
    ]
    readout = build_next_oe_candidates_readout(bottlenecks)

    assert readout["dose_response_available"] is False
    assert [c["reaction"] for c in readout["candidates"]] == ["R_A", "R_B", "R_C"]  # 影子价格降序
    assert "shape" not in readout["candidates"][0]  # 无剂量响应 → 不带效应/形状


def test_with_dose_response_ranks_by_real_effect_not_shadow_price() -> None:
    # R_A 影子价格最大，但 R_B 的 OE 真实效应更大 —— 说明"binding"≠"松开涨得多"。
    bottlenecks = [
        _bottleneck("R_A", "ER 折叠 / DSB", 100.0),
        _bottleneck("R_B", "Golgi", 50.0),
        _bottleneck("R_C", "Golgi", 30.0),  # 不在剂量响应里
    ]
    dose_response = {
        "baseline_objective": 1.0,
        "shapes_by_reaction": {
            "R_A": {"reaction_id": "R_A", "shape": "saturating", "max_relative_gain": 0.05, "half_gain_factor": 2.0},
            "R_B": {"reaction_id": "R_B", "shape": "linear", "max_relative_gain": 0.09, "half_gain_factor": None},
        },
    }

    readout = build_next_oe_candidates_readout(bottlenecks, dose_response)

    assert readout["dose_response_available"] is True
    # 按真实效应排：R_B(0.09) > R_A(0.05) > R_C(无剂量数据=0)
    assert [c["reaction"] for c in readout["candidates"]] == ["R_B", "R_A", "R_C"]
    by_reaction = {c["reaction"]: c for c in readout["candidates"]}
    assert by_reaction["R_A"]["shape"] == "saturating"
    assert by_reaction["R_B"]["effect"] == 0.09
    assert by_reaction["R_C"].get("shape") is None  # 瓶颈里有但没扫剂量响应 → 无形状


def test_top_n_caps_candidates() -> None:
    bottlenecks = [_bottleneck(f"R_{i}", "层", float(100 - i)) for i in range(10)]
    readout = build_next_oe_candidates_readout(bottlenecks, top_n=3)
    assert len(readout["candidates"]) == 3
    assert [c["reaction"] for c in readout["candidates"]] == ["R_0", "R_1", "R_2"]  # 影子价格 top-3
