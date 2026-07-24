from __future__ import annotations

import json

import pandas as pd

from app.services.genome_wide_screen_shortlist import build_shortlist_readout


def _oe(
    gene_id: str,
    ratio: float | None,
    *,
    common_name: str = "",
    candidate_kind: str = "gene",
    secretory_process: str = "ER 折叠 / DSB",
    growth_retention: float = 1.0,
    mapping_confidence: str = "curated",
    target_id: str = "hLF",
) -> dict[str, object]:
    return {
        "target_id": target_id,
        "gene_id": gene_id,
        "common_name": common_name,
        "candidate_kind": candidate_kind,
        "intervention_type": "OE",
        "secretion_ratio_vs_wildtype": ratio,
        "growth_retention_ratio": growth_retention,
        "secretory_process": secretory_process,
        "mapping_confidence": mapping_confidence,
    }


def _ko(gene_id: str, ratio: float | None, *, target_id: str = "hLF") -> dict[str, object]:
    return {
        "target_id": target_id,
        "gene_id": gene_id,
        "common_name": "",
        "candidate_kind": "gene",
        "intervention_type": "KO",
        "secretion_ratio_vs_wildtype": ratio,
        "growth_retention_ratio": 1.0,
        "secretory_process": "代谢 / 其它",
        "mapping_confidence": "curated",
    }


def test_shortlist_ranked_by_relative_improvement_and_excludes_ko_and_non_wins() -> None:
    frame = pd.DataFrame(
        [
            _oe("R_A", 1.08, common_name="PDI1/ERO1"),
            _oe("R_B", 1.02, common_name="OCH1"),
            _oe("R_C", 1.05, common_name="AP-1"),
            _oe("R_D", 1.00, common_name="FLAT"),  # no improvement -> excluded
            _oe("R_E", 0.90, common_name="WORSE"),  # worse than wildtype -> excluded
            _ko("R_KO", 1.20),  # KO -> excluded (this is an OE shortlist)
        ]
    )

    readout = build_shortlist_readout(frame, "hLF")

    names = [row["candidate"] for row in readout["oe_shortlist"]]
    assert names == ["PDI1/ERO1", "AP-1", "OCH1"]
    assert readout["has_strong_oe_lever"] is True
    assert abs(float(readout["top_effect"]) - 0.08) < 1e-9
    assert readout["oe_shortlist"][0]["layer"] == "ER 折叠 / DSB"
    assert readout["oe_shortlist"][0]["confidence"] == "curated"


def test_shortlist_excludes_complex_oe_hypothesis_guesses() -> None:
    frame = pd.DataFrame(
        [
            _oe("R_REAL", 1.05, common_name="REAL_OE"),
            # An untested "whole complex overexpressed at the same ratio" guess with a big ratio
            # must NOT be presented as an ordinary OE win topping the shortlist.
            _oe("R_GUESS", 1.30, common_name="WHOLE_COMPLEX_GUESS", candidate_kind="complex_oe_hypothesis"),
        ]
    )

    readout = build_shortlist_readout(frame, "hLF")

    assert [row["candidate"] for row in readout["oe_shortlist"]] == ["REAL_OE"]


def test_candidate_falls_back_to_gene_id_when_no_common_name() -> None:
    frame = pd.DataFrame([_oe("R_NONAME", 1.03, common_name="")])
    readout = build_shortlist_readout(frame, "hLF")
    assert readout["oe_shortlist"][0]["candidate"] == "R_NONAME"


def test_value_of_information_flags_near_tie_among_top_but_not_separated_leader() -> None:
    frame = pd.DataFrame(
        [
            _oe("R_LEAD", 1.10, common_name="LEADER"),
            _oe("R_TIE1", 1.021, common_name="TIE_A"),
            _oe("R_TIE2", 1.020, common_name="TIE_B"),
        ]
    )

    voi = build_shortlist_readout(frame, "hLF")["value_of_information"]

    assert [row["candidate_id"] for row in voi["ranked_candidates"]] == ["LEADER", "TIE_A", "TIE_B"]
    assert voi["has_actionable_ambiguity"] is True
    tied = {candidate for item in voi["information_items"] for candidate in item["candidates"]}
    assert tied == {"TIE_A", "TIE_B"}  # only the near-tied pair; the separated leader is not flagged


def test_value_of_information_no_ambiguity_when_scores_separated() -> None:
    frame = pd.DataFrame(
        [
            _oe("R1", 1.10, common_name="A"),
            _oe("R2", 1.05, common_name="B"),
            _oe("R3", 1.02, common_name="C"),
        ]
    )
    voi = build_shortlist_readout(frame, "hLF")["value_of_information"]
    assert voi["has_actionable_ambiguity"] is False
    assert voi["information_items"] == []


def test_no_oe_up_candidates_gives_empty_shortlist_and_no_voi() -> None:
    frame = pd.DataFrame([_ko("R_KO", 1.2), _oe("R_FLAT", 0.99, common_name="X")])
    readout = build_shortlist_readout(frame, "hLF")
    assert readout["oe_shortlist"] == []
    assert readout["has_strong_oe_lever"] is False
    assert readout["value_of_information"] == {}


def test_r1_floors_absent_degrades_gracefully_without_dir() -> None:
    frame = pd.DataFrame([_oe("R1", 1.05, common_name="A")])
    readout = build_shortlist_readout(frame, "hLF", r1_readout_dir=None)
    assert readout["why_limited_floors"] == []
    assert readout["r1_available"] is False


def test_r1_floors_loaded_and_sorted_from_cached_json(tmp_path) -> None:
    (tmp_path / "target_bottleneck_lp_attribution_hLF.json").write_text(
        json.dumps(
            {
                "floor_constraints_not_oe_addressable": [
                    {"reaction_id": "weak_floor", "abs_marginal": 10.0},
                    {"reaction_id": "strong_floor", "abs_marginal": 5000.0},
                    {"reaction_id": "mid_floor", "abs_marginal": 200.0},
                ]
            }
        ),
        encoding="utf-8",
    )
    frame = pd.DataFrame([_oe("R1", 1.05, common_name="A")])

    readout = build_shortlist_readout(frame, "hLF", r1_readout_dir=tmp_path)

    ids = [f["reaction_id"] for f in readout["why_limited_floors"]]
    assert ids == ["strong_floor", "mid_floor", "weak_floor"]  # sorted by abs_marginal desc
    assert readout["r1_available"] is True


def test_growth_risky_candidates_flagged() -> None:
    frame = pd.DataFrame(
        [
            _oe("R_SAFE", 1.10, common_name="SAFE", growth_retention=1.0),
            _oe("R_RISK", 1.05, common_name="RISK", growth_retention=0.5),
        ]
    )
    readout = build_shortlist_readout(frame, "hLF")
    assert readout["growth_risky_candidates"] == ["RISK"]


def test_dose_response_absent_degrades_gracefully() -> None:
    frame = pd.DataFrame([_oe("R_A", 1.08, common_name="PDI1/ERO1")])
    readout = build_shortlist_readout(frame, "hLF", dose_response_dir=None)
    assert readout["dose_response_available"] is False
    assert readout["dose_response"] == {}
    assert "shape" not in readout["oe_shortlist"][0]  # no shape attached when no cache


def test_dose_response_shapes_attached_from_cache_by_reaction_id(tmp_path) -> None:
    (tmp_path / "hLF_dose_response.json").write_text(
        json.dumps(
            {
                "target_id": "hLF",
                "tested_factors": [1.25, 1.5, 2.0, 3.0, 5.0, 8.0],
                "baseline_objective": 1.0,
                "shapes_by_reaction": {
                    "R_A": {"reaction_id": "R_A", "shape": "saturating", "max_relative_gain": 0.08, "half_gain_factor": 2.0},
                    "R_B": {"reaction_id": "R_B", "shape": "linear", "max_relative_gain": 0.05, "half_gain_factor": None},
                },
                "warnings": ["w"],
            }
        ),
        encoding="utf-8",
    )
    frame = pd.DataFrame(
        [
            _oe("R_A", 1.08, common_name="PDI1/ERO1"),
            _oe("R_B", 1.05, common_name="AP-1"),
            _oe("R_C", 1.02, common_name="OCH1"),  # not in cache -> no shape
        ]
    )

    readout = build_shortlist_readout(frame, "hLF", dose_response_dir=tmp_path)

    assert readout["dose_response_available"] is True
    by_name = {row["candidate"]: row for row in readout["oe_shortlist"]}
    assert by_name["PDI1/ERO1"]["shape"] == "saturating"
    assert by_name["PDI1/ERO1"]["shape_half_gain_factor"] == 2.0
    assert by_name["AP-1"]["shape"] == "linear"
    assert by_name["OCH1"].get("shape") is None  # reaction absent from cache stays unshaped
    assert readout["dose_response"]["tested_factors"] == [1.25, 1.5, 2.0, 3.0, 5.0, 8.0]


def test_condition_matrix_absent_degrades_gracefully() -> None:
    frame = pd.DataFrame([_oe("R_A", 1.08, common_name="PDI1/ERO1")])
    readout = build_shortlist_readout(frame, "hLF", condition_matrix_dir=None)
    assert readout["condition_matrix_available"] is False
    assert readout["condition_matrix"] == {}
    assert "cross_condition_robustness" not in readout["oe_shortlist"][0]


def test_condition_matrix_robustness_attached_from_cache_by_reaction_id(tmp_path) -> None:
    (tmp_path / "hLF_condition_matrix.json").write_text(
        json.dumps(
            {
                "target_id": "hLF",
                "mu": 0.1,
                "conditions": ["glycerol", "glucose"],
                "per_reaction_across_conditions": {
                    # 两条件形状一致 -> 跨条件稳健
                    "R_A": {"glycerol": {"reaction_id": "R_A", "shape": "saturating"}, "glucose": {"shape": "saturating"}},
                    # 两条件形状不同 -> 条件敏感
                    "R_B": {"glycerol": {"shape": "saturating"}, "glucose": {"shape": "linear"}},
                    # 只有一个条件有数据 -> 仅单条件
                    "R_C": {"glucose": {"shape": "saturating"}},
                },
            }
        ),
        encoding="utf-8",
    )
    frame = pd.DataFrame(
        [
            _oe("R_A", 1.08, common_name="PDI1/ERO1"),
            _oe("R_B", 1.05, common_name="AP-1"),
            _oe("R_C", 1.03, common_name="OCH1"),
            _oe("R_D", 1.02, common_name="NOTINMATRIX"),  # 不在矩阵里 -> 不附标注
        ]
    )

    readout = build_shortlist_readout(frame, "hLF", condition_matrix_dir=tmp_path)

    assert readout["condition_matrix_available"] is True
    assert readout["condition_matrix"]["conditions"] == ["glycerol", "glucose"]
    by_name = {row["candidate"]: row for row in readout["oe_shortlist"]}
    assert by_name["PDI1/ERO1"]["cross_condition_robustness"] == "cross_condition_stable"
    assert by_name["AP-1"]["cross_condition_robustness"] == "cross_condition_sensitive"
    assert by_name["OCH1"]["cross_condition_robustness"] == "cross_condition_single"
    assert "cross_condition_robustness" not in by_name["NOTINMATRIX"]
