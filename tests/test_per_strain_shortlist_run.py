from __future__ import annotations

import csv
from pathlib import Path

from app.services.per_strain_shortlist_run import build_modified_strain_shortlist
from app.services.strain_baseline_service import ingest_tradeoff_csv_into_baseline_cache

from pcsec_pichia.screens.gene_perturbation_map import PROCESS_LABELS

_CALIBER = dict(
    target_id="hLF", carbon_source_id="glucose", media_type=4, mu=0.10,
    enable_ribosome_translation_constraint=False, enable_misfolding_constraint=False,
)


def _write_baseline_csv(path: Path) -> None:
    fields = [
        "target_id", "gene_id", "common_name", "candidate_kind", "intervention_type",
        "secretory_process", "affected_reactions", "secretion_ratio_vs_wildtype",
        "growth_retention_ratio", "mapping_confidence", "support_status",
    ]
    fold = PROCESS_LABELS["disulfide_folding"]
    trans = PROCESS_LABELS["golgi_surface_transport"]
    metab = PROCESS_LABELS["metabolic_or_other"]
    rows = [
        ["hLF", "G_oe_fold", "PDI1", "gene", "OE", fold, "r_pdi", "1.30", "0.98", "high", "supported"],
        ["hLF", "G_oe_trans", "SEC4", "gene", "OE", trans, "r_sec4", "1.20", "0.99", "high", "supported"],
        ["hLF", "G_oe_metab", "PGK1", "gene", "OE", metab, "r_pgk", "1.10", "1.0", "low", "supported"],
        ["hLF", "G_ko_fold", "HAC1", "gene", "KO", fold, "r_hac", "1.25", "0.90", "medium", "supported"],
        ["hLF", "G_ko_trans", "VPS10", "gene", "KO", trans, "r_vps", "1.15", "0.95", "medium", "supported"],
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)


def _ingest_baseline(tmp_path: Path) -> None:
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    _write_baseline_csv(csv_path)
    ingest_tradeoff_csv_into_baseline_cache(
        csv_path=csv_path, source_run="test_baseline", cache_dir=tmp_path,
        target_id="hLF", carbon_source_id="glucose", media_type=4, growth_rate=0.10,
        enable_ribosome_translation_constraint=False, enable_misfolding_constraint=False,
    )


def _fake_analyze(*, ko_reaction_ids, oe_reaction_ids, **kwargs) -> dict:
    # 野生型与改造后瓶颈都落 folding（folding 是耦合层）；transport/glyco 不 binding。
    return {
        "oe_actionable_bottlenecks": [{"reaction_id": "sec_Pdi1p", "secretory_process": "disulfide_folding"}],
        "floor_constraints_not_oe_addressable": [],
        "modified_solve_success": True,
        "applied_modifications": {
            "ko_reaction_ids": list(ko_reaction_ids), "oe_reaction_ids": list(oe_reaction_ids), "oe_factor": 2.0,
        },
    }


def test_reports_needs_baseline_when_cache_missing(tmp_path) -> None:
    out = build_modified_strain_shortlist(
        oe_reaction_ids=("sec_Pdi1p_complex_formation",), cache_dir=tmp_path, _analyze=_fake_analyze, **_CALIBER
    )
    assert out["available"] is False
    assert out["needs_baseline_build"] is True
    assert out["oe_candidates"] == [] and out["ko_candidates"] == []


def test_l1_tags_shortlist_reuse_vs_stale(tmp_path) -> None:
    _ingest_baseline(tmp_path)
    out = build_modified_strain_shortlist(
        oe_reaction_ids=("sec_Pdi1p_complex_formation",), cache_dir=tmp_path, _analyze=_fake_analyze, **_CALIBER
    )
    assert out["available"] is True and out["layer"] == "L1"
    assert out["modified_solve_success"] is True

    oe = {c["gene_id"]: c for c in out["oe_candidates"]}
    ko = {c["gene_id"]: c for c in out["ko_candidates"]}
    # 分泌专属层且与瓶颈无关（transport）→ 可复用；瓶颈层（folding）→ 失效；代谢桶 → 保守失效。
    assert oe["G_oe_trans"]["reuse_status"] == "reusable"
    assert oe["G_oe_fold"]["reuse_status"] == "stale"
    assert oe["G_oe_metab"]["reuse_status"] == "stale"
    assert ko["G_ko_trans"]["reuse_status"] == "reusable"
    assert ko["G_ko_fold"]["reuse_status"] == "stale"

    assert out["oe_reusable_count"] == 1 and out["oe_stale_count"] == 2
    assert out["ko_reusable_count"] == 1 and out["ko_stale_count"] == 1
    assert "folding" in out["affected_modules"] and "metabolic" in out["affected_modules"]
    # 候选保留野生型相对效应（复用值），按效应降序。
    assert oe["G_oe_fold"]["wildtype_effect"] > oe["G_oe_metab"]["wildtype_effect"]
    assert out["oe_candidates"][0]["wildtype_effect"] >= out["oe_candidates"][-1]["wildtype_effect"]


def test_shortlist_drops_non_improving_and_respects_top_n(tmp_path) -> None:
    # 加一行 ratio<=1（无提升）应被剔除；top_n=1 只留最强。
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    _write_baseline_csv(csv_path)
    with csv_path.open("a", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerow(
            ["hLF", "G_oe_flat", "NEUT", "gene", "OE", PROCESS_LABELS["o_glycan_processing"], "r_x", "0.98", "1.0", "low", "supported"]
        )
    ingest_tradeoff_csv_into_baseline_cache(
        csv_path=csv_path, cache_dir=tmp_path, target_id="hLF", carbon_source_id="glucose",
        media_type=4, growth_rate=0.10, enable_ribosome_translation_constraint=False,
        enable_misfolding_constraint=False,
    )
    out = build_modified_strain_shortlist(
        oe_reaction_ids=("sec_Pdi1p_complex_formation",), top_n=1, cache_dir=tmp_path,
        _analyze=_fake_analyze, **_CALIBER
    )
    assert len(out["oe_candidates"]) == 1
    assert out["oe_candidates"][0]["gene_id"] == "G_oe_fold"  # 效应最大（1.30-1）
    assert all(c["gene_id"] != "G_oe_flat" for c in out["oe_candidates"])  # ratio<=1 剔除
