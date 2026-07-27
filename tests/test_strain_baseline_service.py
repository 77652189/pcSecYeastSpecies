from __future__ import annotations

import csv
from pathlib import Path

from app.services.strain_baseline_service import (
    ingest_tradeoff_csv_into_baseline_cache,
    load_strain_baseline_readout,
)

from pcsec_pichia.strain_modifications import StrainModifications

_CALIBER = dict(
    target_id="hLF",
    carbon_source_id="glucose",
    media_type=4,
    growth_rate=0.10,
    enable_ribosome_translation_constraint=False,
    enable_misfolding_constraint=False,
)


def _write_tradeoff_csv(path: Path) -> None:
    fields = [
        "target_id", "gene_id", "common_name", "candidate_kind", "intervention_type",
        "secretory_process", "secretion_ratio_vs_wildtype", "growth_retention_ratio",
        "mapping_confidence", "support_status",
    ]
    rows = [
        ["hLF", "G_OE", "PDI1", "gene", "OE", "折叠 (folding)", "1.25", "0.98", "high", "supported"],
        ["hLF", "G_KO", "PEP4", "gene", "KO", "降解 (degradation)", "1.10", "0.90", "medium", "supported"],
        ["OPN_ALPHA_FULL_PROJECT", "G_X", "X", "gene", "OE", "折叠 (folding)", "2.0", "1.0", "high", "supported"],
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)


def test_readout_reports_miss_when_no_baseline_cached(tmp_path) -> None:
    out = load_strain_baseline_readout(cache_dir=tmp_path, **_CALIBER)
    assert out["available"] is False
    assert out["model_variant_fingerprint"] == "wildtype"
    assert out["oe_rows"] == [] and out["ko_rows"] == []
    assert "后台基线" in out["reason"]
    assert out["caliber"]["carbon_source_id"] == "glucose"


def test_ingest_then_readout_hits_and_partitions(tmp_path) -> None:
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    _write_tradeoff_csv(csv_path)

    ingest = ingest_tradeoff_csv_into_baseline_cache(
        csv_path=csv_path, source_run="overnight_full", cache_dir=tmp_path, **_CALIBER
    )
    assert ingest["candidate_count"] == 2  # 只留 hLF 的 KO/OE，丢弃 OPN 行
    assert ingest["model_variant_fingerprint"] == "wildtype"

    out = load_strain_baseline_readout(cache_dir=tmp_path, **_CALIBER)
    assert out["available"] is True
    assert out["source_run"] == "overnight_full"
    assert {r["gene_id"] for r in out["oe_rows"]} == {"G_OE"}
    assert {r["gene_id"] for r in out["ko_rows"]} == {"G_KO"}
    assert out["oe_rows"][0]["secretion_ratio_vs_wildtype"] == 1.25


def test_readout_is_caliber_sensitive(tmp_path) -> None:
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    _write_tradeoff_csv(csv_path)
    ingest_tradeoff_csv_into_baseline_cache(csv_path=csv_path, cache_dir=tmp_path, **_CALIBER)

    # 同口径命中，换 μ 即未命中（不得跨口径误复用）。
    assert load_strain_baseline_readout(cache_dir=tmp_path, **_CALIBER)["available"] is True
    other = dict(_CALIBER, growth_rate=0.20)
    assert load_strain_baseline_readout(cache_dir=tmp_path, **other)["available"] is False


def test_modified_strain_does_not_hit_wildtype_baseline(tmp_path) -> None:
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    _write_tradeoff_csv(csv_path)
    ingest_tradeoff_csv_into_baseline_cache(csv_path=csv_path, cache_dir=tmp_path, **_CALIBER)

    # 野生型基线在缓存里；带改造去读（不同 model_variant_fingerprint）必须未命中，不能误当野生型。
    mods = StrainModifications(oe_reaction_ids=("sec_Pdi1p_complex_formation",))
    out = load_strain_baseline_readout(cache_dir=tmp_path, modifications=mods, **_CALIBER)
    assert out["available"] is False
    assert out["model_variant_fingerprint"].startswith("strain-")
