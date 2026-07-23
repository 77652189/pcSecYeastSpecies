from __future__ import annotations

import dataclasses

from pcsec_pichia.media import list_carbon_source_formulations
from pcsec_pichia.process_anchors import (
    ProcessGrowthAnchor,
    growth_rate_for,
    list_process_growth_anchors,
    load_process_growth_anchor,
)


def test_process_growth_anchor_fields_stay_mechanism_level_only() -> None:
    """保密护栏：锚点结构只允许机制层字段；不得混入菌株/位点/产量/原始 OD 字段。"""
    field_names = {f.name for f in dataclasses.fields(ProcessGrowthAnchor)}

    assert field_names == {
        "anchor_id",
        "target_family",
        "carbon_source_id",
        "process_role",
        "growth_rate",
        "calibration_status",
        "robustness_note",
        "provenance",
    }
    # 明确禁止的机密字段名不得出现。
    forbidden = {"strain", "strain_id", "genotype", "locus", "titer", "yield", "od", "od_curve", "elisa"}
    assert field_names & forbidden == set()


def test_hlf_growth_rate_anchors_lock_mechanism_level_values() -> None:
    """锁定 hLF 两相 μ 锚点（机制层，源自在手发酵验证）。"""
    growth = load_process_growth_anchor("hlf_glycerol_growth")
    production = load_process_growth_anchor("hlf_glucose_production")

    assert (growth.carbon_source_id, growth.growth_rate) == ("glycerol", 0.10)
    assert (production.carbon_source_id, production.growth_rate) == ("glucose", 0.013)
    assert growth.process_role == "growth_phase"
    assert production.process_role == "production_phase"
    # 未对齐外部实测产量 → 机制层三档停在 internally_calibrated（与 ADR-006 一致）。
    assert growth.calibration_status == "internally_calibrated"
    assert production.calibration_status == "internally_calibrated"

    # 生物学次序：限量补料生产相远慢于甘油生长相。
    assert production.growth_rate < growth.growth_rate

    # 便捷取值一致。
    assert growth_rate_for("hlf_glucose_production") == 0.013


def test_process_anchor_carbon_sources_exist_in_media() -> None:
    """锚点引用的碳源必须是已知 media 碳源，防止拼写漂移。"""
    known = {f.carbon_source_id for f in list_carbon_source_formulations()}
    for anchor in list_process_growth_anchors():
        assert anchor.carbon_source_id in known


def test_process_anchor_provenance_flags_confidential_local_only() -> None:
    """每个锚点都要声明原始数据仅本地私有区、不入 git。"""
    for anchor in list_process_growth_anchors():
        assert "私有区" in anchor.provenance
        assert "不入 git" in anchor.provenance
