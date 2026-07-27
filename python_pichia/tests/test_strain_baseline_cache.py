from __future__ import annotations

from pcsec_pichia.strain_baseline_cache import (
    BASELINE_ROW_FIELDS,
    StrainBaselineCacheKey,
    baseline_to_payload,
    cache_key_digest,
    cache_path,
    cached_baseline_build,
    distill_tradeoff_rows,
    load_cached_baseline,
    store_baseline,
    strain_fingerprint,
)
from pcsec_pichia.strain_modifications import StrainModifications


def _key(**overrides) -> StrainBaselineCacheKey:
    base = dict(
        target_id="hLF",
        carbon_source_id="glucose",
        media_type=4,
        growth_rate=0.10,
        write_ribosome_translation_constraint=False,
        write_misfolding_constraints=False,
    )
    base.update(overrides)
    return StrainBaselineCacheKey(**base)


def _records() -> list[dict[str, object]]:
    return [
        {
            "target_id": "hLF", "gene_id": "G_OE", "common_name": "PDI1", "candidate_kind": "gene",
            "intervention_type": "OE", "secretory_process": "折叠 (folding)",
            "secretion_ratio_vs_wildtype": "1.25", "growth_retention_ratio": "0.98",
            "mapping_confidence": "high", "support_status": "supported",
        },
        {
            "target_id": "hLF", "gene_id": "G_KO", "common_name": "PEP4", "candidate_kind": "gene",
            "intervention_type": "KO", "secretory_process": "降解 (degradation)",
            "secretion_ratio_vs_wildtype": "1.10", "growth_retention_ratio": "0.90",
            "mapping_confidence": "medium", "support_status": "supported",
        },
        {  # 错的 target → 蒸馏时丢弃
            "target_id": "OPN_ALPHA_FULL_PROJECT", "gene_id": "G_OTHER", "intervention_type": "OE",
            "secretion_ratio_vs_wildtype": "9.9", "secretory_process": "折叠 (folding)",
        },
        {  # 非 KO/OE（跳过的候选）→ 蒸馏时丢弃
            "target_id": "hLF", "gene_id": "G_SKIP", "intervention_type": "",
            "secretion_ratio_vs_wildtype": None, "secretory_process": "",
        },
        {  # 数值不可解析 → 安全降级为 None（不炸）
            "target_id": "hLF", "gene_id": "G_BAD", "intervention_type": "OE",
            "secretion_ratio_vs_wildtype": "n/a", "growth_retention_ratio": "",
            "secretory_process": "运输 (transport)",
        },
    ]


def test_strain_fingerprint_wildtype_and_set_semantics() -> None:
    assert strain_fingerprint(None) == "wildtype"
    assert strain_fingerprint(StrainModifications()) == "wildtype"

    mods_a = StrainModifications(ko_reaction_ids=("R_b", "R_a"), oe_reaction_ids=("sec_X_formation",))
    mods_b = StrainModifications(ko_reaction_ids=("R_a", "R_b"), oe_reaction_ids=("sec_X_formation",))
    fp_a = strain_fingerprint(mods_a)
    assert fp_a.startswith("strain-")
    assert fp_a == strain_fingerprint(mods_b)  # 集合语义：KO 顺序不改变指纹
    assert fp_a != "wildtype"

    # oe_factor 是口径的一部分 → 不同剂量必须不同指纹。
    assert strain_fingerprint(StrainModifications(oe_reaction_ids=("sec_X_formation",), oe_factor=2.0)) != \
        strain_fingerprint(StrainModifications(oe_reaction_ids=("sec_X_formation",), oe_factor=4.0))
    # KO 一个反应 vs OE 同名反应，语义不同 → 指纹不同。
    assert strain_fingerprint(StrainModifications(ko_reaction_ids=("R_a",))) != \
        strain_fingerprint(StrainModifications(oe_reaction_ids=("R_a",)))


def test_cache_key_digest_is_deterministic_and_field_sensitive() -> None:
    base = _key()
    assert cache_key_digest(base) == cache_key_digest(_key())  # 确定性

    mutations = (
        _key(target_id="OPN_ALPHA_FULL_PROJECT"),
        _key(carbon_source_id="methanol"),
        _key(media_type=2),
        _key(growth_rate=0.013),
        _key(write_ribosome_translation_constraint=True),
        _key(write_misfolding_constraints=True),
        _key(mode="precise"),
        _key(compatibility_mode="reported"),
        _key(solver_method="highs-ipm"),
        _key(model_variant_fingerprint="strain-abc123"),  # 改造变体不得与野生型同键
        _key(schema_version="old-unknown-classification"),  # 旧分类基线永不命中
    )
    base_digest = cache_key_digest(base)
    for mutated in mutations:
        assert cache_key_digest(mutated) != base_digest


def test_distill_tradeoff_rows_filters_by_target_and_type_and_coerces() -> None:
    rows = distill_tradeoff_rows(_records(), target_id="hLF")
    gene_ids = {row["gene_id"] for row in rows}
    assert gene_ids == {"G_OE", "G_KO", "G_BAD"}  # 丢弃错 target 与非 KO/OE

    by_gene = {row["gene_id"]: row for row in rows}
    assert by_gene["G_OE"]["secretion_ratio_vs_wildtype"] == 1.25  # 数值转 float
    assert by_gene["G_OE"]["secretory_process"] == "折叠 (folding)"
    assert by_gene["G_BAD"]["secretion_ratio_vs_wildtype"] is None  # 不可解析安全降级
    assert by_gene["G_BAD"]["growth_retention_ratio"] is None
    # 只保留约定字段，不夹带筛查表其它列。
    assert set(by_gene["G_KO"].keys()) == set(BASELINE_ROW_FIELDS)


def test_store_and_load_round_trip(tmp_path) -> None:
    key = _key()
    rows = distill_tradeoff_rows(_records(), target_id="hLF")
    store_baseline(key, rows, cache_dir=tmp_path, source_run="overnight_full")

    loaded = load_cached_baseline(key, cache_dir=tmp_path)
    assert loaded is not None
    assert loaded["candidate_count"] == 3
    assert loaded["source_run"] == "overnight_full"
    assert loaded["cache_key_digest"] == cache_key_digest(key)
    assert {row["gene_id"] for row in loaded["rows"]} == {"G_OE", "G_KO", "G_BAD"}

    # 口径不同（μ 换了）→ 未命中（不得复用别的口径的基线）。
    assert load_cached_baseline(_key(growth_rate=0.20), cache_dir=tmp_path) is None


def test_schema_version_bump_isolates_old_unknown_classification_cache(tmp_path) -> None:
    # 旧 unknown 分类基线写在旧 schema_version 键下；当前默认键必须读不到它。
    old_key = _key(schema_version="old-unknown-classification")
    store_baseline(old_key, distill_tradeoff_rows(_records(), "hLF"), cache_dir=tmp_path)
    assert load_cached_baseline(_key(), cache_dir=tmp_path) is None  # 硬保证：不碰旧缓存


def test_digest_mismatch_is_treated_as_miss(tmp_path) -> None:
    key = _key()
    path = cache_path(key, tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('{"cache_key_digest": "deadbeef", "rows": []}', encoding="utf-8")
    assert load_cached_baseline(key, cache_dir=tmp_path) is None


def test_cached_baseline_build_miss_then_hit_and_force(tmp_path) -> None:
    key = _key()
    calls: list[int] = []

    def compute() -> list[dict[str, object]]:
        calls.append(1)
        return distill_tradeoff_rows(_records(), "hLF")

    payload_1, from_cache_1 = cached_baseline_build(key, compute, cache_dir=tmp_path)
    assert from_cache_1 is False and len(calls) == 1
    assert cache_path(key, tmp_path).exists()
    assert payload_1["candidate_count"] == 3

    payload_2, from_cache_2 = cached_baseline_build(key, compute, cache_dir=tmp_path)
    assert from_cache_2 is True and len(calls) == 1  # 命中不得再算
    assert payload_2["cache_key_digest"] == payload_1["cache_key_digest"]

    _, from_cache_3 = cached_baseline_build(key, compute, cache_dir=tmp_path, force=True)
    assert from_cache_3 is False and len(calls) == 2  # force 必须重算


def test_baseline_to_payload_stamps_built_at_and_count() -> None:
    key = _key()
    payload = baseline_to_payload(key, [{"gene_id": "G_OE"}], built_at="2026-07-24T00:00:00+00:00")
    assert payload["built_at"] == "2026-07-24T00:00:00+00:00"
    assert payload["candidate_count"] == 1
    assert payload["cache_key"]["target_id"] == "hLF"
