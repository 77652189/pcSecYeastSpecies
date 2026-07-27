from __future__ import annotations

from app.services.per_strain_layer_reuse import (
    affected_modules,
    bottleneck_modules,
    build_layer_reuse_tags,
    tag_shortlist_reuse,
    to_secretory_module,
)

from pcsec_pichia.screens.gene_perturbation_map import PROCESS_LABELS


def test_to_module_bridges_english_keys_and_display_labels() -> None:
    # 瓶颈侧：classify 英文键
    assert to_secretory_module("disulfide_folding") == "folding"
    assert to_secretory_module("chaperone_folding") == "folding"
    assert to_secretory_module("o_glycan_processing") == "glycosylation"
    assert to_secretory_module("metabolic_or_other") == "metabolic"
    # 候选侧：gene_perturbation_map 展示标签（经 PROCESS_LABELS 逆桥）
    assert to_secretory_module(PROCESS_LABELS["disulfide_folding"]) == "folding"  # "ER 折叠 / DSB"
    assert to_secretory_module(PROCESS_LABELS["erad_misfolding"]) == "folding"  # "错误折叠 / ERAD"
    assert to_secretory_module(PROCESS_LABELS["o_glycan_processing"]) == "glycosylation"
    assert to_secretory_module(PROCESS_LABELS["golgi_surface_transport"]) == "transport"
    assert to_secretory_module(PROCESS_LABELS["metabolic_or_other"]) == "metabolic"
    # 认不出 / 空 → unknown（保守），绝不猜成分泌层
    assert to_secretory_module("") == "unknown"
    assert to_secretory_module(None) == "unknown"
    assert to_secretory_module("完全不认识的标签") == "unknown"


def _result(bottleneck_procs=(), floor_procs=()) -> dict:
    return {
        "oe_actionable_bottlenecks": [{"reaction_id": f"r{i}", "secretory_process": p} for i, p in enumerate(bottleneck_procs)],
        "floor_constraints_not_oe_addressable": [{"reaction_id": f"f{i}", "secretory_process": p} for i, p in enumerate(floor_procs)],
    }


def test_bottleneck_modules_collects_from_bottlenecks_and_floors() -> None:
    result = _result(bottleneck_procs=("disulfide_folding", "metabolic_or_other"), floor_procs=("o_glycan_processing",))
    assert bottleneck_modules(result) == {"folding", "metabolic", "glycosylation"}


def test_affected_modules_union_plus_conservative() -> None:
    wt = _result(bottleneck_procs=("disulfide_folding",))  # 野生型瓶颈=folding
    mod = _result(bottleneck_procs=("disulfide_folding", "o_glycan_processing"))  # 改造后 folding+glyco 顶上
    affected = affected_modules(wt, mod)
    assert "folding" in affected and "glycosylation" in affected
    assert "metabolic" in affected and "unknown" in affected  # 保守模块恒在
    assert "transport" not in affected  # 两状态都不 binding 的分泌层 → 不受影响


def test_tag_marks_uninvolved_secretory_reusable_and_metabolic_conservative() -> None:
    # 瓶颈在 folding；transport 候选无关→可复用；folding 候选→失效；metabolic 候选→保守失效
    wt = _result(bottleneck_procs=("disulfide_folding",))
    mod = _result(bottleneck_procs=("disulfide_folding",))
    affected = affected_modules(wt, mod)
    rows = [
        {"gene_id": "G_transport", "secretory_process": PROCESS_LABELS["golgi_surface_transport"]},
        {"gene_id": "G_folding", "secretory_process": PROCESS_LABELS["disulfide_folding"]},
        {"gene_id": "G_metab", "secretory_process": PROCESS_LABELS["metabolic_or_other"]},
    ]
    tagged = {row["gene_id"]: row for row in tag_shortlist_reuse(rows, affected)}
    assert tagged["G_transport"]["reuse_status"] == "reusable"
    assert tagged["G_transport"]["reuse_module"] == "transport"
    assert tagged["G_folding"]["reuse_status"] == "stale"
    assert tagged["G_metab"]["reuse_status"] == "stale"  # 代谢桶：即便非瓶颈也保守重算


def test_build_layer_reuse_tags_end_to_end_counts_and_caveat() -> None:
    wt = _result(bottleneck_procs=("disulfide_folding",))
    mod = _result(bottleneck_procs=("disulfide_folding", "o_glycan_processing"))
    rows = [
        {"gene_id": "A", "secretory_process": PROCESS_LABELS["golgi_surface_transport"]},  # transport → reusable
        {"gene_id": "B", "secretory_process": PROCESS_LABELS["disulfide_folding"]},  # folding → stale
        {"gene_id": "C", "secretory_process": PROCESS_LABELS["o_glycan_processing"]},  # glyco (改造后瓶颈) → stale
        {"gene_id": "D", "secretory_process": PROCESS_LABELS["metabolic_or_other"]},  # metabolic → stale
    ]
    out = build_layer_reuse_tags(wildtype_result=wt, modified_result=mod, shortlist_rows=rows)
    assert out["reusable_count"] == 1 and out["stale_count"] == 3
    assert out["wildtype_bottleneck_modules"] == ["folding"]
    assert set(out["modified_bottleneck_modules"]) == {"folding", "glycosylation"}
    assert "近似" in out["caveat"] and "重算" in out["caveat"]
