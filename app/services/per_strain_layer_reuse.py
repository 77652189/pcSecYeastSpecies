"""D2 受影响层判定 + 复用打标（服务层 · 纯装配 · 零新增求解 · ADR-004 #1 迭代2）。

分层复用的 L1：给一份**同口径**野生型短名单（D4 基线里的 OE+KO 候选），根据"改造后瓶颈相对野生型
瓶颈的层变化"，给每个候选打 `reuse_status`（可复用 / 已失效），供 D3 决定哪些直接复用、哪些按改造后
重算。本模块**不解任何 LP**——只读 C2 已产的瓶颈层 + 候选层，纯字符串装配。

两处硬约束（D4 时坐实，见 EXECUTION_PLAN 阶段② 迭代2）：

1. **分类两套词表且不同**：LP 瓶颈归因用 `classify_secretory_process`（英文键 `disulfide_folding`…），
   筛查短名单用 `gene_perturbation_map` 展示标签（`ER 折叠 / DSB`…）——不能直接字符串比对。故本模块把
   两边都归并到少数几个**粗分泌模块**（folding/glycosylation/transport/…）做同粒度比较，桥是
   `gene_perturbation_map.PROCESS_LABELS`（英文键 ↔ 展示标签的权威映射，英文键集是 classify 的超集）。

2. **只 5 个粗桶、绝大多数候选落 `metabolic_or_other`**：层级复用只对**分泌专属模块**（折叠/糖基化/
   ER 转运/ERAD/Golgi/翻译/降解）干净有效（稀疏、模块化）；**代谢桶 + 未解析一律保守**——改造后视作
   受影响、不复用。目标改造多在分泌层、代谢候选多为 slack 价值低，故这条限制不致命但必须诚实呈现。
"""

from __future__ import annotations

from typing import Any, Sequence

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.screens.gene_perturbation_map import PROCESS_LABELS  # noqa: E402 - engine constant after bootstrap

# 细分类英文键 → 粗分泌模块。两套词表都先归到这里再比（同粒度）。
_ENGLISH_KEY_TO_MODULE: dict[str, str] = {
    "disulfide_folding": "folding",
    "chaperone_folding": "folding",
    "erad_misfolding": "folding",
    "n_glycan_processing": "glycosylation",
    "o_glycan_processing": "glycosylation",
    "er_translocation": "transport",
    "er_to_golgi_transport": "transport",
    "golgi_surface_transport": "transport",
    "ribosome": "translation",
    "proteasome_degradation": "degradation",
    "secretory_capacity": "secretory_capacity",
    "metabolic_or_other": "metabolic",
    "unknown": "unknown",
}
# 分泌专属模块：稀疏、模块化，与瓶颈结构无关时可干净复用野生型结果。
SECRETORY_SPECIFIC_MODULES = frozenset(
    {"folding", "glycosylation", "transport", "translation", "degradation", "secretory_capacity"}
)
# 保守模块：改造后一律视作受影响（复用不可信），不管瓶颈是否明显涉及。
CONSERVATIVE_MODULES = frozenset({"metabolic", "unknown"})

# 展示标签 → 英文键（PROCESS_LABELS 的逆）；候选侧存的是展示标签，先逆查再归模块。
_LABEL_TO_ENGLISH_KEY: dict[str, str] = {label: key for key, label in PROCESS_LABELS.items()}


def to_secretory_module(process: object) -> str:
    """把任一 `secretory_process` 标注（英文 classify 键 或 展示标签）归并到粗分泌模块。

    识别不了 → `"unknown"`（进而被当作保守模块），绝不猜成某个分泌层而误判可复用。
    """
    if not process:
        return "unknown"
    text = str(process).strip()
    if text in _ENGLISH_KEY_TO_MODULE:  # 已是英文键（瓶颈侧）
        return _ENGLISH_KEY_TO_MODULE[text]
    english_key = _LABEL_TO_ENGLISH_KEY.get(text)  # 是展示标签（候选侧）→ 逆查
    if english_key in _ENGLISH_KEY_TO_MODULE:
        return _ENGLISH_KEY_TO_MODULE[english_key]
    return "unknown"


def bottleneck_modules(analysis_result: dict[str, Any]) -> set[str]:
    """一次 C2 `analyze_next_oe_candidates` 结果里，瓶颈 + floor 涉及的粗分泌模块集。"""
    modules: set[str] = set()
    for field in ("oe_actionable_bottlenecks", "floor_constraints_not_oe_addressable"):
        for entry in analysis_result.get(field) or []:
            if isinstance(entry, dict):
                modules.add(to_secretory_module(entry.get("secretory_process")))
    return modules


def affected_modules(wildtype_result: dict[str, Any], modified_result: dict[str, Any]) -> set[str]:
    """受影响的粗分泌模块。

    = 野生型 ∪ 改造后的瓶颈模块（binding 结构涉及的层，改造可能已扰动其耦合）∪ 保守模块（metabolic/
    unknown 一律纳入）。**分泌专属且两状态都不 binding 的模块**才算不受影响 → 其候选可复用野生型结果。
    """
    return bottleneck_modules(wildtype_result) | bottleneck_modules(modified_result) | set(CONSERVATIVE_MODULES)


def tag_shortlist_reuse(
    shortlist_rows: Sequence[dict[str, Any]],
    affected: set[str],
) -> list[dict[str, Any]]:
    """给野生型短名单每个候选打 `reuse_status`（reusable/stale）+ 模块 + 依据。纯装配。"""
    tagged: list[dict[str, Any]] = []
    for row in shortlist_rows:
        module = to_secretory_module(row.get("secretory_process"))
        if module in CONSERVATIVE_MODULES:
            status = "stale"
            reason = "代谢 / 未解析层——改造后复用不可信，保守按改造后重算"
        elif module in affected:
            status = "stale"
            reason = f"所在分泌模块（{module}）改造前后涉及瓶颈结构，需按改造后重算"
        else:
            status = "reusable"
            reason = f"所在分泌模块（{module}）与瓶颈结构无关，野生型结果可复用"
        tagged.append({**dict(row), "reuse_module": module, "reuse_status": status, "reuse_reason": reason})
    return tagged


REUSE_CAVEAT = (
    "分层复用是近似（LP 全局耦合、退化最优有风险）：只对与瓶颈结构无关的**分泌专属层**复用野生型结果；"
    "**代谢 / 未解析层与瓶颈涉及层一律按改造后重算**。相对信号，非绝对产量。"
)


def build_layer_reuse_tags(
    *,
    wildtype_result: dict[str, Any],
    modified_result: dict[str, Any],
    shortlist_rows: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    """D2 顶层装配：受影响模块集 + 给短名单打标 + 计数 + 诚实 caveat。零求解。"""
    wt_modules = bottleneck_modules(wildtype_result)
    mod_modules = bottleneck_modules(modified_result)
    affected = wt_modules | mod_modules | set(CONSERVATIVE_MODULES)
    tagged = tag_shortlist_reuse(shortlist_rows, affected)
    return {
        "affected_modules": sorted(affected),
        "wildtype_bottleneck_modules": sorted(wt_modules),
        "modified_bottleneck_modules": sorted(mod_modules),
        "tagged_candidates": tagged,
        "reusable_count": sum(1 for row in tagged if row["reuse_status"] == "reusable"),
        "stale_count": sum(1 for row in tagged if row["reuse_status"] == "stale"),
        "caveat": REUSE_CAVEAT,
    }


__all__ = [
    "CONSERVATIVE_MODULES",
    "REUSE_CAVEAT",
    "SECRETORY_SPECIFIC_MODULES",
    "affected_modules",
    "bottleneck_modules",
    "build_layer_reuse_tags",
    "tag_shortlist_reuse",
    "to_secretory_module",
]
