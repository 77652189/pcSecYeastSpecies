"""真实工艺生长速率锚点（机制层，源自在手发酵验证的派生抽象）。

每个锚点是"某目标家族在某工艺相下的真实操作 μ"，供跨条件稳健性 / OE 容量
在真实操作点上求解时**选用**（opt-in）。

保密边界（见 [ADR-006] / [ADR-003]）：原始 OD 曲线、菌株编号、温度/pH 明细为机密，
仅存仓库外本地私有区（`pcSec_wetlab_private/`），**不入 git / 云端**。本模块只保留
可公开的机制层 μ 量级 + 稳健性说明 + provenance，**不含任何菌株 / 位点 / 产量(titer)**。

默认行为不变：默认分泌仿真仍固定 μ=0.10（glucose `corrected_reference`，回归锁定）；
这些锚点是显式选用的操作点，不改变任何默认。
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ProcessGrowthAnchor:
    """真实工艺相的机制层比生长速率锚点。

    仅含可公开字段；不得加入菌株编号、基因位点、产量(titer) 或原始 OD 明细。
    """

    anchor_id: str
    target_family: str  # 产物家族标识（如 "hLF"，已是 repo 内公开目标名）
    carbon_source_id: str  # 对应 media 碳源 id
    process_role: str  # "growth_phase" | "production_phase"
    growth_rate: float  # 机制层 μ 锚点 (h^-1)
    # 机制层标定档：内部一致、量级合理，但未对齐外部实测产量（与 ADR-006 三档一致）。
    calibration_status: str
    robustness_note: str
    provenance: str


_PROCESS_GROWTH_ANCHORS: dict[str, ProcessGrowthAnchor] = {
    "hlf_glycerol_growth": ProcessGrowthAnchor(
        anchor_id="hlf_glycerol_growth",
        target_family="hLF",
        carbon_source_id="glycerol",
        process_role="growth_phase",
        growth_rate=0.10,
        calibration_status="internally_calibrated",
        robustness_note=(
            "甘油生长相：跨温度 20–30℃ / pH 4–7 稳健（全程均值 ~0.10，早期指数峰值 ~0.13，"
            "随甘油耗尽减速）。≈ 模型当前默认固定 μ。"
        ),
        provenance=(
            "在手发酵罐验证派生的机制层 μ 锚点；原始 OD / 菌株 / 温度·pH 明细仅本地私有区，不入 git。"
        ),
    ),
    "hlf_glucose_production": ProcessGrowthAnchor(
        anchor_id="hlf_glucose_production",
        target_family="hLF",
        carbon_source_id="glucose",
        process_role="production_phase",
        growth_rate=0.013,
        calibration_status="internally_calibrated",
        robustness_note=(
            "葡萄糖生产相：跨温度 20–30℃ / pH 4–7 稳健（0.010–0.015）。限量补料生产相，"
            "远低于 glucose μmax；模型中 μ 越低→蛋白预算留给分泌越多。"
        ),
        provenance=(
            "在手发酵罐验证派生的机制层 μ 锚点；原始 OD / 菌株 / 温度·pH 明细仅本地私有区，不入 git。"
        ),
    ),
}


def list_process_growth_anchors() -> tuple[ProcessGrowthAnchor, ...]:
    """列出全部真实工艺生长速率锚点。"""
    return tuple(_PROCESS_GROWTH_ANCHORS.values())


def load_process_growth_anchor(anchor_id: str) -> ProcessGrowthAnchor:
    """按 anchor_id 取锚点；未知 id 抛 KeyError。"""
    return _PROCESS_GROWTH_ANCHORS[anchor_id]


def growth_rate_for(anchor_id: str) -> float:
    """取某工艺操作点的机制层 μ，用于把仿真固定到该真实操作点。"""
    return _PROCESS_GROWTH_ANCHORS[anchor_id].growth_rate


__all__ = [
    "ProcessGrowthAnchor",
    "growth_rate_for",
    "list_process_growth_anchors",
    "load_process_growth_anchor",
]
