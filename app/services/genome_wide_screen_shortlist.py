"""① 候选短名单读出（服务层，复用筛查缓存，零新增 LP 求解）。

把 genome_wide 筛查已算好并缓存的 tradeoff CSV（逐候选 secretion_ratio_vs_wildtype，
已带中文 secretory_process）+ 缓存的 R1 LP 瓶颈读出（per-target 影子价格）合成给
研发/湿实验看的一份读出：

    为什么受限（R1）  +  OE 提升候选短名单（筛查）  +  该测什么（R4 价值-of-information）

效率关键：逐候选的 LP 求解在筛查阶段一次性做完并缓存，本模块**零新增求解**——
R4/排序只是对缓存排序的纯后处理，R1 复用缓存读出。剂量响应形状（R2）需要对短名单
补扫倍数（有界额外求解），本模块不做。

分层：UI（app/ui）禁止 import 引擎，故 R4 的 prioritize/summarize_value_of_information
放在这一服务层调用、返回纯 dict，再由 genome_wide_screen 视图渲染。CLI 版见
python_pichia/tools/build_candidate_shortlist_readout.py（两者逻辑保持一致）。
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.analysis import (  # noqa: E402 - engine import after path bootstrap (services are allowlisted)
    prioritize_value_of_information,
    summarize_value_of_information,
)

# OE 相对提升低于此阈值视作“无实质提升”（与 R2 剂量响应 1e-3 噪声底同量级）。
STRONG_EFFECT_THRESHOLD = 0.01  # 1% 相对提升
GROWTH_RISK_THRESHOLD = 0.9  # 生长保持率低于此值提示有生长代价


def _load_r1_floors(r1_dir: Path | None, target_id: str, top: int = 5) -> list[dict[str, object]]:
    """从缓存 R1 读出取最强的下界 floor（为什么受限）。目录/文件缺失或损坏则返回空。"""
    if r1_dir is None:
        return []
    path = Path(r1_dir) / f"target_bottleneck_lp_attribution_{target_id}.json"
    if not path.exists():
        return []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return []
    floors = data.get("floor_constraints_not_oe_addressable") or []
    rows = [
        {"reaction_id": str(f.get("reaction_id")), "abs_marginal": float(f.get("abs_marginal") or 0.0)}
        for f in floors
        if isinstance(f, dict)
    ]
    rows.sort(key=lambda r: r["abs_marginal"], reverse=True)
    return rows[:top]


def _oe_shortlist(frame: pd.DataFrame, target_id: str, top_n: int) -> list[dict[str, object]]:
    """target 的 OE 提升候选（ratio>1），按相对提升降序，取 top_n。

    剔除 candidate_kind == "complex_oe_hypothesis"：那是“把整个复合体按同一比例整体过表达”的
    未验证猜测（见 analyze_single_target 同样的剔除），混进普通 OE 赢家里会误导。
    """
    sub = frame[(frame["target_id"] == target_id) & (frame["intervention_type"] == "OE")].copy()
    if "candidate_kind" in sub.columns:
        sub = sub[sub["candidate_kind"] != "complex_oe_hypothesis"]
    sub = sub.dropna(subset=["secretion_ratio_vs_wildtype"])
    sub["effect"] = sub["secretion_ratio_vs_wildtype"].astype(float) - 1.0
    sub = sub[sub["effect"] > 0].sort_values("effect", ascending=False).head(top_n)
    rows: list[dict[str, object]] = []
    for _, r in sub.iterrows():
        rows.append(
            {
                "candidate": str(r.get("common_name") or r.get("gene_id")),
                "reaction": str(r.get("gene_id")),
                "layer": str(r.get("secretory_process") or "未解析"),
                "effect": float(r["effect"]),
                "growth_retention": float(r.get("growth_retention_ratio") or 1.0),
                "confidence": str(r.get("mapping_confidence") or ""),
            }
        )
    return rows


def build_shortlist_readout(
    frame: pd.DataFrame,
    target_id: str,
    *,
    top_n: int = 8,
    r1_readout_dir: Path | None = None,
) -> dict[str, object]:
    """合成一个 target 的短名单读出 dict（零新增求解）。

    frame: 已经过 load_gene_tradeoff_csv 归一化的筛查 tradeoff 表。
    r1_readout_dir: 缓存 R1 读出目录（None 或缺文件时“为什么受限”段落优雅降级为空）。
    """
    floors = _load_r1_floors(r1_readout_dir, target_id)
    shortlist = _oe_shortlist(frame, target_id, top_n)

    # R4 价值-of-information：对 OE 短名单排序（分数=相对提升）检近似并列 → 该测什么。零求解。
    # top_k 设为短名单长度，让它扫遍整个短名单的相邻近似并列（不只 top-3）。
    if shortlist:
        voi = summarize_value_of_information(
            prioritize_value_of_information(
                target_id,
                [(row["candidate"], row["effect"]) for row in shortlist],
                top_k=max(2, len(shortlist)),
            )
        )
    else:
        voi = {}

    top_effect = float(shortlist[0]["effect"]) if shortlist else 0.0
    return {
        "target_id": target_id,
        "why_limited_floors": floors,
        "r1_available": bool(floors),
        "oe_shortlist": shortlist,
        "has_strong_oe_lever": top_effect >= STRONG_EFFECT_THRESHOLD,
        "top_effect": top_effect,
        "growth_risky_candidates": [
            row["candidate"] for row in shortlist if float(row["growth_retention"]) < GROWTH_RISK_THRESHOLD
        ],
        "value_of_information": voi,
    }


__all__ = [
    "GROWTH_RISK_THRESHOLD",
    "STRONG_EFFECT_THRESHOLD",
    "build_shortlist_readout",
]
