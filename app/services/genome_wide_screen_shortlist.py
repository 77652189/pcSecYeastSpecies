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
from app.services.gene_name_annotation import safe_gene_display_label

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
        {
            "reaction_id": str(f.get("reaction_id")),
            "abs_marginal": float(f.get("abs_marginal") or 0.0),
            # 分泌层（如 disulfide_folding）——面板拿它显示研究员看得懂的"卡在哪一层"，而非原始反应 id。
            "secretory_process": str(f.get("secretory_process") or ""),
        }
        for f in floors
        if isinstance(f, dict)
    ]
    rows.sort(key=lambda r: r["abs_marginal"], reverse=True)
    return rows[:top]


def _load_dose_response(dose_response_dir: Path | None, target_id: str) -> dict[str, object]:
    """读离线扫描缓存 `{target}_dose_response.json`（reaction_id→形状）。缺失/损坏返回空。

    形状本身是模型+反应的属性、与哪次筛查无关，故按 (target, reaction_id) 缓存、跨 run 可复用。
    生成见 python_pichia/tools/run_shortlist_dose_response.py（有界额外求解，离线/后台跑）。
    """
    if dose_response_dir is None:
        return {}
    path = Path(dose_response_dir) / f"{target_id}_dose_response.json"
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _load_condition_matrix(condition_matrix_dir: Path | None, target_id: str) -> dict[str, object]:
    """读离线 B2 条件矩阵缓存 `{target}_condition_matrix.json`（reaction→各碳源条件形状）。缺失/损坏返回空。

    生成见 python_pichia/tools/run_shortlist_condition_matrix.py（离线/后台跑；干净单碳源、μ=0.10）。
    """
    if condition_matrix_dir is None:
        return {}
    path = Path(condition_matrix_dir) / f"{target_id}_condition_matrix.json"
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _cross_condition_robustness(condition_shapes: dict[str, object]) -> str:
    """按各碳源条件下的剂量响应形状判跨条件稳健性（B4 展示；真·敏感 vs 数值假象留 B3 噪声门）。"""
    present = [str(shape) for shape in condition_shapes.values() if shape]
    if len(present) < 2:
        return "cross_condition_single"
    return "cross_condition_stable" if len(set(present)) == 1 else "cross_condition_sensitive"


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
                # 名字梯队：正式符号 → 标准显示名 → 蛋白描述名 → 策展常用名，位点号始终随行显示
                # （见 gene_name_annotation.safe_gene_display_label）。此前只回退 common_name，
                # 基因候选因此常只显示 locus tag，研究员看不出是什么基因。
                "candidate": safe_gene_display_label(r),
                "reaction": str(r.get("gene_id")),
                "layer": str(r.get("secretory_process") or "未解析"),
                "effect": float(r["effect"]),
                "growth_retention": float(r.get("growth_retention_ratio") or 1.0),
                "confidence": str(r.get("mapping_confidence") or ""),
                # 名字可信度标注：信息量分档、可自查的库 ID、俗名↔位点待复核标记。
                "annotation_tier": str(r.get("annotation_tier") or ""),
                "annotation_accession": str(r.get("annotation_accession") or ""),
                "identity_review": str(r.get("identity_review") or ""),
            }
        )
    return rows


def build_shortlist_readout(
    frame: pd.DataFrame,
    target_id: str,
    *,
    top_n: int = 8,
    r1_readout_dir: Path | None = None,
    dose_response_dir: Path | None = None,
    condition_matrix_dir: Path | None = None,
) -> dict[str, object]:
    """合成一个 target 的短名单读出 dict（本模块零新增求解，只读缓存）。

    frame: 已经过 load_gene_tradeoff_csv 归一化的筛查 tradeoff 表。
    r1_readout_dir: 缓存 R1 读出目录（None 或缺文件时“为什么受限”段落优雅降级为空）。
    dose_response_dir: 离线 R2 剂量响应缓存目录（None 或缺文件时不附形状、面板显示“未扫描”）。
    """
    floors = _load_r1_floors(r1_readout_dir, target_id)
    shortlist = _oe_shortlist(frame, target_id, top_n)

    # R2 剂量响应形状（越加越好/很快到顶）——只读离线扫描缓存并按 reaction_id 附到候选，本模块零求解。
    dose_response = _load_dose_response(dose_response_dir, target_id)
    shapes_by_reaction = dose_response.get("shapes_by_reaction") or {}
    for row in shortlist:
        shape = shapes_by_reaction.get(str(row["reaction"])) if isinstance(shapes_by_reaction, dict) else None
        if isinstance(shape, dict):
            row["shape"] = shape.get("shape")
            row["shape_max_gain"] = shape.get("max_relative_gain")
            row["shape_half_gain_factor"] = shape.get("half_gain_factor")

    # B2/B4 跨条件稳健性——只读离线条件矩阵缓存并按 reaction_id 附到候选，本模块零求解。
    condition_matrix = _load_condition_matrix(condition_matrix_dir, target_id)
    per_reaction_conditions = condition_matrix.get("per_reaction_across_conditions") or {}
    for row in shortlist:
        by_condition = (
            per_reaction_conditions.get(str(row["reaction"])) if isinstance(per_reaction_conditions, dict) else None
        )
        if isinstance(by_condition, dict) and by_condition:
            shapes = {
                str(cond): (shp.get("shape") if isinstance(shp, dict) else None)
                for cond, shp in by_condition.items()
            }
            row["condition_shapes"] = shapes
            row["cross_condition_robustness"] = _cross_condition_robustness(shapes)

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
        "dose_response_available": bool(shapes_by_reaction),
        "dose_response": {
            "tested_factors": dose_response.get("tested_factors"),
            "baseline_objective": dose_response.get("baseline_objective"),
            "reaction_shapes": dose_response.get("reaction_shapes"),
            "warnings": dose_response.get("warnings"),
        }
        if dose_response
        else {},
        "condition_matrix_available": bool(per_reaction_conditions),
        "condition_matrix": {
            "conditions": [str(c) for c in (condition_matrix.get("conditions") or [])],
            "mu": condition_matrix.get("mu"),
        }
        if condition_matrix
        else {},
    }


__all__ = [
    "GROWTH_RISK_THRESHOLD",
    "STRONG_EFFECT_THRESHOLD",
    "build_shortlist_readout",
]
