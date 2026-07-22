"""① 候选短名单读出（复用筛查缓存，零新增 LP 求解）。

把已经算好并缓存的筛查结果（catalog KO/OE tradeoff CSV，逐候选 secretion_ratio_vs_wildtype，
两个 target × KO/OE，已带中文 secretory_process）+ 缓存的 R1 LP 瓶颈读出（per-target 影子价格）
合成一份给研发/湿实验看的读出：

    为什么受限（R1）  +  OE 提升候选短名单（筛查）  +  该测什么（R4 价值-of-information）

效率关键：1000+/244 逐候选的 LP 求解在筛查阶段一次性做完并缓存，本工具**零新增求解**——
R4/排序只是对缓存排序的纯后处理，R1 复用缓存读出。剂量响应形状（R2）需要对短名单补扫倍数
（有界的额外求解），本工具不做，作为后续 --with-dose-response 增强。

用法（从 python_pichia/，src/ 在 PYTHONPATH）:
    python tools/build_candidate_shortlist_readout.py \
        --screen-csv ../local_runs/catalog_reaction_screen_hlf_opn_expanded61_v2/catalog_reaction_tradeoff_rows.csv \
        --r1-readout-dir ../local_runs/r1_readout \
        --output-dir ../local_runs/candidate_shortlist_readout
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from pcsec_pichia.analysis import (  # noqa: E402
    prioritize_value_of_information,
    summarize_value_of_information,
)

# OE 相对提升低于此阈值视作“无实质提升”（与 R2 剂量响应的 1e-3 噪声底一致的量级）。
STRONG_EFFECT_THRESHOLD = 0.01  # 1% 相对提升
GROWTH_RISK_THRESHOLD = 0.9  # 生长保持率低于此值提示有生长代价


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen-csv", type=Path, required=True, help="筛查 tradeoff CSV（复用，已求解）")
    parser.add_argument("--r1-readout-dir", type=Path, required=True, help="缓存的 R1 LP 瓶颈读出目录")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--targets", nargs="+", default=["hLF", "OPN_ALPHA_FULL_PROJECT"])
    parser.add_argument("--top-n", type=int, default=8)
    return parser.parse_args()


def _load_r1_floors(r1_dir: Path, target_id: str, top: int = 5) -> list[dict[str, object]]:
    """从缓存 R1 读出取最强的下界 floor（为什么受限）。缺失则返回空。"""
    path = r1_dir / f"target_bottleneck_lp_attribution_{target_id}.json"
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

    剔除 candidate_kind == "complex_oe_hypothesis"（未验证的“整复合体按同比例过表达”猜测，
    与 analyze_single_target / 服务层 genome_wide_screen_shortlist 保持一致）。
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


def _readout_for_target(frame: pd.DataFrame, r1_dir: Path, target_id: str, top_n: int) -> dict[str, object]:
    floors = _load_r1_floors(r1_dir, target_id)
    shortlist = _oe_shortlist(frame, target_id, top_n)

    # R4 价值-of-information：对 OE 短名单排序（分数=相对提升）检近似并列 → 该测什么。零求解。
    # top_k 设为短名单长度，让它扫遍整个短名单的相邻近似并列（不只 top-3），例如两个 ~2% 的
    # N-糖基化候选模型分不清、该测哪个。
    voi = summarize_value_of_information(
        prioritize_value_of_information(
            target_id,
            [(row["candidate"], row["effect"]) for row in shortlist],
            top_k=max(2, len(shortlist)),
        )
    )

    top_effect = shortlist[0]["effect"] if shortlist else 0.0
    has_strong = top_effect >= STRONG_EFFECT_THRESHOLD
    growth_risky = [row for row in shortlist if row["growth_retention"] < GROWTH_RISK_THRESHOLD]
    return {
        "target_id": target_id,
        "why_limited_floors": floors,
        "oe_shortlist": shortlist,
        "has_strong_oe_lever": has_strong,
        "top_effect": top_effect,
        "growth_risky_candidates": [row["candidate"] for row in growth_risky],
        "value_of_information": voi,
    }


def _fmt_pct(value: float) -> str:
    return f"+{value * 100:.2f}%"


def _markdown(readout: dict[str, object]) -> str:
    target = readout["target_id"]
    floors = readout["why_limited_floors"]
    shortlist = readout["oe_shortlist"]
    voi = readout["value_of_information"]
    lines: list[str] = []
    lines.append(f"# {target}：候选短名单读出（复用筛查缓存，零新增 LP 求解）")
    lines.append("")
    # 一句话结论
    if shortlist and readout["has_strong_oe_lever"]:
        headline = f"OE 提升候选里 **{shortlist[0]['candidate']}**（{shortlist[0]['layer']}）最强（{_fmt_pct(shortlist[0]['effect'])}）"
    else:
        headline = "**没有强 OE 提升杠杆**（最高相对提升低于 1%）——这个靶点大概率不受限于可 OE 的分泌机器上限"
    top_floor = floors[0]["reaction_id"] if floors else "（无缓存 R1 读出）"
    lines.append(f"> 一句话：{target} 的分泌最强约束在 `{top_floor}`（下界/最低要求，OE 动不了）；{headline}。")
    lines.append("")

    lines.append("## 1. 为什么受限（R1 LP 影子价格 · 最强约束层）")
    lines.append("下界=最低要求类约束，承载最大影子价格，是“卡在哪一层”的答案；但 OE 放宽的是上限、对它们无效。")
    if floors:
        lines.append("")
        lines.append("| 反应 | 影子价格(绝对值) |")
        lines.append("|---|---|")
        for f in floors:
            lines.append(f"| `{f['reaction_id']}` | {f['abs_marginal']:.1f} |")
    else:
        lines.append("（无缓存 R1 读出；可用 `run_target_bottleneck_lp_attribution_check.py` 生成）")
    lines.append("")

    lines.append(f"## 2. OE 提升候选短名单（来自筛查缓存，按相对提升排序，top-{len(shortlist)}）")
    lines.append("相对提升 = 该 OE 相对野生型的分泌比值 − 1（固定 2× OE、corrected 培养基下的模型解）。")
    if shortlist:
        lines.append("")
        lines.append("| 排名 | 候选 | 资源层 | 相对提升 | 生长保持 | 证据置信度 |")
        lines.append("|---|---|---|---|---|---|")
        for i, row in enumerate(shortlist, start=1):
            growth = f"{row['growth_retention']:.2f}" + ("（有生长代价）" if row["growth_retention"] < GROWTH_RISK_THRESHOLD else "")
            lines.append(
                f"| {i} | {row['candidate']} | {row['layer']} | {_fmt_pct(row['effect'])} | {growth} | {row['confidence']} |"
            )
    else:
        lines.append("（无 ratio>1 的 OE 候选）")
    lines.append("")

    lines.append("## 3. 该测什么（R4 价值-of-information）")
    lines.append("模型给的是相对排序、不是绝对产量。这里标出顶部名次里模型分不清的候选，并给出最能消解歧义的最小测量。")
    items = voi.get("information_items") or []
    if voi.get("has_actionable_ambiguity") and items:
        lines.append("")
        for it in items:
            cands = "、".join(str(c) for c in (it.get("candidates") or []))
            lines.append(f"- **优先级 {it.get('priority_rank')}**（{cands}）：对这几个做靶点特异的分泌定量湿实验以定序。")
    elif shortlist:
        lines.append("")
        lines.append("- 顶部候选相对提升分离明显，当前排序较可信；优先验证榜首即可。")
    else:
        lines.append("")
        lines.append("- 无可排序的 OE 提升候选，暂无“该测什么”的明确建议。")
    lines.append("")

    lines.append("## 4. 诚实边界")
    lines.append("- 相对信号，**非绝对产量 / mg·L⁻¹**；复用的是筛查在固定 2× OE、corrected 培养基下的模型解。")
    lines.append("- 干预以复合体/反应表达，使用时需人工翻译成“过表达对应基因”。")
    lines.append("- 剂量响应形状（越加越好还是很快到顶）需对短名单补扫倍数（R2），本读出未含。")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    frame = pd.read_csv(args.screen_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for target_id in args.targets:
        readout = _readout_for_target(frame, args.r1_readout_dir, target_id, args.top_n)
        (args.output_dir / f"{target_id}_shortlist_readout.json").write_text(
            json.dumps(readout, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
        )
        md = _markdown(readout)
        (args.output_dir / f"{target_id}_shortlist_readout.md").write_text(md, encoding="utf-8")
        print(f"[{target_id}] wrote readout: strong_oe_lever={readout['has_strong_oe_lever']} "
              f"top={readout['oe_shortlist'][0]['candidate'] if readout['oe_shortlist'] else None} "
              f"voi_items={len(readout['value_of_information'].get('information_items') or [])}")
    print("DONE ->", args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
