"""改造后菌株的"下一步 OE 候选"读出（复用 R1 瓶颈 + R2 有界剂量响应；纯装配、零求解）。

给一个已 solve 的（可能 KO/OE 改造过的）菌株：它的 binding 上限瓶颈复合体
（`oe_actionable_bottlenecks`）就是"下一步 OE 能松哪一层"的线索；再对这些复合体的**有界** OE
剂量响应（R2）量化"松开能涨多少、涨到哪饱和"，按真实效应排序，给出 per-strain 的下一候选。

- 相对信号、复合体级、非绝对产量；不重扫全基因组、不动求解核心。
- 本模块是**纯装配**（零求解）：影子价格瓶颈来自 solve 结果、剂量响应来自调用方（C2 编排）跑好的
  有界 sweep。求解编排不在这里。
"""

from __future__ import annotations

from typing import Any

DEFAULT_TOP_N = 6

NEXT_OE_CANDIDATE_CAVEATS = (
    "复合体级候选：reaction_id 即 OE 目标（放宽该复合体的产能上限）。",
    "瓶颈会转移：松开当前 #1 后下一层约束会顶上来——每改一轮都应重跑。",
    "相对信号、非绝对产量；不含降解 / 表达调控等模型范围外因素。",
)


def _abs_marginal(entry: dict[str, Any]) -> float:
    value = entry.get("abs_marginal")
    if value is None:
        value = abs(float(entry.get("marginal") or 0.0))
    return abs(float(value or 0.0))


def build_next_oe_candidates_readout(
    oe_actionable_bottlenecks: list[dict[str, Any]] | tuple[dict[str, Any], ...],
    dose_response: dict[str, Any] | None = None,
    *,
    top_n: int = DEFAULT_TOP_N,
) -> dict[str, Any]:
    """从改造后菌株的 binding 上限瓶颈 + 有界剂量响应，装配"下一步 OE 候选"读出（零求解）。

    oe_actionable_bottlenecks: solve 结果 `lp_attribution` 里的 binding 上限复合体条目
        （含 `reaction_id` / `secretory_process` / `abs_marginal`|`marginal`）。
    dose_response: 对这些复合体跑的**有界** OE 剂量响应结果（`shapes_by_reaction`: reaction_id→形状）；
        None / 空则只按影子价格排、不带效应。有剂量响应则按**真实效应**排（谁 OE 真涨得多），
        因为影子价格只说"binding"、不说"松开能涨多少"。
    """
    entries = [
        entry
        for entry in (oe_actionable_bottlenecks or [])
        if isinstance(entry, dict) and entry.get("reaction_id")
    ]
    entries = sorted(entries, key=_abs_marginal, reverse=True)[: max(0, int(top_n))]

    shapes = (dose_response or {}).get("shapes_by_reaction") or {}
    dose_available = bool(shapes)

    candidates: list[dict[str, Any]] = []
    for entry in entries:
        reaction_id = str(entry["reaction_id"])
        row: dict[str, Any] = {
            "reaction": reaction_id,
            "layer": str(entry.get("secretory_process") or "未解析"),
            "shadow_price": _abs_marginal(entry),
        }
        shape = shapes.get(reaction_id) if isinstance(shapes, dict) else None
        if isinstance(shape, dict):
            row["shape"] = shape.get("shape")
            row["effect"] = shape.get("max_relative_gain")
            row["half_gain_factor"] = shape.get("half_gain_factor")
        candidates.append(row)

    if dose_available:
        candidates.sort(key=lambda row: float(row.get("effect") or 0.0), reverse=True)

    return {
        "candidates": candidates,
        "dose_response_available": dose_available,
        "baseline_objective": (dose_response or {}).get("baseline_objective"),
        "top_n": int(top_n),
        "caveats": list(NEXT_OE_CANDIDATE_CAVEATS),
    }


__all__ = ["DEFAULT_TOP_N", "NEXT_OE_CANDIDATE_CAVEATS", "build_next_oe_candidates_readout"]
