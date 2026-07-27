"""D3 分层短名单编排（服务层 · ADR-004 #1 迭代2）。

把 D4 复用地基（同口径野生型 OE+KO 短名单）+ C2 改造后/野生型瓶颈 + D2 受影响层打标 编排成
"改造后菌株的下一步候选短名单"：

- **L1（即时 · 本模块）**：读 D4 缓存的野生型短名单 → 跑一次野生型瓶颈 + 一次改造后瓶颈（C2）→ D2 打标
  （可复用 / 已失效）。给出"哪些先前候选对这个改造后菌株还可信、哪些要重算"，**不重算**（快）。
- **L2（按需 · 见 `recompute_stale_candidates`）**：只对**已失效**候选按改造后基线重算——OE 复用 R2 剂量响应
  sweep、KO 复用 D1 `run_knockout_screen`（都带 `strain_modifications`）——与复用值合并重排。

诚实边界（D2）：复用只对与瓶颈结构无关的**分泌专属层**干净；代谢/未解析层保守（一律已失效）。相对信号、
非绝对产量。若该口径下没有 D4 野生型基线 → 诚实报"需先跑后台全基因组基线"（不伪造短名单）。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Sequence

from app import ensure_python_pichia_on_path
from app.services.per_strain_layer_reuse import REUSE_CAVEAT, build_layer_reuse_tags
from app.services.strain_baseline_service import load_strain_baseline_readout

ensure_python_pichia_on_path()

DEFAULT_SHORTLIST_TOP_N = 8


def _select_shortlist(rows: Sequence[dict[str, Any]], top_n: int) -> list[dict[str, Any]]:
    """从基线候选行里取"有实质提升"（secretion_ratio>1）的 top_n，按相对效应降序。

    effect = ratio - 1（与 build_shortlist_readout 同口径）；ratio 缺失/≤1 的剔除（只留帮上忙的杠杆）。
    """
    scored: list[dict[str, Any]] = []
    for row in rows:
        ratio = row.get("secretion_ratio_vs_wildtype")
        if ratio in (None, ""):
            continue
        try:
            effect = float(ratio) - 1.0
        except (TypeError, ValueError):
            continue
        if effect <= 0:
            continue
        scored.append(
            {
                "candidate": str(row.get("common_name") or row.get("gene_id") or ""),
                "gene_id": str(row.get("gene_id") or ""),
                "secretory_process": row.get("secretory_process") or "",
                "affected_reactions": row.get("affected_reactions") or "",
                "intervention_type": str(row.get("intervention_type") or "").upper(),
                "wildtype_effect": effect,
                "growth_retention": float(row.get("growth_retention_ratio") or 1.0),
                "mapping_confidence": row.get("mapping_confidence") or "",
            }
        )
    scored.sort(key=lambda r: r["wildtype_effect"], reverse=True)
    return scored[: max(0, int(top_n))]


def _analyze_bottlenecks(
    *,
    target_id: str,
    ko_reaction_ids: Sequence[str],
    oe_reaction_ids: Sequence[str],
    oe_factor: float,
    mu: float,
    media_type: int,
    carbon_source_id: str,
    compatibility_mode: str,
    enable_ribosome_translation_constraint: bool,
    enable_misfolding_constraint: bool,
    root: Path | None,
) -> dict[str, Any]:
    """薄壳：调 C2 引擎拿一次瓶颈归因（懒 import，方便单测 monkeypatch）。"""
    from pcsec_pichia.next_oe_candidates import analyze_next_oe_candidates

    return analyze_next_oe_candidates(
        target_id=target_id,
        ko_reaction_ids=tuple(ko_reaction_ids),
        oe_reaction_ids=tuple(oe_reaction_ids),
        oe_factor=oe_factor,
        mu=mu,
        media_type=media_type,
        carbon_source_id=carbon_source_id,
        compatibility_mode=compatibility_mode,
        enable_ribosome_translation_constraint=enable_ribosome_translation_constraint,
        enable_misfolding_constraint=enable_misfolding_constraint,
        root=root,
    )


def build_modified_strain_shortlist(
    *,
    target_id: str,
    ko_reaction_ids: Sequence[str] = (),
    oe_reaction_ids: Sequence[str] = (),
    oe_factor: float = 2.0,
    mu: float = 0.10,
    media_type: int = 4,
    carbon_source_id: str = "glucose",
    compatibility_mode: str = "corrected",
    enable_ribosome_translation_constraint: bool = False,
    enable_misfolding_constraint: bool = False,
    top_n: int = DEFAULT_SHORTLIST_TOP_N,
    cache_dir: Path | str | None = None,
    root: Path | None = None,
    _analyze: Callable[..., dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """L1：改造后菌株的下一步候选短名单（复用地基 + 瓶颈打标，不重算）。

    未命中 D4 野生型基线（该口径没跑过后台全量）→ 返回 `available=False, needs_baseline_build=True`，
    面板据此指引先跑后台基线（不伪造短名单）。`_analyze` 仅供测试注入。
    """
    baseline_kwargs = dict(
        target_id=target_id,
        carbon_source_id=carbon_source_id,
        media_type=media_type,
        growth_rate=mu,
        enable_ribosome_translation_constraint=enable_ribosome_translation_constraint,
        enable_misfolding_constraint=enable_misfolding_constraint,
        compatibility_mode=compatibility_mode,
    )
    if cache_dir is not None:
        baseline_kwargs["cache_dir"] = cache_dir
    baseline = load_strain_baseline_readout(**baseline_kwargs)
    if not baseline.get("available"):
        return {
            "available": False,
            "needs_baseline_build": True,
            "caliber": baseline.get("caliber"),
            "reason": baseline.get("reason")
            or "该口径下还没有野生型全基因组后台基线，无法复用；请先跑后台基线构建。",
            "oe_candidates": [],
            "ko_candidates": [],
        }

    analyze = _analyze or _analyze_bottlenecks
    common = dict(
        target_id=target_id, oe_factor=oe_factor, mu=mu, media_type=media_type,
        carbon_source_id=carbon_source_id, compatibility_mode=compatibility_mode,
        enable_ribosome_translation_constraint=enable_ribosome_translation_constraint,
        enable_misfolding_constraint=enable_misfolding_constraint, root=root,
    )
    wildtype_result = analyze(ko_reaction_ids=(), oe_reaction_ids=(), **common)
    modified_result = analyze(ko_reaction_ids=ko_reaction_ids, oe_reaction_ids=oe_reaction_ids, **common)

    oe_shortlist = _select_shortlist(baseline.get("oe_rows") or [], top_n)
    ko_shortlist = _select_shortlist(baseline.get("ko_rows") or [], top_n)

    oe_tags = build_layer_reuse_tags(
        wildtype_result=wildtype_result, modified_result=modified_result, shortlist_rows=oe_shortlist
    )
    ko_tags = build_layer_reuse_tags(
        wildtype_result=wildtype_result, modified_result=modified_result, shortlist_rows=ko_shortlist
    )
    return {
        "available": True,
        "layer": "L1",
        "target_id": target_id,
        "caliber": baseline.get("caliber"),
        "baseline_built_at": baseline.get("built_at"),
        "baseline_source_run": baseline.get("source_run"),
        "modified_solve_success": bool(modified_result.get("modified_solve_success")),
        "affected_modules": oe_tags["affected_modules"],  # OE/KO 同一 wt/mod 瓶颈 → affected 相同
        "wildtype_bottleneck_modules": oe_tags["wildtype_bottleneck_modules"],
        "modified_bottleneck_modules": oe_tags["modified_bottleneck_modules"],
        "oe_candidates": oe_tags["tagged_candidates"],
        "ko_candidates": ko_tags["tagged_candidates"],
        "oe_reusable_count": oe_tags["reusable_count"],
        "oe_stale_count": oe_tags["stale_count"],
        "ko_reusable_count": ko_tags["reusable_count"],
        "ko_stale_count": ko_tags["stale_count"],
        "applied_modifications": modified_result.get("applied_modifications"),
        "caveat": REUSE_CAVEAT,
    }


def _split_reactions(affected_reactions: object) -> list[str]:
    if not affected_reactions:
        return []
    return [r.strip() for r in str(affected_reactions).split(";") if r.strip()]


def _recompute_effects(**kwargs: Any) -> dict[str, Any]:
    """薄壳：调 C2/L2 引擎重算（懒 import，方便单测 monkeypatch）。"""
    from pcsec_pichia.next_oe_candidates import recompute_modified_strain_candidate_effects

    return recompute_modified_strain_candidate_effects(**kwargs)


def _apply_recomputed(
    candidates: Sequence[dict[str, Any]],
    effects: dict[str, float],
    *,
    by: str,
) -> list[dict[str, Any]]:
    """把改造后重算效应并回候选并按有效效应重排（纯装配）。

    已失效且拿到重算值 → `effective_effect` = 改造后效应、`recompute_status="recomputed"`；已失效但无值
    （改造后不可行 / 缺）→ 回退野生型值 + `recompute_status="recompute_failed"`（显式标，不假装重算过）；
    可复用 → `effective_effect` = 野生型值、`recompute_status="reused"`。
    """
    merged: list[dict[str, Any]] = []
    for candidate in candidates:
        row = dict(candidate)
        wildtype_effect = float(row.get("wildtype_effect") or 0.0)
        if row.get("reuse_status") == "stale":
            keys = _split_reactions(row.get("affected_reactions")) if by == "reactions" else [row.get("gene_id")]
            recomputed = [effects[k] for k in keys if k and k in effects]
            if recomputed:
                effect = max(recomputed) if by == "reactions" else recomputed[0]
                row["modified_effect"] = effect
                row["effective_effect"] = effect
                row["recompute_status"] = "recomputed"
            else:
                row["effective_effect"] = wildtype_effect
                row["recompute_status"] = "recompute_failed"
        else:
            row["effective_effect"] = wildtype_effect
            row["recompute_status"] = "reused"
        merged.append(row)
    merged.sort(key=lambda r: r.get("effective_effect") or 0.0, reverse=True)
    return merged


def recompute_stale_candidates(
    l1_readout: dict[str, Any],
    *,
    dose_response_factors: Sequence[float] = (),
    root: Path | None = None,
    _recompute: Callable[..., dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """L2：对 L1 里已失效候选在改造后菌株上重算（OE 复用 R2、KO 复用 D1），与复用值合并重排。

    未命中基线（L1 已诚实报缺）→ 原样透传。`_recompute` 仅供测试注入。
    """
    if not l1_readout.get("available"):
        return dict(l1_readout)
    caliber = l1_readout.get("caliber") or {}
    mods = l1_readout.get("applied_modifications") or {}
    oe_candidates = l1_readout.get("oe_candidates") or []
    ko_candidates = l1_readout.get("ko_candidates") or []

    stale_oe_reactions = [
        reaction
        for candidate in oe_candidates
        if candidate.get("reuse_status") == "stale"
        for reaction in _split_reactions(candidate.get("affected_reactions"))
    ]
    stale_ko_genes = [
        candidate["gene_id"]
        for candidate in ko_candidates
        if candidate.get("reuse_status") == "stale" and candidate.get("gene_id")
    ]

    recompute = _recompute or _recompute_effects
    effects = recompute(
        target_id=l1_readout.get("target_id"),
        ko_reaction_ids=tuple(mods.get("ko_reaction_ids") or ()),
        oe_reaction_ids=tuple(mods.get("oe_reaction_ids") or ()),
        oe_factor=float(mods.get("oe_factor") or 2.0),
        stale_oe_reactions=stale_oe_reactions,
        stale_ko_genes=stale_ko_genes,
        mu=float(caliber.get("growth_rate", 0.10)),
        media_type=int(caliber.get("media_type", 4)),
        carbon_source_id=str(caliber.get("carbon_source_id", "glucose")),
        compatibility_mode=str(caliber.get("compatibility_mode", "corrected")),
        enable_ribosome_translation_constraint=bool(caliber.get("write_ribosome_translation_constraint", False)),
        enable_misfolding_constraint=bool(caliber.get("write_misfolding_constraints", False)),
        dose_response_factors=tuple(dose_response_factors),
        root=root,
    )
    return {
        **l1_readout,
        "layer": "L2",
        "oe_candidates": _apply_recomputed(oe_candidates, effects.get("oe_effects") or {}, by="reactions"),
        "ko_candidates": _apply_recomputed(ko_candidates, effects.get("ko_effects") or {}, by="gene"),
        "recomputed_oe_count": effects.get("recomputed_oe_count", 0),
        "recomputed_ko_count": effects.get("recomputed_ko_count", 0),
        "recompute_warnings": list(effects.get("warnings") or []),
    }


__all__ = [
    "DEFAULT_SHORTLIST_TOP_N",
    "build_modified_strain_shortlist",
    "recompute_stale_candidates",
]
