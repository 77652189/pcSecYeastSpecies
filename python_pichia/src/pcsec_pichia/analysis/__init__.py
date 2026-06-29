from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from pcsec_pichia.secretion_plan import SecretionPlanResult, build_secretion_plan
from pcsec_pichia.targets import TargetSpec, target_features


@dataclass(frozen=True)
class ProteinCostItem:
    category: str
    label: str
    basis: str
    raw_value: float
    relative_score: float
    interpretation: str


@dataclass(frozen=True)
class ProteinCostAnalysisResult:
    target_id: str
    protein_id: str
    route_kind: str
    sequence_features: dict[str, Any]
    ptm_counts: dict[str, int]
    stage_counts: dict[str, int]
    cost_items: tuple[ProteinCostItem, ...]
    total_relative_score: float
    dominant_cost_categories: tuple[str, ...]
    warnings: tuple[str, ...]
    result_status: str = "draft_explanatory"


def analyze_target_protein_cost(
    target: TargetSpec,
    plan: SecretionPlanResult | None = None,
) -> ProteinCostAnalysisResult:
    resolved_plan = plan or build_secretion_plan(target)
    features = _target_feature_payload(target)
    raw_items = _cost_items_from_target_and_plan(target, resolved_plan, features)
    cost_items = _normalise_cost_scores(raw_items)
    return ProteinCostAnalysisResult(
        target_id=target.target_id,
        protein_id=target.protein_id,
        route_kind=resolved_plan.route_kind,
        sequence_features=features,
        ptm_counts=dict(resolved_plan.ptm_counts),
        stage_counts=dict(resolved_plan.stage_counts),
        cost_items=cost_items,
        total_relative_score=round(sum(item.relative_score for item in cost_items), 3),
        dominant_cost_categories=_dominant_categories(cost_items),
        warnings=tuple(_cost_warnings(target)),
    )


def summarize_protein_cost_analysis(result: ProteinCostAnalysisResult) -> dict[str, object]:
    return {
        "target_id": result.target_id,
        "protein_id": result.protein_id,
        "route_kind": result.route_kind,
        "sequence_features": result.sequence_features,
        "ptm_counts": result.ptm_counts,
        "stage_counts": result.stage_counts,
        "cost_items": build_cost_item_table(result),
        "total_relative_score": result.total_relative_score,
        "dominant_cost_categories": list(result.dominant_cost_categories),
        "warnings": list(result.warnings),
        "result_status": result.result_status,
    }


def build_cost_item_table(result: ProteinCostAnalysisResult) -> tuple[dict[str, object], ...]:
    return tuple(asdict(item) for item in result.cost_items)


def _target_feature_payload(target: TargetSpec) -> dict[str, Any]:
    features = target_features(target)
    return {
        "full_sequence_length": int(features.get("full_length") or 0),
        "mature_sequence_length": int(features.get("mature_length") or 0),
        "leader_sequence_length": int(features.get("leader_length") or 0),
        "signal_peptide_length": int(features.get("signal_peptide_length") or 0),
        "protein_mw": float(features.get("protein_mw") or 0.0),
        "cysteine_count": int(features.get("cysteine_count") or 0),
        "inferred_disulfide_pairs_from_cysteines": int(
            features.get("inferred_disulfide_pairs_from_cysteines") or 0
        ),
        "n_glycosylation_motif_positions": list(features.get("n_glycosylation_motif_positions") or []),
        "ser_thr_count": int(features.get("ser_thr_count") or 0),
        "valid_mature_sequence": bool(features.get("valid_mature_sequence")),
    }


def _cost_items_from_target_and_plan(
    target: TargetSpec,
    plan: SecretionPlanResult,
    features: dict[str, Any],
) -> tuple[ProteinCostItem, ...]:
    full_len = float(features["full_sequence_length"])
    mature_len = float(features["mature_sequence_length"])
    leader_len = float(features["leader_sequence_length"])
    signal_len = float(features["signal_peptide_length"])
    stage_counts = plan.stage_counts
    dsb = int(target.disulfide_sites)
    ng = int(target.n_glycosylation_sites)
    og = int(target.o_glycosylation_sites)
    motif_count = len(features["n_glycosylation_motif_positions"])
    ser_thr_count = int(features["ser_thr_count"])
    ptm_burden = dsb + ng + og + int(target.gpi_sites) + int(target.transmembrane)

    return (
        ProteinCostItem(
            "translation",
            "翻译负担",
            "full sequence length",
            full_len,
            0.0,
            "目标蛋白越长，翻译相关氨基酸与核糖体占用越高。",
        ),
        ProteinCostItem(
            "leader_signal_processing",
            "信号肽/leader 处理",
            "leader length + signal peptide length",
            leader_len + signal_len,
            0.0,
            "分泌 leader 和信号肽会增加进入分泌途径前后的处理负担。",
        ),
        ProteinCostItem(
            "er_translocation",
            "ER 转运",
            "through ER route and translocation reactions",
            float(stage_counts.get("translocation", 0) * 20 + (10 if target.through_er else 0)),
            0.0,
            "经 ER 分泌的目标蛋白需要转运与 Sec 通路相关步骤。",
        ),
        ProteinCostItem(
            "folding_dsb",
            "二硫键折叠 DSB",
            "declared DSB count",
            float(dsb * 18),
            0.0,
            "声明的二硫键数量越多，PDI/氧化折叠相关负担越高。",
        ),
        ProteinCostItem(
            "n_glycosylation",
            "N-糖基化 NG",
            "declared NG count and motif count",
            float(ng * 25 + (motif_count * 2 if ng > 0 else 0)),
            0.0,
            "N-糖基化位点和可见 motif 会增加 ER/Golgi 加工解释成本。",
        ),
        ProteinCostItem(
            "o_glycosylation",
            "O-糖基化 OG",
            "declared OG count and Ser/Thr support",
            float(og * 12 + min(ser_thr_count, og)),
            0.0,
            "O-糖基化位点会增加 ER/Golgi O-linked processing 负担。",
        ),
        ProteinCostItem(
            "misfolding_erad",
            "错误折叠/ERAD 风险",
            "length, PTM burden, and misfolding plan rows",
            float(mature_len / 20 + ptm_burden * 3 + stage_counts.get("misfolding", 0) * 10),
            0.0,
            "长序列和复杂 PTM 更容易带来折叠质量控制压力。",
        ),
        ProteinCostItem(
            "transport_secretion",
            "囊泡运输与分泌",
            "ER-Golgi, Golgi, final transport, exchange rows",
            float(
                (
                    stage_counts.get("er_to_golgi", 0)
                    + stage_counts.get("golgi_processing", 0)
                    + stage_counts.get("final_transport", 0)
                    + stage_counts.get("exchange", 0)
                )
                * 10
            ),
            0.0,
            "ER-Golgi 与最终分泌步骤越多，分泌路径解释成本越高。",
        ),
    )


def _normalise_cost_scores(items: tuple[ProteinCostItem, ...]) -> tuple[ProteinCostItem, ...]:
    total = sum(max(0.0, item.raw_value) for item in items)
    if total <= 0:
        return items
    return tuple(
        ProteinCostItem(
            category=item.category,
            label=item.label,
            basis=item.basis,
            raw_value=round(item.raw_value, 6),
            relative_score=round(max(0.0, item.raw_value) / total * 100.0, 3),
            interpretation=item.interpretation,
        )
        for item in items
    )


def _dominant_categories(items: tuple[ProteinCostItem, ...], limit: int = 3) -> tuple[str, ...]:
    ordered = sorted(items, key=lambda item: item.relative_score, reverse=True)
    return tuple(item.category for item in ordered[: max(0, int(limit))] if item.relative_score > 0)


def _cost_warnings(target: TargetSpec) -> list[str]:
    warnings = [
        "当前成本分析是 Python draft explanatory score，不代表真实发酵产量或湿实验成本。",
        "该分析不使用 LP shadow price，也不改变 corrected pipeline 的求解目标或约束。",
    ]
    if target.target_id == "hLF":
        warnings.append("hLF 使用项目定义 710aa 序列；结果不是 old MATLAB hLF fully aligned。")
    if target.source.startswith("request") or target.source.startswith("json:") or target.target_id.upper().startswith("CUSTOM"):
        if target.disulfide_sites == 0 and target.n_glycosylation_sites == 0 and target.o_glycosylation_sites == 0:
            warnings.append("custom target 的 DSB/NG/OG 均为 0；未声明 PTM 不会被自动推断进成本项。")
    return warnings


__all__ = [
    "ProteinCostAnalysisResult",
    "ProteinCostItem",
    "analyze_target_protein_cost",
    "build_cost_item_table",
    "summarize_protein_cost_analysis",
]
