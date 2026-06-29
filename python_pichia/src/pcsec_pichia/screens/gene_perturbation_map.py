from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from pcsec_pichia.screens._prototype_adapter import classify_secretory_process
from pcsec_pichia.screens.candidate_resolution import reactions_for_gene


PROCESS_LABELS = {
    "ribosome": "翻译",
    "proteasome_degradation": "蛋白降解",
    "disulfide_folding": "ER 折叠 / DSB",
    "n_glycan_processing": "N-糖基化 NG",
    "o_glycan_processing": "O-糖基化 OG",
    "chaperone_folding": "ER 折叠 / 分子伴侣",
    "erad_misfolding": "错误折叠 / ERAD",
    "er_translocation": "ER 转运",
    "er_to_golgi_transport": "ER 到 Golgi 转运",
    "golgi_surface_transport": "Golgi 到胞外运输",
    "secretory_capacity": "分泌容量",
    "metabolic_or_other": "代谢或其它反应",
    "unknown": "未解析",
}


@dataclass(frozen=True)
class GeneReactionMapping:
    reaction_id: str | None
    secretory_process: str
    mapping_level: str
    mapping_confidence: str
    complex_id: str | None
    complex_subunit_ids: tuple[str, ...]
    complex_subunit_stoichiometry: tuple[float, ...]
    interpretation: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "reaction_id": self.reaction_id,
            "secretory_process": self.secretory_process,
            "mapping_level": self.mapping_level,
            "mapping_confidence": self.mapping_confidence,
            "complex_id": self.complex_id,
            "complex_subunit_ids": list(self.complex_subunit_ids),
            "complex_subunit_stoichiometry": list(self.complex_subunit_stoichiometry),
            "interpretation": self.interpretation,
        }


@dataclass(frozen=True)
class GenePerturbationMapping:
    gene_id: str
    resolved: bool
    reaction_count: int
    secretory_processes: tuple[str, ...]
    mapping_status: str
    mapping_confidence: str
    warnings: tuple[str, ...]
    reactions: tuple[GeneReactionMapping, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene_id": self.gene_id,
            "resolved": self.resolved,
            "reaction_count": self.reaction_count,
            "secretory_processes": list(self.secretory_processes),
            "mapping_status": self.mapping_status,
            "mapping_confidence": self.mapping_confidence,
            "warnings": list(self.warnings),
            "reactions": [reaction.to_dict() for reaction in self.reactions],
        }


@dataclass(frozen=True)
class GenePerturbationMapResult:
    genes: tuple[GenePerturbationMapping, ...]

    @property
    def rows(self) -> tuple[dict[str, Any], ...]:
        rows: list[dict[str, Any]] = []
        for gene in self.genes:
            for reaction in gene.reactions:
                rows.append(
                    {
                        "gene_id": gene.gene_id,
                        "resolved": gene.resolved,
                        "reaction_count": gene.reaction_count,
                        "gene_secretory_processes": list(gene.secretory_processes),
                        "mapping_status": gene.mapping_status,
                        "gene_mapping_confidence": gene.mapping_confidence,
                        "gene_warnings": list(gene.warnings),
                        **reaction.to_dict(),
                    }
                )
        return tuple(rows)

    def to_dict(self) -> dict[str, Any]:
        return {
            "genes": [gene.to_dict() for gene in self.genes],
            "rows": list(self.rows),
        }


def build_gene_perturbation_map(
    model: Any,
    gene_ids: tuple[str, ...],
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> GenePerturbationMapResult:
    """Build a read-only gene -> reaction -> secretion-process explanation map.

    This helper never runs a solver and never changes the model. It is intended
    for KO/OE preflight explanation, especially for gene-level OE proxy inputs.
    """
    subunit_metadata = complex_subunits or {}
    genes = tuple(
        _build_gene_mapping(model, gene_id, subunit_metadata)
        for gene_id in _dedupe_gene_ids(gene_ids)
    )
    return GenePerturbationMapResult(genes=genes)


def build_reaction_perturbation_mapping(
    reaction_id: str | None,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> GeneReactionMapping:
    """Build the reaction-level explanation used by formal KO/OE result rows."""
    if not reaction_id:
        return GeneReactionMapping(
            reaction_id=None,
            secretory_process="未解析",
            mapping_level="unresolved",
            mapping_confidence="unresolved",
            complex_id=None,
            complex_subunit_ids=(),
            complex_subunit_stoichiometry=(),
            interpretation="未解析到可解释的模型反应。",
        )
    return _build_reaction_mapping(str(reaction_id), complex_subunits or {})


def _build_gene_mapping(
    model: Any,
    gene_id: str,
    complex_subunits: dict[str, list[dict[str, object]]],
) -> GenePerturbationMapping:
    matched_reactions = reactions_for_gene(model, gene_id)
    if gene_id not in getattr(model, "gene_index", {}) or not matched_reactions:
        warning = "基因未在模型 gene_index 中找到，或无法解析到模型反应。"
        return GenePerturbationMapping(
            gene_id=gene_id,
            resolved=False,
            reaction_count=0,
            secretory_processes=(),
            mapping_status="unresolved",
            mapping_confidence="unresolved",
            warnings=(warning,),
            reactions=(
                GeneReactionMapping(
                    reaction_id=None,
                    secretory_process="未解析",
                    mapping_level="unresolved",
                    mapping_confidence="unresolved",
                    complex_id=None,
                    complex_subunit_ids=(),
                    complex_subunit_stoichiometry=(),
                    interpretation=warning,
                ),
            ),
        )

    reaction_mappings = tuple(
        _build_reaction_mapping(reaction_id, complex_subunits)
        for reaction_id in matched_reactions
    )
    process_labels = tuple(
        dict.fromkeys(reaction.secretory_process for reaction in reaction_mappings)
    )
    confidence = _gene_confidence(
        reaction.mapping_confidence for reaction in reaction_mappings
    )
    return GenePerturbationMapping(
        gene_id=gene_id,
        resolved=True,
        reaction_count=len(reaction_mappings),
        secretory_processes=process_labels,
        mapping_status="resolved",
        mapping_confidence=confidence,
        warnings=(),
        reactions=reaction_mappings,
    )


def _build_reaction_mapping(
    reaction_id: str,
    complex_subunits: dict[str, list[dict[str, object]]],
) -> GeneReactionMapping:
    complex_id = _complex_id_for_reaction(reaction_id)
    subunits = complex_subunits.get(complex_id or "", [])
    process_code = classify_secretory_process(reaction_id)
    process_label = _secretory_process_label(process_code)
    level, confidence = _mapping_level_and_confidence(process_code, bool(subunits))
    return GeneReactionMapping(
        reaction_id=reaction_id,
        secretory_process=process_label,
        mapping_level=level,
        mapping_confidence=confidence,
        complex_id=complex_id,
        complex_subunit_ids=tuple(str(item["subunit_id"]) for item in subunits),
        complex_subunit_stoichiometry=tuple(
            float(item["stoichiometry"]) for item in subunits
        ),
        interpretation=_interpretation(level, process_label),
    )


def _mapping_level_and_confidence(
    process_code: str,
    has_complex_subunits: bool,
) -> tuple[str, str]:
    if process_code == "metabolic_or_other":
        return "metabolic_or_other", "low"
    if process_code == "unknown":
        return "reaction_proxy", "low"
    if has_complex_subunits:
        return "complex_subunit", "medium"
    return "direct_gpr", "high"


def _interpretation(mapping_level: str, process_label: str) -> str:
    if mapping_level == "complex_subunit":
        return (
            f"该基因关联分泌复合体反应（{process_label}）；"
            "OE gene 应解释为 reaction-level OE proxy。"
        )
    if mapping_level == "direct_gpr":
        return (
            f"该基因通过模型 GPR 直接关联该分泌反应（{process_label}）；"
            "KO 可按 GPR 解释，OE 仍按反应代理解释。"
        )
    if mapping_level == "metabolic_or_other":
        return "该基因关联代谢或其它反应；可能间接影响分泌，解释置信度较低。"
    if mapping_level == "reaction_proxy":
        return "该基因关联模型反应，但暂未归类到明确分泌环节；解释置信度较低。"
    return "未解析到可解释的模型反应。"


def _complex_id_for_reaction(reaction_id: str) -> str | None:
    if reaction_id.endswith("_formation"):
        return reaction_id.removesuffix("_formation")
    return None


def _secretory_process_label(process_code: str) -> str:
    return PROCESS_LABELS.get(process_code, process_code)


def _gene_confidence(confidences: Iterable[str]) -> str:
    values = tuple(str(value) for value in confidences)
    if "high" in values:
        return "high"
    if "medium" in values:
        return "medium"
    if "low" in values:
        return "low"
    return "unresolved"


def _dedupe_gene_ids(gene_ids: tuple[str, ...]) -> tuple[str, ...]:
    deduped: list[str] = []
    for gene_id in gene_ids:
        item = str(gene_id or "").strip()
        if item and item not in deduped:
            deduped.append(item)
    return tuple(deduped)


__all__ = [
    "GenePerturbationMapping",
    "GenePerturbationMapResult",
    "GeneReactionMapping",
    "build_gene_perturbation_map",
    "build_reaction_perturbation_mapping",
]
