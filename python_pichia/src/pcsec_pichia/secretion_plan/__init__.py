from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pcsec_pichia.secretion_plan._prototype_adapter import (
    add_misfolding_plan,
    coat_other_stoichiometries,
    dsb_stoichiometries,
    golgi_n_stoichiometries,
    golgi_o_stoichiometries,
    mature_stoichiometries,
    ng_stoichiometries,
    og_stoichiometries,
    opn_like_er_golgi_transport_stoichiometries,
    opn_like_og_misfolding_stoichiometries,
    opn_like_target_stoichiometries,
    opn_like_translocation_stoichiometries,
    is_opn_like_supported_target,
    is_soluble_secretory_supported_target,
    soluble_misfolding_stoichiometries,
    soluble_secretory_target_stoichiometries,
    target_reaction_plan,
    transport_to_secretory_stoichiometries,
)
from pcsec_pichia.targets import TargetSpec


@dataclass(frozen=True)
class SecretionPlanResult:
    target_id: str
    protein_id: str
    supported: bool
    route_kind: str
    reaction_count: int
    stage_counts: dict[str, int]
    reaction_ids: tuple[str, ...]
    ptm_counts: dict[str, int]
    formal_pcsec_simulation_status: str
    raw_plan: dict[str, Any]


def build_secretion_plan(target: TargetSpec) -> SecretionPlanResult:
    raw_plan = target_reaction_plan(target)
    reaction_ids = tuple(str(item["reaction_id"]) for item in raw_plan.get("reactions", []))
    return SecretionPlanResult(
        target_id=target.target_id,
        protein_id=target.protein_id,
        supported=is_supported_target(target),
        route_kind=_route_kind(target),
        reaction_count=int(raw_plan.get("reaction_count", 0)),
        stage_counts={str(key): int(value) for key, value in dict(raw_plan.get("stage_counts", {})).items()},
        reaction_ids=reaction_ids,
        ptm_counts={
            "disulfide_sites": int(target.disulfide_sites),
            "n_glycosylation_sites": int(target.n_glycosylation_sites),
            "o_glycosylation_sites": int(target.o_glycosylation_sites),
            "transmembrane": int(target.transmembrane),
            "gpi_sites": int(target.gpi_sites),
            "cotranslation": int(target.cotranslation),
        },
        formal_pcsec_simulation_status=str(raw_plan.get("formal_pcsec_simulation_status", "")),
        raw_plan=raw_plan,
    )


def is_supported_target(target: TargetSpec) -> bool:
    return is_opn_like_supported_target(target) or is_soluble_secretory_supported_target(target)


def summarize_secretion_plan(target: TargetSpec) -> dict[str, object]:
    result = build_secretion_plan(target)
    return {
        "target_id": result.target_id,
        "protein_id": result.protein_id,
        "supported": result.supported,
        "route_kind": result.route_kind,
        "reaction_count": result.reaction_count,
        "stage_counts": result.stage_counts,
        "ptm_counts": result.ptm_counts,
        "formal_pcsec_simulation_status": result.formal_pcsec_simulation_status,
    }


def _route_kind(target: TargetSpec) -> str:
    if is_opn_like_supported_target(target):
        return "opn_like_soluble_secretory"
    if is_soluble_secretory_supported_target(target):
        return "soluble_secretory"
    return "unsupported"


__all__ = [
    "SecretionPlanResult",
    "add_misfolding_plan",
    "build_secretion_plan",
    "coat_other_stoichiometries",
    "dsb_stoichiometries",
    "golgi_n_stoichiometries",
    "golgi_o_stoichiometries",
    "mature_stoichiometries",
    "ng_stoichiometries",
    "og_stoichiometries",
    "opn_like_er_golgi_transport_stoichiometries",
    "opn_like_og_misfolding_stoichiometries",
    "opn_like_target_stoichiometries",
    "opn_like_translocation_stoichiometries",
    "is_supported_target",
    "summarize_secretion_plan",
    "soluble_misfolding_stoichiometries",
    "soluble_secretory_target_stoichiometries",
    "target_reaction_plan",
    "transport_to_secretory_stoichiometries",
]
