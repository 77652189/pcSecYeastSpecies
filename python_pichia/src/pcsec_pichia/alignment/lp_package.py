from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pcsec_pichia.adapters.lp_parser import LpFileSummary, parse_lp_file
from pcsec_pichia.adapters.lp_writer import write_stoichiometric_lp
from pcsec_pichia.core.pichia_model import PichiaModel
from pcsec_pichia.core.target_protein_plan import PLAN_STAGE_ORDER, TargetProteinBuildPlan


@dataclass(frozen=True)
class CurrentTargetLpInputs:
    """Input bundle for building LP alignment packages."""

    target_model: TargetBuildResult
    fixed_mu_model: PichiaModel
    metabolic_enzymes: Any
    secretory_enzymes: Any
    combined_enzymes: Any
    target_enzymedata: Any


@dataclass(frozen=True)
class TargetBuildResult:
    """Result of adding target reactions to a PichiaModel.

    This replaces the many route-specific ``Target*BranchExtensionResult``
    dataclasses from the legacy extender.
    """

    model: PichiaModel
    protein_id: str
    added_reaction_ids: tuple[str, ...]
    added_reaction_count: int
    added_metabolite_count: int
    exchange_reaction_id: str
    status: str = "target_built"


@dataclass(frozen=True)
class CurrentTargetLpAlignmentPackage:
    """Metadata for a generated LP file used to align Python output against a MATLAB baseline."""

    lp_path: Path
    summary_path: Path
    target_id: str
    protein_id: str
    mu: float
    media_type: int
    objective_reaction: str
    production_ratio: float | None
    status: str
    base_reaction_count: int
    base_metabolite_count: int
    target_reaction_count: int
    target_metabolite_count: int
    added_reaction_count: int
    added_metabolite_count: int
    added_reaction_ids: tuple[str, ...]
    reaction_stage_counts: dict[str, int]
    target_reaction_variables: dict[str, str]
    target_variable_min: int | None
    target_variable_max: int | None
    exchange_reaction_id: str
    exchange_variable: str
    target_secretory_reaction_coefficient_count: int
    merged_secretory_reaction_coefficient_count: int
    merged_combined_protein_count: int
    lp_summary: LpFileSummary

    def to_dict(self) -> dict[str, object]:
        return {
            "lp_path": str(self.lp_path),
            "summary_path": str(self.summary_path),
            "target_id": self.target_id,
            "protein_id": self.protein_id,
            "mu": self.mu,
            "media_type": self.media_type,
            "objective_reaction": self.objective_reaction,
            "production_ratio": self.production_ratio,
            "status": self.status,
            "base_reaction_count": self.base_reaction_count,
            "base_metabolite_count": self.base_metabolite_count,
            "target_reaction_count": self.target_reaction_count,
            "target_metabolite_count": self.target_metabolite_count,
            "added_reaction_count": self.added_reaction_count,
            "added_metabolite_count": self.added_metabolite_count,
            "added_reaction_ids": list(self.added_reaction_ids),
            "reaction_stage_counts": self.reaction_stage_counts,
            "target_reaction_variables": self.target_reaction_variables,
            "target_variable_min": self.target_variable_min,
            "target_variable_max": self.target_variable_max,
            "exchange_reaction_id": self.exchange_reaction_id,
            "exchange_variable": self.exchange_variable,
            "target_secretory_reaction_coefficient_count": self.target_secretory_reaction_coefficient_count,
            "merged_secretory_reaction_coefficient_count": self.merged_secretory_reaction_coefficient_count,
            "merged_combined_protein_count": self.merged_combined_protein_count,
            "lp_summary": self.lp_summary.to_dict(),
        }


def build_alignment_package(
    lp_path: Path,
    summary_path: Path,
    target_model: TargetBuildResult,
    inputs: CurrentTargetLpInputs,
    plan: TargetProteinBuildPlan,
    mu: float,
    media_type: int,
    objective_reaction: str | None = None,
    production_ratio: float | None = None,
    unmodeled_er_protein_fraction: float = 0.040,
) -> CurrentTargetLpAlignmentPackage:
    """Build a :class:`CurrentTargetLpAlignmentPackage` from already-prepared inputs.

    This is the single entry point for producing LP alignment packages,
    replacing the 25+ route-specific methods from the legacy engine.
    """
    model_for_lp = inputs.fixed_mu_model
    if production_ratio is not None and production_ratio > 0:
        model_for_lp = model_for_lp.change_rxn_bounds(
            target_model.exchange_reaction_id,
            lower=production_ratio,
            upper=production_ratio,
        )

    write_stoichiometric_lp(
        model_for_lp,
        lp_path,
        objective_reaction=objective_reaction or target_model.exchange_reaction_id,
        metabolic_enzymes=inputs.metabolic_enzymes,
        secretory_enzymes=inputs.secretory_enzymes,
        combined_enzymes=inputs.combined_enzymes,
        mu=mu,
        unmodeled_er_protein_fraction=unmodeled_er_protein_fraction,
        include_mitochondrial_constraint=True,
        include_proteasome_constraint=True,
        include_ribosome_assembly_constraint=True,
    )

    lp_summary = parse_lp_file(lp_path)

    reaction_variables = _target_reaction_variables(
        inputs.target_model.model, plan.reaction_ids
    )
    variable_indices = [
        int(variable[1:]) for variable in reaction_variables.values()
    ]
    base_reaction_count = (
        len(inputs.target_model.model.rxns) - inputs.target_model.added_reaction_count
    )
    base_metabolite_count = (
        len(inputs.target_model.model.mets)
        - inputs.target_model.added_metabolite_count
    )

    package = CurrentTargetLpAlignmentPackage(
        lp_path=lp_path,
        summary_path=summary_path,
        target_id=plan.target_id,
        protein_id=plan.protein_id,
        mu=mu,
        media_type=media_type,
        objective_reaction=objective_reaction or target_model.exchange_reaction_id,
        production_ratio=production_ratio,
        status=target_model.status,
        base_reaction_count=base_reaction_count,
        base_metabolite_count=base_metabolite_count,
        target_reaction_count=len(inputs.target_model.model.rxns),
        target_metabolite_count=len(inputs.target_model.model.mets),
        added_reaction_count=inputs.target_model.added_reaction_count,
        added_metabolite_count=inputs.target_model.added_metabolite_count,
        added_reaction_ids=inputs.target_model.added_reaction_ids,
        reaction_stage_counts=_reaction_stage_counts(plan),
        target_reaction_variables=reaction_variables,
        target_variable_min=min(variable_indices) if variable_indices else None,
        target_variable_max=max(variable_indices) if variable_indices else None,
        exchange_reaction_id=inputs.target_model.exchange_reaction_id,
        exchange_variable=reaction_variables.get(
            inputs.target_model.exchange_reaction_id, ""
        ),
        target_secretory_reaction_coefficient_count=len(
            inputs.target_enzymedata.reaction_coefficients
        ),
        merged_secretory_reaction_coefficient_count=len(
            inputs.secretory_enzymes.reaction_coefficients
        ),
        merged_combined_protein_count=inputs.combined_enzymes.protein_count,
        lp_summary=lp_summary,
    )
    summary_path.write_text(
        json.dumps(package.to_dict(), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return package


def _reaction_stage_counts(plan: TargetProteinBuildPlan) -> dict[str, int]:
    counts = Counter(reaction.stage for reaction in plan.reactions)
    return {stage: counts[stage] for stage in PLAN_STAGE_ORDER if counts[stage]}


def _target_reaction_variables(
    model: PichiaModel, reaction_ids: tuple[str, ...]
) -> dict[str, str]:
    reaction_index = model.reaction_index
    return {
        reaction_id: f"X{reaction_index[reaction_id] + 1}"
        for reaction_id in reaction_ids
    }


def mu_token(mu: float) -> str:
    """Format a float as a filesystem-safe token (e.g. 0.10 -> '0p10')."""
    return f"{mu:.2f}".replace(".", "p")


def scientific_token(value: float) -> str:
    """Format a float as a filesystem-safe scientific token (e.g. 1e-8 -> '1eM08')."""
    text = f"{value:.0e}"
    return text.replace(".", "p").replace("-", "m").replace("+", "")


__all__ = [
    "CurrentTargetLpAlignmentPackage",
    "CurrentTargetLpInputs",
    "TargetBuildResult",
    "build_alignment_package",
    "mu_token",
    "scientific_token",
]
