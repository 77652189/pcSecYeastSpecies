from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import sparse

from pcsec_pichia.constraints._prototype_adapter import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    TargetSpec,
    build_supported_target_model,
    build_target_enzymedata,
    build_pcsec_constraint_matrices,
    metabolic_coupling_rows,
    misfolding_constraint_rows,
    mitochondrial_rows,
    proteasome_rows,
    protein_mass_rows,
    ribosome_assembly_rows,
    ribosome_translation_rows,
    secretory_coupling_rows,
)


@dataclass(frozen=True)
class PcSecConstraintResult:
    target_id: str
    protein_id: str
    supported: bool
    build_status: str
    exchange_reaction_id: str | None
    mu: float
    write_ribosome_translation_constraint: bool
    write_misfolding_constraints: bool
    constraint_counts: dict[str, int]
    eq_shape: tuple[int, int] | None
    ub_shape: tuple[int, int] | None
    matlab_alignment_status: str
    A_eq: sparse.csr_matrix | None = None
    b_eq: np.ndarray | None = None
    A_ub: sparse.csr_matrix | None = None
    b_ub: np.ndarray | None = None


def build_pcsec_constraints(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float = 0.10,
    total_protein_content: float = 0.37,
    unmodeled_er_protein_fraction: float = 0.040,
    mitochondrial_protein_fraction: float = 0.05,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> PcSecConstraintResult:
    """Build pcSec constraint matrices for a supported target without solving LP."""

    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None:
        return PcSecConstraintResult(
            target_id=target.target_id,
            protein_id=target.protein_id,
            supported=False,
            build_status=build.status,
            exchange_reaction_id=build.exchange_reaction_id,
            mu=float(mu),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
            constraint_counts={},
            eq_shape=None,
            ub_shape=None,
            matlab_alignment_status="pending",
        )

    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = combined.with_target(target_enzymedata)
    fixed_model = build.model.with_bounds({"BIOMASS": (mu, mu)})
    A_eq, b_eq, A_ub, b_ub, counts = build_pcsec_constraint_matrices(
        fixed_model,
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu=mu,
        total_protein_content=total_protein_content,
        unmodeled_er_protein_fraction=unmodeled_er_protein_fraction,
        mitochondrial_protein_fraction=mitochondrial_protein_fraction,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    return PcSecConstraintResult(
        target_id=target.target_id,
        protein_id=target.protein_id,
        supported=True,
        build_status=build.status,
        exchange_reaction_id=build.exchange_reaction_id,
        mu=float(mu),
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
        constraint_counts={str(key): int(value) for key, value in counts.items()},
        eq_shape=(int(A_eq.shape[0]), int(A_eq.shape[1])),
        ub_shape=(int(A_ub.shape[0]), int(A_ub.shape[1])),
        matlab_alignment_status="pending",
        A_eq=A_eq,
        b_eq=b_eq,
        A_ub=A_ub,
        b_ub=b_ub,
    )


def summarize_pcsec_constraints(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float = 0.10,
    total_protein_content: float = 0.37,
    unmodeled_er_protein_fraction: float = 0.040,
    mitochondrial_protein_fraction: float = 0.05,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> dict[str, Any]:
    result = build_pcsec_constraints(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        mu=mu,
        total_protein_content=total_protein_content,
        unmodeled_er_protein_fraction=unmodeled_er_protein_fraction,
        mitochondrial_protein_fraction=mitochondrial_protein_fraction,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    return {
        "target_id": result.target_id,
        "protein_id": result.protein_id,
        "supported": result.supported,
        "build_status": result.build_status,
        "exchange_reaction_id": result.exchange_reaction_id,
        "mu": result.mu,
        "write_ribosome_translation_constraint": result.write_ribosome_translation_constraint,
        "write_misfolding_constraints": result.write_misfolding_constraints,
        "constraint_counts": result.constraint_counts,
        "eq_shape": result.eq_shape,
        "ub_shape": result.ub_shape,
        "matlab_alignment_status": result.matlab_alignment_status,
    }


__all__ = [
    "PcSecConstraintResult",
    "build_pcsec_constraints",
    "build_pcsec_constraint_matrices",
    "metabolic_coupling_rows",
    "misfolding_constraint_rows",
    "mitochondrial_rows",
    "proteasome_rows",
    "protein_mass_rows",
    "ribosome_assembly_rows",
    "ribosome_translation_rows",
    "secretory_coupling_rows",
    "summarize_pcsec_constraints",
]
