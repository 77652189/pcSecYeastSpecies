from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from pcsec_pichia.core.pichia_enzymes import CombinedEnzymeData, MetabolicEnzymeData, SecretoryEnzymeData
from pcsec_pichia.core.pichia_model import PichiaModel, collect_compartment_reactions


LpObjectiveSense = Literal["Maximize", "Minimize"]


@dataclass(frozen=True)
class LpWriteSummary:
    path: Path
    objective_reaction: str
    objective_variable: str
    optimization_sense: LpObjectiveSense
    stoichiometric_constraint_count: int
    metabolic_coupling_constraint_count: int
    secretory_coupling_constraint_count: int
    protein_mass_constraint_count: int
    mitochondrial_constraint_count: int
    proteasome_constraint_count: int
    ribosome_assembly_constraint_count: int
    bounds_count: int
    variable_count: int


def write_stoichiometric_lp(
    model: PichiaModel,
    output_path: Path,
    objective_reaction: str,
    optimization_sense: LpObjectiveSense = "Maximize",
    metabolic_enzymes: MetabolicEnzymeData | None = None,
    secretory_enzymes: SecretoryEnzymeData | None = None,
    combined_enzymes: CombinedEnzymeData | None = None,
    mu: float | None = None,
    total_protein_content: float = 0.37,
    unmodeled_er_protein_fraction: float | None = None,
    include_mitochondrial_constraint: bool = False,
    mitochondrial_protein_fraction: float = 0.05,
    include_proteasome_constraint: bool = False,
    include_ribosome_assembly_constraint: bool = False,
) -> LpWriteSummary:
    if objective_reaction not in model.reaction_index:
        raise KeyError(f"Reaction not found: {objective_reaction}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    objective_index = model.reaction_index[objective_reaction] + 1
    csr = model.s_matrix.tocsr()
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{optimization_sense}\n")
        handle.write(f"obj: X{objective_index}\n")
        handle.write("Subject To\n")
        for row_index in range(csr.shape[0]):
            start = csr.indptr[row_index]
            end = csr.indptr[row_index + 1]
            columns = csr.indices[start:end]
            values = csr.data[start:end]
            if len(columns) > 1:
                order = np.argsort(columns)
                columns = columns[order]
                values = values[order]
            expression = _format_linear_expression(values, columns)
            rhs = float(model.b[row_index]) if len(model.b) > row_index else 0.0
            handle.write(f"C{row_index + 1}: {expression} = {_format_number(rhs)}\n")
        metabolic_count = 0
        if metabolic_enzymes is not None:
            if mu is None or mu <= 0:
                raise ValueError("mu must be a positive number when metabolic enzyme coupling is enabled.")
            metabolic_count = _write_metabolic_coupling_constraints(handle, model, metabolic_enzymes, mu)
        secretory_count = 0
        if secretory_enzymes is not None:
            if mu is None or mu <= 0:
                raise ValueError("mu must be a positive number when secretory enzyme coupling is enabled.")
            secretory_count = _write_secretory_coupling_constraints(
                handle,
                model,
                secretory_enzymes,
                mu,
                start_index=metabolic_count + 1,
            )
        protein_mass_count = 0
        if combined_enzymes is not None:
            if mu is None or mu <= 0:
                raise ValueError("mu must be a positive number when protein mass constraints are enabled.")
            protein_mass_count = _write_protein_mass_constraints(
                handle,
                model,
                combined_enzymes,
                mu=mu,
                total_protein_content=total_protein_content,
                unmodeled_er_protein_fraction=(
                    unmodeled_er_protein_fraction
                    if unmodeled_er_protein_fraction is not None
                    else total_protein_content * 0.04
                ),
                start_index=metabolic_count + secretory_count + 1,
            )
        mitochondrial_count = 0
        if include_mitochondrial_constraint:
            if combined_enzymes is None:
                raise ValueError("combined_enzymes is required when mitochondrial constraint is enabled.")
            if mu is None or mu <= 0:
                raise ValueError("mu must be a positive number when mitochondrial constraint is enabled.")
            mitochondrial_count = _write_mitochondrial_constraint(
                handle,
                model,
                combined_enzymes,
                mu=mu,
                mitochondrial_protein_fraction=mitochondrial_protein_fraction,
            )
        proteasome_count = 0
        if include_proteasome_constraint:
            if combined_enzymes is None:
                raise ValueError("combined_enzymes is required when proteasome constraint is enabled.")
            if mu is None or mu <= 0:
                raise ValueError("mu must be a positive number when proteasome constraint is enabled.")
            proteasome_count = _write_proteasome_constraint(
                handle,
                model,
                combined_enzymes,
                mu=mu,
                start_index=metabolic_count + secretory_count + protein_mass_count + 1,
            )
        ribosome_assembly_count = 0
        if include_ribosome_assembly_constraint:
            if combined_enzymes is None:
                raise ValueError("combined_enzymes is required when ribosome assembly constraint is enabled.")
            if mu is None or mu <= 0:
                raise ValueError("mu must be a positive number when ribosome assembly constraint is enabled.")
            ribosome_assembly_count = _write_ribosome_assembly_constraint(
                handle,
                model,
                combined_enzymes,
                mu=mu,
                start_index=metabolic_count + secretory_count + 6,
            )
        handle.write("Bounds\n")
        for index, (lower, upper) in enumerate(zip(model.lb, model.ub), start=1):
            lower_text = _format_bound(float(lower))
            if float(upper) >= 100:
                upper_text = "+infinity"
            else:
                upper_text = _format_bound(float(upper))
            handle.write(f"{lower_text} <= X{index} <= {upper_text}\n")
        handle.write("End\n")
    return LpWriteSummary(
        path=output_path,
        objective_reaction=objective_reaction,
        objective_variable=f"X{objective_index}",
        optimization_sense=optimization_sense,
        stoichiometric_constraint_count=csr.shape[0],
        metabolic_coupling_constraint_count=metabolic_count,
        secretory_coupling_constraint_count=secretory_count,
        protein_mass_constraint_count=protein_mass_count,
        mitochondrial_constraint_count=mitochondrial_count,
        proteasome_constraint_count=proteasome_count,
        ribosome_assembly_constraint_count=ribosome_assembly_count,
        bounds_count=len(model.rxns),
        variable_count=len(model.rxns),
    )


def _write_metabolic_coupling_constraints(handle, model: PichiaModel, metabolic_enzymes: MetabolicEnzymeData, mu: float) -> int:
    reaction_index = model.reaction_index
    for index, (enzyme_id, kcat) in enumerate(zip(metabolic_enzymes.enzymes, metabolic_enzymes.kcat), start=1):
        reaction_id = metabolic_enzymes.reaction_id_for_enzyme(enzyme_id)
        formation_id = metabolic_enzymes.formation_reaction_id_for_enzyme(enzyme_id)
        if reaction_id not in reaction_index:
            raise KeyError(f"Metabolic reaction not found for enzyme {enzyme_id}: {reaction_id}")
        if formation_id not in reaction_index:
            raise KeyError(f"Formation reaction not found for enzyme {enzyme_id}: {formation_id}")
        reaction_var = reaction_index[reaction_id] + 1
        formation_var = reaction_index[formation_id] + 1
        coefficient = float(kcat) / mu
        handle.write(f"CM{index}: X{reaction_var} - {_format_number(coefficient)} X{formation_var} = 0\n")
    return metabolic_enzymes.enzyme_count


def _write_secretory_coupling_constraints(
    handle,
    model: PichiaModel,
    secretory_enzymes: SecretoryEnzymeData,
    mu: float,
    start_index: int,
) -> int:
    reaction_index = model.reaction_index
    written = 0
    for offset, entry in enumerate(secretory_enzymes.unique_complex_entries()):
        constraint_index = start_index + offset
        formation_id = secretory_enzymes.formation_reaction_id_for_complex(entry.complex_id)
        if formation_id not in reaction_index:
            raise KeyError(f"Formation reaction not found for secretory complex {entry.complex_id}: {formation_id}")
        matched_columns: list[int] = []
        matched_coefficients: list[float] = []
        suffix = f"_{entry.complex_id}"
        for zero_based_index, reaction_id in enumerate(model.rxns):
            if not reaction_id.endswith(suffix):
                continue
            if reaction_id not in secretory_enzymes.reaction_coefficients:
                continue
            matched_columns.append(zero_based_index)
            matched_coefficients.append(secretory_enzymes.reaction_coefficients[reaction_id])
        if not matched_columns:
            continue
        expression = _format_linear_expression(matched_coefficients, matched_columns)
        formation_var = reaction_index[formation_id] + 1
        coefficient = _matlab_divide_by_mu(entry.kcat, mu)
        handle.write(f"CM{constraint_index}: {expression} - {_format_number(coefficient)} X{formation_var} = 0\n")
        written += 1
    return written


def _write_protein_mass_constraints(
    handle,
    model: PichiaModel,
    combined_enzymes: CombinedEnzymeData,
    mu: float,
    total_protein_content: float,
    unmodeled_er_protein_fraction: float,
    start_index: int,
) -> int:
    dilution_columns: list[int] = []
    dilution_coefficients: list[float] = []
    for zero_based_index, reaction_id in enumerate(model.rxns):
        if "_dilution" not in reaction_id or "dummy" in reaction_id:
            continue
        molecular_weight = combined_enzymes.molecular_weight_for_dilution_reaction(reaction_id)
        dilution_columns.append(zero_based_index)
        dilution_coefficients.append(molecular_weight / 1000.0)

    dummy_index = _required_reaction_index(model, "dilute_dummy")
    dummy_er_index = _required_reaction_index(model, "dilute_dummyER")
    dilution_columns.extend([dummy_index, dummy_er_index])
    dilution_coefficients.extend([40.0, 40.0])

    modeled_fraction = model.modeled_protein_fraction("BIOMASS", "PROTEIN[c]")
    total_modeled_protein = total_protein_content * modeled_fraction
    expression = _format_linear_expression(dilution_coefficients, dilution_columns)
    handle.write(f"CM{start_index}: {expression} = {_format_number(mu * total_modeled_protein)}\n")
    handle.write(
        f"CM{start_index + 1}: {_format_number(40.0)} X{dummy_er_index + 1} = "
        f"{_format_number(mu * unmodeled_er_protein_fraction)}\n"
    )
    return 2


def _write_mitochondrial_constraint(
    handle,
    model: PichiaModel,
    combined_enzymes: CombinedEnzymeData,
    mu: float,
    mitochondrial_protein_fraction: float,
) -> int:
    columns: list[int] = []
    coefficients: list[float] = []
    for reaction_id in collect_compartment_reactions(model, ("m", "mm")):
        dilution_reaction_id = f"{reaction_id}_complex_dilution"
        if dilution_reaction_id not in model.reaction_index:
            continue
        enzyme_id = dilution_reaction_id.replace("_dilution", "")
        molecular_weight = combined_enzymes.exact_enzyme_mw(enzyme_id)
        columns.append(model.reaction_index[dilution_reaction_id])
        coefficients.append(molecular_weight / 1000.0)
    expression = _format_linear_expression(coefficients, columns)
    handle.write(f"Cmito: {expression} <= {_format_number(mu * mitochondrial_protein_fraction)}\n")
    return 1


def _write_proteasome_constraint(
    handle,
    model: PichiaModel,
    combined_enzymes: CombinedEnzymeData,
    mu: float,
    start_index: int,
) -> int:
    columns: list[int] = []
    coefficients: list[float] = []
    for zero_based_index, reaction_id in enumerate(model.rxns):
        if reaction_id.endswith("_subunit_degradation"):
            protein_id = reaction_id.replace("_subunit_degradation", "").replace("r_", "")
            protein_length = combined_enzymes.exact_protein_length(protein_id)
            coefficient = protein_length / 467.0
        elif reaction_id.endswith("_sp_degradation"):
            coefficient = 25.0 / 467.0
        else:
            continue
        columns.append(zero_based_index)
        coefficients.append(coefficient)

    if not columns:
        raise ValueError("No proteasome degradation reactions were found in the model.")

    formation_index = _required_reaction_index(model, "Mach_proteasome_complex_formation")
    formation_coefficient = combined_enzymes.exact_enzyme_kcat("Mach_proteasome_complex") / mu
    expression = _format_linear_expression(coefficients, columns)
    handle.write(
        f"CM{start_index}: {expression} - {_format_number(formation_coefficient)} "
        f"X{formation_index + 1} = 0\n"
    )
    return 1


def _write_ribosome_assembly_constraint(
    handle,
    model: PichiaModel,
    combined_enzymes: CombinedEnzymeData,
    mu: float,
    start_index: int,
) -> int:
    ribosome_formation_index = _required_reaction_index(model, "Mach_Ribosome_complex_formation")
    assembly_formation_index = _required_reaction_index(model, "Mach_Ribosome_Assembly_Factors_complex_formation")
    assembly_coefficient = combined_enzymes.exact_enzyme_kcat("Mach_Ribosome_Assembly_Factors_complex") / mu
    handle.write(
        f"CM{start_index}: X{ribosome_formation_index + 1} - "
        f"{_format_number(assembly_coefficient)} X{assembly_formation_index + 1} = 0\n"
    )
    return 1


def _required_reaction_index(model: PichiaModel, reaction_id: str) -> int:
    try:
        return model.reaction_index[reaction_id]
    except KeyError as exc:
        raise KeyError(f"Reaction not found: {reaction_id}") from exc


def _format_linear_expression(values, columns) -> str:
    if len(columns) == 0:
        return "0"
    parts: list[str] = []
    for position, (value, column) in enumerate(zip(values, columns)):
        coefficient = float(value)
        variable = f"X{int(column) + 1}"
        if position == 0:
            parts.append(f"{_format_number(coefficient)} {variable}")
        elif coefficient >= 0:
            parts.append(f"+ {_format_number(coefficient)} {variable}")
        else:
            parts.append(f"{_format_number(coefficient)} {variable}")
    return _wrap_expression_terms(parts)


def _wrap_expression_terms(parts: list[str], terms_per_line: int = 150) -> str:
    if terms_per_line <= 0:
        raise ValueError("terms_per_line must be positive.")
    wrapped: list[str] = []
    for index, part in enumerate(parts):
        if index > 0 and index % terms_per_line == 0:
            wrapped.append("\n ")
        elif index > 0:
            wrapped.append(" ")
        wrapped.append(part)
    return "".join(wrapped)


def _format_number(value: float) -> str:
    return f"{value:.15f}"


def _matlab_divide_by_mu(value: float, mu: float) -> float:
    # MATLAB writeLPGlc emits secretory kcat/mu coefficients with %.15f.
    # For the fixed-mu smoke cases, multiplying by reciprocal(mu) matches
    # MATLAB's emitted secretory coefficients more closely than direct division.
    return float(value) * (1.0 / float(mu))


def _format_bound(value: float) -> str:
    return f"{value:.6f}"
