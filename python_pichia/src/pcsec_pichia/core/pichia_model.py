from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Literal

import numpy as np
from scipy import sparse


BoundType = Literal["lower", "upper", "both"]


@dataclass(frozen=True)
class PichiaModelSummary:
    model_id: str
    source_file: Path
    reaction_count: int
    metabolite_count: int
    gene_count: int
    stoichiometric_shape: tuple[int, int]
    rxn_gene_shape: tuple[int, int] | None
    key_reactions: dict[str, bool]


@dataclass(frozen=True)
class PichiaEnzymeDataSummary:
    source_file: Path
    variable_name: str
    enzyme_count: int
    protein_count: int
    has_subunit_matrix: bool
    has_kcat: bool
    fields: tuple[str, ...]


@dataclass(frozen=True)
class PichiaModel:
    model_id: str
    source_file: Path
    rxns: list[str]
    mets: list[str]
    genes: list[str]
    lb: np.ndarray
    ub: np.ndarray
    c: np.ndarray
    b: np.ndarray
    s_matrix: sparse.spmatrix | sparse.sparray
    rules: list[str]
    gr_rules: list[str]
    rxn_gene_mat: sparse.spmatrix | sparse.sparray | None = None

    @property
    def reaction_index(self) -> dict[str, int]:
        return {rxn: index for index, rxn in enumerate(self.rxns)}

    @property
    def metabolite_index(self) -> dict[str, int]:
        return {met: index for index, met in enumerate(self.mets)}

    @property
    def gene_index(self) -> dict[str, int]:
        return {gene: index for index, gene in enumerate(self.genes)}

    def change_rxn_bounds(self, reaction_id: str, lower: float | None = None, upper: float | None = None) -> "PichiaModel":
        try:
            index = self.reaction_index[reaction_id]
        except KeyError as exc:
            raise KeyError(f"Reaction not found: {reaction_id}") from exc
        lb = np.array(self.lb, dtype=float, copy=True)
        ub = np.array(self.ub, dtype=float, copy=True)
        if lower is not None:
            lb[index] = lower
        if upper is not None:
            ub[index] = upper
        return replace(self, lb=lb, ub=ub)

    def with_reaction_bounds(self, changes: dict[str, tuple[float | None, float | None]]) -> "PichiaModel":
        reaction_index = self.reaction_index
        lb = np.array(self.lb, dtype=float, copy=True)
        ub = np.array(self.ub, dtype=float, copy=True)
        for reaction_id, (lower, upper) in changes.items():
            try:
                index = reaction_index[reaction_id]
            except KeyError as exc:
                raise KeyError(f"Reaction not found: {reaction_id}") from exc
            if lower is not None:
                lb[index] = lower
            if upper is not None:
                ub[index] = upper
        return replace(self, lb=lb, ub=ub)

    def set_bound(self, reaction_id: str, value: float, bound_type: BoundType = "both") -> "PichiaModel":
        if bound_type == "lower":
            return self.change_rxn_bounds(reaction_id, lower=value)
        if bound_type == "upper":
            return self.change_rxn_bounds(reaction_id, upper=value)
        return self.change_rxn_bounds(reaction_id, lower=value, upper=value)

    def set_objective(self, reaction_id: str, coefficient: float = 1.0) -> "PichiaModel":
        try:
            index = self.reaction_index[reaction_id]
        except KeyError as exc:
            raise KeyError(f"Reaction not found: {reaction_id}") from exc
        c = np.zeros_like(np.array(self.c, dtype=float), dtype=float)
        c[index] = coefficient
        return replace(self, c=c)

    def add_reaction(
        self,
        reaction_id: str,
        stoichiometry: dict[str, float],
        lower_bound: float = 0.0,
        upper_bound: float = 1000.0,
        objective_coefficient: float = 0.0,
        rule: str = "",
        gr_rule: str = "",
    ) -> "PichiaModel":
        if reaction_id in self.reaction_index:
            raise ValueError(f"Reaction already exists: {reaction_id}")

        mets = list(self.mets)
        met_index = {metabolite_id: index for index, metabolite_id in enumerate(mets)}
        new_metabolite_count = 0
        for metabolite_id in stoichiometry:
            if metabolite_id in met_index:
                continue
            met_index[metabolite_id] = len(mets)
            mets.append(metabolite_id)
            new_metabolite_count += 1

        matrix = self.s_matrix.tocsr()
        if new_metabolite_count:
            matrix = sparse.vstack(
                [matrix, sparse.csr_matrix((new_metabolite_count, matrix.shape[1]))],
                format="csr",
            )
        new_column = sparse.lil_matrix((len(mets), 1), dtype=float)
        for metabolite_id, coefficient in stoichiometry.items():
            if coefficient:
                new_column[met_index[metabolite_id], 0] = float(coefficient)
        s_matrix = sparse.hstack([matrix, new_column.tocsr()], format="csc")

        b = np.array(self.b, dtype=float, copy=True)
        if new_metabolite_count:
            b = np.concatenate([b, np.zeros(new_metabolite_count, dtype=float)])

        rxn_gene_mat = self.rxn_gene_mat
        if rxn_gene_mat is not None:
            rxn_gene_mat = sparse.vstack(
                [rxn_gene_mat.tocsr(), sparse.csr_matrix((1, rxn_gene_mat.shape[1]))],
                format="csr",
            )

        return replace(
            self,
            rxns=[*self.rxns, reaction_id],
            mets=mets,
            lb=np.concatenate([np.array(self.lb, dtype=float, copy=True), np.array([lower_bound], dtype=float)]),
            ub=np.concatenate([np.array(self.ub, dtype=float, copy=True), np.array([upper_bound], dtype=float)]),
            c=np.concatenate([np.array(self.c, dtype=float, copy=True), np.array([objective_coefficient], dtype=float)]),
            b=b,
            s_matrix=s_matrix,
            rules=[*self.rules, rule],
            gr_rules=[*self.gr_rules, gr_rule],
            rxn_gene_mat=rxn_gene_mat,
        )

    def summary(self) -> PichiaModelSummary:
        rxn_gene_shape = tuple(self.rxn_gene_mat.shape) if self.rxn_gene_mat is not None else None
        return PichiaModelSummary(
            model_id=self.model_id,
            source_file=self.source_file,
            reaction_count=len(self.rxns),
            metabolite_count=len(self.mets),
            gene_count=len(self.genes),
            stoichiometric_shape=tuple(self.s_matrix.shape),
            rxn_gene_shape=rxn_gene_shape,
            key_reactions={reaction_id: reaction_id in self.reaction_index for reaction_id in default_key_reactions()},
        )

    def modeled_protein_fraction(self, biomass_reaction_id: str = "BIOMASS", protein_metabolite_id: str = "PROTEIN[c]") -> float:
        try:
            met_index = self.metabolite_index[protein_metabolite_id]
            rxn_index = self.reaction_index[biomass_reaction_id]
        except KeyError as exc:
            raise KeyError(
                f"Could not locate biomass reaction {biomass_reaction_id!r} or protein metabolite {protein_metabolite_id!r}."
            ) from exc
        return 1.0 + float(self.s_matrix[met_index, rxn_index])


def default_key_reactions() -> tuple[str, ...]:
    return (
        "BIOMASS",
        "BIOMASS_glyc",
        "BIOMASS_meoh",
        "Ex_glc_D",
        "Ex_o2",
        "Ex_glyc",
        "Ex_meoh",
    )


def collect_compartment_reactions(model: PichiaModel, compartment_ids: str | tuple[str, ...] | list[str]) -> list[str]:
    if isinstance(compartment_ids, str):
        compartment_ids = (compartment_ids,)
    tokens = tuple(f"[{compartment_id}]" for compartment_id in compartment_ids)
    met_indices = [index for index, metabolite_id in enumerate(model.mets) if any(token in metabolite_id for token in tokens)]
    if not met_indices:
        return []
    submatrix = model.s_matrix[met_indices, :].tocsc()
    has_stoichiometry = np.diff(submatrix.indptr) > 0
    rules = model.rules or [""] * len(model.rxns)
    return [
        reaction_id
        for index, reaction_id in enumerate(model.rxns)
        if has_stoichiometry[index] and index < len(rules) and str(rules[index]).strip()
    ]


def summary_to_dict(summary: PichiaModelSummary | PichiaEnzymeDataSummary) -> dict[str, Any]:
    data = dict(summary.__dict__)
    data["source_file"] = str(data["source_file"])
    return data
