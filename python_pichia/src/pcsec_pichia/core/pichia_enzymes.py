from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class MetabolicEnzymeData:
    source_file: Path
    enzymes: list[str]
    kcat: np.ndarray

    @property
    def enzyme_count(self) -> int:
        return len(self.enzymes)

    def reaction_id_for_enzyme(self, enzyme_id: str) -> str:
        return enzyme_id.replace("_complex", "")

    def formation_reaction_id_for_enzyme(self, enzyme_id: str) -> str:
        return f"{enzyme_id}_formation"


@dataclass(frozen=True)
class SecretoryComplexEntry:
    complex_id: str
    compartment: str
    kcat: float


@dataclass(frozen=True)
class SecretoryEnzymeData:
    source_file: Path
    reaction_coefficient_sources: tuple[Path, ...]
    complexes: list[str]
    compartments: list[str]
    kcat: np.ndarray
    coefficient_refs: list[str]
    reaction_coefficients: dict[str, float]

    @property
    def complex_count(self) -> int:
        return len(self.unique_complex_entries())

    def unique_complex_entries(self) -> list[SecretoryComplexEntry]:
        first_index_by_key: dict[str, int] = {}
        for index, (complex_id, compartment) in enumerate(zip(self.complexes, self.compartments)):
            key = f"{complex_id}{compartment}"
            if key not in first_index_by_key:
                first_index_by_key[key] = index
        return [
            SecretoryComplexEntry(
                complex_id=self.complexes[index],
                compartment=self.compartments[index],
                kcat=float(self.kcat[index]),
            )
            for _key, index in sorted(first_index_by_key.items())
        ]

    def formation_reaction_id_for_complex(self, complex_id: str) -> str:
        return f"{complex_id}_formation"

    def with_reaction_coefficients(self, extra_coefficients: dict[str, float]) -> "SecretoryEnzymeData":
        merged = dict(self.reaction_coefficients)
        merged.update({reaction_id: float(value) for reaction_id, value in extra_coefficients.items()})
        return replace(self, reaction_coefficients=merged)


@dataclass(frozen=True)
class CombinedEnzymeData:
    source_files: tuple[Path, ...]
    enzymes: list[str]
    kcat: np.ndarray
    enzyme_mw: np.ndarray
    proteins: list[str]
    protein_length: np.ndarray
    protein_mw: np.ndarray

    @property
    def enzyme_count(self) -> int:
        return len(self.enzymes)

    @property
    def protein_count(self) -> int:
        return len(self.proteins)

    def molecular_weight_for_dilution_reaction(self, reaction_id: str) -> float:
        component_name = reaction_id.replace("_dilution", "")
        molecular_weight = self._first_enzyme_mw_containing(component_name)
        if molecular_weight is not None:
            return molecular_weight
        molecular_weight = self._first_enzyme_mw_containing(component_name.replace("_complex", ""))
        if molecular_weight is not None:
            return molecular_weight
        protein_name = component_name.split("_misfolding", 1)[0]
        molecular_weight = self._first_protein_mw_containing(protein_name)
        if molecular_weight is None:
            raise KeyError(f"Could not resolve molecular weight for dilution reaction: {reaction_id}")
        return molecular_weight

    def exact_enzyme_mw(self, enzyme_id: str) -> float:
        for known_enzyme_id, molecular_weight in zip(self.enzymes, self.enzyme_mw):
            if known_enzyme_id == enzyme_id:
                return float(molecular_weight)
        raise KeyError(f"Could not resolve molecular weight for enzyme: {enzyme_id}")

    def exact_enzyme_kcat(self, enzyme_id: str) -> float:
        for known_enzyme_id, kcat in zip(self.enzymes, self.kcat):
            if known_enzyme_id == enzyme_id:
                return float(kcat)
        raise KeyError(f"Could not resolve kcat for enzyme: {enzyme_id}")

    def exact_protein_length(self, protein_id: str) -> float:
        for known_protein_id, protein_length in zip(self.proteins, self.protein_length):
            if known_protein_id == protein_id:
                return float(protein_length)
        raise KeyError(f"Could not resolve protein length for protein: {protein_id}")

    def with_target_proteins(self, target_enzymedata: "TargetProteinEnzymeData") -> "CombinedEnzymeData":
        proteins = [*self.proteins, *target_enzymedata.proteins]
        protein_length = np.concatenate([self.protein_length, target_enzymedata.protein_length])
        protein_mw = np.concatenate([self.protein_mw, target_enzymedata.protein_mw])

        unique_proteins: list[str] = []
        unique_protein_length: list[float] = []
        unique_protein_mw: list[float] = []
        seen: set[str] = set()
        for protein_id, length, molecular_weight in zip(proteins, protein_length, protein_mw):
            if protein_id in seen:
                continue
            seen.add(protein_id)
            unique_proteins.append(protein_id)
            unique_protein_length.append(float(length))
            unique_protein_mw.append(float(molecular_weight))

        return replace(
            self,
            proteins=unique_proteins,
            protein_length=np.array(unique_protein_length, dtype=float),
            protein_mw=np.array(unique_protein_mw, dtype=float),
        )

    def _first_enzyme_mw_containing(self, text: str) -> float | None:
        for enzyme_id, molecular_weight in zip(self.enzymes, self.enzyme_mw):
            if text and text in enzyme_id:
                return float(molecular_weight)
        return None

    def _first_protein_mw_containing(self, text: str) -> float | None:
        for protein_id, molecular_weight in zip(self.proteins, self.protein_mw):
            if text and text in protein_id:
                return float(molecular_weight)
        return None


@dataclass(frozen=True)
class TargetProteinEnzymeData:
    proteins: list[str]
    protein_mw: np.ndarray
    protein_length: np.ndarray
    protein_pst: np.ndarray
    protein_pst_info: tuple[str, ...]
    protein_extra_mw_specific: np.ndarray
    protein_extra_mw: np.ndarray
    protein_loc: list[str]
    kdeg: np.ndarray
    reaction_coefficients: dict[str, float] = field(default_factory=dict)

    @property
    def protein_count(self) -> int:
        return len(self.proteins)

    def index_for_protein(self, protein_id: str) -> int:
        try:
            return self.proteins.index(protein_id)
        except ValueError as exc:
            raise KeyError(f"Could not resolve target protein: {protein_id}") from exc

    def pst_value(self, protein_id: str, pst_name: str) -> float:
        if pst_name not in self.protein_pst_info:
            raise KeyError(f"Unknown PTM field: {pst_name}")
        row_index = self.index_for_protein(protein_id)
        column_index = self.protein_pst_info.index(pst_name)
        return float(self.protein_pst[row_index, column_index])
