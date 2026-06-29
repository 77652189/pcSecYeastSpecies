from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.io import loadmat

from pcsec_pichia.core.pichia_enzymes import CombinedEnzymeData, MetabolicEnzymeData, SecretoryEnzymeData
from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.core.pichia_model import (
    PichiaEnzymeDataSummary,
    PichiaModel,
    PichiaModelSummary,
    summary_to_dict,
)


@dataclass(frozen=True)
class MatStructLoader:
    paths: ProjectPaths

    def load_pcsec_pichia_model(self) -> PichiaModel:
        path = self.paths.pichia_model_mat
        struct = _load_single_struct(path, "model")
        return PichiaModel(
            model_id=str(getattr(struct, "modelID", "pcSecPichia")),
            source_file=path,
            rxns=_string_list(getattr(struct, "rxns")),
            mets=_string_list(getattr(struct, "mets")),
            genes=_string_list(getattr(struct, "genes")),
            lb=_numeric_vector(getattr(struct, "lb")),
            ub=_numeric_vector(getattr(struct, "ub")),
            c=_numeric_vector(getattr(struct, "c")),
            b=_numeric_vector(getattr(struct, "b")),
            s_matrix=_sparse_matrix(getattr(struct, "S")),
            rules=_string_list(getattr(struct, "rules")),
            gr_rules=_string_list(getattr(struct, "grRules")),
            rxn_gene_mat=_optional_sparse_matrix(getattr(struct, "rxnGeneMat", None)),
        )

    def load_pcsec_pichia_summary(self) -> PichiaModelSummary:
        return self.load_pcsec_pichia_model().summary()

    def load_metabolic_enzymedata(self) -> MetabolicEnzymeData:
        path = self.paths.pichia_enzymedata_metabolic_mat
        struct = _load_single_struct(path, "enzymedata")
        return MetabolicEnzymeData(
            source_file=path,
            enzymes=_string_list(getattr(struct, "enzyme")),
            kcat=_numeric_vector(getattr(struct, "kcat")).astype(float),
        )

    def load_secretory_enzymedata(self) -> SecretoryEnzymeData:
        sec_path = self.paths.pichia_enzymedata_sec_mat
        met_path = self.paths.pichia_enzymedata_metabolic_mat
        dummy_path = self.paths.pichia_enzymedata_dummy_er_mat
        sec = _load_single_struct(sec_path, "enzymedataSEC")
        met = _load_single_struct(met_path, "enzymedata")
        dummy = _load_single_struct(dummy_path, "enzymedataDummyER")
        coefficient_rxns = _string_list(getattr(met, "rxns")) + _string_list(getattr(dummy, "rxns"))
        coefficient_values = np.concatenate(
            [
                _numeric_vector(getattr(met, "rxnscoef")).astype(float),
                _numeric_vector(getattr(dummy, "rxnscoef")).astype(float),
            ]
        )
        reaction_coefficients: dict[str, float] = {}
        for reaction_id, value in zip(coefficient_rxns, coefficient_values):
            if reaction_id not in reaction_coefficients:
                reaction_coefficients[reaction_id] = float(value)
        return SecretoryEnzymeData(
            source_file=sec_path,
            reaction_coefficient_sources=(met_path, dummy_path),
            complexes=_string_list(getattr(sec, "enzyme")),
            compartments=_string_list(getattr(sec, "comp")),
            kcat=_numeric_vector(getattr(sec, "kcat")).astype(float),
            coefficient_refs=_string_list(getattr(sec, "coefref")),
            reaction_coefficients=reaction_coefficients,
        )

    def load_combined_enzymedata(self) -> CombinedEnzymeData:
        met_path = self.paths.pichia_enzymedata_metabolic_mat
        sec_path = self.paths.pichia_enzymedata_sec_mat
        machine_path = self.paths.pichia_enzymedata_machine_mat
        structs = [
            _load_single_struct(met_path, "enzymedata"),
            _load_single_struct(sec_path, "enzymedataSEC"),
            _load_single_struct(machine_path, "enzymedataMachine"),
        ]
        enzymes: list[str] = []
        kcat_arrays: list[np.ndarray] = []
        enzyme_mw_arrays: list[np.ndarray] = []
        proteins_all: list[str] = []
        protein_length_all: list[float] = []
        protein_mw_all: list[float] = []
        for struct in structs:
            enzymes.extend(_string_list(getattr(struct, "enzyme")))
            kcat_arrays.append(_numeric_vector(getattr(struct, "kcat")).astype(float))
            enzyme_mw_arrays.append(_numeric_vector(getattr(struct, "enzyme_MW")).astype(float))
            proteins = _string_list(getattr(struct, "proteins"))
            protein_length = _numeric_vector(getattr(struct, "proteinLength")).astype(float)
            protein_mw = _numeric_vector(getattr(struct, "proteinMWs")).astype(float)
            proteins_all.extend(proteins)
            protein_length_all.extend(float(value) for value in protein_length)
            protein_mw_all.extend(float(value) for value in protein_mw)

        unique_proteins: list[str] = []
        unique_protein_length: list[float] = []
        unique_protein_mw: list[float] = []
        seen: set[str] = set()
        for protein_id, protein_length, molecular_weight in zip(proteins_all, protein_length_all, protein_mw_all):
            if protein_id in seen:
                continue
            seen.add(protein_id)
            unique_proteins.append(protein_id)
            unique_protein_length.append(protein_length)
            unique_protein_mw.append(molecular_weight)

        return CombinedEnzymeData(
            source_files=(met_path, sec_path, machine_path),
            enzymes=enzymes,
            kcat=np.concatenate(kcat_arrays),
            enzyme_mw=np.concatenate(enzyme_mw_arrays),
            proteins=unique_proteins,
            protein_length=np.array(unique_protein_length, dtype=float),
            protein_mw=np.array(unique_protein_mw, dtype=float),
        )

    def load_pichia_enzymedata_summaries(self) -> list[PichiaEnzymeDataSummary]:
        files = [
            (self.paths.pichia_enzymedata_sec_mat, "enzymedataSEC"),
            (self.paths.pichia_enzymedata_machine_mat, "enzymedataMachine"),
        ]
        return [self.load_enzymedata_summary(path, variable_name) for path, variable_name in files]

    def load_enzymedata_summary(self, path: Path, variable_name: str) -> PichiaEnzymeDataSummary:
        struct = _load_single_struct(path, variable_name)
        fields = tuple(getattr(struct, "_fieldnames", []) or [])
        enzymes = _string_list(getattr(struct, "enzyme", []))
        proteins = _string_list(getattr(struct, "proteins", []))
        return PichiaEnzymeDataSummary(
            source_file=path,
            variable_name=variable_name,
            enzyme_count=len(enzymes),
            protein_count=len(proteins),
            has_subunit_matrix=hasattr(struct, "subunit") and hasattr(struct, "subunit_stoichiometry"),
            has_kcat=hasattr(struct, "kcat"),
            fields=fields,
        )

    def write_model_summary_json(self, output_path: Path) -> Path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        summary = self.load_pcsec_pichia_summary()
        enzyme_summaries = self.load_pichia_enzymedata_summaries()
        payload = {
            "model": summary_to_dict(summary),
            "enzymedata": [summary_to_dict(item) for item in enzyme_summaries],
        }
        output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        return output_path


def _load_single_struct(path: Path, variable_name: str) -> Any:
    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    if variable_name not in data:
        available = ", ".join(sorted(key for key in data if not key.startswith("__")))
        raise KeyError(f"{path} does not contain MATLAB variable {variable_name!r}. Available: {available}")
    return data[variable_name]


def _string_list(value: Any) -> list[str]:
    array = np.asarray(value, dtype=object)
    if array.shape == ():
        text = _string_scalar(array.item())
        return [text] if text else []
    return [_string_scalar(item) for item in array.reshape(-1)]


def _string_scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _string_scalar(value.item())
        if value.dtype.kind in {"U", "S"}:
            return "".join(str(item) for item in value.reshape(-1)).strip()
        return " ".join(_string_scalar(item) for item in value.reshape(-1)).strip()
    return str(value)


def _numeric_vector(value: Any) -> np.ndarray:
    return np.asarray(value).reshape(-1)


def _sparse_matrix(value: Any) -> sparse.spmatrix | sparse.sparray:
    if sparse.issparse(value):
        return value
    return sparse.csc_matrix(value)


def _optional_sparse_matrix(value: Any) -> sparse.spmatrix | sparse.sparray | None:
    if value is None:
        return None
    return _sparse_matrix(value)
