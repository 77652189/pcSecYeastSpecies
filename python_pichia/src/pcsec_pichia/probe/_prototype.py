from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import loadmat
from scipy.optimize import linprog


AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
RULE_TOKEN_PATTERN = re.compile(r"x\((\d+)\)")
DEFAULT_TRADEOFF_MUS = (0.05, 0.10, 0.20)
IMPROVEMENT_REL_TOLERANCE = 1e-6

AA_MW = {
    "A": 71.08,
    "C": 103.14,
    "D": 115.09,
    "E": 129.11,
    "F": 147.17,
    "G": 57.05,
    "H": 137.14,
    "I": 113.16,
    "K": 128.17,
    "L": 113.16,
    "M": 131.20,
    "N": 114.10,
    "P": 97.12,
    "Q": 128.13,
    "R": 156.19,
    "S": 87.08,
    "T": 101.10,
    "V": 99.13,
    "W": 186.21,
    "Y": 163.17,
}


@dataclass(frozen=True)
class CobraModel:
    source_file: str
    rxns: list[str]
    mets: list[str]
    genes: list[str]
    lb: np.ndarray
    ub: np.ndarray
    b: np.ndarray
    s_matrix: sparse.csc_matrix
    rules: list[str]
    gr_rules: list[str]

    @property
    def reaction_index(self) -> dict[str, int]:
        return {reaction_id: index for index, reaction_id in enumerate(self.rxns)}

    @property
    def gene_index(self) -> dict[str, int]:
        return {gene_id: index for index, gene_id in enumerate(self.genes)}

    def with_bounds(self, changes: dict[str, tuple[float | None, float | None]]) -> "CobraModel":
        reaction_index = self.reaction_index
        lb = self.lb.copy()
        ub = self.ub.copy()
        for reaction_id, (lower, upper) in changes.items():
            if reaction_id not in reaction_index:
                continue
            index = reaction_index[reaction_id]
            if lower is not None:
                lb[index] = lower
            if upper is not None:
                ub[index] = upper
        return CobraModel(
            source_file=self.source_file,
            rxns=self.rxns,
            mets=self.mets,
            genes=self.genes,
            lb=lb,
            ub=ub,
            b=self.b,
            s_matrix=self.s_matrix,
            rules=self.rules,
            gr_rules=self.gr_rules,
        )

    def add_reaction(
        self,
        reaction_id: str,
        stoichiometry: dict[str, float],
        lower_bound: float = 0.0,
        upper_bound: float = 1000.0,
        rule: str = "",
        gr_rule: str = "",
    ) -> "CobraModel":
        if reaction_id in self.reaction_index:
            raise ValueError(f"Reaction already exists: {reaction_id}")

        mets = list(self.mets)
        met_index = {metabolite_id: index for index, metabolite_id in enumerate(mets)}
        added_metabolites = 0
        for metabolite_id in stoichiometry:
            if metabolite_id in met_index:
                continue
            met_index[metabolite_id] = len(mets)
            mets.append(metabolite_id)
            added_metabolites += 1

        matrix = self.s_matrix.tocsr()
        if added_metabolites:
            matrix = sparse.vstack([matrix, sparse.csr_matrix((added_metabolites, matrix.shape[1]))], format="csr")
        new_column = sparse.lil_matrix((len(mets), 1), dtype=float)
        for metabolite_id, coefficient in stoichiometry.items():
            if coefficient:
                new_column[met_index[metabolite_id], 0] = float(coefficient)
        new_matrix = sparse.hstack([matrix, new_column.tocsr()], format="csc")

        b = self.b.copy()
        if added_metabolites:
            b = np.concatenate([b, np.zeros(added_metabolites, dtype=float)])

        return CobraModel(
            source_file=self.source_file,
            rxns=[*self.rxns, reaction_id],
            mets=mets,
            genes=self.genes,
            lb=np.concatenate([self.lb.copy(), np.array([lower_bound], dtype=float)]),
            ub=np.concatenate([self.ub.copy(), np.array([upper_bound], dtype=float)]),
            b=b,
            s_matrix=new_matrix,
            rules=[*self.rules, rule],
            gr_rules=[*self.gr_rules, gr_rule],
        )


@dataclass(frozen=True)
class TargetSpec:
    target_id: str
    protein_id: str
    mature_sequence: str
    leader_sequence: str
    signal_peptide_sequence: str
    through_er: bool
    localization: str
    disulfide_sites: int
    n_glycosylation_sites: int
    o_glycosylation_sites: int
    transmembrane: int
    gpi_sites: int
    cotranslation: int
    source: str

    @property
    def full_sequence(self) -> str:
        return self.leader_sequence + self.mature_sequence


@dataclass(frozen=True)
class SolveResult:
    objective: str
    status: str
    success: bool
    objective_value: float | None
    message: str
    fluxes: dict[str, float]
    sensitivity: dict[str, tuple[float, ...]] | None = None


@dataclass(frozen=True)
class TargetModelBuildResult:
    status: str
    supported: bool
    reason: str
    model: CobraModel | None
    exchange_reaction_id: str | None
    added_reaction_count: int
    added_metabolite_count: int


@dataclass(frozen=True)
class MetabolicEnzymeData:
    enzymes: list[str]
    kcat: np.ndarray


@dataclass(frozen=True)
class SecretoryComplexEntry:
    complex_id: str
    compartment: str
    kcat: float


@dataclass(frozen=True)
class SecretoryEnzymeData:
    complexes: list[str]
    compartments: list[str]
    kcat: np.ndarray
    coefficient_refs: list[str]
    reaction_coefficients: dict[str, float]
    complex_subunits: dict[str, list[dict[str, object]]]

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

    def with_reaction_coefficients(self, extra_coefficients: dict[str, float]) -> "SecretoryEnzymeData":
        merged = dict(self.reaction_coefficients)
        merged.update({reaction_id: float(value) for reaction_id, value in extra_coefficients.items()})
        return SecretoryEnzymeData(
            complexes=self.complexes,
            compartments=self.compartments,
            kcat=self.kcat,
            coefficient_refs=self.coefficient_refs,
            reaction_coefficients=merged,
            complex_subunits=self.complex_subunits,
        )

    def with_complex_kcat_multiplier(self, complex_id: str, factor: float) -> "SecretoryEnzymeData":
        kcat = self.kcat.copy()
        for index, known_id in enumerate(self.complexes):
            if known_id == complex_id:
                kcat[index] = kcat[index] * factor
        return SecretoryEnzymeData(
            complexes=self.complexes,
            compartments=self.compartments,
            kcat=kcat,
            coefficient_refs=self.coefficient_refs,
            reaction_coefficients=self.reaction_coefficients,
            complex_subunits=self.complex_subunits,
        )


@dataclass(frozen=True)
class CombinedEnzymeData:
    enzymes: list[str]
    kcat: np.ndarray
    enzyme_mw: np.ndarray
    proteins: list[str]
    protein_length: np.ndarray
    protein_mw: np.ndarray
    protein_kdeg: np.ndarray

    def with_target(self, target: "TargetProteinEnzymeData") -> "CombinedEnzymeData":
        proteins = [*self.proteins, target.protein_id]
        lengths = np.concatenate([self.protein_length, np.array([target.protein_length], dtype=float)])
        mws = np.concatenate([self.protein_mw, np.array([target.protein_mw], dtype=float)])
        kdegs = np.concatenate([self.protein_kdeg, np.array([0.0], dtype=float)])
        seen: set[str] = set()
        unique_proteins: list[str] = []
        unique_lengths: list[float] = []
        unique_mws: list[float] = []
        unique_kdegs: list[float] = []
        for protein_id, length, mw, kdeg in zip(proteins, lengths, mws, kdegs):
            if protein_id in seen:
                continue
            seen.add(protein_id)
            unique_proteins.append(protein_id)
            unique_lengths.append(float(length))
            unique_mws.append(float(mw))
            unique_kdegs.append(float(kdeg))
        return CombinedEnzymeData(
            enzymes=self.enzymes,
            kcat=self.kcat,
            enzyme_mw=self.enzyme_mw,
            proteins=unique_proteins,
            protein_length=np.array(unique_lengths, dtype=float),
            protein_mw=np.array(unique_mws, dtype=float),
            protein_kdeg=np.array(unique_kdegs, dtype=float),
        )

    def molecular_weight_for_dilution_reaction(self, reaction_id: str) -> float:
        component_name = reaction_id.replace("_dilution", "")
        for query in (component_name, component_name.replace("_complex", "")):
            for enzyme_id, mw in zip(self.enzymes, self.enzyme_mw):
                if query and query in enzyme_id:
                    return float(mw)
        protein_name = component_name.split("_misfolding", 1)[0]
        for known_protein_id, mw in zip(self.proteins, self.protein_mw):
            if protein_name and protein_name in known_protein_id:
                return float(mw)
        raise KeyError(f"Could not resolve molecular weight for dilution reaction: {reaction_id}")

    def exact_enzyme_mw(self, enzyme_id: str) -> float:
        for known_id, value in zip(self.enzymes, self.enzyme_mw):
            if known_id == enzyme_id:
                return float(value)
        raise KeyError(f"Could not resolve enzyme MW: {enzyme_id}")

    def exact_enzyme_kcat(self, enzyme_id: str) -> float:
        for known_id, value in zip(self.enzymes, self.kcat):
            if known_id == enzyme_id:
                return float(value)
        raise KeyError(f"Could not resolve enzyme kcat: {enzyme_id}")

    def with_enzyme_kcat_multiplier(self, enzyme_id: str, factor: float) -> "CombinedEnzymeData":
        kcat = self.kcat.copy()
        for index, known_id in enumerate(self.enzymes):
            if known_id == enzyme_id:
                kcat[index] = kcat[index] * factor
        return CombinedEnzymeData(
            enzymes=self.enzymes,
            kcat=kcat,
            enzyme_mw=self.enzyme_mw,
            proteins=self.proteins,
            protein_length=self.protein_length,
            protein_mw=self.protein_mw,
            protein_kdeg=self.protein_kdeg,
        )

    def exact_protein_length(self, protein_id: str) -> float:
        for known_id, value in zip(self.proteins, self.protein_length):
            if known_id == protein_id:
                return float(value)
        raise KeyError(f"Could not resolve protein length: {protein_id}")

    def exact_protein_kdeg(self, protein_id: str) -> float:
        for known_id, value in zip(self.proteins, self.protein_kdeg):
            if known_id == protein_id:
                return float(value)
        raise KeyError(f"Could not resolve protein kdeg: {protein_id}")


@dataclass(frozen=True)
class TargetProteinEnzymeData:
    protein_id: str
    protein_mw: float
    protein_length: float
    disulfides: float
    n_glycans: float
    o_glycans: float
    gpi_sites: float
    reaction_coefficients: dict[str, float]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def output_dir() -> Path:
    return repo_root() / "local_runs" / "pichia_hlf_opn_probe"


def string_list(value: object) -> list[str]:
    array = np.asarray(value, dtype=object)
    if array.shape == ():
        return [string_scalar(array.item())]
    return [string_scalar(item) for item in array.reshape(-1)]


def string_scalar(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return string_scalar(value.item())
        if value.dtype.kind in {"U", "S"}:
            return "".join(str(item) for item in value.reshape(-1)).strip()
        return " ".join(string_scalar(item) for item in value.reshape(-1)).strip()
    return str(value)


def nonempty_string_items(value: object) -> list[str]:
    items = []
    for item in np.asarray(value, dtype=object).reshape(-1):
        text = string_scalar(item).strip()
        if text:
            items.append(text)
    return items


def numeric_items(value: object) -> list[float]:
    values: list[float] = []
    for item in np.asarray(value, dtype=object).reshape(-1):
        try:
            values.append(float(np.asarray(item).reshape(-1)[0]))
        except Exception:
            values.append(0.0)
    return values


def secretory_complex_subunit_map(sec_struct: object, complexes: list[str]) -> dict[str, list[dict[str, object]]]:
    subunits_raw = np.asarray(getattr(sec_struct, "subunit"), dtype=object)
    stoich_raw = np.asarray(getattr(sec_struct, "subunit_stoichiometry"), dtype=object)
    mapping: dict[str, list[dict[str, object]]] = {}
    for index, complex_id in enumerate(complexes):
        subunit_row = subunits_raw[index] if subunits_raw.shape and subunits_raw.shape[0] > index else []
        stoich_row = stoich_raw[index] if stoich_raw.shape and stoich_raw.shape[0] > index else []
        subunit_ids = nonempty_string_items(subunit_row)
        stoich_values = numeric_items(stoich_row)
        entries: list[dict[str, object]] = []
        for sub_index, subunit_id in enumerate(subunit_ids):
            stoich = stoich_values[sub_index] if sub_index < len(stoich_values) else 1.0
            if stoich <= 0:
                continue
            entries.append({"subunit_id": subunit_id, "stoichiometry": stoich})
        mapping[complex_id] = entries
    return mapping


_MODEL_CACHE: dict[str, CobraModel] = {}


def load_pcsec_pichia_model(root: Path) -> CobraModel:
    cache_key = str(root.resolve())
    if cache_key in _MODEL_CACHE:
        return _MODEL_CACHE[cache_key]

    path = root / "Model" / "pcSecPichia.mat"
    raw = loadmat(path, squeeze_me=True, struct_as_record=False)
    model = raw["model"]
    matrix = getattr(model, "S")
    if not sparse.issparse(matrix):
        matrix = sparse.csc_matrix(matrix)
    else:
        matrix = matrix.tocsc()
    rules = string_list(getattr(model, "rules", []))
    gr_rules = string_list(getattr(model, "grRules", []))
    rxns = string_list(getattr(model, "rxns"))
    if len(rules) < len(rxns):
        rules.extend([""] * (len(rxns) - len(rules)))
    if len(gr_rules) < len(rxns):
        gr_rules.extend([""] * (len(rxns) - len(gr_rules)))
    result = CobraModel(
        source_file=str(path),
        rxns=rxns,
        mets=string_list(getattr(model, "mets")),
        genes=string_list(getattr(model, "genes")),
        lb=np.asarray(getattr(model, "lb"), dtype=float).reshape(-1),
        ub=np.asarray(getattr(model, "ub"), dtype=float).reshape(-1),
        b=np.asarray(getattr(model, "b"), dtype=float).reshape(-1),
        s_matrix=matrix,
        rules=rules,
        gr_rules=gr_rules,
    )
    _MODEL_CACHE[cache_key] = result
    return result


def load_single_struct(path: Path, variable_name: str):
    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    if variable_name not in data:
        available = ", ".join(sorted(key for key in data if not key.startswith("__")))
        raise KeyError(f"{path} does not contain {variable_name!r}. Available: {available}")
    return data[variable_name]


def numeric_vector(value: object) -> np.ndarray:
    return np.asarray(value, dtype=float).reshape(-1)


def load_metabolic_enzymedata(root: Path) -> MetabolicEnzymeData:
    path = root / "Enzymedata" / "pcSecPichia" / "enzymedat_PP.mat"
    struct = load_single_struct(path, "enzymedata")
    return MetabolicEnzymeData(
        enzymes=string_list(getattr(struct, "enzyme")),
        kcat=numeric_vector(getattr(struct, "kcat")),
    )


def load_secretory_enzymedata(root: Path) -> SecretoryEnzymeData:
    base = root / "Enzymedata" / "pcSecPichia"
    sec = load_single_struct(base / "enzymedataSEC_PP.mat", "enzymedataSEC")
    met = load_single_struct(base / "enzymedat_PP.mat", "enzymedata")
    dummy = load_single_struct(base / "enzymedataDummyER_PP.mat", "enzymedataDummyER")
    coefficient_rxns = string_list(getattr(met, "rxns")) + string_list(getattr(dummy, "rxns"))
    coefficient_values = np.concatenate([numeric_vector(getattr(met, "rxnscoef")), numeric_vector(getattr(dummy, "rxnscoef"))])
    reaction_coefficients: dict[str, float] = {}
    for reaction_id, value in zip(coefficient_rxns, coefficient_values):
        if reaction_id and reaction_id not in reaction_coefficients:
            reaction_coefficients[reaction_id] = float(value)
    complexes = string_list(getattr(sec, "enzyme"))
    return SecretoryEnzymeData(
        complexes=complexes,
        compartments=string_list(getattr(sec, "comp")),
        kcat=numeric_vector(getattr(sec, "kcat")),
        coefficient_refs=string_list(getattr(sec, "coefref")),
        reaction_coefficients=reaction_coefficients,
        complex_subunits=secretory_complex_subunit_map(sec, complexes),
    )


def load_combined_enzymedata(root: Path) -> CombinedEnzymeData:
    base = root / "Enzymedata" / "pcSecPichia"
    structs = [
        load_single_struct(base / "enzymedat_PP.mat", "enzymedata"),
        load_single_struct(base / "enzymedataSEC_PP.mat", "enzymedataSEC"),
        load_single_struct(base / "enzymedataMachine_PP.mat", "enzymedataMachine"),
    ]
    enzymes: list[str] = []
    kcat_arrays: list[np.ndarray] = []
    enzyme_mw_arrays: list[np.ndarray] = []
    proteins: list[str] = []
    protein_lengths: list[float] = []
    protein_mws: list[float] = []
    protein_kdegs: list[float] = []
    for struct in structs:
        enzymes.extend(string_list(getattr(struct, "enzyme")))
        kcat_arrays.append(numeric_vector(getattr(struct, "kcat")))
        enzyme_mw_arrays.append(numeric_vector(getattr(struct, "enzyme_MW")))
        local_proteins = string_list(getattr(struct, "proteins"))
        local_lengths = numeric_vector(getattr(struct, "proteinLength"))
        local_mws = numeric_vector(getattr(struct, "proteinMWs"))
        local_kdegs = numeric_vector(getattr(struct, "kdeg")) if hasattr(struct, "kdeg") else np.zeros(len(local_proteins), dtype=float)
        proteins.extend(local_proteins)
        protein_lengths.extend(float(value) for value in local_lengths)
        protein_mws.extend(float(value) for value in local_mws)
        protein_kdegs.extend(float(value) for value in local_kdegs)

    unique_proteins: list[str] = []
    unique_lengths: list[float] = []
    unique_mws: list[float] = []
    unique_kdegs: list[float] = []
    seen: set[str] = set()
    for protein_id, length, mw, kdeg in zip(proteins, protein_lengths, protein_mws, protein_kdegs):
        if protein_id in seen:
            continue
        seen.add(protein_id)
        unique_proteins.append(protein_id)
        unique_lengths.append(length)
        unique_mws.append(mw)
        unique_kdegs.append(kdeg)
    return CombinedEnzymeData(
        enzymes=enzymes,
        kcat=np.concatenate(kcat_arrays),
        enzyme_mw=np.concatenate(enzyme_mw_arrays),
        proteins=unique_proteins,
        protein_length=np.array(unique_lengths, dtype=float),
        protein_mw=np.array(unique_mws, dtype=float),
        protein_kdeg=np.array(unique_kdegs, dtype=float),
    )


def build_target_enzymedata(target: TargetSpec, model: CobraModel, secretory: SecretoryEnzymeData) -> TargetProteinEnzymeData:
    features = target_features(target)
    base = TargetProteinEnzymeData(
        protein_id=target.protein_id,
        protein_mw=float(features["protein_mw"]),
        protein_length=float(len(target.full_sequence)),
        disulfides=float(target.disulfide_sites),
        n_glycans=float(target.n_glycosylation_sites),
        o_glycans=float(target.o_glycosylation_sites),
        gpi_sites=float(target.gpi_sites),
        reaction_coefficients={},
    )
    coefficients: dict[str, float] = {}
    for complex_id, coefficient_ref in zip(secretory.complexes, secretory.coefficient_refs):
        suffix = f"_{complex_id}"
        for reaction_id in model.rxns:
            if reaction_id.startswith("dummyER") or not reaction_id.endswith(suffix):
                continue
            if not reaction_id.startswith(f"{target.protein_id}_"):
                continue
            coefficients[reaction_id] = evaluate_coefficient_ref(coefficient_ref, base)
    return TargetProteinEnzymeData(
        protein_id=base.protein_id,
        protein_mw=base.protein_mw,
        protein_length=base.protein_length,
        disulfides=base.disulfides,
        n_glycans=base.n_glycans,
        o_glycans=base.o_glycans,
        gpi_sites=base.gpi_sites,
        reaction_coefficients=coefficients,
    )


def evaluate_coefficient_ref(coefficient_ref: str, target: TargetProteinEnzymeData) -> float:
    coefficient_ref = coefficient_ref.strip()
    if coefficient_ref == "proteinLength":
        return target.protein_length / 467.0
    if coefficient_ref == "ProteinMW":
        return target.protein_mw / 54580.0
    variables = {
        "DSB": target.disulfides,
        "NG": target.n_glycans,
        "OG": target.o_glycans,
        "GPI": target.gpi_sites,
        "proteinLength": target.protein_length,
        "ProteinMW": target.protein_mw,
    }
    try:
        return float(eval(coefficient_ref, {"__builtins__": {}}, variables))
    except Exception as exc:
        raise ValueError(f"Unsupported coefficient reference: {coefficient_ref}") from exc


def set_media_pp_like_matlab(model: CobraModel, media_type: int = 4) -> CobraModel:
    lb = model.lb.copy()
    ub = model.ub.copy()
    for index, reaction_id in enumerate(model.rxns):
        if reaction_id.startswith("Ex_"):
            lb[index] = 0.0
            ub[index] = 1000.0

    def set_lower(reactions: Iterable[str], lower: float) -> None:
        reaction_index = model.reaction_index
        for reaction_id in reactions:
            index = reaction_index.get(reaction_id)
            if index is not None:
                lb[index] = lower

    minimal = ("Ex_nh4", "Ex_o2", "Ex_pi", "Ex_so4", "Ex_fe2", "Ex_h", "Ex_h2o", "Ex_na1", "Ex_k", "Ex_co2")
    vitamins = ("Ex_btn", "Ex_thm", "Ex_4abz", "Ex_pnto_R", "Ex_inost", "Ex_nac", "Ex_ribflv")
    core_aa = (
        "Ex_arg_L",
        "Ex_asp_L",
        "Ex_glu_L",
        "Ex_gly",
        "Ex_his_L",
        "Ex_ile_L",
        "Ex_leu_L",
        "Ex_lys_L",
        "Ex_met_L",
        "Ex_phe_L",
        "Ex_thr_L",
        "Ex_trp_L",
        "Ex_tyr_L",
        "Ex_val_L",
        "Ex_ura",
    )
    all_aa = core_aa + ("Ex_ala_L", "Ex_asn_L", "Ex_cys_L", "Ex_gln_L", "Ex_pro_L", "Ex_ser_L")
    # NOTE: MATLAB setMediaPP.m line 62 calls setLowerBounds(model, minimal, -1000)
    # but does NOT assign the return value back to model. This means the minimal
    # exchange lower bounds are silently NOT applied in MATLAB. We replicate this
    # behavior here for LP-level alignment. Only vitamins and amino acids are
    # actually applied (those calls DO capture the return value).
    # set_lower(minimal, -1000.0)  # <-- commented out to match MATLAB bug
    if media_type in {2, 3, 4, 5}:
        set_lower(vitamins, -2.0)
    if media_type in {3, 4}:
        set_lower(core_aa, -0.08)
    if media_type == 5:
        set_lower(all_aa, -0.08)

    return CobraModel(
        source_file=model.source_file,
        rxns=model.rxns,
        mets=model.mets,
        genes=model.genes,
        lb=lb,
        ub=ub,
        b=model.b,
        s_matrix=model.s_matrix,
        rules=model.rules,
        gr_rules=model.gr_rules,
    )


def prepare_glucose_model(model: CobraModel, media_type: int = 4) -> CobraModel:
    return set_media_pp_like_matlab(model, media_type=media_type).with_bounds(
        {
            "Ex_glc_D": (-1000.0, None),
            "BIOMASS": (None, 1000.0),
            "Ex_glyc": (0.0, None),
            "BIOMASS_glyc": (0.0, 0.0),
            "Ex_meoh": (0.0, None),
            "BIOMASS_meoh": (0.0, 0.0),
            "Ex_o2": (-1000.0, None),
        }
    )


@dataclass(frozen=True)
class AminoAcidStoichiometry:
    amino_acids: tuple[str, ...]
    translation_substrates: dict[str, str]
    translation_products: dict[str, str]
    degradation_products: dict[str, str]
    energy_substrates: tuple[str, str, str]
    energy_products: tuple[str, str, str, str]

    def translation(self, protein_id: str, sequence: str) -> dict[str, float]:
        counts = Counter(sequence)
        counts["M"] += 1
        stoich: dict[str, float] = {}
        for aa in self.amino_acids:
            count = counts.get(aa, 0)
            if not count:
                continue
            accumulate(stoich, self.translation_substrates[aa], -count)
            accumulate(stoich, self.translation_products[aa], count)
        length = len(sequence)
        total_atp = 2 * length + 3
        total_gtp = 2 * length + 3
        total_h2o = total_atp + total_gtp
        accumulate(stoich, self.energy_substrates[0], -total_h2o)
        accumulate(stoich, self.energy_substrates[1], -total_atp)
        accumulate(stoich, self.energy_substrates[2], -total_gtp)
        accumulate(stoich, self.energy_products[0], total_h2o)
        accumulate(stoich, self.energy_products[1], total_atp)
        accumulate(stoich, self.energy_products[2], total_gtp)
        accumulate(stoich, self.energy_products[3], total_h2o)
        accumulate(stoich, f"{protein_id}_peptide[c]", 1.0)
        return stoich

    def signal_peptide_degradation(self, protein_id: str, sequence: str) -> dict[str, float]:
        return self.degradation(sequence, f"{protein_id}_sp[c]", include_initiator_methionine=False)

    def subunit_degradation(self, protein_id: str, sequence: str) -> dict[str, float]:
        return self.degradation(sequence, f"{protein_id}_subunit[c]", include_initiator_methionine=True)

    def degradation(self, sequence: str, terminal_metabolite: str, include_initiator_methionine: bool) -> dict[str, float]:
        counts = Counter(sequence)
        if include_initiator_methionine:
            counts["M"] += 1
        stoich: dict[str, float] = {}
        for aa in self.amino_acids:
            count = counts.get(aa, 0)
            if count:
                accumulate(stoich, self.degradation_products[aa], count)
        energy_count = math.floor(1.3 * len(sequence))
        accumulate(stoich, self.energy_substrates[0], -energy_count)
        accumulate(stoich, self.energy_substrates[1], -energy_count)
        accumulate(stoich, self.energy_products[0], energy_count)
        accumulate(stoich, self.energy_products[1], energy_count)
        accumulate(stoich, self.energy_products[3], energy_count)
        accumulate(stoich, terminal_metabolite, -1.0)
        return stoich


def load_aa_stoichiometry(root: Path) -> AminoAcidStoichiometry:
    path = root / "Data" / "pcSecPichia" / "aa_id_PP.xlsx"
    cytoplasm = pd.read_excel(path, sheet_name="cytoplasm")
    energy = pd.read_excel(path, sheet_name="energy")
    amino_acids: list[str] = []
    translation_substrates: dict[str, str] = {}
    translation_products: dict[str, str] = {}
    degradation_products: dict[str, str] = {}
    for _, row in cytoplasm.iterrows():
        aa = str(row["aa_id"]).strip()
        amino_acids.append(aa)
        translation_substrates[aa] = str(row["subs_id"]).strip()
        translation_products[aa] = str(row["prod_id"]).strip()
        degradation_products[aa] = str(row["deg_prod_id"]).strip()
    return AminoAcidStoichiometry(
        amino_acids=tuple(amino_acids),
        translation_substrates=translation_substrates,
        translation_products=translation_products,
        degradation_products=degradation_products,
        energy_substrates=tuple(str(item).strip() for item in energy["subs_id"].dropna().tolist()[:3]),
        energy_products=tuple(str(item).strip() for item in energy["prod_id"].dropna().tolist()[:4]),
    )


def accumulate(stoich: dict[str, float], metabolite_id: str, coefficient: float) -> None:
    if not metabolite_id or coefficient == 0:
        return
    stoich[metabolite_id] = stoich.get(metabolite_id, 0.0) + float(coefficient)
    if stoich[metabolite_id] == 0:
        del stoich[metabolite_id]


def solve_maximize(model: CobraModel, objective_reaction: str, key_reactions: Iterable[str] = ()) -> SolveResult:
    reaction_index = model.reaction_index
    if objective_reaction not in reaction_index:
        return SolveResult(objective_reaction, "missing_objective", False, None, "Objective reaction not found.", {})
    objective = np.zeros(len(model.rxns), dtype=float)
    objective[reaction_index[objective_reaction]] = -1.0
    bounds = list(zip(model.lb.tolist(), model.ub.tolist()))
    result = linprog(
        c=objective,
        A_eq=model.s_matrix,
        b_eq=model.b,
        bounds=bounds,
        method="highs",
        options={"presolve": True, "disp": False},
    )
    fluxes: dict[str, float] = {}
    if result.success and result.x is not None:
        for reaction_id in {objective_reaction, *key_reactions}:
            index = reaction_index.get(reaction_id)
            if index is not None:
                fluxes[reaction_id] = float(result.x[index])
    return SolveResult(
        objective=objective_reaction,
        status=str(result.status),
        success=bool(result.success),
        objective_value=float(-result.fun) if result.success else None,
        message=str(result.message),
        fluxes=fluxes,
        sensitivity=_linprog_sensitivity(result) if result.success else None,
    )


def solve_pcsec_maximize(
    model: CobraModel,
    objective_reaction: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    key_reactions: Iterable[str] = (),
    total_protein_content: float = 0.37,
    unmodeled_er_protein_fraction: float = 0.040,
    mitochondrial_protein_fraction: float = 0.05,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> tuple[SolveResult, dict[str, int]]:
    reaction_index = model.reaction_index
    if objective_reaction not in reaction_index:
        return (
            SolveResult(objective_reaction, "missing_objective", False, None, "Objective reaction not found.", {}),
            {},
        )
    objective = np.zeros(len(model.rxns), dtype=float)
    objective[reaction_index[objective_reaction]] = -1.0
    A_eq, b_eq, A_ub, b_ub, counts = build_pcsec_constraint_matrices(
        model,
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
        mu=mu,
        total_protein_content=total_protein_content,
        unmodeled_er_protein_fraction=unmodeled_er_protein_fraction,
        mitochondrial_protein_fraction=mitochondrial_protein_fraction,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    result = linprog(
        c=objective,
        A_eq=A_eq,
        b_eq=b_eq,
        A_ub=A_ub,
        b_ub=b_ub,
        bounds=list(zip(model.lb.tolist(), model.ub.tolist())),
        method="highs",
        options={"presolve": True, "disp": False},
    )
    fluxes: dict[str, float] = {}
    if result.success and result.x is not None:
        for reaction_id in {objective_reaction, *key_reactions}:
            index = reaction_index.get(reaction_id)
            if index is not None:
                fluxes[reaction_id] = float(result.x[index])
    return (
        SolveResult(
            objective=objective_reaction,
            status=str(result.status),
            success=bool(result.success),
            objective_value=float(-result.fun) if result.success else None,
            message=str(result.message),
            fluxes=fluxes,
            sensitivity=_linprog_sensitivity(result) if result.success else None,
        ),
        counts,
    )


def _linprog_sensitivity(result) -> dict[str, tuple[float, ...]]:
    return {
        "eq_marginals": _marginal_tuple(getattr(getattr(result, "eqlin", None), "marginals", ())),
        "ub_marginals": _marginal_tuple(getattr(getattr(result, "ineqlin", None), "marginals", ())),
        "lower_marginals": _marginal_tuple(getattr(getattr(result, "lower", None), "marginals", ())),
        "upper_marginals": _marginal_tuple(getattr(getattr(result, "upper", None), "marginals", ())),
    }


def _marginal_tuple(values) -> tuple[float, ...]:
    try:
        return tuple(float(value) for value in values)
    except TypeError:
        return ()


def build_pcsec_constraint_matrices(
    model: CobraModel,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    total_protein_content: float,
    unmodeled_er_protein_fraction: float,
    mitochondrial_protein_fraction: float,
    write_ribosome_translation_constraint: bool,
    write_misfolding_constraints: bool,
) -> tuple[sparse.csr_matrix, np.ndarray, sparse.csr_matrix, np.ndarray, dict[str, int]]:
    rows_eq: list[sparse.csr_matrix] = [model.s_matrix.tocsr()]
    rhs_eq: list[np.ndarray] = [model.b.astype(float)]
    rows_ub: list[sparse.csr_matrix] = []
    rhs_ub: list[float] = []
    counts: dict[str, int] = {"stoichiometric": model.s_matrix.shape[0]}

    metabolic_rows, metabolic_rhs = metabolic_coupling_rows(model, metabolic, mu)
    append_rows(rows_eq, rhs_eq, metabolic_rows, metabolic_rhs)
    counts["metabolic_coupling"] = len(metabolic_rhs)

    secretory_rows, secretory_rhs = secretory_coupling_rows(model, secretory, mu)
    append_rows(rows_eq, rhs_eq, secretory_rows, secretory_rhs)
    counts["secretory_coupling"] = len(secretory_rhs)

    protein_rows, protein_rhs = protein_mass_rows(
        model,
        combined,
        mu=mu,
        total_protein_content=total_protein_content,
        unmodeled_er_protein_fraction=unmodeled_er_protein_fraction,
    )
    append_rows(rows_eq, rhs_eq, protein_rows, protein_rhs)
    counts["protein_mass"] = len(protein_rhs)

    proteasome_constraint_rows, proteasome_rhs = proteasome_rows(model, combined, mu)
    append_rows(rows_eq, rhs_eq, proteasome_constraint_rows, proteasome_rhs)
    counts["proteasome"] = len(proteasome_rhs)

    ribosome_constraint_rows, ribosome_rhs = ribosome_assembly_rows(model, combined, mu)
    append_rows(rows_eq, rhs_eq, ribosome_constraint_rows, ribosome_rhs)
    counts["ribosome_assembly"] = len(ribosome_rhs)

    if write_ribosome_translation_constraint:
        ribosome_translation_rows_, ribosome_translation_rhs = ribosome_translation_rows(model, combined, mu)
        append_rows(rows_eq, rhs_eq, ribosome_translation_rows_, ribosome_translation_rhs)
        counts["ribosome_translation"] = len(ribosome_translation_rhs)
    else:
        counts["ribosome_translation"] = 0

    if write_misfolding_constraints:
        misfolding_rows_, misfolding_rhs = misfolding_constraint_rows(model, combined)
        append_rows(rows_eq, rhs_eq, misfolding_rows_, misfolding_rhs)
        counts["misfolding"] = len(misfolding_rhs)
    else:
        counts["misfolding"] = 0

    mito_constraint_rows, mito_rhs = mitochondrial_rows(model, combined, mu, mitochondrial_protein_fraction)
    if mito_constraint_rows:
        rows_ub.extend(mito_constraint_rows)
        rhs_ub.extend(mito_rhs)
    counts["mitochondrial"] = len(mito_rhs)

    A_eq = sparse.vstack(rows_eq, format="csr")
    b_eq = np.concatenate(rhs_eq)
    if rows_ub:
        A_ub = sparse.vstack(rows_ub, format="csr")
        b_ub = np.array(rhs_ub, dtype=float)
    else:
        A_ub = sparse.csr_matrix((0, len(model.rxns)))
        b_ub = np.array([], dtype=float)
    counts["eq_total"] = A_eq.shape[0]
    counts["ub_total"] = A_ub.shape[0]
    return A_eq, b_eq, A_ub, b_ub, counts


def append_rows(
    row_groups: list[sparse.csr_matrix],
    rhs_groups: list[np.ndarray],
    rows: list[sparse.csr_matrix],
    rhs: list[float],
) -> None:
    if not rows:
        return
    row_groups.append(sparse.vstack(rows, format="csr"))
    rhs_groups.append(np.array(rhs, dtype=float))


def sparse_row(model: CobraModel, terms: dict[int, float]) -> sparse.csr_matrix:
    columns = []
    values = []
    for index, value in terms.items():
        if value:
            columns.append(index)
            values.append(float(value))
    return sparse.csr_matrix((values, ([0] * len(columns), columns)), shape=(1, len(model.rxns)))


def metabolic_coupling_rows(model: CobraModel, metabolic: MetabolicEnzymeData, mu: float) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    rows: list[sparse.csr_matrix] = []
    rhs: list[float] = []
    for enzyme_id, kcat in zip(metabolic.enzymes, metabolic.kcat):
        reaction_id = enzyme_id.replace("_complex", "")
        formation_id = f"{enzyme_id}_formation"
        if reaction_id not in reaction_index or formation_id not in reaction_index:
            continue
        rows.append(sparse_row(model, {reaction_index[reaction_id]: 1.0, reaction_index[formation_id]: -float(kcat) / mu}))
        rhs.append(0.0)
    return rows, rhs


def secretory_coupling_rows(model: CobraModel, secretory: SecretoryEnzymeData, mu: float) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    rows: list[sparse.csr_matrix] = []
    rhs: list[float] = []
    for entry in secretory.unique_complex_entries():
        formation_id = f"{entry.complex_id}_formation"
        if formation_id not in reaction_index:
            continue
        terms: dict[int, float] = {}
        suffix = f"_{entry.complex_id}"
        for index, reaction_id in enumerate(model.rxns):
            if not reaction_id.endswith(suffix):
                continue
            coefficient = secretory.reaction_coefficients.get(reaction_id)
            if coefficient is None:
                continue
            terms[index] = terms.get(index, 0.0) + float(coefficient)
        if not terms:
            continue
        terms[reaction_index[formation_id]] = terms.get(reaction_index[formation_id], 0.0) - float(entry.kcat) / mu
        rows.append(sparse_row(model, terms))
        rhs.append(0.0)
    return rows, rhs


def protein_mass_rows(
    model: CobraModel,
    combined: CombinedEnzymeData,
    mu: float,
    total_protein_content: float,
    unmodeled_er_protein_fraction: float,
) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    terms: dict[int, float] = {}
    for index, reaction_id in enumerate(model.rxns):
        if "_dilution" not in reaction_id or "dummy" in reaction_id:
            continue
        try:
            mw = combined.molecular_weight_for_dilution_reaction(reaction_id)
        except KeyError:
            continue
        terms[index] = mw / 1000.0
    for reaction_id in ("dilute_dummy", "dilute_dummyER"):
        if reaction_id in reaction_index:
            terms[reaction_index[reaction_id]] = terms.get(reaction_index[reaction_id], 0.0) + 40.0

    modeled_fraction = modeled_protein_fraction(model)
    rows = [sparse_row(model, terms)]
    rhs = [mu * total_protein_content * modeled_fraction]
    if "dilute_dummyER" in reaction_index:
        rows.append(sparse_row(model, {reaction_index["dilute_dummyER"]: 40.0}))
        rhs.append(mu * unmodeled_er_protein_fraction)
    return rows, rhs


def mitochondrial_rows(
    model: CobraModel,
    combined: CombinedEnzymeData,
    mu: float,
    mitochondrial_protein_fraction: float,
) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    terms: dict[int, float] = {}
    for reaction_id in collect_compartment_reactions(model, ("m", "mm")):
        dilution_id = f"{reaction_id}_complex_dilution"
        if dilution_id not in reaction_index:
            continue
        enzyme_id = dilution_id.replace("_dilution", "")
        try:
            mw = combined.exact_enzyme_mw(enzyme_id)
        except KeyError:
            continue
        terms[reaction_index[dilution_id]] = mw / 1000.0
    if not terms:
        return [], []
    return [sparse_row(model, terms)], [mu * mitochondrial_protein_fraction]


def proteasome_rows(model: CobraModel, combined: CombinedEnzymeData, mu: float) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    if "Mach_proteasome_complex_formation" not in reaction_index:
        return [], []
    terms: dict[int, float] = {}
    for index, reaction_id in enumerate(model.rxns):
        if reaction_id.endswith("_subunit_degradation"):
            protein_id = reaction_id.replace("_subunit_degradation", "").replace("r_", "")
            try:
                coefficient = combined.exact_protein_length(protein_id) / 467.0
            except KeyError:
                continue
        elif reaction_id.endswith("_sp_degradation"):
            coefficient = 25.0 / 467.0
        else:
            continue
        terms[index] = coefficient
    if not terms:
        return [], []
    formation_index = reaction_index["Mach_proteasome_complex_formation"]
    terms[formation_index] = terms.get(formation_index, 0.0) - combined.exact_enzyme_kcat("Mach_proteasome_complex") / mu
    return [sparse_row(model, terms)], [0.0]


def ribosome_assembly_rows(model: CobraModel, combined: CombinedEnzymeData, mu: float) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    use_id = "Mach_Ribosome_complex_formation"
    syn_id = "Mach_Ribosome_Assembly_Factors_complex_formation"
    if use_id not in reaction_index or syn_id not in reaction_index:
        return [], []
    coefficient = combined.exact_enzyme_kcat("Mach_Ribosome_Assembly_Factors_complex") / mu
    return [sparse_row(model, {reaction_index[use_id]: 1.0, reaction_index[syn_id]: -coefficient})], [0.0]


def ribosome_translation_rows(model: CobraModel, combined: CombinedEnzymeData, mu: float) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    formation_id = "Mach_Ribosome_complex_formation"
    if formation_id not in reaction_index:
        return [], []
    terms: dict[int, float] = {}
    for index, reaction_id in enumerate(model.rxns):
        if not reaction_id.endswith("_translation"):
            continue
        protein_id = reaction_id.replace("_peptide_translation", "").replace("r_", "")
        try:
            terms[index] = combined.exact_protein_length(protein_id)
        except KeyError:
            continue
    for reaction_id in ("translate_dummy", "translate_dummyER"):
        if reaction_id in reaction_index:
            terms[reaction_index[reaction_id]] = terms.get(reaction_index[reaction_id], 0.0) + 423.0
    if "PROTEINS" in reaction_index:
        terms[reaction_index["PROTEINS"]] = terms.get(reaction_index["PROTEINS"], 0.0) + 4.23
    if not terms:
        return [], []
    kcat_ribosome = 2.06e1 * mu / (0.486 + mu) * 3600.0
    terms[reaction_index[formation_id]] = terms.get(reaction_index[formation_id], 0.0) - kcat_ribosome / mu
    return [sparse_row(model, terms)], [0.0]


def misfolding_constraint_rows(model: CobraModel, combined: CombinedEnzymeData) -> tuple[list[sparse.csr_matrix], list[float]]:
    reaction_index = model.reaction_index
    rows: list[sparse.csr_matrix] = []
    rhs: list[float] = []
    for reaction_id in model.rxns:
        if "_misfold_" not in reaction_id or "dummy" in reaction_id:
            continue
        protein_id = reaction_id.split("_misfold_", 1)[0]
        translation_id = f"r_{protein_id}_peptide_translation"
        if translation_id not in reaction_index or reaction_id not in reaction_index:
            continue
        try:
            kdeg = combined.exact_protein_kdeg(protein_id)
        except KeyError:
            continue
        rows.append(sparse_row(model, {reaction_index[reaction_id]: 1.0, reaction_index[translation_id]: -kdeg}))
        rhs.append(0.0)
    return rows, rhs


def collect_compartment_reactions(model: CobraModel, compartment_ids: tuple[str, ...]) -> list[str]:
    tokens = tuple(f"[{compartment_id}]" for compartment_id in compartment_ids)
    met_indices = [index for index, metabolite_id in enumerate(model.mets) if any(token in metabolite_id for token in tokens)]
    if not met_indices:
        return []
    submatrix = model.s_matrix[met_indices, :].tocsc()
    has_stoich = np.diff(submatrix.indptr) > 0
    return [
        reaction_id
        for index, reaction_id in enumerate(model.rxns)
        if has_stoich[index] and index < len(model.rules) and str(model.rules[index]).strip() not in {"", "[]"}
    ]


def modeled_protein_fraction(model: CobraModel) -> float:
    try:
        met_index = model.mets.index("PROTEIN[c]")
        rxn_index = model.reaction_index["BIOMASS"]
    except (KeyError, ValueError):
        return 1.0
    return 1.0 + float(model.s_matrix[met_index, rxn_index])


def inactive_reactions_for_gene(model: CobraModel, gene_id: str) -> list[str]:
    if gene_id not in model.gene_index:
        return []
    gene_number = model.gene_index[gene_id] + 1
    inactive: list[str] = []
    for reaction_id, rule in zip(model.rxns, model.rules):
        normalized = rule.strip()
        if not normalized or normalized == "[]":
            continue
        if evaluate_rule_after_ko(normalized, gene_number) is False:
            inactive.append(reaction_id)
    return inactive


def evaluate_rule_after_ko(rule: str, knocked_gene_number: int) -> bool | None:
    expression = rule.replace("&", " and ").replace("|", " or ")

    def repl(match: re.Match[str]) -> str:
        gene_number = int(match.group(1))
        return "False" if gene_number == knocked_gene_number else "True"

    expression = RULE_TOKEN_PATTERN.sub(repl, expression)
    if "x(" in expression:
        return None
    try:
        return bool(eval(expression, {"__builtins__": {}}, {}))
    except Exception:
        return None


def run_ko_screen(model: CobraModel, baseline: SolveResult, genes: list[str], objective: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for gene_id in genes:
        inactive = inactive_reactions_for_gene(model, gene_id)
        changes = {reaction_id: (0.0, 0.0) for reaction_id in inactive}
        solved = solve_maximize(model.with_bounds(changes), objective, key_reactions=("BIOMASS", "Ex_glc_D"))
        rows.append(
            {
                "gene": gene_id,
                "inactive_reaction_count": len(inactive),
                "inactive_reactions_preview": inactive[:10],
                "status": solved.status,
                "success": solved.success,
                "objective_value": solved.objective_value,
                "delta_vs_baseline": (
                    solved.objective_value - baseline.objective_value
                    if solved.success and baseline.objective_value is not None and solved.objective_value is not None
                    else None
                ),
            }
        )
    return rows


def run_oe_screen(model: CobraModel, baseline: SolveResult, reactions: list[str], objective: str, factor: float = 2.0) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    reaction_index = model.reaction_index
    for reaction_id in reactions:
        index = reaction_index.get(reaction_id)
        if index is None:
            rows.append({"reaction": reaction_id, "status": "missing_reaction", "success": False})
            continue
        old_upper = float(model.ub[index])
        new_upper = old_upper * factor if old_upper > 0 else factor
        solved = solve_maximize(model.with_bounds({reaction_id: (None, new_upper)}), objective, key_reactions=("BIOMASS", "Ex_glc_D"))
        rows.append(
            {
                "reaction": reaction_id,
                "old_upper": old_upper,
                "new_upper": new_upper,
                "status": solved.status,
                "success": solved.success,
                "objective_value": solved.objective_value,
                "delta_vs_baseline": (
                    solved.objective_value - baseline.objective_value
                    if solved.success and baseline.objective_value is not None and solved.objective_value is not None
                    else None
                ),
            }
        )
    return rows


def run_pcsec_ko_screen(
    model: CobraModel,
    baseline: SolveResult,
    genes: list[str],
    objective: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for gene_id in genes:
        inactive = inactive_reactions_for_gene(model, gene_id)
        changes = {reaction_id: (0.0, 0.0) for reaction_id in inactive}
        solved, counts = solve_pcsec_maximize(
            model.with_bounds(changes),
            objective,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            key_reactions=("BIOMASS", "Ex_glc_D", objective),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        rows.append(
            {
                "gene": gene_id,
                "inactive_reaction_count": len(inactive),
                "inactive_reactions_preview": inactive[:10],
                "status": solved.status,
                "success": solved.success,
                "objective_value": solved.objective_value,
                "delta_vs_baseline": (
                    solved.objective_value - baseline.objective_value
                    if solved.success and baseline.objective_value is not None and solved.objective_value is not None
                    else None
                ),
                "constraint_counts": counts,
            }
        )
    return rows


def run_pcsec_reaction_ko_screen(
    model: CobraModel,
    baseline: SolveResult,
    reactions: list[str],
    objective: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    reaction_index = model.reaction_index
    for reaction_id in reactions:
        index = reaction_index.get(reaction_id)
        if index is None:
            rows.append({"reaction": reaction_id, "status": "missing_reaction", "success": False})
            continue
        old_lower = float(model.lb[index])
        old_upper = float(model.ub[index])
        solved, counts = solve_pcsec_maximize(
            model.with_bounds({reaction_id: (0.0, 0.0)}),
            objective,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            key_reactions=("BIOMASS", "Ex_glc_D", objective),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        rows.append(
            {
                "reaction": reaction_id,
                "old_lower": old_lower,
                "old_upper": old_upper,
                "status": solved.status,
                "success": solved.success,
                "objective_value": solved.objective_value,
                "delta_vs_baseline": (
                    solved.objective_value - baseline.objective_value
                    if solved.success and baseline.objective_value is not None and solved.objective_value is not None
                    else None
                ),
                "constraint_counts": counts,
            }
        )
    return rows


def run_pcsec_oe_screen(
    model: CobraModel,
    baseline: SolveResult,
    reactions: list[str],
    objective: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    factor: float = 2.0,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    reaction_index = model.reaction_index
    for reaction_id in reactions:
        index = reaction_index.get(reaction_id)
        if index is None:
            rows.append({"reaction": reaction_id, "status": "missing_reaction", "success": False})
            continue
        old_upper = float(model.ub[index])
        new_upper = old_upper
        complex_id = reaction_id.replace("_formation", "")
        perturbed_secretory = secretory
        perturbed_combined = combined
        capacity_basis = "reaction_upper_bound"
        if complex_id.startswith("sec_"):
            perturbed_secretory = secretory.with_complex_kcat_multiplier(complex_id, factor)
            capacity_basis = "secretory_complex_kcat_multiplier"
        elif complex_id.startswith("Mach_"):
            perturbed_combined = combined.with_enzyme_kcat_multiplier(complex_id, factor)
            capacity_basis = "machine_complex_kcat_multiplier"
        else:
            new_upper = old_upper * factor if old_upper > 0 else factor
        solved, counts = solve_pcsec_maximize(
            model.with_bounds({reaction_id: (None, new_upper)}),
            objective,
            metabolic=metabolic,
            secretory=perturbed_secretory,
            combined=perturbed_combined,
            mu=mu,
            key_reactions=("BIOMASS", "Ex_glc_D", objective),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        rows.append(
            {
                "reaction": reaction_id,
                "capacity_basis": capacity_basis,
                "capacity_multiplier": factor,
                "old_upper": old_upper,
                "new_upper": new_upper,
                "status": solved.status,
                "success": solved.success,
                "objective_value": solved.objective_value,
                "delta_vs_baseline": (
                    solved.objective_value - baseline.objective_value
                    if solved.success and baseline.objective_value is not None and solved.objective_value is not None
                    else None
                ),
                "constraint_counts": counts,
            }
        )
    return rows


def run_pcsec_growth_tradeoff(
    model: CobraModel,
    objective: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu_points: Iterable[float],
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mu in sorted({float(value) for value in mu_points if float(value) > 0}):
        fixed_model = model.with_bounds({"BIOMASS": (mu, mu)})
        solved, counts = solve_pcsec_maximize(
            fixed_model,
            objective,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            key_reactions=("BIOMASS", "Ex_glc_D", "Ex_o2", objective),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        secretion = solved.objective_value if solved.success else None
        rows.append(
            {
                "mu": mu,
                "success": solved.success,
                "status": solved.status,
                "secretion_flux": secretion,
                "secretion_per_biomass": secretion / mu if secretion is not None and mu > 0 else None,
                "message": solved.message,
                "constraint_counts": counts,
            }
        )
    return rows


def classify_secretory_process(reaction_id: str | None) -> str:
    if not reaction_id:
        return "unknown"
    text = reaction_id.lower()
    if "ribosome" in text:
        return "ribosome"
    if "proteasome" in text or "subunit_degradation" in text or "sp_degradation" in text:
        return "proteasome_degradation"
    if any(token in text for token in ("pdi", "ero1", "dsb")):
        return "disulfide_folding"
    if any(token in text for token in ("ostc", "cwh41", "rot2", "mns1", "och1", "mpol", "mnn2", "ng_")):
        return "n_glycan_processing"
    if any(token in text for token in ("pmt", "ktr", "o_linked", "og_")):
        return "o_glycan_processing"
    if any(token in text for token in ("kar2", "bip", "nefs", "rac", "ssa1", "ydj1", "snl1")):
        return "chaperone_folding"
    if any(token in text for token in ("erad", "ubc6", "ubc7", "hrd", "doa10", "yos9", "cdc48", "png1", "rad23")):
        return "erad_misfolding"
    if any(token in text for token in ("sec61", "sec63", "spc", "srp", "get1", "get2", "get3")):
        return "er_translocation"
    if any(token in text for token in ("copii", "sec12", "sar1", "sec23", "sec24", "sec13", "sec31", "ypt1", "uso1")):
        return "er_to_golgi_transport"
    if any(token in text for token in ("arf1", "pep12", "vps", "chc1", "clc1", "apl", "apm", "aps")):
        return "golgi_surface_transport"
    if reaction_id.startswith("sec_") or "_sec_" in reaction_id:
        return "secretory_capacity"
    return "metabolic_or_other"


def infer_gene_symbols_from_reaction(reaction_id: str | None) -> list[str]:
    if not reaction_id:
        return []
    text = reaction_id
    for prefix in ("sec_", "Mach_"):
        if text.startswith(prefix):
            text = text[len(prefix) :]
    for suffix in ("_formation", "_complex"):
        text = text.replace(suffix, "")
    tokens = [token for token in re.split(r"[_\\s]+", text) if token]
    stopwords = {
        "normal",
        "Golgi",
        "ERGL",
        "HDSVI",
        "HDSVII",
        "Post",
        "translation",
        "PSTA",
        "TC",
        "NG",
        "OG",
        "ERNG",
        "GLNG",
        "GLOG",
        "COPII",
        "Ribosome",
        "Assembly",
        "Factors",
    }
    return [token for token in tokens if token not in stopwords][:12]


def relative_change(delta: float | None, baseline_value: float | None) -> float | None:
    if delta is None or baseline_value is None or abs(baseline_value) < 1e-15:
        return None
    return float(delta) / abs(float(baseline_value))


def classify_candidate_effect(success: bool, relative_delta: float | None) -> str:
    if not success:
        return "infeasible_at_fixed_mu"
    if relative_delta is None:
        return "unknown"
    if relative_delta > 0.01:
        return "strong_improvement"
    if relative_delta > IMPROVEMENT_REL_TOLERANCE:
        return "weak_improvement"
    if relative_delta < -0.01:
        return "strong_decrease"
    if relative_delta < -IMPROVEMENT_REL_TOLERANCE:
        return "weak_decrease"
    return "neutral"


def candidate_table_row(
    row: dict[str, object],
    perturbation: str,
    entity_type: str,
    candidate_id: str,
    baseline_value: float | None,
    fixed_mu: float,
    reaction_id: str | None = None,
    gene_id: str | None = None,
    simulation_basis: str = "",
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> dict[str, object]:
    success = bool(row.get("success"))
    objective_value = row.get("objective_value")
    secretion_flux = float(objective_value) if objective_value is not None else None
    delta_value = row.get("delta_vs_baseline")
    delta = float(delta_value) if delta_value is not None else None
    rel_delta = relative_change(delta, baseline_value)
    process = classify_secretory_process(reaction_id)
    effect = classify_candidate_effect(success, rel_delta)
    complex_id = reaction_id.replace("_formation", "") if reaction_id else None
    subunit_rows = (complex_subunits or {}).get(complex_id or "", [])
    subunit_ids = [str(item["subunit_id"]) for item in subunit_rows]
    subunit_stoich = [float(item["stoichiometry"]) for item in subunit_rows]
    return {
        "candidate_id": candidate_id,
        "perturbation": perturbation,
        "entity_type": entity_type,
        "simulation_basis": simulation_basis,
        "gene": gene_id,
        "reaction": reaction_id,
        "inferred_gene_symbols": infer_gene_symbols_from_reaction(reaction_id),
        "complex_subunit_ids": subunit_ids,
        "complex_subunit_stoichiometry": subunit_stoich,
        "process": process,
        "success": success,
        "effect": effect,
        "secretion_flux": secretion_flux,
        "secretion_delta": delta,
        "relative_secretion_delta": rel_delta,
        "secretion_per_biomass": secretion_flux / fixed_mu if secretion_flux is not None and fixed_mu > 0 else None,
        "growth_mu": fixed_mu if success else None,
        "growth_constraint": f"BIOMASS fixed to {fixed_mu}",
        "score": rel_delta if success and rel_delta is not None else None,
        "raw_status": row.get("status"),
        "affected_reaction_count": row.get("inactive_reaction_count"),
        "affected_reactions_preview": row.get("inactive_reactions_preview"),
        "old_lower": row.get("old_lower"),
        "old_upper": row.get("old_upper"),
        "new_upper": row.get("new_upper"),
        "capacity_basis": row.get("capacity_basis"),
        "capacity_multiplier": row.get("capacity_multiplier"),
    }


def build_candidate_table(
    gpr_ko_rows: list[dict[str, object]],
    capacity_ko_rows: list[dict[str, object]],
    oe_rows: list[dict[str, object]],
    baseline_value: float | None,
    fixed_mu: float,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in gpr_ko_rows:
        gene_id = str(row.get("gene"))
        preview = row.get("inactive_reactions_preview") or []
        reaction_id = str(preview[0]) if preview else None
        rows.append(
            candidate_table_row(
                row,
                perturbation="KO",
                entity_type="model_gpr_gene",
                candidate_id=gene_id,
                baseline_value=baseline_value,
                fixed_mu=fixed_mu,
                reaction_id=reaction_id,
                gene_id=gene_id,
                simulation_basis="COBRA GPR knockout: all reactions disabled when the gene rule becomes false.",
                complex_subunits=complex_subunits,
            )
        )
    for row in capacity_ko_rows:
        reaction_id = str(row.get("reaction"))
        rows.append(
            candidate_table_row(
                row,
                perturbation="KO",
                entity_type="secretory_capacity_proxy",
                candidate_id=reaction_id,
                baseline_value=baseline_value,
                fixed_mu=fixed_mu,
                reaction_id=reaction_id,
                simulation_basis="Proxy knockout: the selected secretory machinery reaction is fixed to zero.",
                complex_subunits=complex_subunits,
            )
        )
    for row in oe_rows:
        reaction_id = str(row.get("reaction"))
        rows.append(
            candidate_table_row(
                row,
                perturbation="OE",
                entity_type="secretory_capacity_proxy",
                candidate_id=reaction_id,
                baseline_value=baseline_value,
                fixed_mu=fixed_mu,
                reaction_id=reaction_id,
                simulation_basis="Proxy overexpression: the selected secretory or machine complex capacity is multiplied when possible; non-complex reactions fall back to upper-bound relaxation.",
                complex_subunits=complex_subunits,
            )
        )
    return sorted(rows, key=candidate_sort_key)


def candidate_sort_key(row: dict[str, object]) -> tuple[int, float, str]:
    score = row.get("score")
    numeric_score = float(score) if score is not None else -1e99
    success_rank = 0 if row.get("success") else 1
    return (success_rank, -numeric_score, str(row.get("candidate_id") or ""))


def parse_mu_points(text: str | None) -> list[float]:
    if text is None:
        return list(DEFAULT_TRADEOFF_MUS)
    cleaned = text.strip()
    if not cleaned:
        return []
    values: list[float] = []
    for chunk in re.split(r"[,;\\s]+", cleaned):
        if not chunk:
            continue
        value = float(chunk)
        if value <= 0:
            raise ValueError(f"Tradeoff mu must be positive: {chunk}")
        values.append(value)
    return values


def clean_sequence(sequence: object) -> str:
    cleaned = re.sub(r"[^A-Za-z]", "", str(sequence)).upper().replace("U", "C")
    if cleaned.endswith("X"):
        cleaned = cleaned[:-1]
    return cleaned.replace("B", "").replace("Z", "").replace("J", "").replace("O", "")


def load_targets(root: Path, targets_json: Path | None = None) -> list[TargetSpec]:
    if targets_json is not None:
        payload = json.loads(targets_json.read_text(encoding="utf-8"))
        return [target_from_mapping(item, f"json:{targets_json}") for item in payload.get("targets", [])]
    return [load_opn_default(root), load_hlf_default(root)]


def target_from_mapping(item: dict[str, object], source: str) -> TargetSpec:
    return TargetSpec(
        target_id=str(item["target_id"]),
        protein_id=str(item["protein_id"]),
        mature_sequence=clean_sequence(item["mature_sequence"]),
        leader_sequence=clean_sequence(item.get("leader_sequence", "")),
        signal_peptide_sequence=clean_sequence(item.get("signal_peptide_sequence", "")),
        through_er=bool(item.get("through_er", True)),
        localization=str(item.get("localization", "e")),
        disulfide_sites=int(item.get("disulfide_sites", 0)),
        n_glycosylation_sites=int(item.get("n_glycosylation_sites", 0)),
        o_glycosylation_sites=int(item.get("o_glycosylation_sites", 0)),
        transmembrane=int(item.get("transmembrane", 0)),
        gpi_sites=int(item.get("gpi_sites", 0)),
        cotranslation=int(item.get("cotranslation", 0)),
        source=source,
    )


def load_opn_default(root: Path) -> TargetSpec:
    csv_path = root / "Data" / "pcSecPichia" / "TargetProtein_OPN_candidates.csv"
    frame = pd.read_csv(csv_path)
    row = frame.loc[frame["abbreviation"] == "OPN_ALPHA_FULL_PROJECT"].iloc[0]

    # 用户提供的毕赤酵母工程序列（取代从 CSV 推导的老序列）
    # 信号肽: MRFPSIFTAVLFAASSALA (酿酒酵母 alpha-factor 信号肽, 19aa)
    # 引导肽: APVNTTTEDETAQIPAEAVIGYSDLEGDFDVAVLPFSNSTNNGLLFINTTIASIAAKEEGVSLEKR (66aa, 无 STE13 EAE spacer)
    # 成熟蛋白: IPVKQADSGSSEEKQLYNKYPDAVATWLNPDPSQKQNLLAPQNAVSSEETNDFKQETLPSKSNESHDHMDDMDDEDDDDHVDSQDSIDSNDSDDVDDTDDSHQSDESHHSDESDELVTDFPTDLPATEVFTPVVPTVDTYDGRGDSVVYGLRSKSKKFRRPDIQYPDATDEDITSHMESEELNGAYKAIPVAQDLNAPSDWDSRGKDSYETSQLDDQSAETHSHKQSRLYKRKANDESNEHSDVIDSQELSKVSREFHSHEFHSHEDMLVVDPKSKEEDKHLKFRISHELDSASSEVN (298aa, 与 UniProt P10451 天然人源 OPN 一致)
    signal = "MRFPSIFTAVLFAASSALA"
    pro_leader = "APVNTTTEDETAQIPAEAVIGYSDLEGDFDVAVLPFSNSTNNGLLFINTTIASIAAKEEGVSLEKR"
    mature = "IPVKQADSGSSEEKQLYNKYPDAVATWLNPDPSQKQNLLAPQNAVSSEETNDFKQETLPSKSNESHDHMDDMDDEDDDDHVDSQDSIDSNDSDDVDDTDDSHQSDESHHSDESDELVTDFPTDLPATEVFTPVVPTVDTYDGRGDSVVYGLRSKSKKFRRPDIQYPDATDEDITSHMESEELNGAYKAIPVAQDLNAPSDWDSRGKDSYETSQLDDQSAETHSHKQSRLYKRKANDESNEHSDVIDSQELSKVSREFHSHEFHSHEDMLVVDPKSKEEDKHLKFRISHELDSASSEVN"

    return TargetSpec(
        target_id="OPN_ALPHA_FULL_PROJECT",
        protein_id="OPN_ALPHA_FULL_PROJECT",
        mature_sequence=mature,
        leader_sequence=signal + pro_leader,
        signal_peptide_sequence=signal,
        through_er=True,
        localization="e",
        disulfide_sites=int(row["Disulfide site"]),
        n_glycosylation_sites=int(row["N-glycosylation site"]),
        o_glycosylation_sites=int(row["O-linked glycisylation "]),
        transmembrane=int(row["Transmembrane"]),
        gpi_sites=int(row["GPI site"]),
        cotranslation=int(row["Cotranslation"]),
        source="用户提供: OPN (接头 85aa + 成熟 298aa, 与 UniProt P10451 一致)",
    )


def load_hlf_default(root: Path) -> TargetSpec:
    xlsx_path = root / "Data" / "pcSecPichia" / "TargetProtein.xlsx"
    frame = pd.read_excel(xlsx_path, sheet_name="protein_info")
    row = frame.loc[frame["abbreviation"] == "hLF"].iloc[0]
    signal = "MKLVFLVLLFLGALGLCLA"
    full_sequence = clean_sequence(
        "MKLVFLVLLFLGALGLCLAGRRRSVQWCAVSQPEATKCFQWQRNMRKVRGPPVSCIKRDSPIQCIQAIAENRADAVTLDGGFIYEAGLAPYKLRPVAAEVYGTERQPRTHYYAVAVVKKGGSFQLNELQGLKSCHTGLRRTAGWNVPIGTLRPFLNWTGPPEPIEAAVARFFSASCVPGADKGQFPNLCRLCAGTGENKCAFSSQEPYFSYSGAFKCLRDGAGDVAFIRESTVFEDLSDEAERDEYELLCPDNTRKPVDKFKDCHLARVPSHAVVARSVNGKEDAIWNLLRQAQEKFGKDKSPKFQLFGSPSGQKDLLFKDSAIGFSRVPPRIDSGLYLGSGYFTAIQNLRKSEEEVAARRARVVWCAVGEQELRKCNQWSGLSEGSVTCSSASTTEDCIALVLKGEADAMSLDGGYVYTAGKCGLVPVLAENYKSQQSSDPDPNCVDRPVEGYLAVAVVRRSDTSLTWNSVKGKKSCHTAVDRTAGWNIPMGLLFNQTGSCKFDEYFSQSCAPGSDPRSNLCALCIGDEQGENKCVPNSNERYYGYTGAFRCLAENAGDVAFVKDVTVLQNTDGNNNEAWAKDLKLADFALLCLDGKRKPVTEARSCHLAMAPNHAVVSRMDKVERLKQVLLHQQAKFGRNGSDCPDKFCLFQSETKNLLFNDNTECLARLHGKTTYEKYLGPQYVAGITNLKKCSTSPLLEACEFLRK"
    )
    if not full_sequence.startswith(signal):
        raise ValueError("Project-defined hLF sequence does not start with the expected human native signal peptide.")
    return TargetSpec(
        target_id="hLF",
        protein_id="hLF",
        mature_sequence=full_sequence[len(signal) :],
        leader_sequence=signal,
        signal_peptide_sequence=signal,
        through_er=bool(row["ThroughER"]),
        localization=str(row["Localization"]),
        disulfide_sites=int(row["Disulfide site"]),
        n_glycosylation_sites=int(row["N-glycosylation site"]),
        o_glycosylation_sites=int(row["O-linked glycisylation "]),
        transmembrane=int(row["Transmembrane"]),
        gpi_sites=int(row["GPI site"]),
        cotranslation=int(row.get("Cotranslation", 0)) if "Cotranslation" in row.index and not pd.isna(row.get("Cotranslation")) else 0,
        source="用户提供: hLF native signal (19aa) + mature hLF (691aa); PTM counts referenced from TargetProtein.xlsx",
    )


def longest_common_suffix(sequences: list[str]) -> str:
    if not sequences:
        return ""
    reversed_sequences = [sequence[::-1] for sequence in sequences]
    common_reversed = []
    for chars in zip(*reversed_sequences):
        if len(set(chars)) != 1:
            break
        common_reversed.append(chars[0])
    return "".join(common_reversed)[::-1]


def target_features(target: TargetSpec) -> dict[str, object]:
    sequence = target.full_sequence
    mature = target.mature_sequence
    aa_counts = Counter(sequence)
    n_glycan_motifs = [m.start() + 1 for m in re.finditer(r"N[^P][ST][^P]", mature)]
    return {
        "target_id": target.target_id,
        "protein_id": target.protein_id,
        "source": target.source,
        "full_length": len(sequence),
        "mature_length": len(mature),
        "leader_length": len(target.leader_sequence),
        "signal_peptide_length": len(target.signal_peptide_sequence),
        "valid_mature_sequence": bool(AA_PATTERN.fullmatch(mature)),
        "protein_mw": round(18.0 + sum(AA_MW[aa] * aa_counts[aa] for aa in AA_MW), 3),
        "cysteine_count": mature.count("C"),
        "inferred_disulfide_pairs_from_cysteines": mature.count("C") // 2,
        "declared_disulfide_sites": target.disulfide_sites,
        "n_glycosylation_motif_positions": n_glycan_motifs,
        "declared_n_glycosylation_sites": target.n_glycosylation_sites,
        "ser_thr_count": mature.count("S") + mature.count("T"),
        "declared_o_glycosylation_sites": target.o_glycosylation_sites,
    }


def target_reaction_plan(target: TargetSpec) -> dict[str, object]:
    protein_id = target.protein_id
    reactions: list[dict[str, str]] = []

    def add(reaction_id: str, stage: str, source: str) -> None:
        reactions.append({"reaction_id": reaction_id, "stage": stage, "source_function": source})

    for suffix in (
        "Post_translation_PSTA_sec_RAC_complex",
        "Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex",
        "Post_translation_PSTA_sec_SEC61SEC63C_complex",
        "Post_translation_PSTA_sec_BIP_NEFS_complex",
        "Post_translation_TC_sec_SPC_complex",
        "export_sp_to_c",
    ):
        add(f"{protein_id}_{suffix}", "translocation", "Post_Translation_translocation")
    if target.disulfide_sites > 0:
        add(f"{protein_id}_DSB_sec_BIP_NEFS_complex", "folding", "addDSB")
        add(f"{protein_id}_DSB_PDI_II_sec_PDI1_ERV2_Ero1p_complex", "folding", "addDSB")
    if target.o_glycosylation_sites > 0:
        add(f"{protein_id}_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex", "folding", "addOG")
    if target.n_glycosylation_sites > 0:
        for suffix in (
            "ERNG_NG_sec_OSTC_complex",
            "ERNG_FLI_NG_sec_Cwh41p_complex",
            "ERNG_FLII_NG_sec_Rot2p_complex",
            "ERNG_FLIII_NG_sec_Rot2p_complex",
            "ERNG_FLIV_NG_sec_Rot2p_complex",
        ):
            add(f"{protein_id}_{suffix}", "folding", "addNG")
    add_misfolding_plan(target, add)
    for suffix in (
        "COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex",
        "COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Shl23p_Lst1p_Erv29p_Bet1p_Bos1p_complex",
        "COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex",
        "COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex",
    ):
        add(f"{protein_id}_{suffix}", "er_to_golgi", "coat_other")
    if target.n_glycosylation_sites > 0:
        for suffix in (
            "GLNG_Golgi_N_linked_glycosylation_I_sec_Och1p_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex",
            "GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pA_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pB_complex",
            "GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pC_complex",
        ):
            add(f"{protein_id}_{suffix}", "golgi_processing", "golgiProcessing_N_PP")
    if target.o_glycosylation_sites > 0:
        add(f"{protein_id}_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex", "golgi_processing", "golgiProcessing_O")
    add(f"{protein_id}_Mature", "maturation", "mature")
    add(
        f"{protein_id}_HDSVI_sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex",
        "final_transport",
        "transportFromGolgiToS",
    )
    add(f"{protein_id}_HDSVII_sec_Vps1p_Chc1p_Clc1p_complex", "final_transport", "transportFromGolgiToS")
    add(f"r_{protein_id}_peptide_translation", "translation", "addTranslationRxns")
    if target.signal_peptide_sequence:
        add(f"r_{protein_id}_SP_degradation", "degradation", "addDegradationRxns")
    add(f"r_{protein_id}_subunit_degradation", "degradation", "addDegradationRxns")
    add(f"{protein_id} exchange", "exchange", "addTargetProtein")
    return {
        "target_id": target.target_id,
        "protein_id": protein_id,
        "reaction_count": len(reactions),
        "stage_counts": dict(Counter(item["stage"] for item in reactions)),
        "reactions": reactions,
        "formal_pcsec_simulation_status": formal_status(target),
    }


def add_misfolding_plan(target: TargetSpec, add) -> None:
    protein_id = target.protein_id
    add(f"{protein_id}_misfold_ERAD_sec_Kar2p_complex", "misfolding", "addMisfold")
    add(
        f"{protein_id}_ERAD2A_sec_Pdi1p_complex" if target.disulfide_sites > 0 else f"{protein_id}_ERAD2B",
        "misfolding",
        "addMisfold",
    )
    if target.n_glycosylation_sites > 0:
        add(f"{protein_id}_ERAD3A_sec_Mns1p_complex", "misfolding", "addMisfold")
        add(f"{protein_id}_ERAD4A_sec_Mnl1p_Pdi1p_complex", "misfolding", "addMisfold")
    else:
        add(f"{protein_id}_ERAD3B", "misfolding", "addMisfold")
        add(f"{protein_id}_ERAD4B", "misfolding", "addMisfold")
    add(f"{protein_id}_ERAD5A" if target.gpi_sites > 0 else f"{protein_id}_ERAD5B", "misfolding", "addMisfold")
    add(f"{protein_id}_ERADL_sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex", "misfolding", "addMisfold")
    add(f"{protein_id}_ERADL_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex", "misfolding", "addMisfold")
    if target.n_glycosylation_sites > 0:
        add(f"{protein_id}_ERAD7A_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex", "misfolding", "addMisfold")
    elif target.o_glycosylation_sites > 0:
        add(f"{protein_id}_ERAD7B_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex", "misfolding", "addMisfold")
    else:
        add(f"{protein_id}_ERAD7C_sec_Dsk2p_Rad23p_Uba1p_complex", "misfolding", "addMisfold")
    add(f"{protein_id}_degradation_misfolding_c", "misfolding", "addMisfold")
    if target.disulfide_sites > 0:
        add(f"{protein_id}_cycle_accumulation_sec_pdi1p_ero1p_complex", "misfolding", "addMisfold")
    else:
        add(f"{protein_id}_cycle_accumulation", "misfolding", "addMisfold")
    add(f"{protein_id}_cycle_accumulation_sec_acc_Kar2p_complex", "misfolding", "addMisfold")
    add(f"{protein_id}_dilution_misfolding_er", "misfolding", "addMisfold")


def build_supported_target_model(model: CobraModel, target: TargetSpec, amino_acids: AminoAcidStoichiometry) -> TargetModelBuildResult:
    if is_opn_like_supported_target(target):
        stoichiometries = opn_like_target_stoichiometries(target, amino_acids)
        status_detail = "OPN-like target reactions were appended to S/lb/ub; pcSec reference constraints are solved as a separate smoke step."
    elif is_soluble_secretory_supported_target(target):
        stoichiometries = soluble_secretory_target_stoichiometries(target, amino_acids)
        status_detail = "Soluble extracellular secretory target reactions, including DSB/N-glycan branches, were appended to S/lb/ub; pcSec reference constraints are solved as a separate smoke step."
    else:
        return TargetModelBuildResult(
            status="unsupported_target_branch",
            supported=False,
            reason="Only post-translational soluble extracellular targets without transmembrane/GPI branches have stoichiometry ported in this scratch prototype.",
            model=None,
            exchange_reaction_id=None,
            added_reaction_count=0,
            added_metabolite_count=0,
        )

    original_reaction_count = len(model.rxns)
    original_metabolite_count = len(model.mets)
    target_model = model
    for reaction_id, stoich in stoichiometries:
        target_model = target_model.add_reaction(reaction_id, stoich, lower_bound=0.0, upper_bound=1000.0)
    exchange_reaction_id = f"{target.protein_id} exchange"
    return TargetModelBuildResult(
        status="stoichiometric_target_model_built",
        supported=True,
        reason=status_detail,
        model=target_model,
        exchange_reaction_id=exchange_reaction_id,
        added_reaction_count=len(target_model.rxns) - original_reaction_count,
        added_metabolite_count=len(target_model.mets) - original_metabolite_count,
    )


def is_opn_like_supported_target(target: TargetSpec) -> bool:
    return (
        target.through_er
        and target.localization == "e"
        and target.disulfide_sites == 0
        and target.n_glycosylation_sites == 0
        and target.o_glycosylation_sites > 0
        and target.transmembrane == 0
        and target.gpi_sites == 0
        and target.cotranslation == 0
    )


def is_soluble_secretory_supported_target(target: TargetSpec) -> bool:
    return (
        target.through_er
        and target.localization == "e"
        and target.transmembrane == 0
        and target.gpi_sites == 0
        and target.cotranslation == 0
    )


def opn_like_target_stoichiometries(
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
) -> list[tuple[str, dict[str, float]]]:
    protein_id = target.protein_id
    full_sequence = target.full_sequence
    signal_length = len(target.signal_peptide_sequence)
    subunit_sequence = full_sequence[signal_length:] if signal_length else full_sequence
    entries: list[tuple[str, dict[str, float]]] = []
    entries.extend(opn_like_translocation_stoichiometries(target))
    entries.extend(opn_like_og_misfolding_stoichiometries(target))
    entries.extend(opn_like_er_golgi_transport_stoichiometries(target))
    entries.append((f"r_{protein_id}_peptide_translation", amino_acids.translation(protein_id, full_sequence)))
    if signal_length:
        entries.append((f"r_{protein_id}_SP_degradation", amino_acids.signal_peptide_degradation(protein_id, target.signal_peptide_sequence)))
    entries.append((f"r_{protein_id}_subunit_degradation", amino_acids.subunit_degradation(protein_id, subunit_sequence)))
    entries.append((f"{protein_id} exchange", {f"{protein_id}_folding[{target.localization}]": -1.0}))
    return entries


def soluble_secretory_target_stoichiometries(
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
) -> list[tuple[str, dict[str, float]]]:
    protein_id = target.protein_id
    full_sequence = target.full_sequence
    signal_length = len(target.signal_peptide_sequence)
    subunit_sequence = full_sequence[signal_length:] if signal_length else full_sequence
    entries: list[tuple[str, dict[str, float]]] = []

    entries.extend(opn_like_translocation_stoichiometries(target))
    peptide_name = protein_id

    dsb_entries, peptide_name = dsb_stoichiometries(target, peptide_name)
    entries.extend(dsb_entries)

    og_entries, peptide_name = og_stoichiometries(target, peptide_name)
    entries.extend(og_entries)

    ng_entries, peptide_name = ng_stoichiometries(target, peptide_name)
    entries.extend(ng_entries)

    misfold_entries, peptide_name = soluble_misfolding_stoichiometries(target, peptide_name)
    entries.extend(misfold_entries)

    coat_entries, peptide_name = coat_other_stoichiometries(target, peptide_name)
    entries.extend(coat_entries)

    golgi_n_entries, peptide_name = golgi_n_stoichiometries(target, peptide_name)
    entries.extend(golgi_n_entries)

    golgi_o_entries, peptide_name = golgi_o_stoichiometries(target, peptide_name)
    entries.extend(golgi_o_entries)

    mature_entries, peptide_name = mature_stoichiometries(target, peptide_name)
    entries.extend(mature_entries)

    entries.extend(transport_to_secretory_stoichiometries(target, peptide_name))
    entries.append((f"r_{protein_id}_peptide_translation", amino_acids.translation(protein_id, full_sequence)))
    if signal_length:
        entries.append((f"r_{protein_id}_SP_degradation", amino_acids.signal_peptide_degradation(protein_id, target.signal_peptide_sequence)))
    entries.append((f"r_{protein_id}_subunit_degradation", amino_acids.subunit_degradation(protein_id, subunit_sequence)))
    entries.append((f"{protein_id} exchange", {f"{protein_id}_folding[{target.localization}]": -1.0}))
    return entries


def dsb_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    if target.disulfide_sites <= 0:
        return [], peptide_name
    protein_id = target.protein_id
    length = len(target.full_sequence) / 40.0
    entries = [
        (
            f"{protein_id}_DSB_sec_BIP_NEFS_complex",
            {
                f"{peptide_name}[er]": -1.0,
                "atp[er]": -length,
                "h2o[er]": -length,
                f"{peptide_name}_Kar2ATPcplx[er]": 1.0,
                "adp[er]": length,
                "h[er]": length,
                "pi[er]": length,
            },
        ),
        (
            f"{protein_id}_DSB_PDI_II_sec_PDI1_ERV2_Ero1p_complex",
            {
                f"{peptide_name}_Kar2ATPcplx[er]": -1.0,
                "PDI-ox[er]": -float(target.disulfide_sites),
                f"{peptide_name}_DSB[er]": 1.0,
                "PDI[er]": float(target.disulfide_sites),
            },
        ),
    ]
    return entries, f"{peptide_name}_DSB"


def og_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    if target.o_glycosylation_sites <= 0:
        return [], peptide_name
    protein_id = target.protein_id
    entries = [
        (
            f"{protein_id}_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex",
            {
                f"{peptide_name}[er]": -1.0,
                "dolmanp[er]": -float(target.o_glycosylation_sites),
                f"{peptide_name}_OG_M1[er]": 1.0,
                "dolp[er]": float(target.o_glycosylation_sites),
            },
        )
    ]
    return entries, f"{peptide_name}_OG_M1"


def ng_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    if target.n_glycosylation_sites <= 0:
        return [], peptide_name
    protein_id = target.protein_id
    ng = float(target.n_glycosylation_sites)
    entries = [
        (f"{protein_id}_ERNG_NG_sec_OSTC_complex", {f"{peptide_name}[er]": -1.0, "g3m8mpdol[er]": -ng, f"{peptide_name}_G3M9[er]": 1.0, "dolp[er]": ng}),
        (f"{protein_id}_ERNG_FLI_NG_sec_Cwh41p_complex", {f"{peptide_name}_G3M9[er]": -1.0, "h2o[er]": -ng, f"{peptide_name}_G2M9[er]": 1.0, "glc_D[er]": ng}),
        (f"{protein_id}_ERNG_FLII_NG_sec_Rot2p_complex", {f"{peptide_name}_G2M9[er]": -1.0, "h2o[er]": -ng, f"{peptide_name}_G1M9[er]": 1.0, "glc_D[er]": ng}),
        (f"{protein_id}_ERNG_FLIII_NG_sec_Rot2p_complex", {f"{peptide_name}_G1M9[er]": -1.0, "h2o[er]": -ng, f"{peptide_name}_M9[er]": 1.0, "glc_D[er]": ng}),
        (f"{protein_id}_ERNG_FLIV_NG_sec_Mns1p_complex", {f"{peptide_name}_M9[er]": -1.0, "h2o[er]": -ng, f"{peptide_name}_M8[er]": 1.0, "man[er]": ng}),
    ]
    return entries, f"{peptide_name}_M8"


def opn_like_translocation_stoichiometries(target: TargetSpec) -> list[tuple[str, dict[str, float]]]:
    protein_id = target.protein_id
    length_factor = matlab_round_positive(len(target.full_sequence) / 40.0)
    return [
        (
            f"{protein_id}_Post_translation_PSTA_sec_RAC_complex",
            {f"{protein_id}_peptide[c]": -1.0, f"{protein_id}_translocate_1[c]": 1.0},
        ),
        (
            f"{protein_id}_Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex",
            {
                f"{protein_id}_translocate_1[c]": -1.0,
                "atp[c]": -1.0,
                "h2o[c]": -1.0,
                f"{protein_id}_translocate_2[c]": 1.0,
                "adp[c]": 1.0,
                "h[c]": 1.0,
                "pi[c]": 1.0,
            },
        ),
        (
            f"{protein_id}_Post_translation_PSTA_sec_SEC61SEC63C_complex",
            {f"{protein_id}_translocate_2[c]": -1.0, f"{protein_id}_translocate_3[c]": 1.0},
        ),
        (
            f"{protein_id}_Post_translation_PSTA_sec_BIP_NEFS_complex",
            {
                f"{protein_id}_translocate_3[c]": -1.0,
                "atp[c]": -float(length_factor),
                "h2o[c]": -float(length_factor),
                f"{protein_id}[er]": 1.0,
                "adp[c]": float(length_factor),
                "h[c]": float(length_factor),
                "pi[c]": float(length_factor),
            },
        ),
        (
            f"{protein_id}_Post_translation_TC_sec_SPC_complex",
            {
                f"{protein_id}_translocate_3[c]": -1.0,
                "h2o[c]": -1.0,
                f"{protein_id}[er]": 1.0,
                f"{protein_id}_sp[er]": 1.0,
            },
        ),
        (f"{protein_id}_export_sp_to_c", {f"{protein_id}_sp[er]": -1.0, f"{protein_id}_sp[c]": 1.0}),
    ]


def opn_like_og_misfolding_stoichiometries(target: TargetSpec) -> list[tuple[str, dict[str, float]]]:
    protein_id = target.protein_id
    length_factor = float(matlab_round_positive(len(target.full_sequence) / 40.0))
    accumulation_factor = 10.0 * length_factor
    man = float(target.o_glycosylation_sites)
    peptide_name = protein_id
    entries: list[tuple[str, dict[str, float]]] = []
    if target.o_glycosylation_sites > 0:
        entries.append(
            (
                f"{protein_id}_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex",
                {
                    f"{peptide_name}[er]": -1.0,
                    "dolmanp[er]": -float(target.o_glycosylation_sites),
                    f"{peptide_name}_OG_M1[er]": 1.0,
                    "dolp[er]": float(target.o_glycosylation_sites),
                },
            )
        )
        peptide_name = f"{protein_id}_OG_M1"
    entries.extend(
        [
            (
                f"{protein_id}_misfold_ERAD_sec_Kar2p_complex",
                {
                    f"{peptide_name}[er]": -1.0,
                    "atp[er]": -length_factor,
                    "h2o[er]": -length_factor,
                    f"{peptide_name}_misf[er]": 1.0,
                    "adp[er]": length_factor,
                    "h[er]": length_factor,
                    "pi[er]": length_factor,
                },
            ),
            (f"{protein_id}_ERAD2B", {f"{peptide_name}_misf[er]": -1.0, f"{peptide_name}_misf_G1[er]": 1.0}),
            (f"{protein_id}_ERAD3B", {f"{peptide_name}_misf_G1[er]": -1.0, f"{peptide_name}_misf_G2[er]": 1.0}),
            (f"{protein_id}_ERAD4B", {f"{peptide_name}_misf_G2[er]": -1.0, f"{peptide_name}_misf_G3[er]": 1.0}),
            (f"{protein_id}_ERAD5B", {f"{peptide_name}_misf_G3[er]": -1.0, f"{peptide_name}_misf_G4[er]": 1.0}),
            (
                f"{protein_id}_ERADL_sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex",
                {f"{peptide_name}_misf_G4[er]": -1.0, f"{peptide_name}_misf_G5[er]": 1.0},
            ),
            (
                f"{protein_id}_ERADL_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex",
                {
                    f"{peptide_name}_misf_G5[er]": -1.0,
                    "Ubiquitin_for_Transfer[c]": -8.0,
                    f"{peptide_name}_misf_G6[c]": 1.0,
                    "Ubiquitin[c]": 8.0,
                },
            ),
            (
                f"{protein_id}_ERAD7B_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex",
                {f"{peptide_name}_misf_G6[c]": -1.0, f"{protein_id}_misfolding[c]": 1.0, "man[er]": man},
            ),
            (f"{protein_id}_degradation_misfolding_c", {f"{protein_id}_misfolding[c]": -1.0, f"{protein_id}_subunit[c]": 1.0}),
            (f"{protein_id}_cycle_accumulation", {f"{peptide_name}_misf[er]": -1.0, f"{peptide_name}_misf2[er]": 1.0}),
            (
                f"{protein_id}_cycle_accumulation_sec_acc_Kar2p_complex",
                {
                    f"{peptide_name}_misf2[er]": -1.0,
                    "atp[er]": -accumulation_factor,
                    "h2o[er]": -accumulation_factor,
                    f"{peptide_name}_misfolding_acc[er]": 1.0,
                    "adp[er]": accumulation_factor,
                    "h[er]": accumulation_factor,
                    "pi[er]": accumulation_factor,
                },
            ),
            (f"{protein_id}_dilution_misfolding_er", {f"{peptide_name}_misfolding_acc[er]": -1.0}),
        ]
    )
    return entries


def soluble_misfolding_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    protein_id = target.protein_id
    length_factor = float(matlab_round_positive(len(target.full_sequence) / 40.0))
    ng = float(target.n_glycosylation_sites)
    og = float(target.o_glycosylation_sites)
    dsb = float(target.disulfide_sites)
    if ng > 0:
        base = peptide_name.split("_M8", 1)[0]
        peptide = f"{base}_M9"
        man = 7.0 * ng + og
        nac = 2.0 * ng
    else:
        peptide = peptide_name
        man = og
        nac = 0.0

    entries: list[tuple[str, dict[str, float]]] = [
        (
            f"{protein_id}_misfold_ERAD_sec_Kar2p_complex",
            {
                f"{peptide}[er]": -1.0,
                "atp[er]": -length_factor,
                "h2o[er]": -length_factor,
                f"{peptide}_misf[er]": 1.0,
                "adp[er]": length_factor,
                "h[er]": length_factor,
                "pi[er]": length_factor,
            },
        )
    ]
    if dsb > 0:
        entries.append(
            (
                f"{protein_id}_ERAD2A_sec_Pdi1p_complex",
                {
                    f"{peptide}_misf[er]": -1.0,
                    "gthrd[er]": -2.0 * dsb,
                    f"{peptide}_misf_G1[er]": 1.0,
                    "gthox[er]": dsb,
                    "h[er]": 2.0 * dsb,
                },
            )
        )
    else:
        entries.append((f"{protein_id}_ERAD2B", {f"{peptide}_misf[er]": -1.0, f"{peptide}_misf_G1[er]": 1.0}))

    if ng > 0:
        entries.extend(
            [
                (f"{protein_id}_ERAD3A_sec_Mns1p_complex", {f"{peptide}_misf_G1[er]": -1.0, "h2o[er]": -ng, f"{peptide}_misf_G2[er]": 1.0, "man[er]": ng}),
                (f"{protein_id}_ERAD4A_sec_Mnl1p_Pdi1p_complex", {f"{peptide}_misf_G2[er]": -1.0, "h2o[er]": -ng, f"{peptide}_misf_G3[er]": 1.0, "man[er]": ng}),
            ]
        )
    else:
        entries.extend(
            [
                (f"{protein_id}_ERAD3B", {f"{peptide}_misf_G1[er]": -1.0, f"{peptide}_misf_G2[er]": 1.0}),
                (f"{protein_id}_ERAD4B", {f"{peptide}_misf_G2[er]": -1.0, f"{peptide}_misf_G3[er]": 1.0}),
            ]
        )

    entries.extend(
        [
            (f"{protein_id}_ERAD5B", {f"{peptide}_misf_G3[er]": -1.0, f"{peptide}_misf_G4[er]": 1.0}),
            (
                f"{protein_id}_ERADL_sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex",
                {f"{peptide}_misf_G4[er]": -1.0, f"{peptide}_misf_G5[er]": 1.0},
            ),
            (
                f"{protein_id}_ERADL_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex",
                {
                    f"{peptide}_misf_G5[er]": -1.0,
                    "Ubiquitin_for_Transfer[c]": -8.0,
                    f"{peptide}_misf_G6[c]": 1.0,
                    "Ubiquitin[c]": 8.0,
                },
            ),
        ]
    )
    if ng > 0:
        entries.append(
            (
                f"{protein_id}_ERAD7A_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex",
                {f"{peptide}_misf_G6[c]": -1.0, f"{protein_id}_misfolding[c]": 1.0, "man[er]": man, "acgam1p[c]": nac},
            )
        )
    elif og > 0:
        entries.append(
            (
                f"{protein_id}_ERAD7B_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex",
                {f"{peptide}_misf_G6[c]": -1.0, f"{protein_id}_misfolding[c]": 1.0, "man[er]": man},
            )
        )
    else:
        entries.append((f"{protein_id}_ERAD7C_sec_Dsk2p_Rad23p_Uba1p_complex", {f"{peptide}_misf_G6[c]": -1.0, f"{protein_id}_misfolding[c]": 1.0}))

    entries.append((f"{protein_id}_degradation_misfolding_c", {f"{protein_id}_misfolding[c]": -1.0, f"{protein_id}_subunit[c]": 1.0}))
    if dsb > 0:
        entries.append(
            (
                f"{protein_id}_cycle_accumulation_sec_pdi1p_ero1p_complex",
                {
                    f"{peptide}_misf[er]": -1.0,
                    "gthrd[er]": -20.0 * dsb,
                    "o2[er]": -10.0 * dsb,
                    "h2o2[er]": 10.0 * dsb,
                    "gthox[er]": 10.0 * dsb,
                    f"{peptide}_misf2[er]": 1.0,
                },
            )
        )
    else:
        entries.append((f"{protein_id}_cycle_accumulation", {f"{peptide}_misf[er]": -1.0, f"{peptide}_misf2[er]": 1.0}))
    entries.extend(
        [
            (
                f"{protein_id}_cycle_accumulation_sec_acc_Kar2p_complex",
                {
                    f"{peptide}_misf2[er]": -1.0,
                    "atp[er]": -10.0 * length_factor,
                    "h2o[er]": -10.0 * length_factor,
                    f"{peptide}_misfolding_acc[er]": 1.0,
                    "adp[er]": 10.0 * length_factor,
                    "h[er]": 10.0 * length_factor,
                    "pi[er]": 10.0 * length_factor,
                },
            ),
            (f"{protein_id}_dilution_misfolding_er", {f"{peptide}_misfolding_acc[er]": -1.0}),
        ]
    )
    next_peptide = f"{peptide.split('_M9', 1)[0]}_M8" if ng > 0 else peptide
    return entries, next_peptide


def coat_other_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    protein_id = target.protein_id
    entries = [
        (
            f"{protein_id}_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex",
            {f"{peptide_name}[er]": -1.0, "gtp[er]": -1.0, "h2o[er]": -1.0, f"{peptide_name}_COP_coated[er]": 1.0, "gdp[er]": 1.0, "h[er]": 1.0, "pi[er]": 1.0},
        ),
        (
            f"{protein_id}_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Shl23p_Lst1p_Erv29p_Bet1p_Bos1p_complex",
            {f"{peptide_name}[er]": -1.0, "gtp[er]": -1.0, "h2o[er]": -1.0, f"{peptide_name}_COP_coated[er]": 1.0, "gdp[er]": 1.0, "h[er]": 1.0, "pi[er]": 1.0},
        ),
        (
            f"{protein_id}_COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex",
            {f"{peptide_name}_COP_coated[er]": -1.0, f"{peptide_name}_COP_coated[c]": 1.0},
        ),
        (
            f"{protein_id}_COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex",
            {f"{peptide_name}_COP_coated[c]": -1.0, "gtp[c]": -1.0, "h2o[c]": -1.0, f"{peptide_name}[g]": 1.0, "gdp[c]": 1.0, "h[c]": 1.0, "pi[c]": 1.0},
        ),
    ]
    return entries, peptide_name


def golgi_n_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    if target.n_glycosylation_sites <= 0:
        return [], peptide_name
    protein_id = target.protein_id
    ng = float(target.n_glycosylation_sites)
    entries = [
        (f"{protein_id}_GLNG_Golgi_N_linked_glycosylation_I_sec_Och1p_complex", {f"{peptide_name}[g]": -1.0, "gdpmann[g]": -ng, f"{peptide_name}_GNG_G1[g]": 1.0, "gdp[g]": ng}),
        (f"{protein_id}_GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex", {f"{peptide_name}_GNG_G1[g]": -1.0, "gdpmann[g]": -(9.0 * 0.2 * ng), f"{peptide_name}_GNG_G2[g]": 1.0, "gdp[g]": 9.0 * 0.2 * ng}),
        (f"{protein_id}_GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex", {f"{peptide_name}_GNG_G2[g]": -1.0, "gdpmann[g]": -(30.0 * 0.2 * ng), f"{peptide_name}_GNG_G3[g]": 1.0, "gdp[g]": 30.0 * 0.2 * ng}),
        (f"{protein_id}_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pA_complex", {f"{peptide_name}_GNG_G3[g]": -1.0, "gdpmann[g]": -(0.2 * ng), f"{peptide_name}_GNG_G4[g]": 1.0, "gdp[g]": 0.2 * ng}),
        (f"{protein_id}_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pB_complex", {f"{peptide_name}_GNG_G3[g]": -1.0, "gdpmann[g]": -(0.2 * ng), f"{peptide_name}_GNG_G4[g]": 1.0, "gdp[g]": 0.2 * ng}),
        (f"{protein_id}_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pC_complex", {f"{peptide_name}_GNG_G3[g]": -1.0, "gdpmann[g]": -(0.2 * ng), f"{peptide_name}_GNG_G4[g]": 1.0, "gdp[g]": 0.2 * ng}),
    ]
    return entries, f"{peptide_name}_GNG_G4"


def golgi_o_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    if target.o_glycosylation_sites <= 0:
        return [], peptide_name
    protein_id = target.protein_id
    o_mannose_count = float(2 * target.o_glycosylation_sites)
    entries = [
        (
            f"{protein_id}_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex",
            {f"{peptide_name}[g]": -1.0, "gdpmann[g]": -o_mannose_count, f"{peptide_name}_GOG_G1[g]": 1.0, "gdp[g]": o_mannose_count},
        )
    ]
    return entries, f"{peptide_name}_GOG_G1"


def mature_stoichiometries(target: TargetSpec, peptide_name: str) -> tuple[list[tuple[str, dict[str, float]]], str]:
    protein_id = target.protein_id
    mature_peptide = f"{peptide_name}_mature"
    return [(f"{protein_id}_Mature", {f"{peptide_name}[g]": -1.0, f"{mature_peptide}[g]": 1.0})], mature_peptide


def transport_to_secretory_stoichiometries(target: TargetSpec, peptide_name: str) -> list[tuple[str, dict[str, float]]]:
    protein_id = target.protein_id
    return [
        (
            f"{protein_id}_HDSVI_sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex",
            {f"{peptide_name}[g]": -1.0, "gtp[c]": -1.0, "h2o[c]": -1.0, f"{peptide_name}[ce]": 1.0, "gdp[c]": 1.0, "pi[c]": 1.0, "h[c]": 1.0},
        ),
        (
            f"{protein_id}_HDSVII_sec_Vps1p_Chc1p_Clc1p_complex",
            {f"{peptide_name}[ce]": -1.0, "gtp[c]": -1.0, "h2o[c]": -1.0, f"{protein_id}_folding[e]": 1.0, "gdp[c]": 1.0, "pi[c]": 1.0, "h[c]": 1.0},
        ),
    ]


def opn_like_er_golgi_transport_stoichiometries(target: TargetSpec) -> list[tuple[str, dict[str, float]]]:
    protein_id = target.protein_id
    peptide_name = f"{protein_id}_OG_M1" if target.o_glycosylation_sites > 0 else protein_id
    golgi_peptide = f"{peptide_name}_GOG_G1"
    mature_peptide = f"{golgi_peptide}_mature"
    o_mannose_count = float(2 * target.o_glycosylation_sites)
    return [
        (
            f"{protein_id}_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex",
            {
                f"{peptide_name}[er]": -1.0,
                "gtp[er]": -1.0,
                "h2o[er]": -1.0,
                f"{peptide_name}_COP_coated[er]": 1.0,
                "gdp[er]": 1.0,
                "h[er]": 1.0,
                "pi[er]": 1.0,
            },
        ),
        (
            f"{protein_id}_COPII_normal_ERGL1A_sec_Sec12p_Sar1p_Shl23p_Lst1p_Erv29p_Bet1p_Bos1p_complex",
            {
                f"{peptide_name}[er]": -1.0,
                "gtp[er]": -1.0,
                "h2o[er]": -1.0,
                f"{peptide_name}_COP_coated[er]": 1.0,
                "gdp[er]": 1.0,
                "h[er]": 1.0,
                "pi[er]": 1.0,
            },
        ),
        (
            f"{protein_id}_COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex",
            {f"{peptide_name}_COP_coated[er]": -1.0, f"{peptide_name}_COP_coated[c]": 1.0},
        ),
        (
            f"{protein_id}_COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex",
            {
                f"{peptide_name}_COP_coated[c]": -1.0,
                "gtp[c]": -1.0,
                "h2o[c]": -1.0,
                f"{peptide_name}[g]": 1.0,
                "gdp[c]": 1.0,
                "h[c]": 1.0,
                "pi[c]": 1.0,
            },
        ),
        (
            f"{protein_id}_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex",
            {
                f"{peptide_name}[g]": -1.0,
                "gdpmann[g]": -o_mannose_count,
                f"{golgi_peptide}[g]": 1.0,
                "gdp[g]": o_mannose_count,
            },
        ),
        (f"{protein_id}_Mature", {f"{golgi_peptide}[g]": -1.0, f"{mature_peptide}[g]": 1.0}),
        (
            f"{protein_id}_HDSVI_sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex",
            {
                f"{mature_peptide}[g]": -1.0,
                "gtp[c]": -1.0,
                "h2o[c]": -1.0,
                f"{mature_peptide}[ce]": 1.0,
                "gdp[c]": 1.0,
                "pi[c]": 1.0,
                "h[c]": 1.0,
            },
        ),
        (
            f"{protein_id}_HDSVII_sec_Vps1p_Chc1p_Clc1p_complex",
            {
                f"{mature_peptide}[ce]": -1.0,
                "gtp[c]": -1.0,
                "h2o[c]": -1.0,
                f"{protein_id}_folding[e]": 1.0,
                "gdp[c]": 1.0,
                "pi[c]": 1.0,
                "h[c]": 1.0,
            },
        ),
    ]


def matlab_round_positive(value: float) -> int:
    return int(value + 0.5)


def target_secretion_smoke(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    ko_genes: list[str],
    oe_reactions: list[str],
    fixed_mu: float = 0.10,
    tradeoff_mus: Iterable[float] = DEFAULT_TRADEOFF_MUS,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> dict[str, object]:
    build = build_supported_target_model(model, target, amino_acids)
    payload: dict[str, object] = {
        "build": {
            "status": build.status,
            "supported": build.supported,
            "reason": build.reason,
            "exchange_reaction_id": build.exchange_reaction_id,
            "added_reaction_count": build.added_reaction_count,
            "added_metabolite_count": build.added_metabolite_count,
        },
        "fixed_mu": fixed_mu,
        "formal_pcsec_prediction": False,
        "interpretation": "Stoichiometric smoke is always reported for supported branches; pcSec reference smoke adds writeLPGlc-style protein/secretory coupling when available.",
        "pcsec_constraint_options": {
            "write_ribosome_translation_constraint": write_ribosome_translation_constraint,
            "write_misfolding_constraints": write_misfolding_constraints,
        },
    }
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        return payload

    fixed_model = build.model.with_bounds({"BIOMASS": (fixed_mu, fixed_mu)})
    baseline = solve_maximize(
        fixed_model,
        build.exchange_reaction_id,
        key_reactions=("BIOMASS", "Ex_glc_D", "Ex_o2", build.exchange_reaction_id),
    )
    payload["target_exchange_max"] = asdict(baseline)
    if baseline.success:
        payload["ko_screen_draft"] = run_ko_screen(fixed_model, baseline, ko_genes, build.exchange_reaction_id)
        payload["oe_screen_draft"] = run_oe_screen(fixed_model, baseline, oe_reactions, build.exchange_reaction_id)
    else:
        payload["ko_screen_draft"] = []
        payload["oe_screen_draft"] = []

    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = combined.with_target(target_enzymedata)
    pcsec_result, pcsec_counts = solve_pcsec_maximize(
        fixed_model,
        build.exchange_reaction_id,
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu=fixed_mu,
        key_reactions=("BIOMASS", "Ex_glc_D", "Ex_o2", build.exchange_reaction_id),
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    payload["target_pcsec_reference_max"] = asdict(pcsec_result)
    payload["target_pcsec_constraint_counts"] = pcsec_counts
    payload["target_reaction_coefficient_count"] = len(target_enzymedata.reaction_coefficients)
    payload["pcsec_growth_tradeoff"] = run_pcsec_growth_tradeoff(
        build.model,
        build.exchange_reaction_id,
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu_points=tradeoff_mus,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    if pcsec_result.success:
        pcsec_ko_rows = run_pcsec_ko_screen(
            fixed_model,
            pcsec_result,
            ko_genes,
            build.exchange_reaction_id,
            metabolic=metabolic,
            secretory=target_secretory,
            combined=target_combined,
            mu=fixed_mu,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        pcsec_capacity_ko_rows = run_pcsec_reaction_ko_screen(
            fixed_model,
            pcsec_result,
            oe_reactions,
            build.exchange_reaction_id,
            metabolic=metabolic,
            secretory=target_secretory,
            combined=target_combined,
            mu=fixed_mu,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        pcsec_oe_rows = run_pcsec_oe_screen(
            fixed_model,
            pcsec_result,
            oe_reactions,
            build.exchange_reaction_id,
            metabolic=metabolic,
            secretory=target_secretory,
            combined=target_combined,
            mu=fixed_mu,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        payload["pcsec_ko_screen_draft"] = pcsec_ko_rows
        payload["pcsec_capacity_ko_screen_draft"] = pcsec_capacity_ko_rows
        payload["pcsec_oe_screen_draft"] = pcsec_oe_rows
        payload["candidate_table"] = build_candidate_table(
            pcsec_ko_rows,
            pcsec_capacity_ko_rows,
            pcsec_oe_rows,
            baseline_value=pcsec_result.objective_value,
            fixed_mu=fixed_mu,
            complex_subunits=target_secretory.complex_subunits,
        )
    else:
        payload["pcsec_ko_screen_draft"] = []
        payload["pcsec_capacity_ko_screen_draft"] = []
        payload["pcsec_oe_screen_draft"] = []
        payload["candidate_table"] = []
    return payload


def formal_status(target: TargetSpec) -> dict[str, object]:
    if is_soluble_secretory_supported_target(target):
        return {
            "status": "pcsec_reference_smoke_supported",
            "reason": "Soluble extracellular secretory branch can be appended to S/lb/ub and solved with reference writeLPGlc-style protein/secretory coupling. MATLAB LP/solver alignment is still required before formal production prediction.",
        }
    return {
        "status": "not_formal_yet",
        "reason": "Target requires branches beyond the current soluble extracellular post-translational secretion path.",
    }


def default_ko_genes(model: CobraModel, limit: int) -> list[str]:
    candidates: list[str] = []
    for reaction_id, gr_rule in zip(model.rxns, model.gr_rules):
        if gr_rule and gr_rule != "[]" and reaction_id.endswith("_fwd"):
            candidates.extend(re.findall(r"PAS_[A-Za-z0-9_\-]+", gr_rule))
        if len(candidates) >= limit:
            break
    deduped: list[str] = []
    for gene in candidates:
        if gene in model.gene_index and gene not in deduped:
            deduped.append(gene)
    return deduped[:limit]


def default_oe_reactions(model: CobraModel, limit: int) -> list[str]:
    preferred = [
        "sec_BIP_NEFS_complex_formation",
        "sec_Kar2p_complex_formation",
        "sec_PDI1_ERV2_Ero1p_complex_formation",
        "sec_OSTC_complex_formation",
        "sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex_formation",
        "sec_Cwh41p_complex_formation",
        "sec_Rot2p_complex_formation",
        "sec_Mns1p_complex_formation",
        "sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex_formation",
        "sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex_formation",
        "sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex_formation",
        "sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex_formation",
        "sec_Vps1p_Chc1p_Clc1p_complex_formation",
        "Mach_proteasome_complex_formation",
        "Mach_Ribosome_complex_formation",
    ]
    existing = [reaction_id for reaction_id in preferred if reaction_id in model.reaction_index]
    if len(existing) >= limit:
        return existing[:limit]
    for reaction_id in model.rxns:
        if reaction_id.endswith("_formation") and reaction_id not in existing:
            existing.append(reaction_id)
        if len(existing) >= limit:
            break
    return existing


def format_number(value: object, digits: int = 6) -> str:
    if value is None:
        return "None"
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{numeric:.{digits}g}"


def format_candidate_rows(rows: list[dict[str, object]], perturbation: str, limit: int = 3) -> str:
    selected = [row for row in rows if row.get("perturbation") == perturbation]
    if not selected:
        return "none"
    parts: list[str] = []
    for row in selected[:limit]:
        subunits = row.get("complex_subunit_ids") or []
        subunit_text = ",".join(str(item) for item in subunits[:6]) if subunits else "none"
        parts.append(
            "{candidate} [{effect}, {entity}, process={process}, delta={delta}, rel={rel}, subunits={subunits}]".format(
                candidate=row.get("candidate_id"),
                effect=row.get("effect"),
                entity=row.get("entity_type"),
                process=row.get("process"),
                delta=format_number(row.get("secretion_delta")),
                rel=format_number(row.get("relative_secretion_delta")),
                subunits=subunit_text,
            )
        )
    return "; ".join(parts)


def format_tradeoff_rows(rows: list[dict[str, object]], limit: int = 5) -> str:
    if not rows:
        return "none"
    parts: list[str] = []
    for row in rows[:limit]:
        parts.append(
            "mu={mu}: success={success}, secretion={secretion}, secretion/biomass={ratio}".format(
                mu=format_number(row.get("mu")),
                success=row.get("success"),
                secretion=format_number(row.get("secretion_flux")),
                ratio=format_number(row.get("secretion_per_biomass")),
            )
        )
    return "; ".join(parts)


def write_outputs(report: dict[str, object], output_prefix: str) -> tuple[Path, Path, Path, Path]:
    out = output_dir()
    safe_prefix = re.sub(r"[^A-Za-z0-9_.-]+", "_", output_prefix).strip("._") or "probe"
    summary_path = out / f"{safe_prefix}_summary.json"
    report_path = out / f"{safe_prefix}_REPORT.md"
    candidates_path = out / f"{safe_prefix}_candidates.csv"
    tradeoff_path = out / f"{safe_prefix}_tradeoff.csv"
    summary_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        "# pcSecPichia hLF/OPN Probe Report",
        "",
        f"- Model: `{report['model']['source_file']}`",
        f"- Reactions: `{report['model']['reaction_count']}`",
        f"- Metabolites: `{report['model']['metabolite_count']}`",
        f"- Genes: `{report['model']['gene_count']}`",
        f"- Baseline GEM FBA: `{report['gem_layer']['baseline_growth']['success']}` objective `{report['gem_layer']['baseline_growth']['objective_value']}`",
        "",
        "## Target Inputs",
    ]
    for target in report["targets"]:
        smoke = target["target_secretion_smoke"]
        exchange = smoke.get("target_exchange_max", {})
        exchange_value = exchange.get("objective_value") if isinstance(exchange, dict) else None
        pcsec_exchange = smoke.get("target_pcsec_reference_max", {})
        pcsec_value = pcsec_exchange.get("objective_value") if isinstance(pcsec_exchange, dict) else None
        pcsec_counts = smoke.get("target_pcsec_constraint_counts") or {}
        constraint_options = smoke.get("pcsec_constraint_options") or {}
        pcsec_ko_rows = smoke.get("pcsec_ko_screen_draft") or []
        pcsec_capacity_ko_rows = smoke.get("pcsec_capacity_ko_screen_draft") or []
        pcsec_oe_rows = smoke.get("pcsec_oe_screen_draft") or []
        candidate_rows = smoke.get("candidate_table") or []
        tradeoff_rows = smoke.get("pcsec_growth_tradeoff") or []
        lines.extend(
            [
                "",
                f"### {target['features']['target_id']}",
                f"- Source: `{target['features']['source']}`",
                f"- Mature length: `{target['features']['mature_length']}`",
                f"- Leader length: `{target['features']['leader_length']}`",
                f"- Declared DSB/NG/OG: `{target['features']['declared_disulfide_sites']}` / `{target['features']['declared_n_glycosylation_sites']}` / `{target['features']['declared_o_glycosylation_sites']}`",
                f"- Planned target reactions: `{target['reaction_plan']['reaction_count']}`",
                f"- Formal pcSec status: `{target['reaction_plan']['formal_pcsec_simulation_status']['status']}`",
                f"- Target model build: `{smoke['build']['status']}`",
                f"- Stoichiometric exchange smoke: `{exchange.get('success') if isinstance(exchange, dict) else None}` objective `{exchange_value}`",
                f"- pcSec reference exchange smoke: `{pcsec_exchange.get('success') if isinstance(pcsec_exchange, dict) else None}` objective `{pcsec_value}`",
                f"- pcSec constraints: eq `{pcsec_counts.get('eq_total')}` / ub `{pcsec_counts.get('ub_total')}`",
                f"- Optional constraints: ribosome_translation `{constraint_options.get('write_ribosome_translation_constraint')}` rows `{pcsec_counts.get('ribosome_translation')}`; misfolding `{constraint_options.get('write_misfolding_constraints')}` rows `{pcsec_counts.get('misfolding')}`",
                f"- Target rxnscoef count: `{smoke.get('target_reaction_coefficient_count')}`",
                f"- pcSec KO draft rows: `{len(pcsec_ko_rows)}`; best delta `{best_delta(pcsec_ko_rows)}`",
                f"- pcSec capacity-KO proxy rows: `{len(pcsec_capacity_ko_rows)}`; best delta `{best_delta(pcsec_capacity_ko_rows)}`",
                f"- pcSec OE draft rows: `{len(pcsec_oe_rows)}`; best delta `{best_delta(pcsec_oe_rows)}`",
                f"- Candidate table rows: `{len(candidate_rows)}`",
                f"- Top KO candidates: {format_candidate_rows(candidate_rows, 'KO')}",
                f"- Top OE candidates: {format_candidate_rows(candidate_rows, 'OE')}",
                f"- Growth tradeoff: {format_tradeoff_rows(tradeoff_rows)}",
                f"- Boundary: {target['reaction_plan']['formal_pcsec_simulation_status']['reason']}",
            ]
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "For soluble extracellular targets such as OPN and hLF, target stoichiometry and the reference `writeLPGlc` protein/secretory constraints can now be solved as a smoke test.",
            "This still needs MATLAB LP/solver alignment before the numbers should be treated as formal production predictions.",
            "Branches outside post-translational soluble extracellular secretion remain future work. Ribosome translation and misfolding equality constraints are available as explicit stage3 options and still require MATLAB LP alignment.",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    candidate_rows: list[dict[str, object]] = []
    tradeoff_rows: list[dict[str, object]] = []
    for target in report["targets"]:
        target_id = target["features"]["target_id"]
        smoke = target["target_secretion_smoke"]
        for row in smoke.get("candidate_table") or []:
            candidate_rows.append({"target_id": target_id, **row})
        for row in smoke.get("pcsec_growth_tradeoff") or []:
            tradeoff_rows.append({"target_id": target_id, **row})
    pd.DataFrame(candidate_rows).to_csv(candidates_path, index=False)
    pd.DataFrame(tradeoff_rows).to_csv(tradeoff_path, index=False)
    return summary_path, report_path, candidates_path, tradeoff_path


def best_delta(rows: list[dict[str, object]]) -> float | None:
    values = [row.get("delta_vs_baseline") for row in rows if row.get("delta_vs_baseline") is not None]
    if not values:
        return None
    return float(max(values))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets-json", type=Path)
    parser.add_argument("--ko-limit", type=int, default=3)
    parser.add_argument("--oe-limit", type=int, default=5)
    parser.add_argument("--tradeoff-mu", default="0.05,0.10,0.20")
    parser.add_argument("--enable-ribosome-translation-constraint", action="store_true")
    parser.add_argument("--enable-misfolding-constraint", action="store_true")
    parser.add_argument("--output-prefix", default="probe")
    args = parser.parse_args()

    root = repo_root()
    model = prepare_glucose_model(load_pcsec_pichia_model(root), media_type=4)
    amino_acids = load_aa_stoichiometry(root)
    metabolic = load_metabolic_enzymedata(root)
    secretory = load_secretory_enzymedata(root)
    combined = load_combined_enzymedata(root)
    baseline = solve_maximize(model, "BIOMASS", key_reactions=("Ex_glc_D", "Ex_o2", "BIOMASS"))
    ko_genes = default_ko_genes(model, args.ko_limit)
    oe_reactions = default_oe_reactions(model, args.oe_limit)
    tradeoff_mus = parse_mu_points(args.tradeoff_mu)
    targets = load_targets(root, args.targets_json)

    report = {
        "scope": {
            "scratch_dir": str(output_dir()),
            "write_policy": "created/modified files are restricted to this scratch directory",
            "formal_target_secretion_prediction": False,
        },
        "model": {
            "source_file": model.source_file,
            "reaction_count": len(model.rxns),
            "metabolite_count": len(model.mets),
            "gene_count": len(model.genes),
            "s_shape": list(model.s_matrix.shape),
            "nonempty_gene_rules": sum(1 for item in model.rules if item and item != "[]"),
        },
        "gem_layer": {
            "baseline_growth": asdict(baseline),
            "ko_screen_draft": run_ko_screen(model, baseline, ko_genes, "BIOMASS") if baseline.success else [],
            "oe_screen_draft": run_oe_screen(model, baseline, oe_reactions, "BIOMASS") if baseline.success else [],
            "interpretation": "GEM-layer KO/OE is a solver-path proof. Target secretion KO/OE must wait for formal target pcSec constraints.",
        },
        "targets": [
            {
                "input": asdict(target),
                "features": target_features(target),
                "reaction_plan": target_reaction_plan(target),
                "target_secretion_smoke": target_secretion_smoke(
                    model,
                    target,
                    amino_acids,
                    metabolic,
                    secretory,
                    combined,
                    ko_genes,
                    oe_reactions,
                    tradeoff_mus=tradeoff_mus,
                    write_ribosome_translation_constraint=args.enable_ribosome_translation_constraint,
                    write_misfolding_constraints=args.enable_misfolding_constraint,
                ),
            }
            for target in targets
        ],
        "next_migration_steps": [
            "Align OPN and hLF pcSec reference smoke results against MATLAB local_opn_pichia_glc/writeLPGlc LP and solver output.",
            "Audit hLF DSB and N-glycosylation branch stoichiometry against MATLAB addTargetProtein before allowing formal hLF production predictions.",
            "Align optional writeLPGlc ribosome translation and misfolding equality constraints against MATLAB LP output before enabling them by default.",
            "Expand KO/OE screens beyond smoke-size limits after pcSec OPN/hLF alignment is verified.",
        ],
    }
    summary_path, report_path, candidates_path, tradeoff_path = write_outputs(report, args.output_prefix)
    print(
        json.dumps(
            {
                "summary": str(summary_path),
                "report": str(report_path),
                "candidates": str(candidates_path),
                "tradeoff": str(tradeoff_path),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
