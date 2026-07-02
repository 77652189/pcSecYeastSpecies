from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Literal


BoundRole = Literal[
    "amino_acid",
    "biomass",
    "carbon_source",
    "gas",
    "ion",
    "legacy_compatibility",
    "nitrogen",
    "phosphate",
    "sulfate",
    "vitamin",
    "water_proton_balance",
]
SpecStatus = Literal["active", "planned"]


MINIMAL_INORGANIC_EXCHANGES: tuple[str, ...] = (
    "Ex_nh4",
    "Ex_o2",
    "Ex_pi",
    "Ex_so4",
    "Ex_fe2",
    "Ex_h",
    "Ex_h2o",
    "Ex_na1",
    "Ex_k",
    "Ex_co2",
)

MATLAB_LEGACY_COST_MEDIUM_EXCHANGES: tuple[str, ...] = (
    "Ex_h2o",
    "Ex_nh4",
    "Ex_pi",
    "Ex_h",
    "Ex_so4",
    "Ex_fe2",
    "Ex_co2",
    "Ex_k",
    "Ex_na1",
)

YNB_VITAMINS: tuple[str, ...] = (
    "Ex_btn",
    "Ex_thm",
    "Ex_4abz",
    "Ex_pnto_R",
    "Ex_inost",
    "Ex_nac",
    "Ex_ribflv",
)

CORE_AMINO_ACIDS: tuple[str, ...] = (
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

ALL_AMINO_ACIDS: tuple[str, ...] = CORE_AMINO_ACIDS + (
    "Ex_ala_L",
    "Ex_asn_L",
    "Ex_cys_L",
    "Ex_gln_L",
    "Ex_pro_L",
    "Ex_ser_L",
)


@dataclass(frozen=True)
class MediumBound:
    reaction_id: str
    lower_bound: float | None
    upper_bound: float | None
    role: BoundRole
    source: str
    rationale: str
    confidence: str = "curated"


@dataclass(frozen=True)
class CarbonSourceSpec:
    carbon_source_id: str
    label: str
    status: SpecStatus
    bounds: tuple[MediumBound, ...]
    blocked_reactions: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class CarbonSourceFormulation:
    carbon_source_id: str
    active_uptake_reaction_ids: tuple[str, ...]
    blocked_uptake_reaction_ids: tuple[str, ...]
    candidate_growth_reaction_ids: tuple[str, ...]
    selected_growth_reaction_id: str
    carbon_objective_weights: dict[str, float]
    formulation_status: str
    warnings: tuple[str, ...] = ()
    matlab_alignment_note: str = ""


@dataclass(frozen=True)
class BaseMediumSpec:
    base_medium_id: str
    label: str
    status: SpecStatus
    legacy_media_type: int | None
    supplement_bounds: tuple[MediumBound, ...]
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class CompatibilityOverlaySpec:
    overlay_id: str
    label: str
    status: SpecStatus
    scope: str
    bounds: tuple[MediumBound, ...]
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class MediumConditionSpec:
    condition_id: str
    label: str
    status: SpecStatus
    carbon_source_id: str
    base_medium_id: str
    compatibility_overlay_id: str = "none"
    oxygen_mode: str = "aerobic"
    open_minimal_inorganic_exchanges: bool = True
    block_misfold_dilution: bool = False
    legacy_media_type: int | None = None
    notes: tuple[str, ...] = ()


def list_carbon_source_specs() -> tuple[CarbonSourceSpec, ...]:
    return tuple(_CARBON_SOURCES.values())


def list_carbon_source_formulations() -> tuple[CarbonSourceFormulation, ...]:
    return tuple(_CARBON_SOURCE_FORMULATIONS.values())


def list_base_medium_specs() -> tuple[BaseMediumSpec, ...]:
    return tuple(_BASE_MEDIA.values())


def list_compatibility_overlay_specs() -> tuple[CompatibilityOverlaySpec, ...]:
    return tuple(_COMPATIBILITY_OVERLAYS.values())


def list_medium_condition_specs() -> tuple[MediumConditionSpec, ...]:
    return tuple(_MEDIUM_CONDITIONS.values())


def load_carbon_source_spec(carbon_source_id: str) -> CarbonSourceSpec:
    return _CARBON_SOURCES[carbon_source_id]


def load_carbon_source_formulation(carbon_source_id: str) -> CarbonSourceFormulation:
    return _CARBON_SOURCE_FORMULATIONS[carbon_source_id]


def load_base_medium_spec(base_medium_id: str) -> BaseMediumSpec:
    return _BASE_MEDIA[base_medium_id]


def load_compatibility_overlay_spec(overlay_id: str) -> CompatibilityOverlaySpec:
    return _COMPATIBILITY_OVERLAYS[overlay_id]


def load_medium_condition_spec(condition_id: str) -> MediumConditionSpec:
    return _MEDIUM_CONDITIONS[condition_id]


def compose_medium_condition_bounds(condition_id: str) -> tuple[MediumBound, ...]:
    condition = load_medium_condition_spec(condition_id)
    base = load_base_medium_spec(condition.base_medium_id)
    carbon = load_carbon_source_spec(condition.carbon_source_id)
    overlay = load_compatibility_overlay_spec(condition.compatibility_overlay_id)
    bounds: dict[str, MediumBound] = {}
    for bound in base.supplement_bounds:
        bounds[bound.reaction_id] = bound
    if condition.open_minimal_inorganic_exchanges:
        for reaction_id in MINIMAL_INORGANIC_EXCHANGES:
            bounds[reaction_id] = _minimal_inorganic_bound(reaction_id)
    for bound in carbon.bounds:
        bounds[bound.reaction_id] = bound
    for reaction_id in carbon.blocked_reactions:
        bounds[reaction_id] = MediumBound(
            reaction_id=reaction_id,
            lower_bound=0.0,
            upper_bound=0.0 if reaction_id.startswith("BIOMASS") else None,
            role="carbon_source",
            source=carbon.carbon_source_id,
            rationale=f"blocked by {carbon.carbon_source_id} carbon-source condition",
        )
    for bound in overlay.bounds:
        bounds[bound.reaction_id] = bound
    return tuple(bounds.values())


def summarize_medium_condition(condition_id: str) -> dict[str, Any]:
    condition = load_medium_condition_spec(condition_id)
    bounds = compose_medium_condition_bounds(condition_id)
    formulation = load_carbon_source_formulation(condition.carbon_source_id)
    return {
        **asdict(condition),
        "scientific_status": medium_condition_scientific_status(condition.condition_id),
        "warnings": list(medium_condition_warnings(condition.condition_id)),
        "carbon_source_formulation": asdict(formulation),
        "bound_count": len(bounds),
        "active_uptake_reactions": [
            bound.reaction_id for bound in bounds if bound.lower_bound is not None and bound.lower_bound < 0
        ],
        "closed_uptake_reactions": [
            bound.reaction_id for bound in bounds if bound.lower_bound == 0.0
        ],
        "bounds": [asdict(bound) for bound in bounds],
    }


def medium_condition_scientific_status(condition_id: str) -> str:
    condition = load_medium_condition_spec(condition_id)
    if condition.compatibility_overlay_id == "matlab_legacy_cost":
        return "matlab_legacy_artifact_compatibility"
    by_carbon = {
        "glucose": "default_python_corrected_reference",
        "methanol": "draft_methanol_induction_boundary_requires_calibration",
        "glycerol": "draft_growth_carbon_boundary_not_induction_proxy",
        "glucose_glycerol": "draft_co_carbon_boundary_requires_promoter_context",
        "glycerol_methanol": "draft_mixed_feed_boundary_requires_calibration",
    }
    return by_carbon.get(condition.carbon_source_id, "draft_medium_boundary")


def medium_condition_warnings(condition_id: str) -> tuple[str, ...]:
    condition = load_medium_condition_spec(condition_id)
    warnings = []
    if condition.carbon_source_id == "glucose":
        if condition.compatibility_overlay_id == "matlab_legacy_cost":
            warnings.append(
                "MATLAB legacy medium overlay 仅用于历史 artifact 对齐，不应作为 corrected 默认培养基。"
            )
        return tuple(warnings)
    warnings.append(
        "当前非葡萄糖培养基是 Python draft 边界条件；尚未声明与旧 MATLAB baseline fully aligned。"
    )
    if condition.carbon_source_id == "glycerol":
        warnings.append("甘油条件适合作为生长碳源边界，不应直接解释为 AOX1 甲醇诱导生产条件。")
    elif condition.carbon_source_id == "methanol":
        warnings.append("甲醇条件可作为诱导相关边界，但仍需要目标和工艺特异的校准。")
    elif condition.carbon_source_id == "glucose_glycerol":
        warnings.append(
            "葡萄糖+甘油条件只打开共碳源摄取；当前模型未加入摄取比例、碳源偏好或启动子阻遏约束。"
        )
    elif condition.carbon_source_id == "glycerol_methanol":
        warnings.append("甘油+甲醇可作为 mixed-feed 边界，但仍需要摄取比例和诱导阶段校准。")
    return tuple(warnings)


def diff_medium_conditions(left_condition_id: str, right_condition_id: str) -> dict[str, Any]:
    left = {bound.reaction_id: bound for bound in compose_medium_condition_bounds(left_condition_id)}
    right = {bound.reaction_id: bound for bound in compose_medium_condition_bounds(right_condition_id)}
    reaction_ids = sorted(set(left) | set(right))
    rows = []
    for reaction_id in reaction_ids:
        left_bound = left.get(reaction_id)
        right_bound = right.get(reaction_id)
        if left_bound == right_bound:
            continue
        rows.append(
            {
                "reaction_id": reaction_id,
                "left_lower_bound": left_bound.lower_bound if left_bound else None,
                "right_lower_bound": right_bound.lower_bound if right_bound else None,
                "left_upper_bound": left_bound.upper_bound if left_bound else None,
                "right_upper_bound": right_bound.upper_bound if right_bound else None,
                "left_source": left_bound.source if left_bound else None,
                "right_source": right_bound.source if right_bound else None,
            }
        )
    return {
        "left_condition_id": left_condition_id,
        "right_condition_id": right_condition_id,
        "difference_count": len(rows),
        "differences": rows,
    }


def _minimal_inorganic_bound(reaction_id: str) -> MediumBound:
    role_by_reaction = {
        "Ex_nh4": "nitrogen",
        "Ex_o2": "gas",
        "Ex_pi": "phosphate",
        "Ex_so4": "sulfate",
        "Ex_fe2": "ion",
        "Ex_h": "water_proton_balance",
        "Ex_h2o": "water_proton_balance",
        "Ex_na1": "ion",
        "Ex_k": "ion",
        "Ex_co2": "gas",
    }
    return MediumBound(
        reaction_id=reaction_id,
        lower_bound=-1000.0,
        upper_bound=None,
        role=role_by_reaction.get(reaction_id, "ion"),  # type: ignore[arg-type]
        source="python_corrected_minimal_inorganic",
        rationale="current Python corrected medium opens minimal inorganic exchanges for glucose-reference FBA",
    )


def _bound(
    reaction_id: str,
    lower: float | None,
    upper: float | None,
    role: BoundRole,
    source: str,
    rationale: str,
    confidence: str = "curated",
) -> MediumBound:
    return MediumBound(
        reaction_id=reaction_id,
        lower_bound=lower,
        upper_bound=upper,
        role=role,
        source=source,
        rationale=rationale,
        confidence=confidence,
    )


def _vitamin_bounds(source: str) -> tuple[MediumBound, ...]:
    return tuple(_bound(reaction_id, -2.0, None, "vitamin", source, "YNB vitamin uptake bound") for reaction_id in YNB_VITAMINS)


def _amino_acid_bounds(reaction_ids: tuple[str, ...], source: str) -> tuple[MediumBound, ...]:
    return tuple(
        _bound(reaction_id, -0.08, None, "amino_acid", source, "legacy pcSecPichia amino-acid supplement uptake bound")
        for reaction_id in reaction_ids
    )


def _closed_amino_acid_bounds(reaction_ids: tuple[str, ...], source: str) -> tuple[MediumBound, ...]:
    return tuple(
        _bound(reaction_id, 0.0, None, "amino_acid", source, "amino-acid uptake not included in this supplement layer")
        for reaction_id in reaction_ids
    )


_CARBON_SOURCES: dict[str, CarbonSourceSpec] = {
    "glucose": CarbonSourceSpec(
        carbon_source_id="glucose",
        label="Glucose carbon source",
        status="active",
        bounds=(
            _bound("Ex_glc_D", -1000.0, None, "carbon_source", "current_glucose_reference", "glucose uptake enabled"),
            _bound("Ex_o2", -1000.0, None, "gas", "current_glucose_reference", "aerobic glucose reference condition"),
            _bound("BIOMASS", None, 1000.0, "biomass", "current_glucose_reference", "biomass upper bound opened"),
            _bound("Ex_glyc", 0.0, None, "carbon_source", "current_glucose_reference", "glycerol uptake blocked under glucose reference"),
            _bound("Ex_meoh", 0.0, None, "carbon_source", "current_glucose_reference", "methanol uptake blocked under glucose reference"),
        ),
        blocked_reactions=("BIOMASS_glyc", "BIOMASS_meoh"),
        notes=("Matches the current Python glucose-reference path.",),
    ),
    "glycerol": CarbonSourceSpec(
        carbon_source_id="glycerol",
        label="Glycerol carbon source",
        status="active",
        bounds=(
            _bound("Ex_glyc", -1000.0, None, "carbon_source", "python_carbon_source_boundary", "glycerol uptake enabled", "draft"),
            _bound("Ex_glc_D", 0.0, None, "carbon_source", "python_carbon_source_boundary", "glucose blocked under glycerol condition", "draft"),
            _bound("Ex_meoh", 0.0, None, "carbon_source", "python_carbon_source_boundary", "methanol blocked under glycerol condition", "draft"),
            _bound("Ex_o2", -1000.0, None, "gas", "python_carbon_source_boundary", "aerobic glycerol condition", "draft"),
            _bound("BIOMASS_glyc", 0.0, 1000.0, "biomass", "python_carbon_source_boundary", "glycerol biomass reaction opened", "draft"),
        ),
        blocked_reactions=("BIOMASS", "BIOMASS_meoh"),
        notes=("Boundary configuration is available; biological calibration remains condition-specific.",),
    ),
    "methanol": CarbonSourceSpec(
        carbon_source_id="methanol",
        label="Methanol carbon source",
        status="active",
        bounds=(
            _bound("Ex_meoh", -1000.0, None, "carbon_source", "python_carbon_source_boundary", "methanol uptake enabled", "draft"),
            _bound("Ex_glc_D", 0.0, None, "carbon_source", "python_carbon_source_boundary", "glucose blocked under methanol condition", "draft"),
            _bound("Ex_glyc", 0.0, None, "carbon_source", "python_carbon_source_boundary", "glycerol blocked under methanol condition", "draft"),
            _bound("Ex_o2", -1000.0, None, "gas", "python_carbon_source_boundary", "aerobic methanol condition", "draft"),
            _bound("BIOMASS_meoh", 0.0, 1000.0, "biomass", "python_carbon_source_boundary", "methanol biomass reaction opened", "draft"),
        ),
        blocked_reactions=("BIOMASS", "BIOMASS_glyc"),
        notes=("Boundary configuration is available; induction-specific validation remains condition-specific.",),
    ),
    "glucose_glycerol": CarbonSourceSpec(
        carbon_source_id="glucose_glycerol",
        label="Glucose + glycerol mixed carbon source",
        status="active",
        bounds=(
            _bound("Ex_glc_D", -1000.0, None, "carbon_source", "python_carbon_source_boundary", "glucose uptake enabled in mixed condition", "draft"),
            _bound("Ex_glyc", -1000.0, None, "carbon_source", "python_carbon_source_boundary", "glycerol uptake enabled in mixed condition", "draft"),
            _bound("Ex_meoh", 0.0, None, "carbon_source", "python_carbon_source_boundary", "methanol blocked under glucose/glycerol condition", "draft"),
            _bound("Ex_o2", -1000.0, None, "gas", "python_carbon_source_boundary", "aerobic mixed carbon condition", "draft"),
            _bound("BIOMASS", 0.0, 1000.0, "biomass", "python_carbon_source_boundary", "glucose biomass reaction opened in mixed condition", "draft"),
            _bound("BIOMASS_glyc", 0.0, 1000.0, "biomass", "python_carbon_source_boundary", "glycerol biomass reaction opened in mixed condition", "draft"),
        ),
        blocked_reactions=("BIOMASS_meoh",),
        notes=("Boundary configuration is available; mixed-feed biology still needs dedicated validation.",),
    ),
    "glycerol_methanol": CarbonSourceSpec(
        carbon_source_id="glycerol_methanol",
        label="Glycerol + methanol mixed carbon source",
        status="active",
        bounds=(
            _bound("Ex_glyc", -1000.0, None, "carbon_source", "python_carbon_source_boundary", "glycerol uptake enabled in mixed condition", "draft"),
            _bound("Ex_meoh", -1000.0, None, "carbon_source", "python_carbon_source_boundary", "methanol uptake enabled in mixed condition", "draft"),
            _bound("Ex_glc_D", 0.0, None, "carbon_source", "python_carbon_source_boundary", "glucose blocked under glycerol/methanol condition", "draft"),
            _bound("Ex_o2", -1000.0, None, "gas", "python_carbon_source_boundary", "aerobic mixed carbon condition", "draft"),
            _bound("BIOMASS_glyc", 0.0, 1000.0, "biomass", "python_carbon_source_boundary", "glycerol biomass reaction opened in mixed condition", "draft"),
            _bound("BIOMASS_meoh", 0.0, 1000.0, "biomass", "python_carbon_source_boundary", "methanol biomass reaction opened in mixed condition", "draft"),
        ),
        blocked_reactions=("BIOMASS",),
        notes=("Boundary configuration is available; induction-transition validation remains condition-specific.",),
    ),
}

_CARBON_SOURCE_FORMULATIONS: dict[str, CarbonSourceFormulation] = {
    "glucose": CarbonSourceFormulation(
        carbon_source_id="glucose",
        active_uptake_reaction_ids=("Ex_glc_D",),
        blocked_uptake_reaction_ids=("Ex_glyc", "Ex_meoh"),
        candidate_growth_reaction_ids=("BIOMASS",),
        selected_growth_reaction_id="BIOMASS",
        carbon_objective_weights={"Ex_glc_D": 1.0},
        formulation_status="corrected_reference",
        matlab_alignment_note="Legacy MATLAB writeLPGlc is primarily compatible with this glucose/BIOMASS formulation.",
    ),
    "glycerol": CarbonSourceFormulation(
        carbon_source_id="glycerol",
        active_uptake_reaction_ids=("Ex_glyc",),
        blocked_uptake_reaction_ids=("Ex_glc_D", "Ex_meoh"),
        candidate_growth_reaction_ids=("BIOMASS_glyc",),
        selected_growth_reaction_id="BIOMASS_glyc",
        carbon_objective_weights={"Ex_glyc": 1.0},
        formulation_status="draft_carbon_source_boundary",
        warnings=(
            "Glycerol formulation is a Python carbon-source boundary; it is not an AOX1 methanol-induction proxy.",
        ),
        matlab_alignment_note="Legacy MATLAB writeLPGlc must be patched or replaced because it fixes BIOMASS internally.",
    ),
    "methanol": CarbonSourceFormulation(
        carbon_source_id="methanol",
        active_uptake_reaction_ids=("Ex_meoh",),
        blocked_uptake_reaction_ids=("Ex_glc_D", "Ex_glyc"),
        candidate_growth_reaction_ids=("BIOMASS_meoh",),
        selected_growth_reaction_id="BIOMASS_meoh",
        carbon_objective_weights={"Ex_meoh": 1.0},
        formulation_status="draft_induction_boundary_requires_calibration",
        warnings=(
            "Methanol formulation opens the model boundary but does not by itself calibrate induction, uptake, or promoter regulation.",
        ),
        matlab_alignment_note="Current legacy MATLAB artifact is not a reliable same-condition methanol baseline without a carbon-source-aware writer.",
    ),
    "glucose_glycerol": CarbonSourceFormulation(
        carbon_source_id="glucose_glycerol",
        active_uptake_reaction_ids=("Ex_glc_D", "Ex_glyc"),
        blocked_uptake_reaction_ids=("Ex_meoh",),
        candidate_growth_reaction_ids=("BIOMASS", "BIOMASS_glyc"),
        selected_growth_reaction_id="BIOMASS",
        carbon_objective_weights={"Ex_glc_D": 1.0, "Ex_glyc": 1.0},
        formulation_status="draft_mixed_carbon_boundary",
        warnings=(
            "Mixed-carbon objective weights do not impose feed ratios, uptake preference, or promoter repression.",
            "Fixed-growth probes select BIOMASS and close BIOMASS_glyc unless a future calibrated mixed-growth mode is added.",
        ),
        matlab_alignment_note="Use a patched weighted-carbon LP objective for comparison; raw writeLPGlc remains a glucose-style artifact.",
    ),
    "glycerol_methanol": CarbonSourceFormulation(
        carbon_source_id="glycerol_methanol",
        active_uptake_reaction_ids=("Ex_glyc", "Ex_meoh"),
        blocked_uptake_reaction_ids=("Ex_glc_D",),
        candidate_growth_reaction_ids=("BIOMASS_glyc", "BIOMASS_meoh"),
        selected_growth_reaction_id="BIOMASS_glyc",
        carbon_objective_weights={"Ex_glyc": 1.0, "Ex_meoh": 1.0},
        formulation_status="draft_mixed_carbon_boundary_requires_calibration",
        warnings=(
            "Glycerol+methanol is a boundary probe, not a calibrated transition-feed or induction-stage prediction.",
            "Fixed-growth probes select BIOMASS_glyc and close BIOMASS_meoh unless a future calibrated mixed-growth mode is added.",
        ),
        matlab_alignment_note="Legacy MATLAB writeLPGlc is not a reliable same-condition baseline for this mixed non-glucose formulation.",
    ),
}

_BASE_MEDIA: dict[str, BaseMediumSpec] = {
    "ynb_minimal": BaseMediumSpec(
        base_medium_id="ynb_minimal",
        label="YNB minimal supplements",
        status="active",
        legacy_media_type=2,
        supplement_bounds=_vitamin_bounds("legacy_media_type_2"),
        notes=("Corresponds to legacy media_type=2 supplement layer.",),
    ),
    "ynb_core_aa": BaseMediumSpec(
        base_medium_id="ynb_core_aa",
        label="YNB + core amino-acid supplements",
        status="active",
        legacy_media_type=4,
        supplement_bounds=(
            *_vitamin_bounds("legacy_media_type_4"),
            *_amino_acid_bounds(CORE_AMINO_ACIDS, "legacy_media_type_4"),
            *_closed_amino_acid_bounds(
                tuple(reaction_id for reaction_id in ALL_AMINO_ACIDS if reaction_id not in CORE_AMINO_ACIDS),
                "legacy_media_type_4",
            ),
        ),
        notes=("Corresponds to the current default legacy media_type=4 supplement layer.",),
    ),
    "ynb_all_aa": BaseMediumSpec(
        base_medium_id="ynb_all_aa",
        label="YNB + all amino-acid supplements",
        status="active",
        legacy_media_type=5,
        supplement_bounds=(*_vitamin_bounds("legacy_media_type_5"), *_amino_acid_bounds(ALL_AMINO_ACIDS, "legacy_media_type_5")),
        notes=("Corresponds to legacy media_type=5 supplement layer.",),
    ),
}

_COMPATIBILITY_OVERLAYS: dict[str, CompatibilityOverlaySpec] = {
    "none": CompatibilityOverlaySpec(
        overlay_id="none",
        label="No compatibility overlay",
        status="active",
        scope="default_simulation",
        bounds=(),
    ),
    "matlab_legacy_cost": CompatibilityOverlaySpec(
        overlay_id="matlab_legacy_cost",
        label="MATLAB legacy protein-cost medium overlay",
        status="active",
        scope="protein_cost_slope_only",
        bounds=tuple(
            _bound(
                reaction_id,
                0.0,
                None,
                "legacy_compatibility",
                "MATLAB SimulateProteinCost legacy artifact",
                "matches MATLAB protein-cost artifact exchange lower bound",
                "legacy",
            )
            for reaction_id in MATLAB_LEGACY_COST_MEDIUM_EXCHANGES
        ),
        notes=("Do not apply to corrected default secretion simulation.",),
    ),
}


def _medium_condition(
    carbon_source_id: str,
    base_medium_id: str,
    *,
    notes: tuple[str, ...],
) -> MediumConditionSpec:
    carbon = load_carbon_source_spec(carbon_source_id)
    base = load_base_medium_spec(base_medium_id)
    return MediumConditionSpec(
        condition_id=f"{carbon_source_id}_{base_medium_id}_corrected",
        label=f"{carbon.label} + {base.label}, Python corrected",
        status="active",
        carbon_source_id=carbon_source_id,
        base_medium_id=base_medium_id,
        legacy_media_type=base.legacy_media_type,
        notes=notes,
    )


_MEDIUM_CONDITIONS: dict[str, MediumConditionSpec] = {
    "glucose_ynb_core_aa_corrected": MediumConditionSpec(
        condition_id="glucose_ynb_core_aa_corrected",
        label="Glucose + YNB core amino acids, Python corrected",
        status="active",
        carbon_source_id="glucose",
        base_medium_id="ynb_core_aa",
        compatibility_overlay_id="none",
        legacy_media_type=4,
        notes=("Current formal Python default condition.",),
    ),
    "glucose_ynb_core_aa_matlab_legacy_cost": MediumConditionSpec(
        condition_id="glucose_ynb_core_aa_matlab_legacy_cost",
        label="Glucose + YNB core amino acids, MATLAB protein-cost legacy overlay",
        status="active",
        carbon_source_id="glucose",
        base_medium_id="ynb_core_aa",
        compatibility_overlay_id="matlab_legacy_cost",
        legacy_media_type=4,
        notes=("Use only for MATLAB SimulateProteinCost artifact comparison.",),
    ),
    "glycerol_ynb_core_aa_corrected": MediumConditionSpec(
        condition_id="glycerol_ynb_core_aa_corrected",
        label="Glycerol + YNB core amino acids, Python corrected",
        status="active",
        carbon_source_id="glycerol",
        base_medium_id="ynb_core_aa",
        legacy_media_type=4,
        notes=("Boundary support is active; full biological calibration is not yet claimed.",),
    ),
    "methanol_ynb_core_aa_corrected": MediumConditionSpec(
        condition_id="methanol_ynb_core_aa_corrected",
        label="Methanol + YNB core amino acids, Python corrected",
        status="active",
        carbon_source_id="methanol",
        base_medium_id="ynb_core_aa",
        legacy_media_type=4,
        notes=("Boundary support is active; methanol induction biology still needs dedicated validation.",),
    ),
    "glycerol_methanol_ynb_core_aa_corrected": MediumConditionSpec(
        condition_id="glycerol_methanol_ynb_core_aa_corrected",
        label="Glycerol + methanol + YNB core amino acids, Python corrected",
        status="active",
        carbon_source_id="glycerol_methanol",
        base_medium_id="ynb_core_aa",
        legacy_media_type=4,
        notes=("Boundary support is active; mixed-feed/transition biology still needs dedicated validation.",),
    ),
    "glucose_glycerol_ynb_core_aa_corrected": MediumConditionSpec(
        condition_id="glucose_glycerol_ynb_core_aa_corrected",
        label="Glucose + glycerol + YNB core amino acids, Python corrected",
        status="active",
        carbon_source_id="glucose_glycerol",
        base_medium_id="ynb_core_aa",
        legacy_media_type=4,
        notes=("Boundary support is active; mixed-feed biology still needs dedicated validation.",),
    ),
}

_MEDIUM_CONDITIONS.update(
    {
        f"{carbon_source_id}_{base_medium_id}_corrected": _medium_condition(
            carbon_source_id,
            base_medium_id,
            notes=(
                "Boundary support is active; full biological calibration is not yet claimed.",
            ),
        )
        for carbon_source_id in _CARBON_SOURCES
        for base_medium_id in ("ynb_minimal", "ynb_all_aa")
    }
)


__all__ = [
    "ALL_AMINO_ACIDS",
    "CORE_AMINO_ACIDS",
    "MATLAB_LEGACY_COST_MEDIUM_EXCHANGES",
    "MINIMAL_INORGANIC_EXCHANGES",
    "YNB_VITAMINS",
    "BaseMediumSpec",
    "CarbonSourceFormulation",
    "CarbonSourceSpec",
    "CompatibilityOverlaySpec",
    "MediumBound",
    "MediumConditionSpec",
    "compose_medium_condition_bounds",
    "diff_medium_conditions",
    "list_base_medium_specs",
    "list_carbon_source_formulations",
    "list_carbon_source_specs",
    "list_compatibility_overlay_specs",
    "list_medium_condition_specs",
    "load_base_medium_spec",
    "load_carbon_source_formulation",
    "load_carbon_source_spec",
    "load_compatibility_overlay_spec",
    "load_medium_condition_spec",
    "medium_condition_scientific_status",
    "medium_condition_warnings",
    "summarize_medium_condition",
]
