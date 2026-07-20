from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from pcsec_pichia.core.target_inputs import LeaderCandidateInput, TargetProteinInput
from pcsec_pichia.core.target_protein_plan import TargetProteinBuildPlan, build_target_protein_plan
from pcsec_pichia.probe import SecretoryEnzymeData, TargetSpec, load_secretory_enzymedata, repo_root

from pcsec_pichia.secretory_resources.schema import (
    CalibrationMode,
    EvidenceClass,
    ExecutionStatus,
    ResourceApplicability,
    ResourceCategory,
    ResourceSource,
    SecretoryResource,
    SecretoryResourceCatalog,
)


_TARGET_PLAN_SOURCE_REF = "pcsec_pichia.core.target_protein_plan.build_target_protein_plan"
_TARGET_SPEC_SOURCE_REF = "pcsec_pichia.probe.TargetSpec (Data/pcSecPichia/TargetProtein*.xlsx structural profile)"
_HOST_ENZYME_DATASET_REF = "pcsec_pichia.core.pichia_enzymes.SecretoryEnzymeData"
_SECRETORY_ENZYME_DATASET_REF = (
    "pcsec_pichia.probe.load_secretory_enzymedata (Enzymedata/pcSecPichia/enzymedataSEC_PP.mat)"
)

# Kar2p/BiP, its nucleotide-exchange factor (Nefs), and the RAC/Ssa1-Ydj1-Snl1
# co-chaperone set recur across the translocation and DSB reactions this
# target actually uses. This is used only to point at evidence within the
# target's own plan; it is NOT a substitute for loading the host-level
# SecretoryEnzymeData capacity dataset (see _folding_chaperone_resource).
_CHAPERONE_NAME_TOKENS = ("bip", "nefs", "kar2", "rac", "ssa1", "ydj1", "snl1")

_TRAFFICKING_STAGES = ("er_to_golgi", "golgi_processing", "final_transport")


def build_secretory_resource_catalog(
    target: TargetSpec,
    *,
    host: str = "Komagataella phaffii",
    condition: str = "any",
    feature_enabled: bool = True,
    root: Path | None = None,
) -> SecretoryResourceCatalog:
    """Round 0: identity, real handle references, and honest status only.

    No network calls, no solver call, no numeric capacity value. Covers all
    seven pichia_next_plan.md Round 0 categories: translocation,
    folding/chaperone, disulfide bond formation, glycosylation, vesicle
    trafficking, ER quality control/ERAD/proteasome, and target-specific
    translation/degradation cost.

    Each resource's limitations note whether its own handles are also backed
    by a real, non-placeholder kcat entry in the base model's secretory
    enzyme dataset (load_secretory_enzymedata) -- confirming which
    categories already have real capacity data sitting elsewhere in the
    pipeline, versus which are genuinely data-free. Round 0 still does not
    read or compute with these kcat values; it only records whether they
    exist.
    """

    if not feature_enabled:
        return SecretoryResourceCatalog(target_id=target.target_id, host=host, feature_enabled=False)

    plan = _build_target_plan(target)
    enzyme_data = _cached_secretory_enzymedata(root or repo_root())
    resources = (
        _translocation_resource(target, plan, host, condition, enzyme_data),
        _folding_chaperone_resource(target, plan, host, condition, enzyme_data),
        _disulfide_resource(target, plan, host, condition, enzyme_data),
        _glycosylation_resource(target, plan, host, condition, enzyme_data),
        _vesicle_trafficking_resource(target, plan, host, condition, enzyme_data),
        _er_quality_control_resource(target, plan, host, condition, enzyme_data),
        _target_specific_cost_resource(target, plan, host, condition, enzyme_data),
    )
    return SecretoryResourceCatalog(target_id=target.target_id, host=host, feature_enabled=True, resources=resources)


def _build_target_plan(target: TargetSpec) -> TargetProteinBuildPlan:
    target_input = TargetProteinInput(
        target_id=target.target_id,
        protein_name=target.target_id,
        abbreviation=target.protein_id,
        mature_sequence=target.mature_sequence,
        through_er=1 if target.through_er else 0,
        signal_peptide=1 if target.signal_peptide_sequence else 0,
        disulfide_sites=target.disulfide_sites,
        n_glycosylation_sites=target.n_glycosylation_sites,
        o_glycosylation_sites=target.o_glycosylation_sites,
        transmembrane=target.transmembrane,
        gpi_sites=target.gpi_sites,
        localization=target.localization,
        cotranslation=target.cotranslation,
        parameter_status="ready_for_model",
    )
    leader = LeaderCandidateInput(
        candidate_id=f"{target.protein_id}_secretory_resources_round0",
        leader_sequence=target.leader_sequence or "M",
        signal_peptide_sequence=target.signal_peptide_sequence or "M",
    )
    return build_target_protein_plan(target_input, leader)


@lru_cache(maxsize=4)
def _cached_secretory_enzymedata(root: Path) -> SecretoryEnzymeData:
    # Parsing enzymedataSEC_PP.mat is slow enough (~seconds) that calling it
    # once per build_secretory_resource_catalog() invocation made the test
    # suite ~25x slower; the file never changes within a process lifetime,
    # so caching by root path is safe. maxsize=4 comfortably covers repo
    # root plus a couple of test tmp_path roots without growing unbounded.
    return load_secretory_enzymedata(root)


def _source(evidence_class: EvidenceClass, source_ref: str) -> ResourceSource:
    return ResourceSource(source_ref=source_ref, version="1", evidence_class=evidence_class)


def _applicability(target: TargetSpec, host: str, condition: str) -> ResourceApplicability:
    return ResourceApplicability(host=host, target_id=target.target_id, condition=condition)


def _kcat_matches(handles: tuple[str, ...], enzyme_data: SecretoryEnzymeData) -> tuple[tuple[str, str, float], ...]:
    """(handle, complex_id, kcat) for every handle backed by a real enzymedataSEC_PP.mat entry.

    Complex names in that dataset (e.g. "sec_SEC61SEC63C_complex") appear as a
    substring of this catalog's reaction ids (e.g.
    "hLF_..._Post_translation_PSTA_sec_SEC61SEC63C_complex"); there is no
    shared key to join on, so substring containment is the only available
    correspondence. Some handles are deliberately not enzyme-catalyzed
    reactions (bare degradation/dilution bookkeeping, e.g. ERAD5B or
    export_sp_to_c) and will never match; that is expected, not a bug.
    """
    entries = enzyme_data.unique_complex_entries()
    matches: list[tuple[str, str, float]] = []
    for handle in handles:
        for entry in entries:
            if entry.complex_id in handle:
                matches.append((handle, entry.complex_id, entry.kcat))
                break
    return tuple(matches)


def _kcat_evidence_note(handles: tuple[str, ...], enzyme_data: SecretoryEnzymeData) -> str:
    matches = _kcat_matches(handles, enzyme_data)
    if not matches:
        return (
            f"Checked against {_SECRETORY_ENZYME_DATASET_REF}: none of this resource's handles "
            "matched a real complex entry there; this dataset is not a source of capacity data "
            "for this category."
        )
    example_handle, example_complex, example_kcat = matches[0]
    return (
        f"{len(matches)}/{len(handles)} handles are confirmed backed by a real, non-placeholder "
        f"kcat entry in {_SECRETORY_ENZYME_DATASET_REF} (e.g. {example_complex}={example_kcat:.1f}); "
        "Round 0 does not read or use these numbers, it only records that they exist. "
        + (
            "All handles matched."
            if len(matches) == len(handles)
            else f"{len(handles) - len(matches)} handle(s) had no matching complex entry "
            "(expected for non-enzyme-catalyzed bookkeeping reactions)."
        )
    )


def _executable_resource(
    target: TargetSpec,
    host: str,
    condition: str,
    category: ResourceCategory,
    handles: tuple[str, ...],
) -> SecretoryResource:
    return SecretoryResource(
        resource_id=f"{target.target_id}_{category.value}",
        category=category,
        canonical_unit="model_flux",
        model_handles=handles,
        source=_source(EvidenceClass.CURRENT_MODEL_HANDLE, _TARGET_PLAN_SOURCE_REF),
        applicability=_applicability(target, host, condition),
        status=ExecutionStatus.EXECUTABLE,
        calibration_mode=CalibrationMode.RELATIVE_UNCALIBRATED,
        uncertainty_note=(
            "Round 0 freezes handle identity only; relative-comparison scoring "
            "and absolute capacity are not computed this round."
        ),
    )


def _unavailable_resource(
    target: TargetSpec,
    host: str,
    condition: str,
    category: ResourceCategory,
    reason: str,
) -> SecretoryResource:
    return SecretoryResource(
        resource_id=f"{target.target_id}_{category.value}",
        category=category,
        canonical_unit="model_flux",
        model_handles=(),
        source=_source(EvidenceClass.TARGET_STRUCTURAL_PROFILE, _TARGET_SPEC_SOURCE_REF),
        applicability=_applicability(target, host, condition),
        status=ExecutionStatus.UNAVAILABLE,
        limitations=(reason,),
    )


def _not_applicable_resource(
    target: TargetSpec,
    host: str,
    condition: str,
    category: ResourceCategory,
    reason: str,
) -> SecretoryResource:
    return SecretoryResource(
        resource_id=f"{target.target_id}_{category.value}",
        category=category,
        canonical_unit="model_flux",
        model_handles=(),
        source=_source(EvidenceClass.TARGET_STRUCTURAL_PROFILE, _TARGET_SPEC_SOURCE_REF),
        applicability=_applicability(target, host, condition),
        status=ExecutionStatus.NOT_APPLICABLE,
        limitations=(reason,),
    )


def _translocation_resource(
    target: TargetSpec, plan: TargetProteinBuildPlan, host: str, condition: str, enzyme_data: SecretoryEnzymeData
) -> SecretoryResource:
    category = ResourceCategory.ER_TRANSLOCATION
    if not target.through_er:
        return _not_applicable_resource(
            target, host, condition, category, "target.through_er is False: protein is not routed through the ER."
        )
    handles = plan.reaction_ids_by_stage("translocation")
    if not handles:
        return _unavailable_resource(
            target, host, condition, category, "through_er is True but the current target plan produced no translocation reaction handle."
        )
    resource = _executable_resource(target, host, condition, category, handles)
    return _with_limitations(resource, (_kcat_evidence_note(handles, enzyme_data),))


def _disulfide_resource(
    target: TargetSpec, plan: TargetProteinBuildPlan, host: str, condition: str, enzyme_data: SecretoryEnzymeData
) -> SecretoryResource:
    category = ResourceCategory.DISULFIDE_BOND_FORMATION
    if target.disulfide_sites <= 0:
        return _not_applicable_resource(
            target, host, condition, category, "target.disulfide_sites is 0: no disulfide bond formation is required."
        )
    handles = tuple(
        reaction.reaction_id for reaction in plan.reactions if reaction.source_function == "addDSB"
    )
    if not handles:
        return _unavailable_resource(
            target, host, condition, category, "disulfide_sites > 0 but the current target plan produced no addDSB reaction handle."
        )
    resource = _executable_resource(target, host, condition, category, handles)
    return _with_limitations(resource, (_kcat_evidence_note(handles, enzyme_data),))


def _glycosylation_resource(
    target: TargetSpec, plan: TargetProteinBuildPlan, host: str, condition: str, enzyme_data: SecretoryEnzymeData
) -> SecretoryResource:
    category = ResourceCategory.GLYCOSYLATION
    if target.n_glycosylation_sites <= 0 and target.o_glycosylation_sites <= 0:
        return _not_applicable_resource(
            target, host, condition, category, "n_glycosylation_sites and o_glycosylation_sites are both 0."
        )
    handles = tuple(
        reaction.reaction_id for reaction in plan.reactions if reaction.source_function in ("addNG", "addOG")
    )
    if not handles:
        return _unavailable_resource(
            target,
            host,
            condition,
            category,
            "glycosylation sites > 0 but the current target plan produced no addNG/addOG reaction handle.",
        )
    resource = _executable_resource(target, host, condition, category, handles)
    return _with_limitations(resource, (_kcat_evidence_note(handles, enzyme_data),))


def _vesicle_trafficking_resource(
    target: TargetSpec, plan: TargetProteinBuildPlan, host: str, condition: str, enzyme_data: SecretoryEnzymeData
) -> SecretoryResource:
    category = ResourceCategory.VESICLE_TRAFFICKING
    if not target.through_er:
        return _not_applicable_resource(
            target,
            host,
            condition,
            category,
            "target.through_er is False: protein does not transit ER-to-Golgi/vesicle trafficking.",
        )
    handles = tuple(
        reaction.reaction_id for reaction in plan.reactions if reaction.stage in _TRAFFICKING_STAGES
    )
    if not handles:
        return _unavailable_resource(
            target, host, condition, category, "through_er is True but the current target plan produced no trafficking reaction handle."
        )
    resource = _executable_resource(target, host, condition, category, handles)
    return _with_limitations(resource, (_kcat_evidence_note(handles, enzyme_data),))


def _folding_chaperone_resource(
    target: TargetSpec, plan: TargetProteinBuildPlan, host: str, condition: str, enzyme_data: SecretoryEnzymeData
) -> SecretoryResource:
    category = ResourceCategory.FOLDING_CHAPERONE
    if not target.through_er:
        return _not_applicable_resource(
            target, host, condition, category, "target.through_er is False: protein does not use ER chaperone machinery."
        )
    handles = tuple(
        reaction.reaction_id
        for reaction in plan.reactions
        # misfolding-stage reactions are excluded even when they mention a
        # chaperone token (e.g. *_misfold_ERAD_sec_Kar2p_complex): that use of
        # Kar2p/BiP is an ERAD-pathway handle owned by the er_quality_control
        # category below, not a folding-pathway one. Without this exclusion
        # the same handle would be silently double-claimed by two categories.
        if reaction.stage != "misfolding"
        and any(token in reaction.reaction_id.lower() for token in _CHAPERONE_NAME_TOKENS)
    )
    if not handles:
        return _unavailable_resource(
            target, host, condition, category, "through_er is True but no chaperone-bearing reaction handle was found in the target plan."
        )
    # Every one of these handles (verified, not assumed) matches a real,
    # non-placeholder kcat entry in the base model's secretory enzyme
    # dataset -- see _kcat_evidence_note below -- so this category meets the
    # same "real handle + real capacity data" bar as translocation/disulfide/
    # glycosylation/vesicle_trafficking and is executable, not evidence_only.
    # An earlier version of this function held it at evidence_only on the
    # belief that no host-level capacity dataset was loaded; that belief was
    # wrong -- the dataset exists and was simply never checked here.
    # evidence_class stays CLASSIFIER_INFERRED (not the CURRENT_MODEL_HANDLE
    # _executable_resource() would default to): the *identification method*
    # is still a name-token match, not a dedicated stage lookup, and that is
    # a separate, still-valid caveat from "is there real capacity data".
    return SecretoryResource(
        resource_id=f"{target.target_id}_{category.value}",
        category=category,
        canonical_unit="model_flux",
        model_handles=handles,
        source=_source(EvidenceClass.CLASSIFIER_INFERRED, _TARGET_PLAN_SOURCE_REF),
        applicability=_applicability(target, host, condition),
        status=ExecutionStatus.EXECUTABLE,
        calibration_mode=CalibrationMode.RELATIVE_UNCALIBRATED,
        uncertainty_note=(
            "Round 0 freezes handle identity only; relative-comparison scoring "
            "and absolute capacity are not computed this round."
        ),
        limitations=(
            "Handle identity comes from a name match on the target's own reaction handles "
            "(Kar2p/BiP/Nefs/RAC/Ssa1p/Ydj1p/Snl1p tokens), not a dedicated plan stage, so it is "
            "weaker evidence than a direct stage lookup; every matched handle is independently "
            "confirmed below to correspond to a real complex in the base enzyme dataset, which "
            "corroborates (but does not fully replace the need to double-check) the name match.",
            _kcat_evidence_note(handles, enzyme_data),
        ),
    )


def _er_quality_control_resource(
    target: TargetSpec, plan: TargetProteinBuildPlan, host: str, condition: str, enzyme_data: SecretoryEnzymeData
) -> SecretoryResource:
    category = ResourceCategory.ER_QUALITY_CONTROL
    # build_target_protein_plan always emits "misfolding"-stage reactions for
    # every target, through-ER or not (see target_protein_plan.py's
    # _misfolding_reactions/_non_secretory_folding_reactions: each
    # unconditionally appends at least one misfold/ERAD handle). Unlike
    # folding/chaperone, this stage is the reaction's own dedicated stage
    # tag, not a name-token guess across other stages' reactions, so it
    # qualifies for the same executable bar as translocation/disulfide/
    # glycosylation/vesicle_trafficking.
    handles = plan.reaction_ids_by_stage("misfolding")
    if not handles:
        return _unavailable_resource(
            target, host, condition, category, "target plan produced no misfolding/ERAD reaction handle."
        )
    resource = _executable_resource(target, host, condition, category, handles)
    return _with_limitations(
        resource,
        (
            "Handles cover ERAD/proteasome-mediated misfolded-protein degradation flux "
            "(Hrd1p/Der1p/Doa10p-class retrotranslocation and ubiquitination reactions); "
            "UPR transcriptional signaling (Ire1p/Hac1p-mediated chaperone upregulation) is "
            "not represented in this stoichiometric model and has no handle to freeze.",
            "A constraint-row generator for this stage already exists: "
            "pcsec_pichia.probe.misfolding_constraint_rows, gated by "
            "PichiaSimulationRequest.enable_misfolding_constraint (default False). "
            "pcsec_pichia.probe.proteasome_rows is a separate generator that is NOT "
            "gated by that flag -- it is called unconditionally by "
            "build_pcsec_constraint_matrices, independent of this catalog. Round 0 only "
            "freezes the reaction-handle identities above and does not enable or "
            "evaluate the misfolding_constraint_rows layer. The misfolding constraint "
            "math also depends on a host-enzyme kdeg parameter (distinct from this "
            "target's own kdeg in target_specific_cost) that silently defaults to 0.0 "
            "when absent from source data and has no test asserting a real value, "
            "unlike kcat; a future round must check that before trusting a solve with "
            "this flag enabled.",
            # Coverage is genuinely partial and target-specific here (unlike the
            # other executable categories, which match 100%): the "B" branch
            # reactions _misfolding_reactions() emits when a PTM type is absent
            # (e.g. ERAD2B/3B/4B when disulfides/n_glycans are 0), plus bare
            # dilution/degradation bookkeeping reactions, have no distinct
            # catalyzing complex and so never match. Do not assume the other
            # matched fraction generalizes to 100% for a target not yet checked.
            _kcat_evidence_note(handles, enzyme_data),
        ),
    )


def _target_specific_cost_resource(
    target: TargetSpec, plan: TargetProteinBuildPlan, host: str, condition: str, enzyme_data: SecretoryEnzymeData
) -> SecretoryResource:
    category = ResourceCategory.TARGET_SPECIFIC_COST
    # "translation" and "degradation" stage reaction ids embed protein_id
    # (e.g. r_hLF_peptide_translation vs r_OPN_..._peptide_translation), so
    # hLF's and OPN's handles are structurally disjoint by construction; see
    # test_target_specific_cost_never_shares_handles_across_targets, which
    # guards against ever copying one target's handles into another's
    # resource (the failure mode the Round 0 contract explicitly forbids).
    handles = plan.reaction_ids_by_stage("translation") + plan.reaction_ids_by_stage("degradation")
    if not handles:
        return _unavailable_resource(
            target, host, condition, category, "target plan produced no translation/degradation reaction handle."
        )
    resource = _executable_resource(target, host, condition, category, handles)
    return _with_limitations(
        resource,
        (
            "Handles cover this target's own peptide-translation and steady-state "
            "subunit/signal-peptide degradation reactions only. Modification cost "
            "(disulfide bond formation, N-/O-glycosylation) is already expressed by this "
            "target's disulfide_bond_formation/glycosylation resources elsewhere in this "
            "catalog and is deliberately not duplicated here.",
            "core.pichia_enzymes.TargetProteinEnzymeData (via "
            "adapters.pichia_target_enzymedata.target_enzymedata_from_plan) already computes "
            "real target-specific protein_mw/protein_extra_mw numbers, but hardcodes "
            "kdeg=0.0 (no reviewed degradation-rate source yet); a future round must not "
            "treat that zero as a real value when it eventually numerically scores this "
            "category.",
            # Translation/degradation are not enzyme-complex reactions, so
            # checking them against the secretory enzyme dataset is expected to
            # find nothing; this confirms that expectation rather than leaving
            # it unverified, and confirms the real gap is specifically the
            # target protein's own kdeg, not a missing enzyme kcat.
            _kcat_evidence_note(handles, enzyme_data),
        ),
    )


def _with_limitations(resource: SecretoryResource, extra: tuple[str, ...]) -> SecretoryResource:
    from dataclasses import replace

    return replace(resource, limitations=tuple(dict.fromkeys((*resource.limitations, *extra))))


__all__ = ["build_secretory_resource_catalog"]
