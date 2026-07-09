from __future__ import annotations

from pcsec_pichia.external_refs import (
    ExternalGprCandidateEvidence,
    ExternalReferenceProvenance,
    load_external_reference_cache,
    map_external_gpr_candidates_to_model,
    write_external_gpr_mapping_outputs,
)


def test_map_external_gpr_candidates_confirms_only_reaction_and_gene_mapped_rules() -> None:
    mapped = map_external_gpr_candidates_to_model(
        (_candidate("r_ext", "YBR160W"),),
        current_model_reaction_ids=("R_PIC_1",),
        current_model_gene_ids=("PAS_chr1-1_0001",),
        reaction_crosswalk={"r_ext": "R_PIC_1"},
        gene_crosswalk={"YBR160W": "PAS_chr1-1_0001"},
    )

    assert len(mapped) == 1
    candidate = mapped[0]
    assert candidate.candidate_status == "model_gpr_confirmed"
    assert candidate.gpr_transfer_status == "model_gpr_confirmed"
    assert candidate.mapped_pichia_reaction_id == "R_PIC_1"
    assert candidate.mapped_pichia_gene_ids == ("PAS_chr1-1_0001",)
    assert candidate.reaction_mapping_status == "model_reaction_mapped"
    assert candidate.gene_mapping_status == "model_gene_mapped"
    assert candidate.confidence == "mapped_external_gpr"


def test_map_external_gpr_candidates_keeps_unmapped_reaction_candidate_only() -> None:
    mapped = map_external_gpr_candidates_to_model(
        (_candidate("r_ext", "YBR160W"),),
        current_model_reaction_ids=("R_OTHER",),
        current_model_gene_ids=("PAS_chr1-1_0001",),
        reaction_crosswalk={},
        gene_crosswalk={"YBR160W": "PAS_chr1-1_0001"},
    )

    assert mapped[0].candidate_status == "reaction_mapping_required"
    assert mapped[0].mapped_pichia_reaction_id is None
    assert "external reaction is not mapped" in mapped[0].blocking_reasons[0]


def test_map_external_gpr_candidates_keeps_unmapped_gene_candidate_only() -> None:
    mapped = map_external_gpr_candidates_to_model(
        (_candidate("r_ext", "YBR160W"),),
        current_model_reaction_ids=("R_PIC_1",),
        current_model_gene_ids=("PAS_chr1-1_0001",),
        reaction_crosswalk={"r_ext": "R_PIC_1"},
        gene_crosswalk={},
    )

    assert mapped[0].candidate_status == "gene_mapping_required"
    assert mapped[0].mapped_pichia_reaction_id == "R_PIC_1"
    assert mapped[0].mapped_pichia_gene_ids == ()
    assert "external gene rule is not mapped" in mapped[0].blocking_reasons[0]


def test_map_external_gpr_candidates_marks_crosswalk_outside_current_model() -> None:
    mapped = map_external_gpr_candidates_to_model(
        (_candidate("r_ext", "YBR160W"),),
        current_model_reaction_ids=("R_PIC_1",),
        current_model_gene_ids=("PAS_chr1-1_0001",),
        reaction_crosswalk={"r_ext": "R_NOT_IN_MODEL"},
        gene_crosswalk={"YBR160W": "PAS_chr1-1_0001"},
    )

    assert mapped[0].candidate_status == "not_in_current_model"
    assert "mapped reaction is not present" in mapped[0].blocking_reasons[0]


def test_map_external_gpr_candidates_preserves_source_rule_missing_status() -> None:
    mapped = map_external_gpr_candidates_to_model(
        (_candidate("r_ext", None),),
        current_model_reaction_ids=("R_PIC_1",),
        current_model_gene_ids=("PAS_chr1-1_0001",),
        reaction_crosswalk={"r_ext": "R_PIC_1"},
        gene_crosswalk={"YBR160W": "PAS_chr1-1_0001"},
    )

    assert mapped[0].candidate_status == "source_rule_missing"
    assert mapped[0].confidence == "manual_review_required"


def test_map_external_gpr_candidates_flags_conflicting_mapped_rules() -> None:
    mapped = map_external_gpr_candidates_to_model(
        (_candidate("r_ext", "YBR160W"), _candidate("r_ext", "YDR123C")),
        current_model_reaction_ids=("R_PIC_1",),
        current_model_gene_ids=("PAS_chr1-1_0001", "PAS_chr1-1_0002"),
        reaction_crosswalk={"r_ext": "R_PIC_1"},
        gene_crosswalk={"YBR160W": "PAS_chr1-1_0001", "YDR123C": "PAS_chr1-1_0002"},
    )

    assert {candidate.candidate_status for candidate in mapped} == {"conflicting_gpr_sources"}
    assert all(candidate.confidence == "manual_review_required" for candidate in mapped)


def test_map_external_gpr_candidates_does_not_upgrade_manual_review_candidate() -> None:
    candidate = _candidate("r_ext", "YBR160W", candidate_status="manual_review_required")

    mapped = map_external_gpr_candidates_to_model(
        (candidate,),
        current_model_reaction_ids=("R_PIC_1",),
        current_model_gene_ids=("PAS_chr1-1_0001",),
        reaction_crosswalk={"r_ext": "R_PIC_1"},
        gene_crosswalk={"YBR160W": "PAS_chr1-1_0001"},
    )

    assert mapped[0].candidate_status == "manual_review_required"
    assert mapped[0].confidence == "manual_review_required"


def test_write_external_gpr_mapping_outputs_roundtrips_candidate_cache(tmp_path) -> None:
    mapped = map_external_gpr_candidates_to_model(
        (_candidate("r_ext", "YBR160W"),),
        current_model_reaction_ids=("R_PIC_1",),
        current_model_gene_ids=("PAS_chr1-1_0001",),
        reaction_crosswalk={"r_ext": "R_PIC_1"},
        gene_crosswalk={"YBR160W": "PAS_chr1-1_0001"},
    )

    outputs = write_external_gpr_mapping_outputs(mapped, tmp_path)

    assert outputs.candidate_count == 1
    assert outputs.candidates_path.name == "external_gpr_candidate_evidence.jsonl"
    loaded = load_external_reference_cache(outputs.candidates_path)
    assert loaded == mapped
    report = outputs.report_path.read_text(encoding="utf-8")
    assert "External GPR Mapping Report" in report
    assert "model_gpr_confirmed" in report


def _candidate(
    reaction_id: str,
    gene_rule: str | None,
    *,
    candidate_status: str | None = None,
) -> ExternalGprCandidateEvidence:
    status = candidate_status or ("source_rule_missing" if gene_rule is None else "external_gpr_candidate")
    return ExternalGprCandidateEvidence(
        provenance=ExternalReferenceProvenance(
            source_database="biomodels",
            source_version="Kp.1.0",
            source_url="https://example.test/kp",
            source_query="Kp.1.0",
            retrieved_at="2026-07-09T00:00:00Z",
            raw_record_sha256="f" * 64,
        ),
        external_model_id="Kp.1.0",
        external_reaction_id=reaction_id,
        external_gene_rule=gene_rule,
        candidate_status=status,
        gpr_transfer_status=status,
    )
