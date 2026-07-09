from __future__ import annotations

from pcsec_pichia.external_refs import (
    ExternalFetchConfig,
    ExternalGeneFunctionEvidence,
    ExternalReactionAssociation,
    ExternalReferenceProvenance,
    build_external_gpr_candidates,
    classify_gpr_transfer_status,
    fetch_external_model_reaction_associations,
    load_external_reference_cache,
    parse_sbml_gpr_associations,
    write_external_reference_cache_bundle,
)


def test_parse_sbml_gpr_associations_extracts_reaction_gene_rule(tmp_path) -> None:
    sbml_path = tmp_path / "toy.xml"
    sbml_path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2">
  <model id="toy_model">
    <fbc:listOfGeneProducts>
      <fbc:geneProduct fbc:id="G_YBR160W" fbc:label="YBR160W"/>
      <fbc:geneProduct fbc:id="G_YDR123C" fbc:label="YDR123C"/>
    </fbc:listOfGeneProducts>
    <listOfReactions>
      <reaction id="r_1234" name="example reaction">
        <fbc:geneProductAssociation>
          <fbc:and>
            <fbc:geneProductRef fbc:geneProduct="G_YBR160W"/>
            <fbc:geneProductRef fbc:geneProduct="G_YDR123C"/>
          </fbc:and>
        </fbc:geneProductAssociation>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
""",
        encoding="utf-8",
    )

    associations = parse_sbml_gpr_associations(
        sbml_path,
        source_database="yeast-gem",
        source_model_id="toy_model",
    )

    assert len(associations) == 1
    assert associations[0].external_reaction_id == "r_1234"
    assert associations[0].external_reaction_name == "example reaction"
    assert associations[0].external_gene_ids == ("YBR160W", "YDR123C")
    assert associations[0].gene_rule == "(YBR160W and YDR123C)"
    assert associations[0].association_status == "external_gpr_candidate"


def test_fetch_external_model_reaction_associations_reads_offline_cache_without_network(tmp_path) -> None:
    write_external_reference_cache_bundle(
        (
            _association(
                source_database="yeast-gem",
                external_reaction_id="r_1234",
                external_gene_ids=("YBR160W",),
                gene_rule="YBR160W",
            ),
            _association(
                source_database="bigg",
                external_model_id="other-model",
                external_reaction_id="R_OTHER",
                external_gene_ids=("OTHER",),
                gene_rule="OTHER",
            ),
        ),
        tmp_path,
    )

    associations = fetch_external_model_reaction_associations(
        source_database="yeast-gem",
        model_id="yeast-GEM",
        gene_or_reaction_query="YBR160W",
        config=ExternalFetchConfig(offline_cache_dir=str(tmp_path)),
    )

    assert len(associations) == 1
    assert associations[0].external_reaction_id == "r_1234"
    assert associations[0].gene_rule == "YBR160W"


def test_build_external_gpr_candidates_confirms_only_mapped_current_model_gene_and_reaction() -> None:
    candidates = build_external_gpr_candidates(
        pichia_gene_id="PAS_chr1-1_0001",
        query_gene_id="YBR160W",
        gene_function_evidence=(_gene_function("PAS_chr1-1_0001"),),
        reaction_associations=(
            _association(
                external_reaction_id="r_1234",
                external_gene_ids=("YBR160W",),
                gene_rule="YBR160W",
                mapped_pichia_gene_ids=("PAS_chr1-1_0001",),
            ),
        ),
        current_model_reaction_ids=("R_PIC_1234",),
        reaction_crosswalk={"r_1234": "R_PIC_1234"},
    )

    assert len(candidates) == 1
    candidate = candidates[0]
    assert candidate.pichia_gene_id == "PAS_chr1-1_0001"
    assert candidate.query_gene_id == "YBR160W"
    assert candidate.candidate_status == "model_gpr_confirmed"
    assert candidate.gpr_transfer_status == "model_gpr_confirmed"
    assert candidate.mapped_pichia_reaction_id == "R_PIC_1234"
    assert candidate.mapped_pichia_gene_ids == ("PAS_chr1-1_0001",)
    assert candidate.supporting_gene_evidence
    assert not candidate.blocking_reasons


def test_external_gpr_candidate_mapping_fields_roundtrip_through_cache(tmp_path) -> None:
    candidate = build_external_gpr_candidates(
        pichia_gene_id="PAS_chr1-1_0001",
        query_gene_id="YBR160W",
        gene_function_evidence=(_gene_function("PAS_chr1-1_0001"),),
        reaction_associations=(
            _association(
                external_reaction_id="r_1234",
                external_gene_ids=("YBR160W",),
                gene_rule="YBR160W",
                mapped_pichia_gene_ids=("PAS_chr1-1_0001",),
            ),
        ),
        current_model_reaction_ids=("R_PIC_1234",),
        reaction_crosswalk={"r_1234": "R_PIC_1234"},
    )[0]

    write_external_reference_cache_bundle((candidate,), tmp_path)
    loaded = load_external_reference_cache(tmp_path / "external_reference_records.jsonl")

    assert loaded == (candidate,)
    assert loaded[0].record_type == "gpr_candidate"
    assert loaded[0].candidate_status == "model_gpr_confirmed"


def test_build_external_gpr_candidates_requires_reaction_and_gene_mapping() -> None:
    reaction_missing = build_external_gpr_candidates(
        pichia_gene_id="PAS_chr1-1_0001",
        query_gene_id="YBR160W",
        gene_function_evidence=(),
        reaction_associations=(
            _association(
                external_reaction_id="r_missing",
                external_gene_ids=("YBR160W",),
                gene_rule="YBR160W",
            ),
        ),
        current_model_reaction_ids=("R_PIC_1234",),
        reaction_crosswalk={},
    )[0]
    gene_missing = build_external_gpr_candidates(
        pichia_gene_id="PAS_chr1-1_0001",
        query_gene_id="YBR160W",
        gene_function_evidence=(),
        reaction_associations=(
            _association(
                external_reaction_id="r_1234",
                external_gene_ids=("YBR160W",),
                gene_rule="YBR160W",
            ),
        ),
        current_model_reaction_ids=("R_PIC_1234",),
        reaction_crosswalk={"r_1234": "R_PIC_1234"},
    )[0]

    assert reaction_missing.candidate_status == "reaction_mapping_required"
    assert "external reaction is not mapped to a current Pichia model reaction" in reaction_missing.blocking_reasons
    assert gene_missing.candidate_status == "gene_mapping_required"
    assert "external gene rule is not mapped to a current Pichia model gene" in gene_missing.blocking_reasons


def test_build_external_gpr_candidates_flags_conflicting_source_rules() -> None:
    candidates = build_external_gpr_candidates(
        pichia_gene_id="PAS_chr1-1_0001",
        query_gene_id="YBR160W",
        gene_function_evidence=(),
        reaction_associations=(
            _association(
                external_reaction_id="r_1234",
                external_gene_ids=("YBR160W",),
                gene_rule="YBR160W",
                mapped_pichia_gene_ids=("PAS_chr1-1_0001",),
            ),
            _association(
                external_reaction_id="r_1234",
                external_gene_ids=("YDR123C",),
                gene_rule="YDR123C",
                mapped_pichia_gene_ids=("PAS_chr1-1_0001",),
            ),
        ),
        current_model_reaction_ids=("R_PIC_1234",),
        reaction_crosswalk={"r_1234": "R_PIC_1234"},
    )

    assert {candidate.candidate_status for candidate in candidates} == {"conflicting_gpr_sources"}
    assert all("conflicting external GPR rules" in candidate.blocking_reasons for candidate in candidates)


def test_classify_gpr_transfer_status_keeps_source_rules_candidate_only_until_mapped() -> None:
    status, reasons = classify_gpr_transfer_status(
        gene_mapping_status="external_gene_rule_only",
        reaction_mapping_status="external_reaction_only",
        source_gene_rule="YBR160W",
        mapped_model_reaction_id=None,
        in_current_model_gene_index=False,
    )

    assert status == "reaction_mapping_required"
    assert reasons == ("external reaction is not mapped to a current Pichia model reaction",)
    assert status != "experiment_calibrated"


def _provenance(query: str = "YBR160W", *, source_database: str = "yeast-gem") -> ExternalReferenceProvenance:
    return ExternalReferenceProvenance(
        source_database=source_database,
        source_version="test",
        source_url=f"https://example.test/{query}",
        source_query=query,
        retrieved_at="2026-07-09T00:00:00Z",
        raw_record_sha256="d" * 64,
    )


def _association(
    *,
    source_database: str = "yeast-gem",
    external_model_id: str = "yeast-GEM",
    external_reaction_id: str,
    external_gene_ids: tuple[str, ...],
    gene_rule: str,
    mapped_pichia_gene_ids: tuple[str, ...] = (),
) -> ExternalReactionAssociation:
    return ExternalReactionAssociation(
        provenance=_provenance(external_reaction_id, source_database=source_database),
        external_model_id=external_model_id,
        external_reaction_id=external_reaction_id,
        external_reaction_name="Example reaction",
        external_gene_ids=external_gene_ids,
        gene_rule=gene_rule,
        mapped_pichia_gene_ids=mapped_pichia_gene_ids,
        association_status="external_gpr_candidate",
    )


def _gene_function(gene_id: str) -> ExternalGeneFunctionEvidence:
    return ExternalGeneFunctionEvidence(
        provenance=_provenance(gene_id),
        gene_id=gene_id,
        protein_name="Mapped protein",
        go_terms=("GO:0006886",),
        evidence_scope="reviewed_structured_annotation",
    )
