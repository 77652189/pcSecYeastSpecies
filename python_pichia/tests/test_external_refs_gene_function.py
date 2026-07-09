from __future__ import annotations

from pcsec_pichia.external_refs import (
    attach_ko_oe_external_gene_evidence,
    build_gene_function_evidence,
    build_ko_oe_external_gene_evidence,
    classify_gene_function_confidence,
)
from pcsec_pichia.external_refs.schema import ExternalReferenceProvenance, ExternalReferenceRecord


def test_build_gene_function_evidence_extracts_annotation_fields() -> None:
    records = (
        _record(
            "uniprot",
            "P12345",
            gene_id="PAS_chr1-1_0001",
            gene_name="SEC1",
            protein_name="Secretion protein Sec1",
            function_description="Vesicle trafficking annotation",
            ec_numbers=("2.7.1.1",),
            go_terms=("GO:0006886",),
            pathways=("secretory pathway",),
            orthology=("KOG0001",),
            reviewed=True,
        ),
    )

    evidence = build_gene_function_evidence(
        internal_gene_id="PAS_chr1-1_0001",
        external_records=records,
    )

    assert len(evidence) == 1
    assert evidence[0].gene_id == "PAS_chr1-1_0001"
    assert evidence[0].protein_name == "Secretion protein Sec1"
    assert evidence[0].ec_numbers == ("2.7.1.1",)
    assert evidence[0].go_terms == ("GO:0006886",)
    assert evidence[0].pathways == ("secretory pathway",)
    assert evidence[0].orthology == ("KOG0001",)
    assert evidence[0].reviewed is True
    assert evidence[0].evidence_scope == "reviewed_structured_annotation"


def test_gene_function_confidence_never_claims_phenotype_or_experiment_calibration() -> None:
    evidence = build_gene_function_evidence(
        internal_gene_id="PAS_chr1-1_0001",
        external_records=(
            _record(
                "uniprot",
                "P12345",
                gene_id="PAS_chr1-1_0001",
                gene_name="SEC1",
                protein_name="Secretion protein Sec1",
                go_terms=("GO:0006886",),
                reviewed=True,
            ),
        ),
    )[0]

    confidence, warnings = classify_gene_function_confidence(evidence)

    assert confidence == "reviewed_structured_annotation"
    assert warnings == ("external annotation is not phenotype evidence and must not calibrate recommendation_tier",)
    assert "experiment_calibrated" not in confidence


def test_build_ko_oe_external_gene_evidence_summarizes_name_and_function_without_gpr_claim() -> None:
    evidence = build_ko_oe_external_gene_evidence(
        pichia_gene_id="PAS_chr1-1_0001",
        standard_name="SEC1",
        external_records=(
            _record(
                "uniprot",
                "P12345",
                gene_id="PAS_chr1-1_0001",
                gene_name="SEC1",
                protein_name="Secretion protein Sec1",
                reviewed=True,
            ),
        ),
        model_executable_gene_id="PAS_chr1-1_0001",
        model_gpr_executable=True,
    )

    payload = evidence.to_dict()

    assert evidence.external_name_status == "external_match_confirmed"
    assert len(evidence.function_evidence) == 1
    assert evidence.gpr_candidates == ()
    assert evidence.model_gpr_executable is True
    assert payload["function_evidence"][0]["evidence_scope"] == "reviewed_name_annotation"


def test_attach_ko_oe_external_gene_evidence_adds_fields_without_changing_recommendation_tier() -> None:
    rows = (
        {
            "gene_id": "PAS_chr1-1_0001",
            "canonical_gene_id": "PAS_chr1-1_0001",
            "input_gene_id": "SEC1",
            "common_name": "SEC1",
            "intervention_type": "KO",
            "recommendation_tier": "model_executable",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
        },
    )
    merged = attach_ko_oe_external_gene_evidence(
        rows,
        (
            _record(
                "uniprot",
                "P12345",
                gene_id="PAS_chr1-1_0001",
                gene_name="SEC1",
                protein_name="Secretion protein Sec1",
                function_description="Annotation-only secretion function",
                go_terms=("GO:0006886",),
                reviewed=True,
            ),
        ),
    )

    row = merged[0]

    assert row["recommendation_tier"] == "model_executable"
    assert row["external_gene_function_sources"] == ("uniprot",)
    assert row["external_gene_function_confidence"] == ("reviewed_structured_annotation",)
    assert row["external_gene_function_evidence"][0]["function_description"] == (
        "Annotation-only secretion function"
    )
    assert row["ko_oe_external_gene_evidence"]["external_name_status"] == "external_match_confirmed"
    assert "phenotype_evidence" not in row["ko_oe_external_gene_evidence"]


def test_attach_ko_oe_external_gene_evidence_keeps_unmatched_rows_manual_review_only() -> None:
    rows = (
        {
            "gene_id": "PAS_chr1-1_9999",
            "input_gene_id": "NO_MATCH",
            "intervention_type": "OE_gene_proxy",
            "recommendation_tier": "manual_review_required",
        },
    )
    merged = attach_ko_oe_external_gene_evidence(
        rows,
        (
            _record(
                "uniprot",
                "P12345",
                gene_id="PAS_chr1-1_0001",
                gene_name="SEC1",
                protein_name="Secretion protein Sec1",
            ),
        ),
    )

    row = merged[0]

    assert row["recommendation_tier"] == "manual_review_required"
    assert row["external_gene_function_evidence"] == ()
    assert row["ko_oe_external_gene_evidence"]["external_name_status"] == "external_reference_missing"
    assert row["ko_oe_external_gene_evidence"]["manual_review_reasons"]


def test_attach_ko_oe_external_gene_evidence_accepts_name_audit_style_gene_fields() -> None:
    rows = (
        {
            "internal_gene_id": "PAS_chr1-1_0001",
            "internal_common_name": "SEC1",
            "recommendation_tier": "evidence_supported",
        },
    )
    merged = attach_ko_oe_external_gene_evidence(
        rows,
        (
            _record(
                "uniprot",
                "P12345",
                gene_id="PAS_chr1-1_0001",
                gene_name="SEC1",
                protein_name="Secretion protein Sec1",
                go_terms=("GO:0006886",),
                reviewed=True,
            ),
        ),
    )

    row = merged[0]

    assert row["recommendation_tier"] == "evidence_supported"
    assert row["external_gene_function_sources"] == ("uniprot",)
    assert row["ko_oe_external_gene_evidence"]["pichia_gene_id"] == "PAS_chr1-1_0001"
    assert row["ko_oe_external_gene_evidence"]["external_name_status"] == "external_match_confirmed"


def _record(
    source: str,
    accession: str,
    *,
    gene_id: str,
    gene_name: str,
    protein_name: str | None = None,
    function_description: str | None = None,
    ec_numbers: tuple[str, ...] = (),
    go_terms: tuple[str, ...] = (),
    pathways: tuple[str, ...] = (),
    orthology: tuple[str, ...] = (),
    reviewed: bool | None = None,
) -> ExternalReferenceRecord:
    return ExternalReferenceRecord(
        provenance=ExternalReferenceProvenance(
            source_database=source,
            source_version="test",
            source_url=f"https://example.test/{accession}",
            source_query=accession,
            retrieved_at="2026-07-09T00:00:00Z",
            raw_record_sha256="c" * 64,
        ),
        taxon_id="4922",
        organism="Komagataella phaffii",
        primary_accession=accession,
        gene_id=gene_id,
        gene_name=gene_name,
        locus_tag=gene_id,
        aliases=(gene_name,),
        protein_name=protein_name,
        function_description=function_description,
        ec_numbers=ec_numbers,
        go_terms=go_terms,
        pathways=pathways,
        orthology=orthology,
        reviewed=reviewed,
    )
