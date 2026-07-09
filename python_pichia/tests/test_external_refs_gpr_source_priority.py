from __future__ import annotations

import json

from pcsec_pichia.external_refs import (
    CONFLICTING_GPR_SOURCES,
    ExternalGprCandidateEvidence,
    ExternalReferenceProvenance,
    rank_external_gpr_sources,
    write_gpr_source_priority_outputs,
)


def test_rank_external_gpr_sources_prefers_current_model_then_pichia_specific_models() -> None:
    ranked = rank_external_gpr_sources(
        (
            _candidate("gpruler", "GPRuler", "R1", "AUTO1"),
            _candidate("yeast-gem", "Yeast9", "R1", "YBR160W"),
            _candidate("publication", "iPichia", "R1", "PAS_gene"),
            _candidate("pichia-current", "current Pichia GEM", "R1", "PAS_gene", status="model_gpr_confirmed"),
            _candidate("biomodels", "Kp.1.0", "R1", "KP_gene"),
            _candidate("uniprot", "UniProt", "R1", "annotation"),
        )
    )

    assert [row.priority_tier for row in ranked] == [
        "current_model_gpr",
        "pichia_specific_model_gpr",
        "pichia_literature_model_gpr",
        "homology_supported_yeast_gpr",
        "annotation_only",
        "automatic_rule_candidate",
    ]
    assert ranked[0].priority_rank == 0
    assert ranked[-1].manual_review_required is True


def test_rank_external_gpr_sources_flags_conflicting_rules_without_merging() -> None:
    ranked = rank_external_gpr_sources(
        (
            _candidate("publication", "iPichia", "R_PIC_1", "PAS_gene_A", mapped_reaction="R_PIC_1"),
            _candidate("biomodels", "Kp.1.0", "R_KP_1", "PAS_gene_B", mapped_reaction="R_PIC_1"),
        )
    )

    assert {row.conflict_status for row in ranked} == {CONFLICTING_GPR_SOURCES}
    assert all(row.manual_review_required for row in ranked)
    assert all("conflicting external GPR rules" in row.warnings for row in ranked)


def test_write_gpr_source_priority_outputs(tmp_path) -> None:
    outputs = write_gpr_source_priority_outputs(
        rank_external_gpr_sources((_candidate("publication", "iPichia", "R1", "PAS_gene"),)),
        tmp_path,
    )

    priority_payload = json.loads(outputs.priority_path.read_text(encoding="utf-8"))
    assert priority_payload["record_count"] == 1
    assert priority_payload["records"][0]["priority_tier"] == "pichia_specific_model_gpr"
    assert outputs.conflicts_path.exists()
    assert "GPR Source Priority Report" in outputs.report_path.read_text(encoding="utf-8")


def _candidate(
    source_database: str,
    external_model_id: str,
    reaction_id: str,
    gene_rule: str,
    *,
    status: str = "external_gpr_candidate",
    mapped_reaction: str | None = None,
) -> ExternalGprCandidateEvidence:
    return ExternalGprCandidateEvidence(
        provenance=ExternalReferenceProvenance(
            source_database=source_database,
            source_version=external_model_id,
            source_url=f"https://example.test/{external_model_id}",
            source_query=external_model_id,
            retrieved_at="2026-07-09T00:00:00Z",
            raw_record_sha256="e" * 64,
        ),
        external_model_id=external_model_id,
        external_reaction_id=reaction_id,
        external_gene_rule=gene_rule,
        candidate_status=status,
        mapped_pichia_reaction_id=mapped_reaction,
        gpr_transfer_status=status,
    )
