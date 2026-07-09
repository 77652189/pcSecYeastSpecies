from __future__ import annotations

import ast
import json
from dataclasses import asdict
from pathlib import Path

from pcsec_pichia.homology.cache_schema import ExternalNameReference, NameAuditRow, RuleTransferAuditRow
from pcsec_pichia.homology.crosswalk import write_name_audit_cache, write_rule_transfer_audit_cache
from pcsec_pichia.external_refs import (
    ExternalGeneFunctionEvidence,
    ExternalGprCandidateEvidence,
    ExternalReferenceProvenance,
    write_external_reference_cache_bundle,
)

from app.services.pichia_homology_audit_service import (
    export_homology_audit_rows,
    homology_audit_cache_status,
    load_homology_audit_browser_data,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_missing_cache_returns_empty_state(tmp_path: Path) -> None:
    payload = load_homology_audit_browser_data(cache_root=tmp_path / "missing")

    assert payload["cache_status"]["cache_available"] is False
    assert payload["cache_status"]["row_count"] == 0
    assert payload["name_audit_rows"] == []
    assert payload["rule_transfer_audit_rows"] == []
    assert "build_pichia_homology_cache.py" in payload["cache_status"]["recommended_build_command"]
    assert payload["cache_status"]["external_cache_available"] is False
    assert "build_pichia_external_name_reference_cache.py" in payload["cache_status"][
        "recommended_external_build_command"
    ]


def test_service_loads_and_filters_name_and_rule_transfer_audits(tmp_path: Path) -> None:
    run_dir = _write_cache_run(tmp_path / "run1")

    payload = load_homology_audit_browser_data(
        cache_root=run_dir,
        query="kar",
        name_consistency_status="name_confirmed_by_rbh",
        rule_transfer_status="rule_transfer_ready",
        is_rbh=True,
        in_model_gene_index=True,
        min_identity=50,
    )

    assert payload["cache_status"]["cache_available"] is True
    assert payload["cache_status"]["row_count"] == 2
    assert len(payload["name_audit_rows"]) == 1
    assert len(payload["rule_transfer_audit_rows"]) == 1
    assert payload["name_audit_rows"][0]["name_consistency_status"] == "name_confirmed_by_rbh"
    assert payload["name_audit_rows"][0]["external_crosscheck_status"] == "external_match_confirmed"
    assert payload["name_audit_rows"][0]["external_crosscheck_sources"] == ["UniProt:2026_01:C4R"]
    assert payload["rule_transfer_audit_rows"][0]["rule_transfer_status"] == "rule_transfer_ready"


def test_service_selects_latest_valid_cache_run(tmp_path: Path) -> None:
    base = tmp_path / "cache"
    old_run = _write_cache_run(base / "old")
    new_run = _write_cache_run(base / "new")
    old_time = new_run.stat().st_mtime - 100
    for path in old_run.iterdir():
        path.touch()
    import os

    os.utime(old_run, (old_time, old_time))

    status = homology_audit_cache_status(cache_root=base)

    assert status["cache_available"] is True
    assert status["cache_root"] == str(new_run)


def test_service_reports_external_reference_cache_status(tmp_path: Path) -> None:
    run_dir = _write_cache_run(tmp_path / "run1", with_external=True, with_external_reference=True)

    status = homology_audit_cache_status(cache_root=run_dir)
    payload = load_homology_audit_browser_data(cache_root=run_dir)

    assert status["external_cache_available"] is True
    assert status["external_cache_path"] == str(run_dir / "external_name_references.jsonl")
    assert status["external_reference_count"] == 2
    assert status["external_sources"] == ["NCBI", "UniProt"]
    assert status["external_source_counts"] == {"NCBI": 1, "UniProt": 1}
    assert "build_pichia_external_name_reference_cache.py" in status["recommended_external_build_command"]
    assert status["external_reference_cache"]["cache_available"] is True
    assert status["external_reference_cache"]["record_type_counts"] == {"gene_function": 1, "gpr_candidate": 1}
    assert len(payload["external_reference_rows"]) == 2
    assert {row["evidence_kind"] for row in payload["external_reference_rows"]} == {"gene_function", "gpr_candidate"}


def test_export_returns_utf8_with_header() -> None:
    rows = [
        {
            "internal_common_name": "KAR2 / BiP",
            "query_symbol": "KAR2",
            "name_consistency_status": "name_confirmed_by_rbh",
        }
    ]

    exported = export_homology_audit_rows(rows, file_format="csv")

    text = exported.decode("utf-8")
    assert text.splitlines()[0] == "internal_common_name,query_symbol,name_consistency_status"
    assert "KAR2 / BiP" in text


def test_homology_audit_service_does_not_import_blast_runtime() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_homology_audit_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    imports: set[str] = set()
    called_names: set[str] = set()
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module)
        elif isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                called_names.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                called_names.add(node.func.attr)

    assert "pcsec_pichia.homology.blast_runner" not in imports
    assert "pcsec_pichia.homology.external_fetch" not in imports
    assert not {"run_blastp", "make_blast_db"} & called_names
    assert not {"fetch_external_name_references", "fetch_uniprot_name_reference"} & called_names


def _write_cache_run(path: Path, *, with_external: bool = False, with_external_reference: bool = False) -> Path:
    path.mkdir(parents=True)
    name_rows = (
        NameAuditRow(
            internal_gene_id="PAS_chr2-1_0140",
            internal_common_name="KAR2 / BiP",
            internal_sequence_id="YJL034W",
            external_accession="C4R",
            external_gene_name="KAR2",
            external_locus_tag="PAS_chr2-1_0140",
            external_aliases=("KAR2",),
            identity_pct=75.0,
            query_coverage=95.0,
            subject_coverage=95.0,
            evalue=1e-100,
            is_rbh=True,
            in_model_gene_index=True,
            name_consistency_status="name_confirmed_by_rbh",
            review_status="model_ready_rbh_high_confidence",
            external_crosscheck_status="external_match_confirmed",
            external_crosscheck_sources=("UniProt:2026_01:C4R",),
            external_crosscheck_warnings=(),
            warnings=(),
        ),
        NameAuditRow(
            internal_gene_id="",
            internal_common_name="PDI1",
            internal_sequence_id="YCL043C",
            external_accession="C4P",
            external_gene_name="PDI1",
            external_locus_tag="PAS_chr4_0844",
            external_aliases=("PDI1",),
            identity_pct=45.0,
            query_coverage=90.0,
            subject_coverage=90.0,
            evalue=1e-80,
            is_rbh=True,
            in_model_gene_index=False,
            name_consistency_status="name_confirmed_by_rbh",
            review_status="rbh_not_in_model",
            warnings=("homology evidence only",),
        ),
    )
    rule_rows = (
        RuleTransferAuditRow(
            internal_common_name="KAR2 / BiP",
            query_symbol="KAR2",
            sce_orf="YJL034W",
            pichia_gene_id="PAS_chr2-1_0140",
            pichia_model_gene_id="PAS_chr2-1_0140",
            is_rbh=True,
            in_model_gene_index=True,
            identity_pct=75.0,
            query_coverage=95.0,
            subject_coverage=95.0,
            evalue=1e-100,
            homology_review_status="model_ready_rbh_high_confidence",
            rule_transfer_status="rule_transfer_ready",
            warnings=(),
        ),
        RuleTransferAuditRow(
            internal_common_name="PDI1",
            query_symbol="PDI1",
            sce_orf="YCL043C",
            pichia_gene_id="PAS_chr4_0844",
            pichia_model_gene_id="",
            is_rbh=True,
            in_model_gene_index=False,
            identity_pct=45.0,
            query_coverage=90.0,
            subject_coverage=90.0,
            evalue=1e-80,
            homology_review_status="rbh_not_in_model",
            rule_transfer_status="rule_transfer_supported_not_model_operable",
            warnings=("homology evidence only",),
        ),
    )
    write_name_audit_cache(name_rows, path / "sce_to_pichia_name_audit.jsonl", path / "sce_to_pichia_name_audit.tsv")
    write_rule_transfer_audit_cache(
        rule_rows,
        path / "sce_to_pichia_rule_transfer_audit.jsonl",
        path / "sce_to_pichia_rule_transfer_audit.tsv",
    )
    (path / "homology_audit_summary.json").write_text(
        json.dumps(
            {
                "blast_status": "completed",
                "homology_row_count": 2,
                "name_audit_row_count": 2,
                "rule_transfer_row_count": 2,
            }
        ),
        encoding="utf-8",
    )
    if with_external:
        _write_external_reference_cache(path)
    if with_external_reference:
        _write_unified_external_reference_cache(path)
    return path


def _write_external_reference_cache(path: Path) -> None:
    references = (
        ExternalNameReference(
            source_database="UniProt",
            source_version="2026_01",
            taxon="Komagataella phaffii",
            accession="C4R",
            gene_name="KAR2",
            locus_tag="PAS_chr2-1_0140",
            aliases=("BiP",),
            retrieved_at="2026-07-08T00:00:00+00:00",
        ),
        ExternalNameReference(
            source_database="NCBI",
            source_version="gene",
            taxon="Komagataella phaffii",
            accession="12345",
            gene_name="PDI1",
            locus_tag="PAS_chr4_0844",
            aliases=(),
            retrieved_at="2026-07-08T00:00:00+00:00",
        ),
    )
    rows = []
    for reference in references:
        payload = asdict(reference)
        payload["aliases"] = list(reference.aliases)
        payload["warnings"] = list(reference.warnings)
        rows.append(json.dumps(payload, sort_keys=True))
    (path / "external_name_references.jsonl").write_text("\n".join(rows) + "\n", encoding="utf-8")


def _write_unified_external_reference_cache(path: Path) -> None:
    records = (
        ExternalGeneFunctionEvidence(
            provenance=_external_ref_provenance("PAS_chr2-1_0140", source_database="uniprot"),
            gene_id="PAS_chr2-1_0140",
            protein_name="Kar2 protein",
            function_description="Annotation-only folding function",
            go_terms=("GO:0006457",),
            evidence_scope="reviewed_structured_annotation",
        ),
        ExternalGprCandidateEvidence(
            provenance=_external_ref_provenance("YJL034W", source_database="yeast-gem"),
            external_model_id="yeast-GEM",
            external_reaction_id="r_1234",
            external_gene_rule="YJL034W",
            candidate_status="gene_mapping_required",
            pichia_gene_id="PAS_chr2-1_0140",
            query_gene_id="YJL034W",
            gpr_transfer_status="gene_mapping_required",
            blocking_reasons=("external gene rule is not mapped to a current Pichia model gene",),
        ),
    )
    write_external_reference_cache_bundle(records, path)


def _external_ref_provenance(query: str, *, source_database: str) -> ExternalReferenceProvenance:
    return ExternalReferenceProvenance(
        source_database=source_database,
        source_version="test",
        source_url=f"https://example.test/{query}",
        source_query=query,
        retrieved_at="2026-07-09T00:00:00Z",
        raw_record_sha256="f" * 64,
    )
