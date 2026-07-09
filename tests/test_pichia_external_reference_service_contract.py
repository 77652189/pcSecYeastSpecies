from __future__ import annotations

import ast
from pathlib import Path

from app.services.pichia_external_reference_service import (
    export_external_reference_rows,
    load_external_reference_browser_rows,
    load_external_reference_status,
    submit_external_reference_refresh,
)
from pcsec_pichia.external_refs import (
    EXTERNAL_GPR_MAPPING_CANDIDATES_FILENAME,
    ExternalGeneFunctionEvidence,
    ExternalGprCandidateEvidence,
    ExternalModelInventoryRecord,
    ExternalReactionAssociation,
    ExternalReferenceProvenance,
    ExternalReferenceRecord,
    GprSourcePriorityRecord,
    write_external_reference_cache_bundle,
    write_external_model_inventory,
    write_gpr_source_priority_outputs,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_missing_external_reference_cache_returns_refresh_command(tmp_path: Path) -> None:
    status = load_external_reference_status(tmp_path / "missing")

    assert status["cache_available"] is False
    assert status["record_count"] == 0
    records_path = Path(status["records_path"])
    assert records_path.name == "external_reference_records.jsonl"
    assert records_path.parent.name == "missing"
    assert "external_reference_records.jsonl" in status["missing_files"]
    assert "build_pichia_external_reference_cache.py" in status["recommended_refresh_command"]


def test_service_loads_external_reference_status_and_browser_rows(tmp_path: Path) -> None:
    cache_dir = _write_external_cache(tmp_path / "external")

    status = load_external_reference_status(cache_dir)
    rows = load_external_reference_browser_rows(cache_dir)

    assert status["cache_available"] is True
    assert status["record_count"] == 4
    assert status["source_counts"] == {"uniprot": 2, "yeast-gem": 2}
    assert status["record_type_counts"] == {
        "external_reference": 1,
        "gene_function": 1,
        "gpr_candidate": 1,
        "reaction_association": 1,
    }
    assert status["retrieved_at_range"]["first"] == "2026-07-09T00:00:00Z"
    assert status["retrieved_at_range"]["last"] == "2026-07-09T00:00:00Z"
    assert {row["evidence_kind"] for row in rows} == {
        "external_reference",
        "gene_function",
        "gpr_candidate",
        "reaction_association",
    }
    function_row = next(row for row in rows if row["evidence_kind"] == "gene_function")
    assert function_row["function_description"] == "Annotation-only secretion function"
    assert function_row["manual_review_reasons"] == []
    gpr_row = next(row for row in rows if row["evidence_kind"] == "gpr_candidate")
    assert gpr_row["gpr_transfer_status"] == "gene_mapping_required"
    assert gpr_row["manual_review_reasons"] == [
        "external gene rule is not mapped to a current Pichia model gene"
    ]


def test_service_surfaces_external_model_gpr_summary_from_mapping_cache(tmp_path: Path) -> None:
    cache_dir = _write_external_model_gpr_cache(tmp_path / "external_model_gpr")

    status = load_external_reference_status(cache_dir)
    rows = load_external_reference_browser_rows(cache_dir)

    assert status["cache_available"] is True
    assert status["external_model_sources"] == ["Kp.1.0"]
    assert status["external_gpr_candidate_count"] == 1
    assert status["best_external_gpr_source"] == "biomodels:Kp.1.0"
    assert status["gpr_source_priority"]["best_priority_tier"] == "pichia_literature_model_gpr"
    assert status["external_gpr_mapping_status"] == {"gene_mapping_required": 1}
    assert status["external_model_inventory"]["record_count"] == 1
    assert "conflicting external GPR rules" in status["external_gpr_conflict_warnings"]

    gpr_row = rows[0]
    assert gpr_row["evidence_kind"] == "gpr_candidate"
    assert gpr_row["external_model_sources"] == ["Kp.1.0"]
    assert gpr_row["gpr_source_priority"] == "pichia_literature_model_gpr"
    assert gpr_row["external_gpr_candidate_count"] == 1
    assert gpr_row["external_gpr_mapping_status"] == "gene_mapping_required"
    assert "conflicting external GPR rules" in gpr_row["external_gpr_conflict_warnings"]


def test_service_filters_external_rows_by_query_and_kind(tmp_path: Path) -> None:
    cache_dir = _write_external_cache(tmp_path / "external")

    rows = load_external_reference_browser_rows(
        cache_dir,
        query="YBR160W",
        evidence_kind="gpr_candidate",
        manual_review_only=True,
    )

    assert len(rows) == 1
    assert rows[0]["query_gene_id"] == "YBR160W"
    assert rows[0]["evidence_kind"] == "gpr_candidate"


def test_submit_external_reference_refresh_returns_manual_command_without_network(tmp_path: Path) -> None:
    result = submit_external_reference_refresh(
        homology_run_dir=tmp_path / "homology_run",
        sources=("uniprot", "sgd"),
        limit=10,
    )

    assert result["submitted"] is False
    assert result["network_performed"] is False
    assert "--sources uniprot,sgd" in result["command"]
    assert "--limit 10" in result["command"]
    assert "local_runs\\pichia_external_reference_cache" in result["command"]


def test_export_external_reference_rows_returns_utf8_table() -> None:
    rows = [
        {
            "evidence_kind": "gpr_candidate",
            "source_database": "yeast-gem",
            "query_gene_id": "YBR160W",
            "manual_review_reasons": ["external gene rule is not mapped"],
        }
    ]

    exported = export_external_reference_rows(rows, file_format="csv").decode("utf-8")

    assert exported.splitlines()[0].startswith("evidence_kind,source_database")
    assert "YBR160W" in exported
    assert "external gene rule is not mapped" in exported


def test_external_reference_service_does_not_import_ui_or_network_clients() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_external_reference_service.py"
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

    assert not any(module.startswith("app.ui") for module in imports)
    assert not {"fetch_external_references", "default_http_get", "request_json"} & called_names


def _write_external_cache(path: Path) -> Path:
    records = (
        ExternalReferenceRecord(
            provenance=_provenance("C4R"),
            taxon_id="4922",
            organism="Komagataella phaffii",
            primary_accession="C4R",
            gene_id="PAS_chr1-1_0001",
            gene_name="SEC1",
            locus_tag="PAS_chr1-1_0001",
            aliases=("SEC1",),
            protein_name="Secretion protein Sec1",
        ),
        ExternalGeneFunctionEvidence(
            provenance=_provenance("PAS_chr1-1_0001"),
            gene_id="PAS_chr1-1_0001",
            protein_name="Secretion protein Sec1",
            function_description="Annotation-only secretion function",
            go_terms=("GO:0006886",),
            evidence_scope="reviewed_structured_annotation",
        ),
        ExternalReactionAssociation(
            provenance=_provenance("r_1234", source_database="yeast-gem"),
            external_model_id="yeast-GEM",
            external_reaction_id="r_1234",
            external_gene_ids=("YBR160W",),
            gene_rule="YBR160W",
        ),
        ExternalGprCandidateEvidence(
            provenance=_provenance("YBR160W", source_database="yeast-gem"),
            external_model_id="yeast-GEM",
            external_reaction_id="r_1234",
            external_gene_rule="YBR160W",
            candidate_status="gene_mapping_required",
            pichia_gene_id="PAS_chr1-1_0001",
            query_gene_id="YBR160W",
            mapped_pichia_reaction_id="R_PIC_1234",
            gene_mapping_status="external_gene_rule_only",
            reaction_mapping_status="model_reaction_mapped",
            gpr_transfer_status="gene_mapping_required",
            confidence="manual_review_required",
            blocking_reasons=("external gene rule is not mapped to a current Pichia model gene",),
        ),
    )
    write_external_reference_cache_bundle(records, path)
    return path


def _write_external_model_gpr_cache(path: Path) -> Path:
    candidate = ExternalGprCandidateEvidence(
        provenance=_provenance("PAS_chr1-1_0001", source_database="biomodels"),
        external_model_id="Kp.1.0",
        external_reaction_id="R_KP_SEC",
        external_gene_rule="KP_GENE",
        candidate_status="gene_mapping_required",
        pichia_gene_id="PAS_chr1-1_0001",
        query_gene_id="KP_GENE",
        mapped_pichia_reaction_id="R_PIC_SEC",
        gene_mapping_status="external_gene_rule_only",
        reaction_mapping_status="model_reaction_mapped",
        gpr_transfer_status="gene_mapping_required",
        confidence="manual_review_required",
        blocking_reasons=("external gene rule is not mapped to a current Pichia model gene",),
        mapping_warnings=("conflicting external GPR rules",),
    )
    write_external_reference_cache_bundle(
        (candidate,),
        path,
        records_filename=EXTERNAL_GPR_MAPPING_CANDIDATES_FILENAME,
    )
    write_external_model_inventory(
        (
            ExternalModelInventoryRecord(
                model_id="Kp.1.0",
                model_name="Kp.1.0 genome-scale model",
                organism="Komagataella phaffii",
                source_database_or_repository="BioModels",
                source_url="https://example.test/kp10",
                publication_url="https://doi.org/10.1002/bit.26380",
                license="CC BY 4.0",
                available_artifact_types=("SBML",),
                download_status="downloadable",
                local_path="",
                checksum_sha256="",
                has_gpr=True,
                has_gene_ids=True,
                has_reaction_ids=True,
                has_sbml=True,
                notes="fixture",
            ),
        ),
        path,
    )
    write_gpr_source_priority_outputs(
        (
            GprSourcePriorityRecord(
                candidate_cache_key=candidate.cache_key,
                source_database="biomodels",
                external_model_id="Kp.1.0",
                external_reaction_id="R_KP_SEC",
                mapped_pichia_reaction_id="R_PIC_SEC",
                external_gene_rule="KP_GENE",
                priority_rank=2,
                priority_tier="pichia_literature_model_gpr",
                conflict_status="conflicting_gpr_sources",
                manual_review_required=True,
                warnings=("conflicting external GPR rules",),
            ),
        ),
        path,
    )
    return path


def _provenance(query: str, *, source_database: str = "uniprot") -> ExternalReferenceProvenance:
    return ExternalReferenceProvenance(
        source_database=source_database,
        source_version="test",
        source_url=f"https://example.test/{query}",
        source_query=query,
        retrieved_at="2026-07-09T00:00:00Z",
        raw_record_sha256="e" * 64,
    )
