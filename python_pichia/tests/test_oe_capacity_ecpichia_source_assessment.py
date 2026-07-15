from __future__ import annotations

import hashlib
import json
import zipfile
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from pcsec_pichia.external_refs.capacity_sources import (
    ECPICHIA_G6PDH2_TABLE_PROFILE,
    ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE,
    RetrievalMode,
    cache_ecpichia_supplement_table_source,
    cache_ecpichia_supplement_source,
    cache_uniprot_identity_source,
)
from pcsec_pichia.external_refs.clients import ExternalHttpResponse
from pcsec_pichia.oe_capacity import external_candidate_audit as audit_module
from pcsec_pichia.oe_capacity.external_candidate_audit import (
    ExternalCapacityAuditRequest,
    run_external_capacity_candidate_audit,
)
from pcsec_pichia.oe_capacity.external_candidate_io import (
    parse_ecpichia_g6pdh2_table_evidence,
    parse_ecpichia_g6pdh2_source_assessment,
)
from pcsec_pichia.oe_capacity.external_candidate_evaluation import (
    evaluate_ecpichia_g6pdh2_provenance,
)
from pcsec_pichia.oe_capacity.external_candidate_schema import CapacityModelBinding
from pcsec_pichia.oe_capacity.schema import OECapacityValidationError


REPO_ROOT = Path(__file__).resolve().parents[2]


def _fixture_yaml(*, concentration: str = "0.752073171936811") -> str:
    return f'''models:
    - !!omap
      - id: "G6PDH2"
      - name: "glucose 6-phosphate dehydrogenase"
      - metabolites: !!omap
          - prot_C4R099: -0.00200309027777778
      - lower_bound: 0
      - upper_bound: 1000
      - gene_reaction_rule: "PAS_chr2-1_0308"
      - eccodes: "1.1.1.49"
    - !!omap
      - id: "usage_prot_C4R099"
      - metabolites: !!omap
          - prot_C4R099: -1
      - lower_bound: -0.752073171936811
      - upper_bound: 0
      - gene_reaction_rule: "PAS_chr2-1_0308"
    - !!omap
      - id: "prot_pool_exchange"
      - metabolites: !!omap
          - prot_pool: -1
      - lower_bound: -219.25
      - upper_bound: 0
kinetics:
  - !!omap
    - id: "G6PDH2"
    - kcat: 8000
    - source: "brenda"
    - eccodes: "1.1.1.49"
    - enzymes: !!omap
        - C4R099: 1
proteins:
  - !!omap
    - genes: "PAS_chr2-1_0308"
    - enzymes: "C4R099"
    - mw: 57689
    - sequence: "MTEST"
    - concs: {concentration}
'''


def _profile_for(path: Path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return replace(
        ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE,
        artifact_sha256=digest,
    )


def _table_profile_for(path: Path):
    return replace(
        ECPICHIA_G6PDH2_TABLE_PROFILE,
        artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )


def _write_table_docx(
    path: Path,
    *,
    headers: tuple[str, ...] = (
        "Rxns",
        "kcat, 1/sec",
        "source",
        "ec #",
        "genes",
        "enzymes",
        "mw, g/mole",
        "concs, g/L",
    ),
    row: tuple[str, ...] = (
        "G6PDH2",
        "8000",
        "brenda",
        "1.1.1.49",
        "PAS_chr1-3_0138",
        "C4QVC5",
        "46856",
        "NaN",
    ),
) -> None:
    def xml_row(values: tuple[str, ...]) -> str:
        return "<w:tr>" + "".join(
            f"<w:tc><w:p><w:r><w:t>{value}</w:t></w:r></w:p></w:tc>"
            for value in values
        ) + "</w:tr>"

    document = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        f"<w:body><w:tbl>{xml_row(headers)}{xml_row(row)}</w:tbl></w:body>"
        "</w:document>"
    )
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("word/document.xml", document)


def _cache_identity(output_dir: Path) -> None:
    payload = json.dumps(
        {
            "results": [
                {
                    "primaryAccession": "C4R099",
                    "uniProtkbId": "A0A1E4RTV1_PICPA",
                    "entryType": "UniProtKB unreviewed (TrEMBL)",
                    "genes": [
                        {
                            "geneName": {"value": "G6PDH2"},
                            "orderedLocusNames": [{"value": "PAS_chr2-1_0308"}],
                        }
                    ],
                    "organism": {
                        "taxonId": 644223,
                        "scientificName": "Komagataella phaffii",
                    },
                    "proteinDescription": {
                        "recommendedName": {
                            "fullName": {
                                "value": "Glucose-6-phosphate 1-dehydrogenase"
                            }
                        }
                    },
                }
            ]
        }
    )
    cache_uniprot_identity_source(
        "PAS_chr2-1_0308",
        output_dir,
        http_get=lambda url, config: ExternalHttpResponse(
            status_code=200,
            text=payload,
            url=url,
            headers={"X-UniProt-Release": "2026_02"},
        ),
    )


def test_ecpichia_parser_retains_raw_values_conflicts_and_replays(tmp_path: Path) -> None:
    source_path = tmp_path / "Supplementary 8.yml"
    source_path.write_text(_fixture_yaml(), encoding="utf-8")
    profile = _profile_for(source_path)
    cache_dir = tmp_path / "cache"

    source = cache_ecpichia_supplement_source(
        profile,
        cache_dir,
        source_file=source_path,
    )
    evidence = parse_ecpichia_g6pdh2_source_assessment(source)

    assert evidence.gene_id == "PAS_chr2-1_0308"
    assert evidence.enzyme_id == "C4R099"
    assert evidence.kcat_per_s == 8000
    assert evidence.molecular_weight_g_per_mol == 57689
    assert evidence.reported_concentration == pytest.approx(0.752073171936811)
    assert evidence.reaction_protein_coefficient == pytest.approx(
        -(57689 / (8000 * 3600))
    )
    assert (
        "source_unit_missing_and_supplement_header_requires_review"
        in evidence.conflicts
    )
    assert "formation_flux_conversion_missing" in evidence.missing_information
    assert source.terms_reviewed is False

    replayed = cache_ecpichia_supplement_source(
        profile,
        cache_dir,
        retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
    )
    assert replayed.retrieval_mode is RetrievalMode.OFFLINE_REPLAY
    assert parse_ecpichia_g6pdh2_source_assessment(replayed) == replace(
        evidence,
        source=replayed,
    )


def test_ecpichia_parser_rejects_incomplete_or_incoherent_raw_values(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "Supplementary 8.yml"
    source_path.write_text(_fixture_yaml(concentration="0.5"), encoding="utf-8")
    source = cache_ecpichia_supplement_source(
        _profile_for(source_path),
        tmp_path / "cache",
        source_file=source_path,
    )

    evidence = parse_ecpichia_g6pdh2_source_assessment(source)
    assert evidence.reported_concentration == 0.5


def test_ecpichia_table_parser_replays_and_rejects_tampering(tmp_path: Path) -> None:
    source_path = tmp_path / "Supplementary 11.docx"
    _write_table_docx(source_path)
    profile = _table_profile_for(source_path)
    cache_dir = tmp_path / "table-cache"
    source = cache_ecpichia_supplement_table_source(
        profile, cache_dir, source_file=source_path
    )
    evidence = parse_ecpichia_g6pdh2_table_evidence(source)
    assert evidence.gene_id == "PAS_chr1-3_0138"
    assert evidence.enzyme_id == "C4QVC5"
    assert evidence.reported_concentration_text == "NaN"
    replayed = cache_ecpichia_supplement_table_source(
        profile, cache_dir, retrieval_mode=RetrievalMode.OFFLINE_REPLAY
    )
    assert parse_ecpichia_g6pdh2_table_evidence(replayed).reaction_id == "G6PDH2"
    Path(replayed.cache_path).write_bytes(Path(replayed.cache_path).read_bytes() + b"x")
    with pytest.raises(OECapacityValidationError, match="sha256 mismatch"):
        cache_ecpichia_supplement_table_source(
            profile, cache_dir, retrieval_mode=RetrievalMode.OFFLINE_REPLAY
        )


def test_ecpichia_table_parser_rejects_wrong_headers(tmp_path: Path) -> None:
    source_path = tmp_path / "Supplementary 11.docx"
    _write_table_docx(
        source_path,
        headers=("wrong", "kcat", "source", "ec", "genes", "enzymes", "mw", "concs"),
    )
    source = cache_ecpichia_supplement_table_source(
        _table_profile_for(source_path), tmp_path / "cache", source_file=source_path
    )
    with pytest.raises(OECapacityValidationError, match="headers"):
        parse_ecpichia_g6pdh2_table_evidence(source)


def test_ecpichia_table_parser_requires_g6pdh2_row(tmp_path: Path) -> None:
    source_path = tmp_path / "Supplementary 11.docx"
    _write_table_docx(
        source_path,
        row=("OTHER", "8000", "brenda", "1.1.1.49", "gene", "enzyme", "1", "NaN"),
    )
    source = cache_ecpichia_supplement_table_source(
        _table_profile_for(source_path), tmp_path / "cache", source_file=source_path
    )
    with pytest.raises(OECapacityValidationError, match="exactly one complete G6PDH2"):
        parse_ecpichia_g6pdh2_table_evidence(source)


def test_ecpichia_closure_retains_formula_and_binding_conflicts(tmp_path: Path) -> None:
    yaml_path = tmp_path / "Supplementary 8.yml"
    yaml_path.write_text(_fixture_yaml(), encoding="utf-8")
    yaml_source = cache_ecpichia_supplement_source(
        _profile_for(yaml_path), tmp_path / "yaml-cache", source_file=yaml_path
    )
    table_path = tmp_path / "Supplementary 11.docx"
    _write_table_docx(table_path)
    table_source = cache_ecpichia_supplement_table_source(
        _table_profile_for(table_path), tmp_path / "table-cache", source_file=table_path
    )
    binding = CapacityModelBinding(
        target_id="hLF",
        context_id="glucose_mu_0.1",
        mapping_id="hLF:g6pdh2",
        model_fingerprint="fingerprint:hLF",
        gene_id="PAS_chr2-1_0308",
        enzyme_id="G6PDH2_no_1_fwd_complex",
        reaction_id="G6PDH2_no_1_fwd",
        formation_or_dilution_reaction_id="G6PDH2_no_1_fwd_complex_formation",
        mapping_evidence=("current_model_mapping",),
    )
    closure = evaluate_ecpichia_g6pdh2_provenance(
        parse_ecpichia_g6pdh2_source_assessment(yaml_source),
        parse_ecpichia_g6pdh2_table_evidence(table_source),
        model_bindings=(binding,),
        formal_context_id="glucose_mu_0.1",
    )
    assert closure.completion_outcome == "architecture_decision_required"
    assert len(closure.source_artifacts) == 2
    assert all(len(item["raw_sha256"]) == 64 for item in closure.source_artifacts)
    assert closure.coefficient_trace["matches"] is True
    assert closure.coefficient_trace["gecko_reference"]["source_artifact_verified"] is False
    assert (
        closure.coefficient_trace["kcat_provenance_assessment"][
            "source_artifact_verified"
        ]
        is False
    )
    assert closure.nominal_capacity is None
    assert closure.promotion_preview_available is False
    assert "supplement_table_yaml_gene_conflict" in closure.source_conflicts
    assert "supplement_table_yaml_enzyme_conflict" in closure.source_conflicts
    assert "concentration_unit_conflict_g_per_L_vs_mg_per_gDCW" in closure.missing_information
    assert (
        "kcat_research_artifact_not_bound_to_source_cache"
        in closure.missing_information
    )
    assert closure.model_bindings[0].context_id == "glucose_mu_0.1"


def test_ecpichia_closure_records_formula_mismatch_instead_of_promoting(
    tmp_path: Path,
) -> None:
    yaml_path = tmp_path / "Supplementary 8.yml"
    yaml_path.write_text(
        _fixture_yaml().replace("prot_C4R099: -0.00200309027777778", "prot_C4R099: -0.003"),
        encoding="utf-8",
    )
    yaml_source = cache_ecpichia_supplement_source(
        _profile_for(yaml_path), tmp_path / "yaml-cache", source_file=yaml_path
    )
    table_path = tmp_path / "Supplementary 11.docx"
    _write_table_docx(table_path)
    table_source = cache_ecpichia_supplement_table_source(
        _table_profile_for(table_path), tmp_path / "table-cache", source_file=table_path
    )
    closure = evaluate_ecpichia_g6pdh2_provenance(
        parse_ecpichia_g6pdh2_source_assessment(yaml_source),
        parse_ecpichia_g6pdh2_table_evidence(table_source),
        model_bindings=(),
        formal_context_id="glucose_mu_0.1",
    )
    assert closure.coefficient_trace["matches"] is False
    assert "yaml_coefficient_mw_kcat_formula_mismatch" in closure.source_conflicts
    assert closure.nominal_capacity is None


def test_ecpichia_profile_and_offline_metadata_are_not_forgeable(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "Supplementary 8.yml"
    source_path.write_text(_fixture_yaml(), encoding="utf-8")
    profile = _profile_for(source_path)
    with pytest.raises(OECapacityValidationError, match="license_id"):
        cache_ecpichia_supplement_source(
            replace(profile, license_id=""),
            tmp_path / "invalid-license",
            source_file=source_path,
        )

    cache_dir = tmp_path / "cache"
    cache_ecpichia_supplement_source(
        profile,
        cache_dir,
        source_file=source_path,
    )
    metadata_path = cache_dir / "ecpichia-supplementary-8.source.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["query"] = json.dumps({"gene_id": "forged"})
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(OECapacityValidationError, match="reviewed query"):
        cache_ecpichia_supplement_source(
            profile,
            cache_dir,
            retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
        )


def test_ecpichia_assessment_never_creates_capacity_or_promotion_ready_candidate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_path = tmp_path / "Supplementary 8.yml"
    source_path.write_text(_fixture_yaml(), encoding="utf-8")
    profile = _profile_for(source_path)
    monkeypatch.setattr(
        audit_module,
        "ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE",
        profile,
    )
    monkeypatch.setattr(
        audit_module,
        "prepare_external_candidate_runtime",
        lambda *args, **kwargs: SimpleNamespace(gene_capacity_catalog={}),
    )
    monkeypatch.setattr(
        audit_module,
        "build_capacity_model_binding",
        lambda catalog, **kwargs: CapacityModelBinding(
            target_id=kwargs["target_id"],
            context_id=kwargs["context_id"],
            mapping_id=f"{kwargs['target_id']}:g6pdh2",
            model_fingerprint=f"fingerprint:{kwargs['target_id']}",
            gene_id="PAS_chr2-1_0308",
            enzyme_id="G6PDH2_no_1_fwd_complex",
            reaction_id="G6PDH2_no_1_fwd",
            formation_or_dilution_reaction_id=(
                "G6PDH2_no_1_fwd_complex_formation"
            ),
            mapping_evidence=tuple(kwargs["mapping_evidence"]),
            external_gene_id=kwargs["external_gene_id"],
            external_protein_id=kwargs["external_protein_id"],
        ),
    )
    identity_cache = tmp_path / "identity"
    _cache_identity(identity_cache)
    ecpichia_cache = tmp_path / "ecpichia-cache"
    cache_ecpichia_supplement_source(
        profile,
        ecpichia_cache,
        source_file=source_path,
    )

    outputs = run_external_capacity_candidate_audit(
        ExternalCapacityAuditRequest(
            repo_root=REPO_ROOT,
            output_dir=tmp_path / "local_runs" / "oe_capacity" / "round6a" / "a0b",
            offline_replay=True,
            identity_cache_dir=identity_cache,
            ecpichia_cache_dir=ecpichia_cache,
            ecpichia_assessment=True,
        )
    )
    assert outputs.candidate_count == 0
    assert outputs.promotion_ready_count == 0
    assert outputs.summary()["formal_promotion_performed"] is False
    audit = json.loads(outputs.audit_json_path.read_text(encoding="utf-8"))
    assessment = next(
        item
        for item in audit["source_assessments"]
        if item["assessment_id"] == "ecpichia-supplementary-8-g6pdh2"
    )
    assert assessment["raw_values"]["reported_concentration"] == pytest.approx(
        0.752073171936811
    )
    assert assessment["promotion_ready"] is False
    assert assessment["capacity_value_available"] is False
    assert (
        "source_unit_missing_and_supplement_header_requires_review"
        in assessment["conflicts"]
    )
    assert (
        "supplement_table_yaml_binding_requires_reconciliation"
        in assessment["conflicts"]
    )
    assert audit["formal_promotion_performed"] is False


def test_a0c_audit_writes_architecture_decision_gap_without_capacity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    yaml_path = tmp_path / "Supplementary 8.yml"
    yaml_path.write_text(_fixture_yaml(), encoding="utf-8")
    table_path = tmp_path / "Supplementary 11.docx"
    _write_table_docx(table_path)
    yaml_profile = _profile_for(yaml_path)
    table_profile = _table_profile_for(table_path)
    monkeypatch.setattr(
        audit_module, "ECPICHIA_G6PDH2_SUPPLEMENT_PROFILE", yaml_profile
    )
    monkeypatch.setattr(
        audit_module, "ECPICHIA_G6PDH2_TABLE_PROFILE", table_profile
    )
    monkeypatch.setattr(
        audit_module,
        "prepare_external_candidate_runtime",
        lambda *args, **kwargs: SimpleNamespace(gene_capacity_catalog={}),
    )
    monkeypatch.setattr(
        audit_module,
        "build_capacity_model_binding",
        lambda catalog, **kwargs: CapacityModelBinding(
            target_id=kwargs["target_id"],
            context_id=kwargs["context_id"],
            mapping_id=f"{kwargs['target_id']}:g6pdh2",
            model_fingerprint=f"fingerprint:{kwargs['target_id']}",
            gene_id="PAS_chr2-1_0308",
            enzyme_id="G6PDH2_no_1_fwd_complex",
            reaction_id="G6PDH2_no_1_fwd",
            formation_or_dilution_reaction_id="G6PDH2_no_1_fwd_complex_formation",
            mapping_evidence=tuple(kwargs["mapping_evidence"]),
            external_gene_id=kwargs["external_gene_id"],
            external_protein_id=kwargs["external_protein_id"],
            external_enzyme_id=kwargs.get("external_enzyme_id", ""),
        ),
    )
    identity_cache = tmp_path / "identity"
    _cache_identity(identity_cache)
    yaml_cache = tmp_path / "yaml-cache"
    table_cache = tmp_path / "table-cache"
    cache_ecpichia_supplement_source(
        yaml_profile, yaml_cache, source_file=yaml_path
    )
    cache_ecpichia_supplement_table_source(
        table_profile, table_cache, source_file=table_path
    )
    outputs = run_external_capacity_candidate_audit(
        ExternalCapacityAuditRequest(
            repo_root=REPO_ROOT,
            output_dir=tmp_path / "local_runs" / "oe_capacity" / "round6a" / "a0c",
            offline_replay=True,
            identity_cache_dir=identity_cache,
            ecpichia_cache_dir=yaml_cache,
            ecpichia_table_cache_dir=table_cache,
            ecpichia_provenance_closure=True,
        )
    )
    assert outputs.candidate_count == 0
    assert outputs.promotion_ready_count == 0
    assert outputs.completion_outcome == "architecture_decision_required"
    gap = json.loads(outputs.provenance_gap_json_path.read_text(encoding="utf-8"))
    assert gap["nominal_capacity"] is None
    assert gap["promotion_preview_available"] is False
    assert {item["context_id"] for item in gap["model_bindings"]} == {
        "glucose_mu_0.1"
    }
    assert outputs.provenance_gap_markdown_path.is_file()
    audit = json.loads(outputs.audit_json_path.read_text(encoding="utf-8"))
    assert {item["context_id"] for item in audit["formal_model_bindings"]} == {
        "glucose_mu_0.1"
    }


def test_cli_imports_only_public_audit_api() -> None:
    cli = (REPO_ROOT / "python_pichia" / "tools" / "run_oe_capacity_external_candidate_audit.py").read_text(
        encoding="utf-8"
    )
    assert "from pcsec_pichia.oe_capacity.external_candidate_audit import" in cli
    assert "external_candidate_io import" not in cli
    assert "external_refs.capacity_sources import" not in cli
    assert "--ecpichia-provenance-closure" in cli
    assert "--ecpichia-table-file" in cli
