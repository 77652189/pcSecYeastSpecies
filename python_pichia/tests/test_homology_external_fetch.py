from __future__ import annotations

import importlib.util
import importlib
import json
from pathlib import Path
import sys
import urllib.request

from pcsec_pichia.homology.cache_schema import ExternalNameReference, HomologyCrosswalkRow, NameAuditRow
from pcsec_pichia.homology.crosswalk import build_name_audit_rows, load_external_name_reference_cache
from pcsec_pichia.homology.external_fetch import (
    ExternalFetchConfig,
    ExternalFetchResult,
    HttpResponse,
    fetch_ncbi_name_reference,
    fetch_sgd_name_reference,
    fetch_uniprot_name_reference,
)


def test_uniprot_fake_response_builds_external_name_reference() -> None:
    def fake_get(url: str, config: ExternalFetchConfig) -> HttpResponse:
        assert "rest.uniprot.org/uniprotkb/search" in url
        assert "KAR2" in url
        return HttpResponse(
            200,
            json.dumps(
                {
                    "results": [
                        {
                            "primaryAccession": "C4R8K4",
                            "uniProtkbId": "KAR2_KOMPG",
                            "organism": {"scientificName": "Komagataella phaffii"},
                            "genes": [
                                {
                                    "geneName": {"value": "KAR2"},
                                    "orderedLocusNames": [{"value": "PAS_chr2-1_0140"}],
                                    "synonyms": [{"value": "BiP"}],
                                }
                            ],
                        }
                    ]
                }
            ),
            {"X-UniProt-Release": "2026_01"},
        )

    result = fetch_uniprot_name_reference("KAR2", ExternalFetchConfig(), http_get=fake_get, sleep=lambda _: None)

    assert result.success is True
    ref = result.references[0]
    assert ref.source_database == "UniProt"
    assert ref.source_version == "2026_01"
    assert ref.taxon == "Komagataella phaffii"
    assert ref.accession == "C4R8K4"
    assert ref.gene_name == "KAR2"
    assert ref.locus_tag == "PAS_chr2-1_0140"
    assert ref.aliases == ("BiP",)
    assert ref.retrieved_at


def test_ncbi_fake_gene_response_builds_external_name_reference_without_api_key() -> None:
    seen_urls: list[str] = []

    def fake_get(url: str, config: ExternalFetchConfig) -> HttpResponse:
        seen_urls.append(url)
        assert "api_key=" not in url
        if "esearch.fcgi" in url:
            assert "tool=pcSecYeastSpecies" in url
            assert "email=curator%40example.org" in url
            return HttpResponse(200, json.dumps({"esearchresult": {"idlist": ["12345"]}}))
        return HttpResponse(
            200,
            json.dumps(
                {
                    "result": {
                        "uids": ["12345"],
                        "12345": {
                            "name": "KAR2",
                            "nomenclaturesymbol": "KAR2",
                            "description": "BiP chaperone",
                            "otheraliases": "BiP, GRP78",
                            "organism": {"scientificname": "Komagataella phaffii"},
                        },
                    }
                }
            ),
        )

    config = ExternalFetchConfig(ncbi_email="curator@example.org", retry_count=0)
    result = fetch_ncbi_name_reference("PAS_chr2-1_0140", config, http_get=fake_get, sleep=lambda _: None)

    assert len(seen_urls) == 2
    assert result.success is True
    ref = result.references[0]
    assert ref.source_database == "NCBI"
    assert ref.source_version == "gene"
    assert ref.accession == "12345"
    assert ref.gene_name == "KAR2"
    assert ref.taxon == "Komagataella phaffii"
    assert "BiP" in ref.aliases
    assert "GRP78" in ref.aliases


def test_sgd_fake_response_builds_external_name_reference() -> None:
    def fake_get(url: str, config: ExternalFetchConfig) -> HttpResponse:
        assert url.endswith("/KAR2")
        return HttpResponse(
            200,
            json.dumps(
                {
                    "sgdid": "S000003120",
                    "format_name": "KAR2",
                    "systematic_name": "YJL034W",
                    "aliases": ["BiP", "GRP78"],
                    "organism": "Saccharomyces cerevisiae",
                }
            ),
            {"Last-Modified": "Wed, 08 Jul 2026 00:00:00 GMT"},
        )

    result = fetch_sgd_name_reference("KAR2", ExternalFetchConfig(), http_get=fake_get, sleep=lambda _: None)

    assert result.success is True
    ref = result.references[0]
    assert ref.source_database == "SGD"
    assert ref.source_version == "Wed, 08 Jul 2026 00:00:00 GMT"
    assert ref.accession == "S000003120"
    assert ref.gene_name == "KAR2"
    assert ref.locus_tag == "YJL034W"
    assert ref.aliases == ("BiP", "GRP78")


def test_fetch_timeout_records_warning_and_retries_only_configured_count() -> None:
    calls: list[str] = []
    sleeps: list[float] = []

    def fake_get(url: str, config: ExternalFetchConfig) -> HttpResponse:
        calls.append(url)
        raise TimeoutError("network timeout")

    result = fetch_uniprot_name_reference(
        "KAR2",
        ExternalFetchConfig(retry_count=2, delay_seconds=0.01),
        http_get=fake_get,
        sleep=sleeps.append,
    )

    assert result.success is False
    assert len(calls) == 3
    assert sleeps == [0.01, 0.01]
    assert result.references == ()
    assert "TimeoutError" in result.warnings[-1]
    assert result.error_summary


def test_external_fetch_module_does_not_connect_on_import(monkeypatch) -> None:
    def fail_urlopen(*args, **kwargs):  # pragma: no cover - should never run
        raise AssertionError("import attempted network access")

    monkeypatch.setattr(urllib.request, "urlopen", fail_urlopen)
    import pcsec_pichia.homology.external_fetch as external_fetch

    importlib.reload(external_fetch)


def test_external_reference_builder_collects_deduped_name_audit_queries() -> None:
    builder = _load_builder_module()
    rows = (
        _name_row("KAR2 / BiP", "YJL034W", "C4R8K4", "KAR2", "PAS_chr2-1_0140"),
        _name_row("KAR2 duplicate", "YJL034W", "C4R8K4", "KAR2", "PAS_chr2-1_0140"),
    )

    queries = builder.collect_external_reference_queries(rows)

    assert [(query.match_key, query.query) for query in queries] == [
        ("external_accession", "C4R8K4"),
        ("external_locus_tag", "PAS_chr2-1_0140"),
        ("external_gene_name", "KAR2"),
        ("internal_common_name", "KAR2 / BiP"),
        ("internal_sequence_id", "YJL034W"),
        ("internal_common_name", "KAR2 duplicate"),
    ]


def test_external_reference_builder_writes_loadable_jsonl_and_partial_failure_summary(tmp_path: Path) -> None:
    builder = _load_builder_module()
    rows = (
        _name_row("KAR2 / BiP", "YJL034W", "C4R8K4", "KAR2", "PAS_chr2-1_0140"),
        _name_row("HRD1", "YOL013C", "C4H", "HRD1", "PAS_chr4_0156"),
    )

    def fake_fetch(query: str, config: ExternalFetchConfig) -> tuple[ExternalFetchResult, ...]:
        if query == "C4R8K4":
            return (
                ExternalFetchResult(
                    "UniProt",
                    query,
                    True,
                    (
                        ExternalNameReference(
                            "UniProt",
                            "2026_01",
                            "Komagataella phaffii",
                            "C4R8K4",
                            "KAR2",
                            "PAS_chr2-1_0140",
                            ("BiP",),
                            "2026-07-08T00:00:00+00:00",
                        ),
                    ),
                ),
            )
        return (ExternalFetchResult("UniProt", query, False, warnings=("simulated failure",)),)

    references, summary = builder.build_external_reference_cache(
        rows,
        config=ExternalFetchConfig(enabled_sources=("uniprot",), retry_count=0),
        sources=("uniprot",),
        limit=2,
        fetch_many=fake_fetch,
    )
    output = tmp_path / "external_name_references.jsonl"
    builder.write_external_reference_cache(references, output)

    loaded = load_external_name_reference_cache(output)

    assert summary["query_count"] == 2
    assert summary["success_count"] == 1
    assert summary["failure_count"] == 1
    assert summary["source_counts"] == {"UniProt": 1}
    assert any("simulated failure" in warning for warning in summary["warnings"])
    assert loaded[0].source_database == "UniProt"
    assert loaded[0].warnings == ("query=C4R8K4; match_key=external_accession",)


def test_external_reference_builder_output_drives_name_audit_crosscheck_status() -> None:
    crosswalk = (
        HomologyCrosswalkRow(
            internal_common_name="KAR2 / BiP",
            query_symbol="KAR2",
            sce_orf="YJL034W",
            pichia_gene_id="PAS_chr2-1_0140",
            pichia_model_gene_id="PAS_chr2-1_0140",
            is_rbh=True,
            identity_pct=75.0,
            evalue=1e-100,
            query_coverage=95.0,
            subject_coverage=95.0,
            in_model_gene_index=True,
            review_status="model_ready_rbh_high_confidence",
            external_accession="C4R8K4",
            external_gene_name="KAR2",
            external_locus_tag="PAS_chr2-1_0140",
        ),
    )
    references = (
        ExternalNameReference(
            "UniProt",
            "2026_01",
            "Komagataella phaffii",
            "C4R8K4",
            "KAR2",
            "PAS_chr2-1_0140",
            ("BiP",),
            "2026-07-08T00:00:00+00:00",
        ),
    )

    name_rows = build_name_audit_rows(crosswalk, external_references=references)

    assert name_rows[0].external_crosscheck_status == "external_match_confirmed"


def _load_builder_module():
    path = Path(__file__).resolve().parents[2] / "scripts" / "build_pichia_external_name_reference_cache.py"
    spec = importlib.util.spec_from_file_location("build_pichia_external_name_reference_cache", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _name_row(
    common_name: str,
    sce_orf: str,
    accession: str,
    gene_name: str,
    locus_tag: str,
) -> NameAuditRow:
    return NameAuditRow(
        internal_gene_id=locus_tag,
        internal_common_name=common_name,
        internal_sequence_id=sce_orf,
        external_accession=accession,
        external_gene_name=gene_name,
        external_locus_tag=locus_tag,
        external_aliases=(),
        identity_pct=75.0,
        query_coverage=95.0,
        subject_coverage=95.0,
        evalue=1e-100,
        is_rbh=True,
        in_model_gene_index=True,
        name_consistency_status="name_confirmed_by_rbh",
        review_status="model_ready_rbh_high_confidence",
    )
