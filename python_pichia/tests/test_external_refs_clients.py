from __future__ import annotations

import json
from dataclasses import dataclass

from pcsec_pichia.external_refs import (
    ExternalFetchConfig,
    ExternalFetchResult,
    ExternalHttpResponse,
    ExternalReferenceQuery,
    build_external_reference_cache,
    fetch_external_references,
    load_external_reference_cache,
    load_external_reference_manifest,
)
from pcsec_pichia.external_refs.ncbi import fetch_ncbi_gene_reference
from pcsec_pichia.external_refs.sgd import fetch_sgd_reference
from pcsec_pichia.external_refs.uniprot import fetch_uniprot_reference


def test_uniprot_client_builds_schema_record_with_provenance() -> None:
    query = _query("external_accession", "P12345", preferred_sources=("uniprot",))

    result = fetch_uniprot_reference(
        query,
        ExternalFetchConfig(retry_attempts=1, min_interval_seconds=0),
        http_get=lambda url, config: _response(url, _uniprot_payload()),
        sleep=lambda _seconds: None,
    )

    assert result.success is True
    assert result.http_status == 200
    assert result.raw_record_sha256
    assert result.records[0].provenance.source_database == "uniprot"
    assert result.records[0].provenance.source_url.startswith("https://rest.uniprot.org/")
    assert result.records[0].primary_accession == "P12345"
    assert result.records[0].gene_name == "SEC1"
    assert result.records[0].locus_tag == "PAS_chr1-1_0001"
    assert result.records[0].taxon_id == "4922"
    assert result.records[0].reviewed is True


def test_sgd_client_builds_sce_locus_reference_record() -> None:
    query = _query("sce_homolog", "YDR164C", preferred_sources=("sgd",))

    result = fetch_sgd_reference(
        query,
        ExternalFetchConfig(retry_attempts=1, min_interval_seconds=0),
        http_get=lambda url, config: _response(url, _sgd_payload(), headers={"last-modified": "2026-07-09"}),
        sleep=lambda _seconds: None,
    )

    assert result.success is True
    assert result.records[0].provenance.source_database == "sgd"
    assert result.records[0].provenance.source_version == "2026-07-09"
    assert result.records[0].organism == "Saccharomyces cerevisiae"
    assert result.records[0].taxon_id == "559292"
    assert result.records[0].gene_name == "SEC1"
    assert result.records[0].locus_tag == "YDR164C"


def test_ncbi_client_uses_env_api_key_and_two_step_fetch(monkeypatch) -> None:
    monkeypatch.setenv("NCBI_API_KEY", "test-key")
    seen_urls: list[str] = []

    def fake_get(url: str, config: ExternalFetchConfig) -> ExternalHttpResponse:
        seen_urls.append(url)
        if "esearch.fcgi" in url:
            return _response(url, {"esearchresult": {"idlist": ["123"]}})
        return _response(url, _ncbi_summary_payload())

    result = fetch_ncbi_gene_reference(
        _query("model_gene", "PAS_chr1-1_0001", preferred_sources=("ncbi",)),
        ExternalFetchConfig(retry_attempts=1, min_interval_seconds=0),
        http_get=fake_get,
        sleep=lambda _seconds: None,
    )

    assert result.success is True
    assert len(seen_urls) == 2
    assert all("api_key=test-key" in url for url in seen_urls)
    assert result.records[0].provenance.source_database == "ncbi"
    assert result.records[0].primary_accession == "SEC1"
    assert result.records[0].taxon_id == "4922"


def test_client_retry_records_failure_without_network() -> None:
    attempts: list[str] = []

    def fake_get(url: str, config: ExternalFetchConfig) -> ExternalHttpResponse:
        attempts.append(url)
        return ExternalHttpResponse(status_code=503, text="service unavailable", url=url)

    result = fetch_uniprot_reference(
        _query("external_accession", "P00000", preferred_sources=("uniprot",)),
        ExternalFetchConfig(retry_attempts=2, min_interval_seconds=0),
        http_get=fake_get,
        sleep=lambda _seconds: None,
    )

    assert result.failed is True
    assert result.failure is not None
    assert result.failure.http_status == 503
    assert result.failure.raw_record_sha256
    assert result.attempts == 2
    assert len(attempts) == 2


def test_fetch_external_references_respects_preferred_sources() -> None:
    query = _query("sce_homolog", "YDR164C", preferred_sources=("sgd",))
    clients = (_FakeClient("uniprot"), _FakeClient("sgd"), _FakeClient("ncbi"))

    results = fetch_external_references(
        (query,),
        clients,
        ExternalFetchConfig(sources=("uniprot", "sgd", "ncbi"), min_interval_seconds=0),
        sleep=lambda _seconds: None,
    )

    assert [result.source_database for result in results] == ["sgd"]
    assert clients[0].calls == 0
    assert clients[1].calls == 1
    assert clients[2].calls == 0


def test_fetch_external_references_applies_rate_limit_between_dispatches() -> None:
    query = _query("external_accession", "P12345", preferred_sources=())
    sleeps: list[float] = []

    results = fetch_external_references(
        (query,),
        (_FakeClient("uniprot"), _FakeClient("sgd")),
        ExternalFetchConfig(sources=("uniprot", "sgd"), min_interval_seconds=0.5),
        sleep=lambda seconds: sleeps.append(seconds),
    )

    assert len(results) == 2
    assert sleeps == [0.5]


def test_build_external_reference_cache_writes_manifest_records_and_failures(tmp_path) -> None:
    queries = (
        _query("external_accession", "P12345", preferred_sources=("uniprot",)),
        _query("sce_homolog", "MISSING", preferred_sources=("sgd",)),
    )

    manifest = build_external_reference_cache(
        queries,
        tmp_path,
        config=ExternalFetchConfig(sources=("uniprot", "sgd"), retry_attempts=1, min_interval_seconds=0),
        http_get=lambda url, config: (
            _response(url, _uniprot_payload())
            if "uniprot" in url
            else ExternalHttpResponse(status_code=404, text="{}", url=url)
        ),
        sleep=lambda _seconds: None,
    )

    records = load_external_reference_cache(tmp_path / "external_reference_records.jsonl")
    stored_manifest = load_external_reference_manifest(tmp_path / "external_reference_manifest.json")
    failed_lines = (tmp_path / "failed_queries.jsonl").read_text(encoding="utf-8").splitlines()

    assert manifest == stored_manifest
    assert manifest.query_count == 2
    assert manifest.record_count == 1
    assert manifest.failed_query_count == 1
    assert manifest.input_cache_fingerprint
    assert records[0].primary_accession == "P12345"
    assert len(failed_lines) == 1
    assert json.loads(failed_lines[0])["source_database"] == "sgd"
    assert (tmp_path / "external_reference_summary.md").exists()


def test_build_external_reference_cache_counts_failed_queries_not_failed_sources(tmp_path) -> None:
    query = _query("pichia_gene", "MISSING", preferred_sources=())

    manifest = build_external_reference_cache(
        (query,),
        tmp_path,
        config=ExternalFetchConfig(sources=("uniprot", "sgd"), retry_attempts=1, min_interval_seconds=0),
        http_get=lambda url, config: ExternalHttpResponse(status_code=404, text="{}", url=url),
        sleep=lambda _seconds: None,
    )
    failed_lines = (tmp_path / "failed_queries.jsonl").read_text(encoding="utf-8").splitlines()

    assert manifest.query_count == 1
    assert manifest.failed_query_count == 1
    assert len(failed_lines) == 2


@dataclass
class _FakeClient:
    source_database: str
    calls: int = 0

    def fetch(self, query, config, *, http_get=None, sleep=lambda _seconds: None):
        self.calls += 1
        return ExternalFetchResult(source_database=self.source_database, query=query, success=False)


def _query(query_type: str, query_value: str, *, preferred_sources: tuple[str, ...]) -> ExternalReferenceQuery:
    return ExternalReferenceQuery(
        query_type=query_type,
        query_value=query_value,
        source_context="test",
        source_id="test",
        source_row_id=query_value,
        preferred_sources=preferred_sources,
    )


def _response(url: str, payload: object, headers: dict[str, str] | None = None) -> ExternalHttpResponse:
    return ExternalHttpResponse(
        status_code=200,
        text=json.dumps(payload),
        url=url,
        headers=headers or {},
    )


def _uniprot_payload() -> dict[str, object]:
    return {
        "results": [
            {
                "primaryAccession": "P12345",
                "uniProtkbId": "SEC1_KPH",
                "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "organism": {"taxonId": 4922, "scientificName": "Komagataella phaffii"},
                "genes": [
                    {
                        "geneName": {"value": "SEC1"},
                        "orderedLocusNames": [{"value": "PAS_chr1-1_0001"}],
                        "synonyms": [{"value": "secretory protein 1"}],
                    }
                ],
                "proteinDescription": {
                    "recommendedName": {"fullName": {"value": "Secretion protein Sec1"}}
                },
            }
        ]
    }


def _sgd_payload() -> dict[str, object]:
    return {
        "sgdid": "S000002571",
        "format_name": "SEC1",
        "systematic_name": "YDR164C",
        "aliases": [{"display_name": "SRO7"}],
        "name_description": "Protein involved in secretion",
    }


def _ncbi_summary_payload() -> dict[str, object]:
    return {
        "result": {
            "uids": ["123"],
            "123": {
                "uid": "123",
                "caption": "SEC1",
                "name": "SEC1",
                "description": "Secretion protein",
                "otheraliases": "secretory protein 1",
                "organism": {"taxid": 4922, "scientificname": "Komagataella phaffii"},
            },
        }
    }
