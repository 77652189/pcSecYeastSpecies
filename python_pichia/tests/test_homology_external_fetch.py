from __future__ import annotations

import importlib
import json
import urllib.request

from pcsec_pichia.homology.external_fetch import (
    ExternalFetchConfig,
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
