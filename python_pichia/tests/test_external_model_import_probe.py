from __future__ import annotations

import json
from pathlib import Path

from pcsec_pichia.external_refs import (
    ARTIFACT_MANIFEST_FILENAME,
    ExternalModelArtifactRequest,
    ExternalModelImportProbeRequest,
    cache_external_model_artifacts,
    cobrapy_import_available,
    load_import_probe_requests_from_artifact_cache,
    probe_cobrapy_model_import,
    probe_external_model_imports,
)
from pcsec_pichia.external_refs import model_import_probe


def test_cobrapy_import_probe_returns_unavailable_without_optional_dependency(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(model_import_probe, "_import_cobra", lambda: None)
    request = ExternalModelImportProbeRequest(
        model_id="toy_model",
        artifact_path=str(tmp_path / "toy.xml"),
        artifact_type="SBML",
    )

    result = probe_cobrapy_model_import(request)
    outputs = probe_external_model_imports((request,), tmp_path / "probe")

    assert result.backend == "cobrapy"
    assert result.backend_available is False
    assert result.import_status == "unavailable"
    assert result.manual_review_required is True
    assert "cobrapy_unavailable" in result.warnings
    assert outputs.unavailable_count == 1
    payload = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    assert payload["results"][0]["backend_available"] is False


def test_cobrapy_import_available_handles_broken_import(monkeypatch) -> None:
    def broken_import(name: str):
        raise ImportError(f"broken optional dependency: {name}")

    monkeypatch.setattr(model_import_probe.importlib, "import_module", broken_import)

    assert cobrapy_import_available() is False


def test_import_probe_requests_load_downloaded_artifact_cache_rows(tmp_path) -> None:
    outputs = cache_external_model_artifacts(
        (
            ExternalModelArtifactRequest(
                model_id="toy_model",
                artifact_url="https://example.test/toy.xml",
                artifact_type="SBML",
                filename="toy.xml",
                source_page_url="https://example.test/model",
            ),
        ),
        tmp_path / "artifacts",
        fetcher=lambda _url, _timeout: b"<sbml></sbml>",
    )

    requests = load_import_probe_requests_from_artifact_cache(outputs.manifest_path)
    requests_from_dir = load_import_probe_requests_from_artifact_cache(outputs.manifest_path.parent)

    assert (outputs.manifest_path.parent / ARTIFACT_MANIFEST_FILENAME).exists()
    assert requests == requests_from_dir
    assert len(requests) == 1
    assert requests[0].model_id == "toy_model"
    assert requests[0].artifact_type == "SBML"
    assert Path(requests[0].artifact_path).exists()


def test_import_probe_reports_fake_cobrapy_model_diagnostics(monkeypatch, tmp_path) -> None:
    artifact_path = tmp_path / "toy.xml"
    artifact_path.write_text("<sbml></sbml>", encoding="utf-8")
    request = ExternalModelImportProbeRequest(
        model_id="toy_model",
        artifact_path=str(artifact_path),
        artifact_type="SBML",
        source_page_url="https://example.test/model",
    )

    monkeypatch.setattr(model_import_probe, "_compare_libsbml_import", lambda *_args: ("aligned", ()))
    outputs = probe_external_model_imports((request,), tmp_path / "probe", cobra_module=_FakeCobra())
    result = outputs.results[0]
    manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    report = outputs.report_path.read_text(encoding="utf-8")

    assert result.import_status == "imported"
    assert result.backend_available is True
    assert result.reaction_count == 3
    assert result.metabolite_count == 2
    assert result.gene_count == 2
    assert result.gpr_count == 2
    assert result.objective_reaction == "R_TARGET"
    assert result.libsbml_comparison_status == "aligned"
    assert result.id_sanitization_warnings == ()
    assert result.manual_review_required is False
    assert manifest["imported_count"] == 1
    assert "recommendation tiers" in report
    assert "mg/L" not in report
    assert "absolute secretion yield" in report


def test_import_probe_requires_review_when_core_model_semantics_are_missing(monkeypatch, tmp_path) -> None:
    artifact_path = tmp_path / "toy.xml"
    artifact_path.write_text("<sbml></sbml>", encoding="utf-8")
    monkeypatch.setattr(model_import_probe, "_compare_libsbml_import", lambda *_args: ("aligned", ()))

    result = probe_cobrapy_model_import(
        ExternalModelImportProbeRequest(model_id="toy_model", artifact_path=str(artifact_path), artifact_type="SBML"),
        cobra_module=_FakeCobraEmptySemantics(),
    )

    assert result.import_status == "imported"
    assert result.manual_review_required is True
    assert result.objective_reaction == ""
    assert "no_metabolites_detected" in result.warnings
    assert "no_gpr_rules_detected" in result.warnings
    assert "objective_reaction_not_detected" in result.warnings


class _FakeCobra:
    io = None

    def __init__(self) -> None:
        self.io = _FakeCobraIo()


class _FakeCobraIo:
    def read_sbml_model(self, path: str) -> object:
        assert Path(path).exists()
        return _FakeModel()


class _FakeCobraEmptySemantics:
    class io:
        @staticmethod
        def read_sbml_model(_path: str) -> object:
            return type(
                "EmptySemanticModel",
                (),
                {
                    "reactions": (_FakeReaction("R1", "", 0.0),),
                    "metabolites": (),
                    "genes": (),
                },
            )()


class _FakeReaction:
    def __init__(self, reaction_id: str, gene_reaction_rule: str, objective_coefficient: float) -> None:
        self.id = reaction_id
        self.gene_reaction_rule = gene_reaction_rule
        self.objective_coefficient = objective_coefficient


class _FakeModel:
    reactions = (
        _FakeReaction("R_SOURCE", "", 0.0),
        _FakeReaction("R_HELPER", "PAS_chr1-1", 0.0),
        _FakeReaction("R_TARGET", "PAS_chr2-2 and PAS_chr3-3", 1.0),
    )
    metabolites = ("m1", "m2")
    genes = ("PAS_chr1-1", "PAS_chr2-2")
