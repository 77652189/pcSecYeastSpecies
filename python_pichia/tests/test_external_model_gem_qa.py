from __future__ import annotations

import json
from pathlib import Path

from pcsec_pichia.external_refs import (
    ExternalModelArtifactRequest,
    ExternalModelGemQaRequest,
    cache_external_model_artifacts,
    gem_qa_requests_from_artifact_cache,
    run_external_model_gem_qa,
)
from pcsec_pichia.external_refs import gem_qa, model_import_probe


def test_gem_qa_returns_unavailable_without_cobrapy(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(model_import_probe, "_import_cobra", lambda: None)
    request = ExternalModelGemQaRequest(
        model_id="toy_model",
        artifact_path=str(tmp_path / "toy.xml"),
        artifact_type="SBML",
    )

    outputs = run_external_model_gem_qa((request,), tmp_path / "qa")
    result = outputs.results[0]

    assert result.backend_available is False
    assert result.qa_status == "unavailable"
    assert result.import_status == "unavailable"
    assert result.recommendation_tier_effect == "none"
    assert result.manual_review_reasons == ("unavailable",)
    assert outputs.unavailable_count == 1


def test_gem_qa_runs_basic_cobrapy_checks_and_keeps_metadata_boundary(tmp_path) -> None:
    artifact_path = tmp_path / "toy.xml"
    artifact_path.write_text("<sbml></sbml>", encoding="utf-8")
    request = ExternalModelGemQaRequest(
        model_id="toy_model",
        artifact_path=str(artifact_path),
        artifact_type="SBML",
    )

    outputs = run_external_model_gem_qa((request,), tmp_path / "qa", cobra_module=_FakeCobra(_review_model))
    result = outputs.results[0]
    manifest = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    report = outputs.report_path.read_text(encoding="utf-8")

    assert result.qa_backend == "cobrapy-basic"
    assert result.backend_available is True
    assert result.qa_status == "review_required"
    assert result.import_status == "imported"
    assert result.gpr_coverage == 0.5
    assert result.blocked_reaction_count == 1
    assert result.dead_end_metabolite_count == 2
    assert result.manual_review_reasons == (
        "blocked_reactions_present",
        "dead_end_metabolites_present",
    )
    assert result.recommendation_tier_effect == "none"
    assert manifest["review_required_count"] == 1
    assert "recommendation_tier" in report
    assert "mg/L" not in report
    assert "experimental success rate" in report


def test_gem_qa_keeps_basic_result_when_memote_is_unavailable(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(gem_qa, "_import_memote", lambda: None)
    artifact_path = tmp_path / "toy.xml"
    artifact_path.write_text("<sbml></sbml>", encoding="utf-8")
    request = ExternalModelGemQaRequest(
        model_id="toy_model",
        artifact_path=str(artifact_path),
        artifact_type="SBML",
        run_memote=True,
    )

    outputs = run_external_model_gem_qa((request,), tmp_path / "qa", cobra_module=_FakeCobra(_clean_model))
    result = outputs.results[0]

    assert result.qa_status == "passed_basic"
    assert result.memote_available is False
    assert result.memote_status == "unavailable"
    assert "memote_unavailable" in result.warnings
    assert outputs.passed_basic_count == 1


def test_gem_qa_records_optional_memote_score_without_recommendation_effect(tmp_path) -> None:
    artifact_path = tmp_path / "toy.xml"
    artifact_path.write_text("<sbml></sbml>", encoding="utf-8")
    request = ExternalModelGemQaRequest(
        model_id="toy_model",
        artifact_path=str(artifact_path),
        artifact_type="SBML",
        run_memote=True,
    )

    outputs = run_external_model_gem_qa(
        (request,),
        tmp_path / "qa",
        cobra_module=_FakeCobra(_clean_model),
        memote_module=_FakeMemote(),
    )
    result = outputs.results[0]

    assert result.qa_status == "passed_basic"
    assert result.memote_available is True
    assert result.memote_status == "scored"
    assert result.memote_score == 0.82
    assert result.annotation_score == 0.7
    assert result.stoichiometric_consistency_status == "passed"
    assert result.recommendation_tier_effect == "none"


def test_gem_qa_requests_can_be_loaded_from_artifact_cache(tmp_path) -> None:
    artifact_outputs = cache_external_model_artifacts(
        (
            ExternalModelArtifactRequest(
                model_id="toy_model",
                artifact_url="https://example.test/toy.xml",
                artifact_type="SBML",
                filename="toy.xml",
            ),
        ),
        tmp_path / "artifacts",
        fetcher=lambda _url, _timeout: b"<sbml></sbml>",
    )

    requests = gem_qa_requests_from_artifact_cache(artifact_outputs.manifest_path, run_memote=True)

    assert len(requests) == 1
    assert requests[0].model_id == "toy_model"
    assert requests[0].run_memote is True
    assert Path(requests[0].artifact_path).exists()


class _FakeCobra:
    io = None
    flux_analysis = None

    def __init__(self, model_factory) -> None:
        self.io = _FakeCobraIo(model_factory)
        self.flux_analysis = _FakeFluxAnalysis()


class _FakeCobraIo:
    def __init__(self, model_factory) -> None:
        self._model_factory = model_factory

    def read_sbml_model(self, path: str) -> object:
        assert Path(path).exists()
        return self._model_factory()


class _FakeFluxAnalysis:
    @staticmethod
    def find_blocked_reactions(model: object) -> tuple[str, ...]:
        return tuple(
            reaction.id
            for reaction in model.reactions
            if reaction.lower_bound == 0.0 and reaction.upper_bound == 0.0
        )


class _FakeMemote:
    @staticmethod
    def test_model(_model: object, *, results: bool) -> tuple[int, dict[str, object]]:
        assert results is True
        return 0, {"raw": "result"}

    @staticmethod
    def snapshot_report(_result: object, *, html: bool) -> str:
        assert html is False
        return json.dumps(
            {
                "score": {
                    "total_score": 0.82,
                    "sections": [{"section": "annotation", "score": 0.7}],
                },
                "tests": {"test_stoichiometric_consistency": {"result": "passed"}},
            }
        )


def _review_model() -> object:
    r1 = _FakeReaction("R1", "gene1", 0.0, 1000.0)
    r2 = _FakeReaction("R2", "", 0.0, 0.0)
    r3 = _FakeReaction("R3", "gene2", -1000.0, 1000.0)
    r4 = _FakeReaction("R4", "", 0.0, 1000.0, objective=1.0)
    return _FakeModel(
        reactions=(r1, r2, r3, r4),
        metabolites=(
            _FakeMetabolite((r1, r3)),
            _FakeMetabolite((r2,)),
            _FakeMetabolite(()),
        ),
        genes=("gene1", "gene2"),
    )


def _clean_model() -> object:
    r1 = _FakeReaction("R1", "gene1", 0.0, 1000.0)
    r2 = _FakeReaction("R2", "gene2", 0.0, 1000.0, objective=1.0)
    return _FakeModel(
        reactions=(r1, r2),
        metabolites=(_FakeMetabolite((r1, r2)), _FakeMetabolite((r1, r2))),
        genes=("gene1", "gene2"),
    )


class _FakeModel:
    def __init__(self, reactions: tuple, metabolites: tuple, genes: tuple[str, ...]) -> None:
        self.reactions = reactions
        self.metabolites = metabolites
        self.genes = genes


class _FakeReaction:
    def __init__(
        self,
        reaction_id: str,
        gene_reaction_rule: str,
        lower_bound: float,
        upper_bound: float,
        *,
        objective: float = 0.0,
    ) -> None:
        self.id = reaction_id
        self.gene_reaction_rule = gene_reaction_rule
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.objective_coefficient = objective


class _FakeMetabolite:
    def __init__(self, reactions: tuple) -> None:
        self.reactions = reactions
