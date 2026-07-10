from __future__ import annotations

import csv
import importlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

from pcsec_pichia.external_refs.model_import_probe import (
    ExternalModelImportProbeRequest,
    load_import_probe_requests_from_artifact_cache,
    probe_cobrapy_model_import,
)
from pcsec_pichia.external_refs.schema import utc_now_iso


GEM_QA_MANIFEST_FILENAME = "external_model_gem_qa_manifest.json"
GEM_QA_SUMMARY_TSV_FILENAME = "external_model_gem_qa_summary.tsv"
GEM_QA_REPORT_FILENAME = "external_model_gem_qa_report.md"


@dataclass(frozen=True)
class ExternalModelGemQaRequest:
    model_id: str
    artifact_path: str
    artifact_type: str = ""
    source_page_url: str = ""
    run_memote: bool = False
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        if not self.model_id.strip():
            raise ValueError("model_id must be non-empty.")


@dataclass(frozen=True)
class ExternalModelGemQaResult:
    model_id: str
    artifact_path: str
    qa_backend: str
    backend_available: bool
    qa_status: str
    import_status: str
    memote_available: bool
    memote_status: str
    memote_score: float | None
    stoichiometric_consistency_status: str
    annotation_score: float | None
    gpr_coverage: float | None
    blocked_reaction_count: int | None
    dead_end_metabolite_count: int | None
    manual_review_reasons: tuple[str, ...]
    recommendation_tier_effect: str = "none"
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return _json_ready(asdict(self))


@dataclass(frozen=True)
class ExternalModelGemQaOutputs:
    manifest_path: Path
    summary_tsv_path: Path
    report_path: Path
    generated_at: str
    result_count: int
    passed_basic_count: int
    review_required_count: int
    unavailable_count: int
    results: tuple[ExternalModelGemQaResult, ...]


def gem_qa_requests_from_artifact_cache(
    path: Path,
    *,
    run_memote: bool = False,
) -> tuple[ExternalModelGemQaRequest, ...]:
    return tuple(
        ExternalModelGemQaRequest(
            model_id=request.model_id,
            artifact_path=request.artifact_path,
            artifact_type=request.artifact_type,
            source_page_url=request.source_page_url,
            run_memote=run_memote,
            warnings=request.warnings,
        )
        for request in load_import_probe_requests_from_artifact_cache(path)
    )


def run_external_model_gem_qa(
    requests: Iterable[ExternalModelGemQaRequest],
    output_dir: Path,
    *,
    cobra_module: Any | None = None,
    memote_module: Any | None = None,
) -> ExternalModelGemQaOutputs:
    resolved_requests = tuple(requests)
    results = tuple(
        run_external_model_gem_qa_one(
            request,
            cobra_module=cobra_module,
            memote_module=memote_module,
        )
        for request in resolved_requests
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    generated_at = utc_now_iso()
    manifest_path = output_dir / GEM_QA_MANIFEST_FILENAME
    summary_tsv_path = output_dir / GEM_QA_SUMMARY_TSV_FILENAME
    report_path = output_dir / GEM_QA_REPORT_FILENAME
    manifest_path.write_text(
        json.dumps(_outputs_payload(generated_at, results), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_summary_tsv(results, summary_tsv_path)
    report_path.write_text(render_external_model_gem_qa_report(results, generated_at=generated_at), encoding="utf-8")
    return ExternalModelGemQaOutputs(
        manifest_path=manifest_path,
        summary_tsv_path=summary_tsv_path,
        report_path=report_path,
        generated_at=generated_at,
        result_count=len(results),
        passed_basic_count=sum(1 for result in results if result.qa_status == "passed_basic"),
        review_required_count=sum(1 for result in results if result.qa_status == "review_required"),
        unavailable_count=sum(1 for result in results if result.qa_status == "unavailable"),
        results=results,
    )


def run_external_model_gem_qa_one(
    request: ExternalModelGemQaRequest,
    *,
    cobra_module: Any | None = None,
    memote_module: Any | None = None,
) -> ExternalModelGemQaResult:
    request.validate()
    probe_request = ExternalModelImportProbeRequest(
        model_id=request.model_id,
        artifact_path=request.artifact_path,
        artifact_type=request.artifact_type,
        source_page_url=request.source_page_url,
        warnings=request.warnings,
    )
    probe = probe_cobrapy_model_import(probe_request, cobra_module=cobra_module)
    if probe.import_status != "imported":
        return ExternalModelGemQaResult(
            model_id=request.model_id,
            artifact_path=request.artifact_path,
            qa_backend="cobrapy-basic",
            backend_available=probe.backend_available,
            qa_status="unavailable" if not probe.backend_available else "review_required",
            import_status=probe.import_status,
            memote_available=False,
            memote_status="not_run",
            memote_score=None,
            stoichiometric_consistency_status="not_run",
            annotation_score=None,
            gpr_coverage=None,
            blocked_reaction_count=None,
            dead_end_metabolite_count=None,
            manual_review_reasons=(probe.import_status,),
            warnings=probe.warnings,
        )

    cobra = cobra_module if cobra_module is not None else _import_cobra()
    model = _load_model_for_basic_qa(cobra, request)
    basic = _basic_qa_metrics(model, probe.gpr_count or 0, probe.reaction_count or 0)
    memote = _memote_summary(model, request.run_memote, memote_module=memote_module)
    warnings = tuple(dict.fromkeys((*request.warnings, *memote["warnings"])))
    manual_review_reasons = _manual_review_reasons(basic)
    qa_status = "passed_basic" if not manual_review_reasons else "review_required"
    return ExternalModelGemQaResult(
        model_id=request.model_id,
        artifact_path=request.artifact_path,
        qa_backend="cobrapy-basic+optional-memote" if request.run_memote else "cobrapy-basic",
        backend_available=True,
        qa_status=qa_status,
        import_status=probe.import_status,
        memote_available=bool(memote["available"]),
        memote_status=str(memote["status"]),
        memote_score=memote["score"],
        stoichiometric_consistency_status=str(memote["stoichiometric_consistency_status"]),
        annotation_score=memote["annotation_score"],
        gpr_coverage=basic["gpr_coverage"],
        blocked_reaction_count=basic["blocked_reaction_count"],
        dead_end_metabolite_count=basic["dead_end_metabolite_count"],
        manual_review_reasons=manual_review_reasons,
        warnings=warnings,
    )


def render_external_model_gem_qa_report(
    results: Iterable[ExternalModelGemQaResult],
    *,
    generated_at: str | None = None,
) -> str:
    resolved = tuple(results)
    counts: dict[str, int] = {}
    for result in resolved:
        counts[result.qa_status] = counts.get(result.qa_status, 0) + 1
    lines = [
        "# External GEM QA Report",
        "",
        "This report records external model QA metadata only; it does not alter pcSec recommendations or the default solver path.",
        "",
        "## Summary",
        "",
        f"- generated_at: {generated_at or utc_now_iso()}",
        f"- result_count: {len(resolved)}",
    ]
    for status, count in sorted(counts.items()):
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Results",
            "",
            "| model_id | qa_status | backend | memote | score | gpr_coverage | blocked_rxns | dead_end_mets | review_reasons |",
            "|---|---|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for result in resolved:
        lines.append(
            f"| {result.model_id} | {result.qa_status} | {result.qa_backend} | {result.memote_status} | "
            f"{_fmt(result.memote_score)} | {_fmt(result.gpr_coverage)} | "
            f"{_fmt(result.blocked_reaction_count)} | {_fmt(result.dead_end_metabolite_count)} | "
            f"{'; '.join(result.manual_review_reasons)} |"
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- GEM QA and MEMOTE results remain external model metadata.",
            "- QA results do not change `recommendation_tier`.",
            "- This report does not claim absolute secretion titer, absolute secretion yield, or experimental success rate.",
        ]
    )
    return "\n".join(lines) + "\n"


def _import_cobra() -> Any | None:
    try:
        return importlib.import_module("cobra")
    except Exception:
        return None


def _import_memote() -> Any | None:
    try:
        return importlib.import_module("memote")
    except Exception:
        return None


def _load_model_for_basic_qa(cobra: Any, request: ExternalModelGemQaRequest) -> Any:
    if cobra is None:
        raise RuntimeError("COBRApy model should already be available after successful import probe.")
    artifact_path = Path(request.artifact_path)
    io = getattr(cobra, "io", None)
    artifact_hint = request.artifact_type.lower()
    suffix = artifact_path.suffix.lower()
    if "sbml" in artifact_hint or suffix in {".xml", ".sbml"}:
        return io.read_sbml_model(str(artifact_path))
    if "json" in artifact_hint or suffix == ".json":
        return io.load_json_model(str(artifact_path))
    if "matlab" in artifact_hint or suffix == ".mat":
        return io.load_matlab_model(str(artifact_path))
    raise ValueError(f"unsupported artifact type: {request.artifact_type or suffix or 'unknown'}")


def _basic_qa_metrics(model: Any, gpr_count: int, reaction_count: int) -> dict[str, Any]:
    reactions = tuple(getattr(model, "reactions", ()) or ())
    metabolites = tuple(getattr(model, "metabolites", ()) or ())
    resolved_reaction_count = reaction_count or len(reactions)
    blocked_reactions = tuple(
        reaction
        for reaction in reactions
        if _float_attr(reaction, "lower_bound") == 0.0 and _float_attr(reaction, "upper_bound") == 0.0
    )
    dead_end_metabolites = tuple(
        metabolite
        for metabolite in metabolites
        if len(getattr(metabolite, "reactions", ()) or ()) <= 1
    )
    return {
        "gpr_coverage": None if resolved_reaction_count <= 0 else float(gpr_count) / float(resolved_reaction_count),
        "blocked_reaction_count": len(blocked_reactions),
        "dead_end_metabolite_count": len(dead_end_metabolites),
    }


def _memote_summary(
    model: Any,
    run_memote: bool,
    *,
    memote_module: Any | None,
) -> dict[str, Any]:
    if not run_memote:
        return _memote_payload(False, "not_requested", None, "not_run", None, ())
    memote = memote_module if memote_module is not None else _import_memote()
    if memote is None:
        return _memote_payload(False, "unavailable", None, "not_run", None, ("memote_unavailable",))
    score_model = getattr(memote, "score_model", None)
    if score_model is None:
        return _memote_payload(True, "available_not_run", None, "not_run", None, ("memote_score_model_api_unavailable",))
    summary = score_model(model)
    return _memote_payload(
        True,
        "scored",
        _optional_float(_mapping_get(summary, "memote_score", "score", "total_score")),
        str(_mapping_get(summary, "stoichiometric_consistency_status", "consistency_status") or "not_reported"),
        _optional_float(_mapping_get(summary, "annotation_score")),
        (),
    )


def _memote_payload(
    available: bool,
    status: str,
    score: float | None,
    stoichiometric_consistency_status: str,
    annotation_score: float | None,
    warnings: tuple[str, ...],
) -> dict[str, Any]:
    return {
        "available": available,
        "status": status,
        "score": score,
        "stoichiometric_consistency_status": stoichiometric_consistency_status,
        "annotation_score": annotation_score,
        "warnings": warnings,
    }


def _manual_review_reasons(basic: dict[str, Any]) -> tuple[str, ...]:
    reasons: list[str] = []
    gpr_coverage = basic["gpr_coverage"]
    if gpr_coverage is None or gpr_coverage <= 0.0:
        reasons.append("no_gpr_rules_detected")
    if int(basic["blocked_reaction_count"]) > 0:
        reasons.append("blocked_reactions_present")
    if int(basic["dead_end_metabolite_count"]) > 0:
        reasons.append("dead_end_metabolites_present")
    return tuple(reasons)


def _outputs_payload(
    generated_at: str,
    results: tuple[ExternalModelGemQaResult, ...],
) -> dict[str, Any]:
    return {
        "generated_at": generated_at,
        "result_count": len(results),
        "passed_basic_count": sum(1 for result in results if result.qa_status == "passed_basic"),
        "review_required_count": sum(1 for result in results if result.qa_status == "review_required"),
        "unavailable_count": sum(1 for result in results if result.qa_status == "unavailable"),
        "results": [result.to_dict() for result in results],
    }


def _write_summary_tsv(results: tuple[ExternalModelGemQaResult, ...], path: Path) -> None:
    fieldnames = tuple(ExternalModelGemQaResult.__dataclass_fields__)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({key: _tsv_value(value) for key, value in result.to_dict().items()})


def _float_attr(value: Any, attr: str) -> float | None:
    try:
        return float(getattr(value, attr))
    except (TypeError, ValueError):
        return None


def _mapping_get(value: Any, *keys: str) -> Any:
    if not isinstance(value, dict):
        return None
    for key in keys:
        if key in value:
            return value[key]
    return None


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    return float(value)


def _json_ready(value: object) -> object:
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


def _tsv_value(value: object) -> str:
    if isinstance(value, list):
        return "|".join(str(item) for item in value)
    return "" if value is None else str(value)


def _fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


__all__ = [
    "GEM_QA_MANIFEST_FILENAME",
    "GEM_QA_REPORT_FILENAME",
    "GEM_QA_SUMMARY_TSV_FILENAME",
    "ExternalModelGemQaOutputs",
    "ExternalModelGemQaRequest",
    "ExternalModelGemQaResult",
    "gem_qa_requests_from_artifact_cache",
    "render_external_model_gem_qa_report",
    "run_external_model_gem_qa",
    "run_external_model_gem_qa_one",
]
