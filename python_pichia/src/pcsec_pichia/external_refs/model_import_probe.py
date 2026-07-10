from __future__ import annotations

import csv
import importlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from pcsec_pichia.external_refs.model_artifacts import ARTIFACT_MANIFEST_FILENAME
from pcsec_pichia.external_refs.schema import utc_now_iso


IMPORT_PROBE_MANIFEST_FILENAME = "external_model_import_probe_manifest.json"
IMPORT_PROBE_SUMMARY_TSV_FILENAME = "external_model_import_probe_summary.tsv"
IMPORT_PROBE_REPORT_FILENAME = "external_model_import_probe_report.md"


@dataclass(frozen=True)
class ExternalModelImportProbeRequest:
    model_id: str
    artifact_path: str
    artifact_type: str = ""
    source_page_url: str = ""
    download_status: str = ""
    notes: str = ""
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        if not self.model_id.strip():
            raise ValueError("model_id must be non-empty.")


@dataclass(frozen=True)
class ExternalModelImportProbeResult:
    model_id: str
    artifact_path: str
    artifact_type: str
    backend: str
    backend_available: bool
    import_status: str
    reaction_count: int | None = None
    metabolite_count: int | None = None
    gene_count: int | None = None
    gpr_count: int | None = None
    objective_reaction: str = ""
    libsbml_comparison_status: str = "not_run"
    manual_review_required: bool = True
    source_page_url: str = ""
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return _json_ready(asdict(self))


@dataclass(frozen=True)
class ExternalModelImportProbeOutputs:
    manifest_path: Path
    summary_tsv_path: Path
    report_path: Path
    generated_at: str
    result_count: int
    imported_count: int
    unavailable_count: int
    manual_review_required_count: int
    results: tuple[ExternalModelImportProbeResult, ...]


def cobrapy_import_available() -> bool:
    return _import_cobra() is not None


def load_import_probe_requests_from_artifact_cache(
    path: Path,
) -> tuple[ExternalModelImportProbeRequest, ...]:
    """Read downloaded artifact rows from an artifact-cache manifest or directory."""

    manifest_path = path / ARTIFACT_MANIFEST_FILENAME if path.is_dir() else path
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    results = payload.get("results", ())
    if not isinstance(results, list):
        raise ValueError("artifact manifest results must be a list.")
    requests: list[ExternalModelImportProbeRequest] = []
    for row_number, row in enumerate(results, start=1):
        if not isinstance(row, dict):
            raise ValueError(f"artifact manifest row {row_number} must be an object.")
        requests.append(
            ExternalModelImportProbeRequest(
                model_id=str(row.get("model_id", "")),
                artifact_path=str(row.get("local_path", "")),
                artifact_type=str(row.get("artifact_type", "")),
                source_page_url=str(row.get("source_page_url", "")),
                download_status=str(row.get("download_status", "")),
                warnings=_tuple_from_json(row.get("warnings", ())),
            )
        )
    return tuple(requests)


def probe_cobrapy_model_import(
    request: ExternalModelImportProbeRequest,
    *,
    cobra_module: Any | None = None,
) -> ExternalModelImportProbeResult:
    request.validate()
    cobra = cobra_module if cobra_module is not None else _import_cobra()
    if cobra is None:
        return _unavailable_result(request)

    artifact_path = Path(request.artifact_path) if request.artifact_path else None
    if artifact_path is None or not artifact_path.exists():
        return ExternalModelImportProbeResult(
            model_id=request.model_id,
            artifact_path=request.artifact_path,
            artifact_type=request.artifact_type,
            backend="cobrapy",
            backend_available=True,
            import_status="missing_artifact",
            manual_review_required=True,
            source_page_url=request.source_page_url,
            warnings=tuple(dict.fromkeys((*request.warnings, "artifact_file_missing"))),
        )

    try:
        model = _read_cobrapy_model(cobra, artifact_path, request.artifact_type)
    except ValueError as exc:
        return ExternalModelImportProbeResult(
            model_id=request.model_id,
            artifact_path=str(artifact_path),
            artifact_type=request.artifact_type,
            backend="cobrapy",
            backend_available=True,
            import_status="unsupported_artifact_type",
            manual_review_required=True,
            source_page_url=request.source_page_url,
            warnings=tuple(dict.fromkeys((*request.warnings, str(exc)))),
        )
    except Exception as exc:  # pragma: no cover - parser failures depend on optional COBRApy stack.
        return ExternalModelImportProbeResult(
            model_id=request.model_id,
            artifact_path=str(artifact_path),
            artifact_type=request.artifact_type,
            backend="cobrapy",
            backend_available=True,
            import_status="import_failed",
            manual_review_required=True,
            source_page_url=request.source_page_url,
            warnings=tuple(dict.fromkeys((*request.warnings, f"{type(exc).__name__}: {exc}"))),
        )

    return ExternalModelImportProbeResult(
        model_id=request.model_id,
        artifact_path=str(artifact_path),
        artifact_type=request.artifact_type,
        backend="cobrapy",
        backend_available=True,
        import_status="imported",
        reaction_count=len(getattr(model, "reactions", ()) or ()),
        metabolite_count=len(getattr(model, "metabolites", ()) or ()),
        gene_count=len(getattr(model, "genes", ()) or ()),
        gpr_count=_gpr_count(model),
        objective_reaction=_objective_reaction(model),
        libsbml_comparison_status="not_run",
        manual_review_required=False,
        source_page_url=request.source_page_url,
        warnings=request.warnings,
    )


def probe_external_model_imports(
    requests: Iterable[ExternalModelImportProbeRequest],
    output_dir: Path,
    *,
    cobra_module: Any | None = None,
) -> ExternalModelImportProbeOutputs:
    resolved_requests = tuple(requests)
    results = tuple(
        probe_cobrapy_model_import(request, cobra_module=cobra_module)
        for request in resolved_requests
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    generated_at = utc_now_iso()
    payload = _outputs_payload(generated_at, results)
    manifest_path = output_dir / IMPORT_PROBE_MANIFEST_FILENAME
    summary_tsv_path = output_dir / IMPORT_PROBE_SUMMARY_TSV_FILENAME
    report_path = output_dir / IMPORT_PROBE_REPORT_FILENAME
    manifest_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_summary_tsv(results, summary_tsv_path)
    report_path.write_text(render_external_model_import_probe_report(results, generated_at=generated_at), encoding="utf-8")
    return ExternalModelImportProbeOutputs(
        manifest_path=manifest_path,
        summary_tsv_path=summary_tsv_path,
        report_path=report_path,
        generated_at=generated_at,
        result_count=len(results),
        imported_count=sum(1 for result in results if result.import_status == "imported"),
        unavailable_count=sum(1 for result in results if result.import_status == "unavailable"),
        manual_review_required_count=sum(1 for result in results if result.manual_review_required),
        results=results,
    )


def render_external_model_import_probe_report(
    results: Iterable[ExternalModelImportProbeResult],
    *,
    generated_at: str | None = None,
) -> str:
    resolved = tuple(results)
    counts: dict[str, int] = {}
    for result in resolved:
        counts[result.import_status] = counts.get(result.import_status, 0) + 1
    lines = [
        "# External GEM COBRApy Import Probe",
        "",
        "This probe records external GEM import diagnostics only; it does not alter the pcSec model, recommendation tiers, or default solver path.",
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
            "| model_id | status | backend_available | reactions | metabolites | genes | gpr_rules | manual_review | warnings |",
            "|---|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for result in resolved:
        lines.append(
            f"| {result.model_id} | {result.import_status} | {result.backend_available} | "
            f"{_fmt(result.reaction_count)} | {_fmt(result.metabolite_count)} | {_fmt(result.gene_count)} | "
            f"{_fmt(result.gpr_count)} | {result.manual_review_required} | {'; '.join(result.warnings)} |"
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- COBRApy import success is external model metadata only.",
            "- Import diagnostics do not imply model GPR executability in the current Pichia GEM.",
            "- This report does not claim absolute secretion titer, absolute secretion yield, or experimental success rate.",
        ]
    )
    return "\n".join(lines) + "\n"


def _import_cobra() -> Any | None:
    try:
        return importlib.import_module("cobra")
    except Exception:
        return None


def _unavailable_result(request: ExternalModelImportProbeRequest) -> ExternalModelImportProbeResult:
    return ExternalModelImportProbeResult(
        model_id=request.model_id,
        artifact_path=request.artifact_path,
        artifact_type=request.artifact_type,
        backend="cobrapy",
        backend_available=False,
        import_status="unavailable",
        manual_review_required=True,
        source_page_url=request.source_page_url,
        warnings=tuple(dict.fromkeys((*request.warnings, "cobrapy_unavailable"))),
    )


def _read_cobrapy_model(cobra: Any, path: Path, artifact_type: str) -> Any:
    io = getattr(cobra, "io", None)
    if io is None:
        raise ValueError("cobra.io is unavailable.")
    artifact_hint = artifact_type.lower()
    suffix = path.suffix.lower()
    if "sbml" in artifact_hint or suffix in {".xml", ".sbml"}:
        return io.read_sbml_model(str(path))
    if "json" in artifact_hint or suffix == ".json":
        return io.load_json_model(str(path))
    if "matlab" in artifact_hint or suffix == ".mat":
        return io.load_matlab_model(str(path))
    raise ValueError(f"unsupported artifact type for COBRApy import probe: {artifact_type or suffix or 'unknown'}")


def _gpr_count(model: Any) -> int:
    count = 0
    for reaction in getattr(model, "reactions", ()) or ():
        rule = str(getattr(reaction, "gene_reaction_rule", "") or "").strip()
        if rule:
            count += 1
    return count


def _objective_reaction(model: Any) -> str:
    for reaction in getattr(model, "reactions", ()) or ():
        coefficient = getattr(reaction, "objective_coefficient", 0.0)
        try:
            if float(coefficient):
                return str(getattr(reaction, "id", reaction))
        except (TypeError, ValueError):
            continue
    objective = getattr(model, "objective", None)
    return "" if objective is None else str(objective)


def _outputs_payload(
    generated_at: str,
    results: tuple[ExternalModelImportProbeResult, ...],
) -> dict[str, Any]:
    return {
        "generated_at": generated_at,
        "result_count": len(results),
        "imported_count": sum(1 for result in results if result.import_status == "imported"),
        "unavailable_count": sum(1 for result in results if result.import_status == "unavailable"),
        "manual_review_required_count": sum(1 for result in results if result.manual_review_required),
        "results": [result.to_dict() for result in results],
    }


def _write_summary_tsv(results: tuple[ExternalModelImportProbeResult, ...], path: Path) -> None:
    fieldnames = tuple(ExternalModelImportProbeResult.__dataclass_fields__)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({key: _tsv_value(value) for key, value in result.to_dict().items()})


def _tuple_from_json(value: object) -> tuple[str, ...]:
    if isinstance(value, tuple):
        return tuple(str(item) for item in value)
    if isinstance(value, list):
        return tuple(str(item) for item in value)
    if value in (None, ""):
        return ()
    raise ValueError("warnings must be encoded as a JSON list.")


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
    return "" if value is None else str(value)


__all__ = [
    "IMPORT_PROBE_MANIFEST_FILENAME",
    "IMPORT_PROBE_REPORT_FILENAME",
    "IMPORT_PROBE_SUMMARY_TSV_FILENAME",
    "ExternalModelImportProbeOutputs",
    "ExternalModelImportProbeRequest",
    "ExternalModelImportProbeResult",
    "cobrapy_import_available",
    "load_import_probe_requests_from_artifact_cache",
    "probe_cobrapy_model_import",
    "probe_external_model_imports",
    "render_external_model_import_probe_report",
]
