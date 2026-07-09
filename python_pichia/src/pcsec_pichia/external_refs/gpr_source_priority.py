from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from pcsec_pichia.external_refs.gpr_candidates import CONFLICTING_GPR_SOURCES, MODEL_GPR_CONFIRMED
from pcsec_pichia.external_refs.schema import ExternalGprCandidateEvidence


GPR_SOURCE_PRIORITY_FILENAME = "gpr_source_priority.json"
GPR_SOURCE_CONFLICTS_FILENAME = "gpr_source_conflicts.jsonl"
GPR_SOURCE_PRIORITY_REPORT_FILENAME = "gpr_source_priority_report.md"


@dataclass(frozen=True)
class GprSourcePriorityRecord:
    candidate_cache_key: str
    source_database: str
    external_model_id: str
    external_reaction_id: str
    mapped_pichia_reaction_id: str | None
    external_gene_rule: str | None
    priority_rank: int
    priority_tier: str
    conflict_status: str
    manual_review_required: bool
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class GprSourcePriorityOutputs:
    priority_path: Path
    conflicts_path: Path
    report_path: Path
    record_count: int
    conflict_count: int


def rank_external_gpr_sources(
    candidates: Iterable[ExternalGprCandidateEvidence],
) -> tuple[GprSourcePriorityRecord, ...]:
    rows = tuple(_priority_record(candidate) for candidate in candidates)
    conflict_keys = _conflict_keys(rows)
    resolved: list[GprSourcePriorityRecord] = []
    for row in rows:
        key = _reaction_conflict_key(row)
        if key in conflict_keys:
            warnings = tuple(dict.fromkeys((*row.warnings, "conflicting external GPR rules")))
            resolved.append(
                GprSourcePriorityRecord(
                    candidate_cache_key=row.candidate_cache_key,
                    source_database=row.source_database,
                    external_model_id=row.external_model_id,
                    external_reaction_id=row.external_reaction_id,
                    mapped_pichia_reaction_id=row.mapped_pichia_reaction_id,
                    external_gene_rule=row.external_gene_rule,
                    priority_rank=row.priority_rank,
                    priority_tier=row.priority_tier,
                    conflict_status=CONFLICTING_GPR_SOURCES,
                    manual_review_required=True,
                    warnings=warnings,
                )
            )
        else:
            resolved.append(row)
    return tuple(sorted(resolved, key=lambda row: (row.priority_rank, row.source_database, row.external_model_id)))


def write_gpr_source_priority_outputs(
    records: Iterable[GprSourcePriorityRecord],
    output_dir: Path,
) -> GprSourcePriorityOutputs:
    resolved = tuple(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    priority_path = output_dir / GPR_SOURCE_PRIORITY_FILENAME
    conflicts_path = output_dir / GPR_SOURCE_CONFLICTS_FILENAME
    report_path = output_dir / GPR_SOURCE_PRIORITY_REPORT_FILENAME
    priority_payload = {
        "record_count": len(resolved),
        "records": [_json_ready(asdict(record)) for record in resolved],
    }
    priority_path.write_text(json.dumps(priority_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    conflicts = tuple(record for record in resolved if record.conflict_status == CONFLICTING_GPR_SOURCES)
    with conflicts_path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in conflicts:
            handle.write(json.dumps(_json_ready(asdict(record)), ensure_ascii=False, sort_keys=True))
            handle.write("\n")
    report_path.write_text(render_gpr_source_priority_report(resolved), encoding="utf-8")
    return GprSourcePriorityOutputs(
        priority_path=priority_path,
        conflicts_path=conflicts_path,
        report_path=report_path,
        record_count=len(resolved),
        conflict_count=len(conflicts),
    )


def render_gpr_source_priority_report(records: tuple[GprSourcePriorityRecord, ...]) -> str:
    lines = [
        "# GPR Source Priority Report",
        "",
        "External source ranks are evidence priorities only; conflicts require manual review and are not merged into the current Pichia GEM.",
        "",
        "## Records",
        "",
        "| rank | tier | source | model | reaction | conflict | warnings |",
        "|---:|---|---|---|---|---|---|",
    ]
    for record in records:
        lines.append(
            f"| {record.priority_rank} | {record.priority_tier} | {record.source_database} | "
            f"{record.external_model_id} | {record.mapped_pichia_reaction_id or record.external_reaction_id} | "
            f"{record.conflict_status} | {'; '.join(record.warnings)} |"
        )
    return "\n".join(lines) + "\n"


def _priority_record(candidate: ExternalGprCandidateEvidence) -> GprSourcePriorityRecord:
    rank, tier, warnings = _priority_for_candidate(candidate)
    manual_review = rank >= 4 or bool(candidate.blocking_reasons or candidate.mapping_warnings)
    return GprSourcePriorityRecord(
        candidate_cache_key=candidate.cache_key,
        source_database=candidate.provenance.source_database,
        external_model_id=candidate.external_model_id,
        external_reaction_id=candidate.external_reaction_id,
        mapped_pichia_reaction_id=candidate.mapped_pichia_reaction_id,
        external_gene_rule=candidate.external_gene_rule,
        priority_rank=rank,
        priority_tier=tier,
        conflict_status="none",
        manual_review_required=manual_review,
        warnings=tuple(dict.fromkeys((*candidate.mapping_warnings, *warnings))),
    )


def _priority_for_candidate(candidate: ExternalGprCandidateEvidence) -> tuple[int, str, tuple[str, ...]]:
    model = _normalize(candidate.external_model_id)
    source = _normalize(candidate.provenance.source_database)
    if candidate.candidate_status == MODEL_GPR_CONFIRMED or model in {"currentpichiagem", "currentmodel"}:
        return 0, "current_model_gpr", ()
    if model in {"ipichia", "ecpichia"}:
        return 1, "pichia_specific_model_gpr", ()
    if model in {"kp10", "kp1.0", "kp.1.0", "iaukm"}:
        return 2, "pichia_literature_model_gpr", ()
    if model in {"yeast8", "yeast9", "yeast8yeast9", "yeastgem"} or source in {"yeastgem", "yeast-gem"}:
        return 3, "homology_supported_yeast_gpr", ("cross_species_mapping_required",)
    if source in {"uniprot", "kegg", "ncbi", "sgd"}:
        return 4, "annotation_only", ("annotation_is_not_model_gpr",)
    if source == "gpruler" or model == "gpruler":
        return 5, "automatic_rule_candidate", ("automatic_rule_requires_manual_review",)
    return 6, "manual_review_required", ("unknown_gpr_source_priority",)


def _conflict_keys(rows: tuple[GprSourcePriorityRecord, ...]) -> set[str]:
    rules_by_reaction: dict[str, set[str]] = {}
    for row in rows:
        key = _reaction_conflict_key(row)
        rules_by_reaction.setdefault(key, set()).add(_normalize_rule(row.external_gene_rule))
    return {
        key
        for key, rules in rules_by_reaction.items()
        if len({rule for rule in rules if rule}) > 1
    }


def _reaction_conflict_key(row: GprSourcePriorityRecord) -> str:
    return _normalize(row.mapped_pichia_reaction_id or row.external_reaction_id)


def _normalize_rule(value: object) -> str:
    text = str(value or "")
    for sep in ("(", ")", "/", ",", ";", "|", "+"):
        text = text.replace(sep, " ")
    return " ".join(_normalize(part) for part in text.split() if part.lower() not in {"and", "or"})


def _normalize(value: object) -> str:
    return str(value or "").strip().lower().replace(" ", "").replace("_", "").replace("-", "")


def _json_ready(value: object) -> object:
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "GPR_SOURCE_CONFLICTS_FILENAME",
    "GPR_SOURCE_PRIORITY_FILENAME",
    "GPR_SOURCE_PRIORITY_REPORT_FILENAME",
    "GprSourcePriorityOutputs",
    "GprSourcePriorityRecord",
    "rank_external_gpr_sources",
    "render_gpr_source_priority_report",
    "write_gpr_source_priority_outputs",
]
