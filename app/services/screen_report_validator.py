from __future__ import annotations

import re
from typing import Any


FORBIDDEN_CLAIMS = (
    "mg/L",
    "绝对产量",
    "实验已验证",
    "湿实验已验证",
)
OE_PROXY_FORBIDDEN = (
    "gene-level overexpression",
    "gene-level OE",
    "完整基因过表达",
    "基因级过表达",
)


SUPPORTED_SCHEMA_VERSION = 1
TARGET_BUCKETS = (
    "recommended_ko",
    "recommended_oe",
    "manual_review",
    "not_recommended_or_risky",
)
REQUIRED_TARGET_FIELDS = ("executive_summary", *TARGET_BUCKETS, "evidence_boundaries")


def validate_screen_report_json(fact_pack: dict[str, Any], report: dict[str, Any]) -> dict[str, Any]:
    evidence_by_id = {str(item.get("evidence_id")): item for item in fact_pack.get("evidence_items") or []}
    blocking: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []

    if not isinstance(report, dict):
        return {"verdict": "fail", "blocking_issues": [{"type": "schema", "message": "Report is not a JSON object."}], "warnings": []}
    schema_version = report.get("schema_version")
    if not isinstance(schema_version, int) or isinstance(schema_version, bool) or schema_version != SUPPORTED_SCHEMA_VERSION:
        blocking.append({"type": "schema_version", "message": "Report schema_version must be integer 1."})
    targets = report.get("targets")
    if not isinstance(targets, dict):
        blocking.append({"type": "schema", "message": "Missing targets object."})
        targets = {}
    global_warnings = report.get("global_warnings")
    if not _is_string_list(global_warnings):
        blocking.append({"type": "schema", "location": "global_warnings", "message": "global_warnings must be a list of strings."})

    known_gene_ids, known_reaction_ids, known_source_runs, known_numeric_tokens = _known_tokens(fact_pack)
    all_text = _json_text(report)
    for phrase in FORBIDDEN_CLAIMS:
        if phrase.lower() in all_text.lower():
            blocking.append({"type": "forbidden_claim", "message": f"Forbidden unsupported claim: {phrase}"})
    _check_unknown_tokens(all_text, known_gene_ids, r"\b(?:PAS|AT)[A-Za-z0-9_\-]+\b", "gene_id", blocking)
    _check_unknown_tokens(all_text, known_reaction_ids, r"\bsec_[A-Za-z0-9_]+|\b[A-Z0-9]+_no_[0-9]+_(?:fwd|rvs)\b", "reaction_id", blocking)
    for source_run in re.findall(r"\b(?:ui|overnight|phase|pilot|catalog|gene)[A-Za-z0-9_\-]*\b", all_text):
        if source_run in known_source_runs or source_run in {"gene", "catalog"}:
            continue
        if source_run.startswith(("ui_", "overnight_", "phase", "pilot_")):
            blocking.append({"type": "unsupported_source_run", "message": f"Unknown source_run: {source_run}"})
    for numeric_token in re.findall(r"\b\d+\.\d{3,}\b", all_text):
        if numeric_token not in known_numeric_tokens:
            warnings.append({"type": "numeric_traceability", "message": f"Number not found as exact fact-pack numeric token: {numeric_token}"})

    for target_key in ("hLF", "OPN"):
        if target_key not in targets:
            blocking.append({"type": "missing_target", "location": target_key, "message": f"Missing required target section: {target_key}."})
            continue
        target_report = targets.get(target_key)
        if not isinstance(target_report, dict):
            blocking.append({"type": "schema", "location": target_key, "message": "Target section must be an object."})
            continue
        _validate_target_section(target_key, target_report, evidence_by_id, blocking)
        _validate_minimum_bucket_coverage(
            target_key,
            (fact_pack.get("targets") or {}).get(target_key) or {},
            target_report,
            blocking,
        )

    verdict = "fail" if blocking else "pass"
    return {"verdict": verdict, "blocking_issues": blocking, "warnings": warnings}


def _validate_target_section(
    target_key: str,
    section: dict[str, Any],
    evidence_by_id: dict[str, dict[str, Any]],
    blocking: list[dict[str, Any]],
) -> None:
    missing_fields = [field for field in REQUIRED_TARGET_FIELDS if field not in section]
    for field in missing_fields:
        blocking.append({"type": "missing_field", "location": f"{target_key}.{field}", "message": f"Missing required target field: {field}."})
    if not isinstance(section.get("executive_summary"), str):
        blocking.append({"type": "schema", "location": f"{target_key}.executive_summary", "message": "executive_summary must be a string."})
    if not _is_string_list(section.get("evidence_boundaries")):
        blocking.append({"type": "schema", "location": f"{target_key}.evidence_boundaries", "message": "evidence_boundaries must be a list of strings."})
    _validate_target_prose(target_key, section, evidence_by_id, blocking)

    for bucket in TARGET_BUCKETS:
        rows = section.get(bucket)
        if not isinstance(rows, list):
            blocking.append({"type": "schema", "location": f"{target_key}.{bucket}", "message": "Recommendation bucket must be a list."})
            continue
        for index, row in enumerate(rows):
            location = f"{target_key}.{bucket}[{index}]"
            if not isinstance(row, dict):
                blocking.append({"type": "schema", "location": location, "message": "Recommendation row must be an object."})
                continue
            evidence_id = str(row.get("evidence_id") or "")
            evidence = evidence_by_id.get(evidence_id)
            if not evidence:
                blocking.append({"type": "unsupported_evidence_id", "location": location, "message": f"Unknown evidence_id: {evidence_id}"})
                continue
            if _target_key(str(evidence.get("target_id") or "")) != target_key:
                blocking.append(
                    {
                        "type": "target_mixup",
                        "location": location,
                        "message": f"{evidence_id} belongs to {evidence.get('target_id')}, not {target_key}.",
                        "evidence_id": evidence_id,
                    }
                )
            _validate_bucket_semantics(location, bucket, evidence, blocking)
            for required in ("claim", "rationale", "risk", "next_step"):
                if not str(row.get(required) or "").strip():
                    blocking.append({"type": "missing_field", "location": location, "message": f"Missing {required}.", "evidence_id": evidence_id})
            row_text = _row_prose(row)
            _validate_row_tokens_match_evidence(location, row_text, evidence, blocking)
            if evidence.get("oe_reaction_proxy") and any(term.lower() in row_text.lower() for term in OE_PROXY_FORBIDDEN):
                blocking.append(
                    {
                        "type": "misleading_boundary",
                        "location": location,
                        "message": "OE proxy evidence was described as gene-level overexpression.",
                        "evidence_id": evidence_id,
                    }
                )
            if "homology" in str(evidence.get("recommendation_tier") or "").lower() and "实验验证" in row_text:
                blocking.append(
                    {
                        "type": "misleading_boundary",
                        "location": location,
                        "message": "Homology/annotation evidence cannot be described as experimental validation.",
                        "evidence_id": evidence_id,
                    }
                )


def _validate_row_tokens_match_evidence(
    location: str,
    row_text: str,
    evidence: dict[str, Any],
    blocking: list[dict[str, Any]],
) -> None:
    evidence_id = str(evidence.get("evidence_id") or "")
    allowed_gene_ids = {
        str(evidence.get(key))
        for key in ("gene_id", "canonical_gene_id")
        if str(evidence.get(key) or "").strip()
    }
    allowed_reaction_ids = {
        part
        for part in re.split(r"[;,]\s*", str(evidence.get("reaction_id") or ""))
        if part
    }
    allowed_source_runs = {str(evidence.get("source_run"))} if evidence.get("source_run") else set()
    allowed_numeric_tokens = _numeric_tokens_for_evidence(evidence)

    for token in set(re.findall(r"\b(?:PAS|AT)[A-Za-z0-9_\-]+\b", row_text)):
        if token not in allowed_gene_ids:
            blocking.append(
                {
                    "type": "evidence_token_mismatch",
                    "location": location,
                    "message": f"gene_id {token} does not belong to cited evidence_id {evidence_id}.",
                    "evidence_id": evidence_id,
                }
            )
    for token in set(re.findall(r"\bsec_[A-Za-z0-9_]+|\b[A-Z0-9]+_no_[0-9]+_(?:fwd|rvs)\b", row_text)):
        if token not in allowed_reaction_ids:
            blocking.append(
                {
                    "type": "evidence_token_mismatch",
                    "location": location,
                    "message": f"reaction_id {token} does not belong to cited evidence_id {evidence_id}.",
                    "evidence_id": evidence_id,
                }
            )
    for token in set(re.findall(r"\b(?:ui|overnight|phase|pilot|catalog|gene)[A-Za-z0-9_\-]*\b", row_text)):
        if token in {"gene", "catalog"}:
            continue
        if token.startswith(("ui_", "overnight_", "phase", "pilot_")) and token not in allowed_source_runs:
            blocking.append(
                {
                    "type": "evidence_token_mismatch",
                    "location": location,
                    "message": f"source_run {token} does not belong to cited evidence_id {evidence_id}.",
                    "evidence_id": evidence_id,
                }
            )
    for token in _numeric_tokens_in_prose(row_text):
        if token not in allowed_numeric_tokens:
            blocking.append(
                {
                    "type": "unsupported_numeric_value",
                    "location": location,
                    "message": f"numeric value {token} does not belong to cited evidence_id {evidence_id}.",
                    "evidence_id": evidence_id,
                }
            )


def _validate_bucket_semantics(
    location: str,
    bucket: str,
    evidence: dict[str, Any],
    blocking: list[dict[str, Any]],
) -> None:
    evidence_id = str(evidence.get("evidence_id") or "")
    intervention = str(evidence.get("intervention_type") or "").upper()
    tier = str(evidence.get("recommendation_tier") or "").lower()
    if bucket == "recommended_ko" and intervention != "KO":
        blocking.append({"type": "intervention_mismatch", "location": location, "message": "recommended_ko must cite KO evidence.", "evidence_id": evidence_id})
    if bucket == "recommended_oe" and intervention != "OE":
        blocking.append({"type": "intervention_mismatch", "location": location, "message": "recommended_oe must cite OE evidence.", "evidence_id": evidence_id})
    if bucket in {"recommended_ko", "recommended_oe"} and tier in {
        "not_recommended_growth_risk",
        "manual_review_required",
    }:
        blocking.append(
            {
                "type": "recommendation_tier_mismatch",
                "location": location,
                "message": f"Evidence tier {tier} cannot be placed in a recommended bucket.",
                "evidence_id": evidence_id,
            }
        )


def _validate_target_prose(
    target_key: str,
    section: dict[str, Any],
    evidence_by_id: dict[str, dict[str, Any]],
    blocking: list[dict[str, Any]],
) -> None:
    target_evidence = tuple(
        evidence
        for evidence in evidence_by_id.values()
        if _target_key(str(evidence.get("target_id") or "")) == target_key
    )
    allowed_gene_ids = {
        str(evidence.get(key))
        for evidence in target_evidence
        for key in ("gene_id", "canonical_gene_id")
        if str(evidence.get(key) or "").strip()
    }
    allowed_reaction_ids = {
        part
        for evidence in target_evidence
        for part in re.split(r"[;,]\s*", str(evidence.get("reaction_id") or ""))
        if part
    }
    allowed_source_runs = {
        str(evidence.get("source_run"))
        for evidence in target_evidence
        if evidence.get("source_run")
    }
    allowed_numeric_tokens = {
        token
        for evidence in target_evidence
        for token in _numeric_tokens_for_evidence(evidence)
    }
    prose_fields = [("executive_summary", section.get("executive_summary"))]
    boundaries = section.get("evidence_boundaries")
    if not isinstance(boundaries, list):
        boundaries = []
    prose_fields.extend(
        (f"evidence_boundaries[{index}]", value)
        for index, value in enumerate(boundaries)
        if isinstance(value, str)
    )
    for field, value in prose_fields:
        text = str(value or "")
        location = f"{target_key}.{field}"
        _validate_scoped_prose_tokens(
            location,
            text,
            allowed_gene_ids,
            allowed_reaction_ids,
            allowed_source_runs,
            allowed_numeric_tokens,
            blocking,
        )


def _validate_scoped_prose_tokens(
    location: str,
    text: str,
    allowed_gene_ids: set[str],
    allowed_reaction_ids: set[str],
    allowed_source_runs: set[str],
    allowed_numeric_tokens: set[str],
    blocking: list[dict[str, Any]],
) -> None:
    for token in set(re.findall(r"\b(?:PAS|AT)[A-Za-z0-9_\-]+\b", text)):
        if token not in allowed_gene_ids:
            blocking.append({"type": "target_token_mismatch", "location": location, "message": f"gene_id {token} does not belong to this target."})
    for token in set(re.findall(r"\bsec_[A-Za-z0-9_]+|\b[A-Z0-9]+_no_[0-9]+_(?:fwd|rvs)\b", text)):
        if token not in allowed_reaction_ids:
            blocking.append({"type": "target_token_mismatch", "location": location, "message": f"reaction_id {token} does not belong to this target."})
    for token in set(re.findall(r"\b(?:ui|overnight|phase|pilot|catalog|gene)[A-Za-z0-9_\-]*\b", text)):
        if token in {"gene", "catalog"}:
            continue
        if token.startswith(("ui_", "overnight_", "phase", "pilot_")) and token not in allowed_source_runs:
            blocking.append({"type": "target_token_mismatch", "location": location, "message": f"source_run {token} does not belong to this target."})
    for token in _numeric_tokens_in_prose(text):
        if token not in allowed_numeric_tokens:
            blocking.append({"type": "unsupported_numeric_value", "location": location, "message": f"numeric value {token} is not traceable to this target."})


def _validate_minimum_bucket_coverage(
    target_key: str,
    target_fact_pack: dict[str, Any],
    target_report: dict[str, Any],
    blocking: list[dict[str, Any]],
) -> None:
    category_to_bucket = {
        "useful_ko_candidates": "recommended_ko",
        "useful_oe_candidates": "recommended_oe",
        "growth_risk_candidates": "not_recommended_or_risky",
        "manual_review_candidates": "manual_review",
    }
    for category, bucket in category_to_bucket.items():
        source_rows = target_fact_pack.get(category)
        report_rows = target_report.get(bucket)
        if not isinstance(source_rows, list) or not isinstance(report_rows, list):
            continue
        source_ids = {
            str(row.get("evidence_id"))
            for row in source_rows
            if isinstance(row, dict)
            and row.get("evidence_id")
            and _eligible_for_coverage_category(category, row)
        }
        report_ids = {
            str(row.get("evidence_id"))
            for row in report_rows
            if isinstance(row, dict) and row.get("evidence_id")
        }
        if source_ids and source_ids.isdisjoint(report_ids):
            blocking.append(
                {
                    "type": "omission",
                    "location": f"{target_key}.{bucket}",
                    "message": f"Report bucket does not cite any candidate from fact-pack category {category}.",
                }
            )


def _eligible_for_coverage_category(category: str, row: dict[str, Any]) -> bool:
    if category not in {"useful_ko_candidates", "useful_oe_candidates"}:
        return True
    tier = str(row.get("recommendation_tier") or "").lower()
    return tier not in {"not_recommended_growth_risk", "manual_review_required"}


def _row_prose(row: dict[str, Any]) -> str:
    return "\n".join(str(row.get(field) or "") for field in ("claim", "rationale", "risk", "next_step"))


def _numeric_tokens_in_prose(text: str) -> set[str]:
    return set(
        re.findall(
            r"(?<![A-Za-z0-9_])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?:\s*%)?(?![A-Za-z0-9_])",
            text,
        )
    )


def _is_string_list(value: Any) -> bool:
    return isinstance(value, list) and all(isinstance(item, str) for item in value)


def _known_tokens(fact_pack: dict[str, Any]) -> tuple[set[str], set[str], set[str], set[str]]:
    gene_ids: set[str] = set()
    reaction_ids: set[str] = set()
    source_runs: set[str] = set()
    numeric_tokens: set[str] = set()
    for source in fact_pack.get("source_runs") or []:
        if source.get("run_name"):
            source_runs.add(str(source["run_name"]))
    for item in fact_pack.get("evidence_items") or []:
        for key in ("gene_id", "canonical_gene_id"):
            if item.get(key):
                gene_ids.add(str(item[key]))
        if item.get("reaction_id"):
            reaction_ids.update(part for part in re.split(r"[;,]\s*", str(item["reaction_id"])) if part)
        if item.get("source_run"):
            source_runs.add(str(item["source_run"]))
        for value in (item.get("numeric_fields") or {}).values():
            numeric_tokens.add(str(value))
            try:
                numeric_tokens.add(f"{float(value):.6g}")
            except (TypeError, ValueError):
                pass
    return gene_ids, reaction_ids, source_runs, numeric_tokens


def _numeric_tokens_for_evidence(evidence: dict[str, Any]) -> set[str]:
    tokens: set[str] = set()
    for value in (evidence.get("numeric_fields") or {}).values():
        tokens.add(str(value))
        try:
            tokens.add(f"{float(value):.6g}")
        except (TypeError, ValueError):
            pass
    return tokens


def _check_unknown_tokens(
    text: str,
    known: set[str],
    pattern: str,
    token_type: str,
    blocking: list[dict[str, Any]],
) -> None:
    for token in set(re.findall(pattern, text)):
        if token and token not in known:
            blocking.append({"type": f"unsupported_{token_type}", "message": f"Unknown {token_type}: {token}"})


def _json_text(value: Any) -> str:
    import json

    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _target_key(target_id: str) -> str | None:
    upper = target_id.upper()
    if "HLF" in upper:
        return "hLF"
    if "OPN" in upper:
        return "OPN"
    return None


__all__ = ["validate_screen_report_json"]
