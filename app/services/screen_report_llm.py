from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Protocol


WRITER_SYSTEM_PROMPT = """你是毕赤酵母目标蛋白分泌工程研发报告 writer。
你只能使用用户给出的 fact pack。每条建议必须引用 fact pack 里存在的 evidence_id。
不得新增 fact pack 外的 gene_id、reaction_id、run 或数值。
不得给出 mg/L、绝对产量或实验成功率预测。
不得把 homology_auxiliary/annotation 说成实验验证。
不得把 OE reaction-level proxy 说成完整 gene-level OE。
输出必须是 JSON，schema_version=1，并按 hLF 和 OPN 分区。
"""

JUDGE_SYSTEM_PROMPT = """你是研发报告 judge。你只检查 writer 报告是否忠于 fact pack 和程序 validator 结果。
重点查 unsupported claim、omission、target_mixup、misleading_boundary、readability、prioritization。
输出必须是 JSON，verdict 只能是 pass 或 fail。
"""

_SCREEN_REPORT_ENV_KEYS = {
    "OPENAI_API_KEY",
    "OPENAI_BASE_URL",
    "OPENAI_API_BASE",
    "SCREEN_REPORT_LLM_PROVIDER",
    "REPORT_LLM_PROVIDER",
    "SCREEN_REPORT_LLM_MODEL",
    "SCREEN_REPORT_LLM_BASE_URL",
}


class JsonLlmClient(Protocol):
    def complete_json(self, *, system_prompt: str, payload: dict[str, Any]) -> dict[str, Any]: ...


class OpenAIJsonLlmClient:
    def __init__(self, model: str | None = None, base_url: str | None = None) -> None:
        _load_project_env()
        self.model = model or os.environ.get("SCREEN_REPORT_LLM_MODEL", "gpt-4o-mini")
        self.base_url = (
            base_url
            or os.environ.get("SCREEN_REPORT_LLM_BASE_URL")
            or os.environ.get("OPENAI_BASE_URL")
            or os.environ.get("OPENAI_API_BASE")
            or ""
        )

    def complete_json(self, *, system_prompt: str, payload: dict[str, Any]) -> dict[str, Any]:
        _load_project_env()
        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY 未设置；无法生成 LLM 研发建议报告，但 fact pack 仍可生成和下载。")
        from openai import OpenAI

        client_kwargs: dict[str, str] = {"api_key": api_key}
        if self.base_url:
            client_kwargs["base_url"] = self.base_url
        client = OpenAI(**client_kwargs)
        response = client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": json.dumps(payload, ensure_ascii=False, indent=2)},
            ],
            response_format={"type": "json_object"},
        )
        content = response.choices[0].message.content or "{}"
        parsed = json.loads(content)
        if not isinstance(parsed, dict):
            raise ValueError("LLM returned JSON that is not an object.")
        return parsed


def get_default_screen_report_llm_client() -> JsonLlmClient:
    _load_project_env()
    provider = os.environ.get("SCREEN_REPORT_LLM_PROVIDER", os.environ.get("REPORT_LLM_PROVIDER", "openai")).strip().lower()
    if provider == "openai":
        return OpenAIJsonLlmClient()
    raise ValueError(f"Unknown SCREEN_REPORT_LLM_PROVIDER: {provider!r}. Supported: 'openai'.")


def write_screen_report_draft(
    client: JsonLlmClient,
    fact_pack: dict[str, Any],
    *,
    feedback: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    payload = {
        "fact_pack": _fact_pack_for_writer(fact_pack),
        "required_output_schema": _writer_schema_hint(),
        "feedback_to_fix": feedback or [],
    }
    return client.complete_json(system_prompt=WRITER_SYSTEM_PROMPT, payload=payload)


def judge_screen_report(
    client: JsonLlmClient,
    fact_pack: dict[str, Any],
    report_json: dict[str, Any],
    validator_result: dict[str, Any],
) -> dict[str, Any]:
    payload = {
        "fact_pack": _fact_pack_for_judge(fact_pack, report_json),
        "fact_pack_summary": _fact_pack_summary_for_judge(fact_pack),
        "report_json": report_json,
        "programmatic_validator_result": validator_result,
        "required_output_schema": {
            "verdict": "pass|fail",
            "blocking_issues": [{"type": "unsupported_claim|omission|target_mixup|misleading_boundary|readability|prioritization", "location": "", "message": "", "evidence_id": ""}],
            "required_fixes": [],
        },
    }
    result = client.complete_json(system_prompt=JUDGE_SYSTEM_PROMPT, payload=payload)
    verdict = str(result.get("verdict") or "").lower()
    blocking_issues = result.get("blocking_issues")
    if not isinstance(blocking_issues, list):
        blocking_issues = [{"type": "schema", "message": "Judge blocking_issues must be a list."}]
    else:
        invalid_issue_count = sum(1 for item in blocking_issues if not isinstance(item, dict))
        blocking_issues = [item for item in blocking_issues if isinstance(item, dict)]
        if invalid_issue_count:
            blocking_issues.append(
                {
                    "type": "schema",
                    "message": f"Judge blocking_issues contained {invalid_issue_count} non-object item(s).",
                }
            )
    required_fixes = result.get("required_fixes")
    if not isinstance(required_fixes, list):
        blocking_issues.append({"type": "schema", "message": "Judge required_fixes must be a list."})
        required_fixes = []
    if verdict not in {"pass", "fail"}:
        blocking_issues.append({"type": "schema", "message": "Judge verdict was not pass/fail."})
        verdict = "fail"
    if verdict == "pass" and (blocking_issues or required_fixes):
        if required_fixes and not blocking_issues:
            blocking_issues.append(
                {"type": "schema", "message": "Judge returned pass with non-empty required_fixes."}
            )
        verdict = "fail"
    result["verdict"] = verdict
    result["blocking_issues"] = blocking_issues
    result["required_fixes"] = required_fixes
    return result


def _writer_schema_hint() -> dict[str, Any]:
    empty_target = {
        "executive_summary": "",
        "recommended_ko": [{"evidence_id": "", "claim": "", "rationale": "", "risk": "", "next_step": ""}],
        "recommended_oe": [],
        "manual_review": [],
        "not_recommended_or_risky": [],
        "evidence_boundaries": [],
    }
    return {"schema_version": 1, "targets": {"hLF": empty_target, "OPN": empty_target}, "global_warnings": []}


def _fact_pack_summary_for_judge(fact_pack: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": fact_pack.get("schema_version"),
        "source_runs": fact_pack.get("source_runs"),
        "target_counts": {
            target_key: {
                "candidate_count": target.get("candidate_count"),
                "useful_ko_count": len(target.get("useful_ko_candidates") or []),
                "useful_oe_count": len(target.get("useful_oe_candidates") or []),
                "manual_review_count": len(target.get("manual_review_candidates") or []),
            }
            for target_key, target in (fact_pack.get("targets") or {}).items()
        },
        "evidence_ids": [item.get("evidence_id") for item in fact_pack.get("evidence_items") or []],
    }


def _fact_pack_for_writer(fact_pack: dict[str, Any]) -> dict[str, Any]:
    """Keep LLM input focused on reportable evidence; validator still sees the full fact pack."""

    return {
        "schema_version": fact_pack.get("schema_version"),
        "generated_at": fact_pack.get("generated_at"),
        "source_runs": fact_pack.get("source_runs"),
        "targets": fact_pack.get("targets"),
        "warnings": fact_pack.get("warnings"),
    }


def _fact_pack_for_judge(fact_pack: dict[str, Any], report_json: dict[str, Any]) -> dict[str, Any]:
    cited_ids = _cited_evidence_ids(report_json)
    evidence_by_id = {
        str(item.get("evidence_id")): item
        for item in fact_pack.get("evidence_items") or []
        if item.get("evidence_id")
    }
    return {
        **_fact_pack_for_writer(fact_pack),
        "cited_evidence_items": [
            evidence_by_id[evidence_id]
            for evidence_id in cited_ids
            if evidence_id in evidence_by_id
        ],
    }


def _cited_evidence_ids(report_json: dict[str, Any]) -> list[str]:
    cited: list[str] = []
    targets = report_json.get("targets")
    if not isinstance(targets, dict):
        return cited
    for section in targets.values():
        if not isinstance(section, dict):
            continue
        for bucket in ("recommended_ko", "recommended_oe", "manual_review", "not_recommended_or_risky"):
            rows = section.get(bucket) or []
            if not isinstance(rows, list):
                continue
            for row in rows:
                if isinstance(row, dict) and row.get("evidence_id"):
                    cited.append(str(row["evidence_id"]))
    return list(dict.fromkeys(cited))


_PROJECT_ENV_LOADED = False


def _load_project_env(env_path: Path | None = None) -> None:
    global _PROJECT_ENV_LOADED
    use_cache = env_path is None
    if use_cache and _PROJECT_ENV_LOADED:
        return
    env_path = env_path or (Path(__file__).resolve().parents[2] / ".env")
    if not env_path.exists():
        if use_cache:
            _PROJECT_ENV_LOADED = True
        return
    try:
        lines = env_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        if use_cache:
            _PROJECT_ENV_LOADED = True
        return
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        key = key.strip().lstrip("\ufeff")
        if key not in _SCREEN_REPORT_ENV_KEYS or key in os.environ:
            continue
        os.environ[key] = _clean_env_value(value)
    if use_cache:
        _PROJECT_ENV_LOADED = True


def _clean_env_value(value: str) -> str:
    cleaned = value.strip()
    if len(cleaned) >= 2 and cleaned[0] == cleaned[-1] and cleaned[0] in {'"', "'"}:
        return cleaned[1:-1]
    return cleaned


__all__ = [
    "JsonLlmClient",
    "OpenAIJsonLlmClient",
    "get_default_screen_report_llm_client",
    "judge_screen_report",
    "write_screen_report_draft",
]
