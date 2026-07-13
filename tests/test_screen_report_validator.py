from __future__ import annotations

from app.services.screen_report_validator import validate_screen_report_json


def _fact_pack() -> dict[str, object]:
    return {
        "source_runs": [{"run_name": "fixture_run"}],
        "targets": {"hLF": {}, "OPN": {}},
        "evidence_items": [
            {
                "evidence_id": "hLF-KO-0001",
                "target_id": "hLF",
                "gene_id": "PAS_HLF_KO",
                "reaction_id": "R1_no_1_fwd",
                "source_run": "fixture_run",
                "intervention_type": "KO",
                "oe_reaction_proxy": False,
                "recommendation_tier": "strong_model_candidate",
                "numeric_fields": {"secretion_ratio_vs_wildtype": 1.2},
            },
            {
                "evidence_id": "OPN-OE-0001",
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "gene_id": "PAS_OPN_OE",
                "canonical_gene_id": "PAS_OPN_OE",
                "reaction_id": "sec_SPC_complex_formation",
                "source_run": "overnight_OPN_full",
                "intervention_type": "OE",
                "oe_reaction_proxy": True,
                "recommendation_tier": "promising_but_oe_proxy",
                "numeric_fields": {"secretion_ratio_vs_wildtype": 1.3},
            },
        ],
    }


def _report(evidence_id: str = "hLF-KO-0001", claim: str = "PAS_HLF_KO 使分泌相对提升") -> dict[str, object]:
    return {
        "schema_version": 1,
        "targets": {
            "hLF": {
                "executive_summary": "summary",
                "recommended_ko": [
                    {"evidence_id": evidence_id, "claim": claim, "rationale": "来自 fact pack", "risk": "需复核", "next_step": "小试验证"}
                ],
                "recommended_oe": [],
                "manual_review": [],
                "not_recommended_or_risky": [],
                "evidence_boundaries": [],
            },
            "OPN": {"executive_summary": "", "recommended_ko": [], "recommended_oe": [], "manual_review": [], "not_recommended_or_risky": [], "evidence_boundaries": []},
        },
        "global_warnings": [],
    }


def test_validator_passes_traceable_report() -> None:
    result = validate_screen_report_json(_fact_pack(), _report())

    assert result["verdict"] == "pass"


def test_validator_requires_schema_version() -> None:
    report = _report()
    del report["schema_version"]

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "schema_version" for issue in result["blocking_issues"])


def test_validator_blocks_string_schema_version() -> None:
    report = _report()
    report["schema_version"] = "1"

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "schema_version" for issue in result["blocking_issues"])


def test_validator_blocks_unsupported_integer_schema_version() -> None:
    report = _report()
    report["schema_version"] = 2

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "schema_version" for issue in result["blocking_issues"])


def test_validator_accepts_supported_schema_version() -> None:
    report = _report()
    report["schema_version"] = 1

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "pass"
    assert not any(issue["type"] == "schema_version" for issue in result["blocking_issues"])


def test_validator_requires_both_target_sections_and_complete_schema() -> None:
    report = _report()
    del report["targets"]["OPN"]
    del report["targets"]["hLF"]["evidence_boundaries"]

    result = validate_screen_report_json(_fact_pack(), report)

    issue_types = {issue["type"] for issue in result["blocking_issues"]}
    assert result["verdict"] == "fail"
    assert "missing_target" in issue_types
    assert "missing_field" in issue_types


def test_validator_requires_string_list_global_warnings() -> None:
    report = _report()
    report["global_warnings"] = "not-a-list"

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue.get("location") == "global_warnings" for issue in result["blocking_issues"])


def test_validator_blocks_fabricated_evidence_id_and_gene_id() -> None:
    result = validate_screen_report_json(_fact_pack(), _report("hLF-KO-9999", "PAS_FAKE_GENE 很好"))

    assert result["verdict"] == "fail"
    issue_types = {issue["type"] for issue in result["blocking_issues"]}
    assert "unsupported_evidence_id" in issue_types
    assert "unsupported_gene_id" in issue_types


def test_validator_blocks_mgl_and_oe_proxy_as_gene_level_oe() -> None:
    report = _report()
    report["targets"]["OPN"]["recommended_oe"] = [
        {
            "evidence_id": "OPN-OE-0001",
            "claim": "这是 gene-level OE，可带来 10 mg/L 绝对产量提升",
            "rationale": "来自 fact pack",
            "risk": "低",
            "next_step": "直接做",
        }
    ]

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    issue_types = {issue["type"] for issue in result["blocking_issues"]}
    assert "forbidden_claim" in issue_types
    assert "misleading_boundary" in issue_types


def test_validator_blocks_tokens_from_the_wrong_cited_evidence_item() -> None:
    report = _report(
        "hLF-KO-0001",
        "PAS_OPN_OE from overnight_OPN_full shows secretion ratio 1.300 for hLF.",
    )

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    issue_types = {issue["type"] for issue in result["blocking_issues"]}
    assert "evidence_token_mismatch" in issue_types


def test_validator_blocks_ko_evidence_in_oe_recommendations() -> None:
    report = _report()
    report["targets"]["hLF"]["recommended_oe"] = report["targets"]["hLF"].pop("recommended_ko")

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "intervention_mismatch" for issue in result["blocking_issues"])


def test_validator_blocks_known_gene_borrowed_by_evidence_without_gene() -> None:
    fact_pack = _fact_pack()
    fact_pack["evidence_items"].append(
        {
            "evidence_id": "hLF-KO-0002",
            "target_id": "hLF",
            "intervention_type": "KO",
            "recommendation_tier": "manual_review_required",
            "numeric_fields": {},
        }
    )
    report = _report("hLF-KO-0002", "PAS_OPN_OE 可用于 hLF")
    report["targets"]["hLF"]["recommended_ko"] = []
    report["targets"]["hLF"]["manual_review"] = [
        {"evidence_id": "hLF-KO-0002", "claim": "PAS_OPN_OE 可用于 hLF", "rationale": "待核对", "risk": "未知", "next_step": "人工复核"}
    ]

    result = validate_screen_report_json(fact_pack, report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "evidence_token_mismatch" for issue in result["blocking_issues"])


def test_validator_blocks_fabricated_integer_or_percentage() -> None:
    report = _report(claim="PAS_HLF_KO 预计提升 50% 并进入第 2 轮")

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "unsupported_numeric_value" for issue in result["blocking_issues"])


def test_validator_blocks_target_mixup_in_executive_summary() -> None:
    report = _report()
    report["targets"]["hLF"]["executive_summary"] = "PAS_OPN_OE 是 hLF 的主要候选"

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "target_token_mismatch" for issue in result["blocking_issues"])


def test_validator_blocks_fabricated_number_in_executive_summary() -> None:
    report = _report()
    report["targets"]["hLF"]["executive_summary"] = "预计总体提升 50%"

    result = validate_screen_report_json(_fact_pack(), report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "unsupported_numeric_value" for issue in result["blocking_issues"])


def test_validator_blocks_empty_bucket_when_fact_pack_has_useful_candidates() -> None:
    fact_pack = _fact_pack()
    fact_pack["targets"]["hLF"] = {
        "useful_ko_candidates": [{"evidence_id": "hLF-KO-0001"}],
        "useful_oe_candidates": [],
        "growth_risk_candidates": [],
        "manual_review_candidates": [],
    }
    report = _report()
    report["targets"]["hLF"]["recommended_ko"] = []

    result = validate_screen_report_json(fact_pack, report)

    assert result["verdict"] == "fail"
    assert any(issue["type"] == "omission" for issue in result["blocking_issues"])


def test_validator_does_not_require_growth_risk_candidate_in_recommended_bucket() -> None:
    fact_pack = _fact_pack()
    fact_pack["targets"]["hLF"] = {
        "useful_ko_candidates": [
            {"evidence_id": "hLF-KO-0001", "recommendation_tier": "not_recommended_growth_risk"}
        ],
        "useful_oe_candidates": [],
        "growth_risk_candidates": [{"evidence_id": "hLF-KO-0001"}],
        "manual_review_candidates": [],
    }
    report = _report()
    report["targets"]["hLF"]["recommended_ko"] = []
    report["targets"]["hLF"]["not_recommended_or_risky"] = [
        {
            "evidence_id": "hLF-KO-0001",
            "claim": "PAS_HLF_KO 存在生长风险",
            "rationale": "来自 fact pack",
            "risk": "生长风险",
            "next_step": "不进入推荐列表",
        }
    ]

    result = validate_screen_report_json(fact_pack, report)

    assert not any(
        issue["type"] == "omission" and issue.get("location") == "hLF.recommended_ko"
        for issue in result["blocking_issues"]
    )
