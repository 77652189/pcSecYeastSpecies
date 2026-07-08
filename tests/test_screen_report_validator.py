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
