from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pcsec_pichia.oe_capacity import (
    ConfidenceLevel,
    EvidenceSourceType,
    GPRRole,
    OEExecutionStatus,
)


FIXTURE_PATH = (
    Path(__file__).parent / "fixtures" / "oe_capacity" / "round0_cases.json"
)
REQUIRED_CASE_TYPES = {
    "single_gene",
    "isoenzyme",
    "complex_subunit",
    "mixed",
    "missing_parameter",
    "external_evidence_only",
}
REQUIRED_CASE_FIELDS = {
    "case_id",
    "case_type",
    "gene_id",
    "reaction_id",
    "enzyme_id",
    "gpr_rule",
    "gpr_role",
    "source",
    "execution_status",
    "missing_information",
    "warnings",
}


def _load_cases() -> list[dict[str, Any]]:
    payload = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    return payload["cases"]


def test_round0_fixture_has_unique_ids_all_required_cases_and_contract_fields() -> None:
    cases = _load_cases()
    case_ids = [case["case_id"] for case in cases]

    assert len(case_ids) == len(set(case_ids))
    assert {case["case_type"] for case in cases} == REQUIRED_CASE_TYPES

    valid_roles = {role.value for role in GPRRole}
    valid_statuses = {status.value for status in OEExecutionStatus}
    valid_source_types = {source.value for source in EvidenceSourceType}
    valid_confidences = {confidence.value for confidence in ConfidenceLevel}

    for case in cases:
        assert REQUIRED_CASE_FIELDS <= case.keys()
        assert all(case[field] for field in ("gene_id", "reaction_id", "enzyme_id"))
        assert case["gpr_role"] in valid_roles
        assert case["execution_status"] in valid_statuses
        assert isinstance(case["missing_information"], list)
        assert isinstance(case["warnings"], list)

        source = case["source"]
        assert {"type", "ref", "version", "confidence"} <= source.keys()
        assert source["type"] in valid_source_types
        assert source["confidence"] in valid_confidences
        assert source["ref"]
        assert source["version"]


def test_external_only_fixture_cannot_claim_executable_status() -> None:
    external_case = next(
        case for case in _load_cases() if case["case_type"] == "external_evidence_only"
    )

    assert external_case["execution_status"] == OEExecutionStatus.EXTERNAL_EVIDENCE_ONLY.value
    assert external_case["execution_status"] != OEExecutionStatus.GENE_LEVEL_EXECUTABLE.value
    assert external_case["source"]["type"] in {
        EvidenceSourceType.EXTERNAL_PICHIA_MODEL.value,
        EvidenceSourceType.PICHIA_LITERATURE.value,
        EvidenceSourceType.HOMOLOGY_TRANSFER.value,
        EvidenceSourceType.SMOKE_FIXTURE.value,
    }


def test_complex_without_stoichiometry_is_complex_limited() -> None:
    complex_case = next(
        case for case in _load_cases() if case["case_type"] == "complex_subunit"
    )

    assert complex_case["gpr_role"] == GPRRole.COMPLEX_SUBUNIT.value
    assert complex_case["subunit_stoichiometry"] == []
    assert complex_case["execution_status"] == OEExecutionStatus.COMPLEX_LIMITED.value
    assert "subunit_stoichiometry" in complex_case["missing_information"]


def test_missing_parameter_fixture_names_each_missing_parameter() -> None:
    missing_case = next(
        case for case in _load_cases() if case["case_type"] == "missing_parameter"
    )

    assert missing_case["execution_status"] == OEExecutionStatus.PARTIAL_MAPPING.value
    assert set(missing_case["missing_information"]) >= {
        "kcat",
        "baseline_enzyme_amount",
    }
