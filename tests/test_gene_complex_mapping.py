"""E3 基因↔复合体映射层（ADR-007）：数据契约、门禁、双向查询、无数据时的优雅降级。

最关键的两条门禁（不可放宽）：
- 未复核条目不得进入可执行路径（ADR-001）；
- **亚基化学计量未知时不得声称单基因 OE 能提升复合体容量**——否则就是虚构容量提升。
"""

from __future__ import annotations

import json
from typing import Any

from pcsec_pichia.services import gene_complex_mapping as mapping


def _payload(**overrides: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "pichia_gene_id": "PAS_chr4_0844",
        "complex_reaction_id": "sec_Pdi1p_complex_formation",
        "subunit_role": mapping.SUBUNIT_ROLE_REQUIRED,
        "stoichiometry_status": mapping.STOICHIOMETRY_KNOWN,
        "review_status": mapping.REVIEW_REVIEWED,
        "evidence_source": "curated_literature",
    }
    row.update(overrides)
    return row


def test_fully_specified_reviewed_mapping_may_enter_executable_single_gene_oe() -> None:
    rows, problems = mapping.validate_gene_complex_mapping_payloads([_payload()])

    assert problems == ()
    assert rows[0].may_enter_executable_single_gene_oe is True


def test_unknown_stoichiometry_must_not_claim_single_gene_oe_capacity() -> None:
    """核心边界：化学计量未知 → 不得据此声称单基因 OE 能提升复合体容量。"""
    rows, _ = mapping.validate_gene_complex_mapping_payloads(
        [_payload(stoichiometry_status=mapping.STOICHIOMETRY_UNKNOWN)]
    )

    assert rows[0].is_reviewed is True, "条目本身有效，只是不能用于容量声称"
    assert rows[0].may_enter_executable_single_gene_oe is False


def test_unreviewed_or_rejected_mapping_never_enters_executable_path() -> None:
    for status in (mapping.REVIEW_PENDING, mapping.REVIEW_REJECTED):
        rows, _ = mapping.validate_gene_complex_mapping_payloads([_payload(review_status=status)])
        assert rows[0].may_enter_executable_single_gene_oe is False, status
        assert rows[0].is_reviewed is False, status


def test_non_required_subunit_roles_do_not_unlock_capacity_claim() -> None:
    for role in (mapping.SUBUNIT_ROLE_REPLACEABLE, mapping.SUBUNIT_ROLE_AUXILIARY):
        rows, _ = mapping.validate_gene_complex_mapping_payloads([_payload(subunit_role=role)])
        assert rows[0].may_enter_executable_single_gene_oe is False, role


def test_contract_violations_are_dropped_not_silently_downgraded() -> None:
    """来源不明 / 枚举非法的条目必须丢弃并报告——放行等于把猜测伪装成证据。"""
    payloads = [
        _payload(pichia_gene_id=""),
        _payload(complex_reaction_id=""),
        _payload(subunit_role="谁知道"),
        _payload(stoichiometry_status="maybe"),
        _payload(review_status="looks_fine"),
        _payload(evidence_source=""),
        "not a dict",
    ]

    rows, problems = mapping.validate_gene_complex_mapping_payloads(payloads)

    assert rows == ()
    assert len(problems) == len(payloads)
    assert any("evidence_source" in problem for problem in problems)


def test_duplicate_gene_complex_pairs_are_dropped() -> None:
    rows, problems = mapping.validate_gene_complex_mapping_payloads([_payload(), _payload()])

    assert len(rows) == 1
    assert any("重复" in problem for problem in problems)


def test_both_query_directions_and_reviewed_only_gate() -> None:
    rows, _ = mapping.validate_gene_complex_mapping_payloads(
        [
            _payload(pichia_gene_id="PAS_A"),
            _payload(pichia_gene_id="PAS_B", review_status=mapping.REVIEW_PENDING),
            _payload(pichia_gene_id="PAS_A", complex_reaction_id="sec_Other_complex_formation"),
        ]
    )

    # 反向：这个复合体对应哪几个基因（实验室要动的）
    reviewed = mapping.genes_for_complex(rows, "sec_Pdi1p_complex_formation")
    everything = mapping.genes_for_complex(rows, "sec_Pdi1p_complex_formation", reviewed_only=False)
    assert [row.pichia_gene_id for row in reviewed] == ["PAS_A"]
    assert [row.pichia_gene_id for row in everything] == ["PAS_A", "PAS_B"]

    # 正向：这个基因参与哪些复合体（让基因筛查够得到分泌层）
    complexes = mapping.complexes_for_gene(rows, "PAS_A")
    assert {row.complex_reaction_id for row in complexes} == {
        "sec_Pdi1p_complex_formation",
        "sec_Other_complex_formation",
    }

    assert mapping.genes_for_complex(rows, "") == ()
    assert mapping.complexes_for_gene(rows, "") == ()


def test_draft_from_candidates_never_self_activates() -> None:
    """草稿的意义是把"从零查"降成"打勾"，但**绝不能自己生效**：
    一律最保守三档 → may_enter_executable_single_gene_oe 恒 False。"""
    drafts = mapping.build_draft_mappings_from_candidates(
        [
            {
                "gene_id": "PAS_chr4_0844",
                "source_common_name": "PDI1（单独）",
                "review_reactions": ["sec_Pdi1p_complex_formation"],
                "homology_review_status": "rbh_not_in_model",
            }
        ]
    )

    assert len(drafts) == 1
    draft = drafts[0]
    assert draft.review_status == mapping.REVIEW_PENDING
    assert draft.subunit_role == mapping.SUBUNIT_ROLE_AUXILIARY
    assert draft.stoichiometry_status == mapping.STOICHIOMETRY_UNKNOWN
    assert draft.may_enter_executable_single_gene_oe is False
    assert draft.is_reviewed is False
    assert "同源" in draft.note and "PDI1" in draft.evidence_citation


def test_draft_skips_entries_without_a_resolved_gene_and_dedupes() -> None:
    """复合体/家族名（OST 复合体、KTR）没有单一基因，软件不猜——留给人工查。"""
    drafts = mapping.build_draft_mappings_from_candidates(
        [
            {"gene_id": "", "source_common_name": "OST 复合体", "review_reactions": ["sec_OSTC_complex_formation"]},
            {"gene_id": "PAS_A", "review_reactions": ["sec_X_complex_formation"]},
            {"gene_id": "PAS_A", "review_reactions": ["sec_X_complex_formation"]},
        ]
    )

    assert [(row.pichia_gene_id, row.complex_reaction_id) for row in drafts] == [
        ("PAS_A", "sec_X_complex_formation")
    ]


def test_draft_round_trips_through_the_contract_validator() -> None:
    """起草 → 序列化 → 再校验必须无损，否则导出的文件读不回来。"""
    drafts = mapping.build_draft_mappings_from_candidates(
        [{"gene_id": "PAS_A", "review_reactions": ["sec_X_complex_formation", "sec_Y_complex_formation"]}]
    )
    payload = mapping.serialize_gene_complex_mappings(drafts)
    reloaded, problems = mapping.validate_gene_complex_mapping_payloads(payload["mappings"])

    assert problems == ()
    assert len(reloaded) == len(drafts) == 2
    assert payload["schema_version"] == mapping.GENE_COMPLEX_MAPPING_SCHEMA_VERSION


def test_missing_curation_file_degrades_gracefully(tmp_path) -> None:
    """策展数据未到位是**预期状态**（待拍板），不能抛异常。"""
    rows, notes = mapping.load_gene_complex_mapping_file(tmp_path / "nope.json")

    assert rows == ()
    assert any("未找到" in note for note in notes)


def test_schema_version_mismatch_and_corrupt_file_are_ignored_with_a_note(tmp_path) -> None:
    stale = tmp_path / "stale.json"
    stale.write_text(json.dumps({"schema_version": 999, "mappings": [_payload()]}), encoding="utf-8")
    rows, notes = mapping.load_gene_complex_mapping_file(stale)
    assert rows == ()
    assert any("schema_version" in note for note in notes)

    broken = tmp_path / "broken.json"
    broken.write_text("{not json", encoding="utf-8")
    rows, notes = mapping.load_gene_complex_mapping_file(broken)
    assert rows == ()
    assert notes


def test_round_trip_from_a_real_curation_file(tmp_path) -> None:
    path = tmp_path / "gene_complex_mapping.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": mapping.GENE_COMPLEX_MAPPING_SCHEMA_VERSION,
                "mappings": [_payload(), _payload(pichia_gene_id="PAS_ERO1")],
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    rows, problems = mapping.load_gene_complex_mapping_file(path)
    summary = mapping.summarize_gene_complex_mapping_rows(rows)

    assert problems == ()
    assert summary["mapping_count"] == 2
    assert summary["reviewed_count"] == 2
    assert summary["distinct_complexes"] == 1
    assert summary["single_gene_oe_eligible_count"] == 2
