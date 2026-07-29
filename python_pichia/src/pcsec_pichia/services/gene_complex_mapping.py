"""E3 分泌机器「基因 ↔ 复合体」映射层（ADR-007）。

**这一层解决什么**：模型的 GPR 只覆盖代谢反应；分泌机器是 2793 个复合体形成反应、**零基因关联**。
于是有两个方向的问题答不了：

1. 正向（覆盖面）：基于基因的全基因组筛查触达不到分泌机器复合体；
2. 反向（实验可执行性）：模型说"过表达 `sec_Pdi1p_complex_formation` 有效"，但实验室要过表达
   **哪几个基因**没有依据。

本模块只做数据契约 + 校验 + 门禁 + 查询，**不含策展内容**：映射条目由人工策展提供（生物学判断），
数据缺席时全链路优雅降级为空（与 ADR-005 的 RNA-seq 做法一致）。

**边界（不可放宽）**：
- 本层只扩大"模型够得到的干预点范围"，**不提升绝对预测准确度**；绝对容量按 ADR-002 恒 unavailable。
- 未复核条目不得进入可执行路径（沿用 ADR-001 门禁）。
- **亚基化学计量未知时，不得据此声称单基因 OE 能提升复合体容量**——这是既有约定
  "复合体亚基默认不做单基因 OE"的延续，防止虚构容量提升。
- 不改写 `Model/`、`Enzymedata/` 等受保护科学资产的 GPR，也不写回参考模型；本模块只读。
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


GENE_COMPLEX_MAPPING_SCHEMA_VERSION = 1

# 该基因在复合体里扮演什么角色。决定能不能把"过表达这个基因"当作"提升该复合体容量"。
SUBUNIT_ROLE_REQUIRED = "required_subunit"
SUBUNIT_ROLE_REPLACEABLE = "replaceable_isoenzyme"
SUBUNIT_ROLE_AUXILIARY = "auxiliary"
SUBUNIT_ROLES = frozenset({SUBUNIT_ROLE_REQUIRED, SUBUNIT_ROLE_REPLACEABLE, SUBUNIT_ROLE_AUXILIARY})

STOICHIOMETRY_KNOWN = "known"
STOICHIOMETRY_UNKNOWN = "unknown"
STOICHIOMETRY_STATUSES = frozenset({STOICHIOMETRY_KNOWN, STOICHIOMETRY_UNKNOWN})

REVIEW_REVIEWED = "reviewed"
REVIEW_PENDING = "pending_review"
REVIEW_REJECTED = "rejected"
REVIEW_STATUSES = frozenset({REVIEW_REVIEWED, REVIEW_PENDING, REVIEW_REJECTED})


@dataclass(frozen=True)
class GeneComplexMapping:
    """一条策展映射：某个 Pichia 基因参与某个模型复合体形成反应。"""

    pichia_gene_id: str
    complex_reaction_id: str
    subunit_role: str
    stoichiometry_status: str
    review_status: str
    evidence_source: str
    evidence_version: str = ""
    evidence_retrieved_at: str = ""
    evidence_license: str = ""
    evidence_citation: str = ""
    note: str = ""

    @property
    def is_reviewed(self) -> bool:
        return self.review_status == REVIEW_REVIEWED

    @property
    def may_enter_executable_single_gene_oe(self) -> bool:
        """能否把这条映射用于"单基因 OE 提升该复合体容量"的可执行路径。

        三个条件同时满足才行：已人工复核 + 亚基化学计量已知 + 该基因是必需亚基。
        这只是**允许考虑**，不等于断言提升幅度——幅度仍由模型求解给出，且仍是相对信号。
        """
        return (
            self.is_reviewed
            and self.stoichiometry_status == STOICHIOMETRY_KNOWN
            and self.subunit_role == SUBUNIT_ROLE_REQUIRED
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def validate_gene_complex_mapping_payloads(
    payloads: list[dict[str, Any]] | tuple[dict[str, Any], ...],
) -> tuple[tuple[GeneComplexMapping, ...], tuple[str, ...]]:
    """按 ADR-007 数据契约校验策展条目。返回 (通过的条目, 问题描述)。

    契约不满足的条目**被丢弃而不是降级放行**——放行等于把来源不明的映射伪装成证据。
    """
    valid: list[GeneComplexMapping] = []
    problems: list[str] = []
    seen: set[tuple[str, str]] = set()

    for index, payload in enumerate(payloads):
        if not isinstance(payload, dict):
            problems.append(f"第 {index} 条不是对象，已丢弃。")
            continue
        gene_id = str(payload.get("pichia_gene_id") or "").strip()
        reaction_id = str(payload.get("complex_reaction_id") or "").strip()
        label = f"第 {index} 条（{gene_id or '?'} -> {reaction_id or '?'}）"
        if not gene_id or not reaction_id:
            problems.append(f"{label} 缺 pichia_gene_id 或 complex_reaction_id，已丢弃。")
            continue
        key = (gene_id, reaction_id)
        if key in seen:
            problems.append(f"{label} 与前面重复，已丢弃。")
            continue

        subunit_role = str(payload.get("subunit_role") or "").strip()
        stoichiometry_status = str(payload.get("stoichiometry_status") or "").strip()
        review_status = str(payload.get("review_status") or "").strip()
        evidence_source = str(payload.get("evidence_source") or "").strip()
        if subunit_role not in SUBUNIT_ROLES:
            problems.append(f"{label} subunit_role={subunit_role!r} 不在 {sorted(SUBUNIT_ROLES)}，已丢弃。")
            continue
        if stoichiometry_status not in STOICHIOMETRY_STATUSES:
            problems.append(
                f"{label} stoichiometry_status={stoichiometry_status!r} 不在 {sorted(STOICHIOMETRY_STATUSES)}，已丢弃。"
            )
            continue
        if review_status not in REVIEW_STATUSES:
            problems.append(f"{label} review_status={review_status!r} 不在 {sorted(REVIEW_STATUSES)}，已丢弃。")
            continue
        if not evidence_source:
            problems.append(f"{label} 缺 evidence_source（来源必须可追溯），已丢弃。")
            continue

        seen.add(key)
        valid.append(
            GeneComplexMapping(
                pichia_gene_id=gene_id,
                complex_reaction_id=reaction_id,
                subunit_role=subunit_role,
                stoichiometry_status=stoichiometry_status,
                review_status=review_status,
                evidence_source=evidence_source,
                evidence_version=str(payload.get("evidence_version") or ""),
                evidence_retrieved_at=str(payload.get("evidence_retrieved_at") or ""),
                evidence_license=str(payload.get("evidence_license") or ""),
                evidence_citation=str(payload.get("evidence_citation") or ""),
                note=str(payload.get("note") or ""),
            )
        )
    return tuple(valid), tuple(problems)


def load_gene_complex_mapping_file(path: Path) -> tuple[tuple[GeneComplexMapping, ...], tuple[str, ...]]:
    """读策展映射文件（只读）。文件缺失 / 损坏 / schema 不符时返回空 + 说明，不抛异常。"""
    if not path.exists():
        return (), (f"未找到策展映射文件：{path}（尚无策展数据，本层优雅降级为空）。",)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, UnicodeDecodeError) as exc:
        return (), (f"策展映射文件无法解析：{exc}",)
    if not isinstance(payload, dict):
        return (), ("策展映射文件顶层必须是对象。",)
    version = payload.get("schema_version")
    if version != GENE_COMPLEX_MAPPING_SCHEMA_VERSION:
        return (), (
            f"策展映射 schema_version={version!r}，当前期望 {GENE_COMPLEX_MAPPING_SCHEMA_VERSION}；已忽略。",
        )
    rows = payload.get("mappings")
    if not isinstance(rows, list):
        return (), ("策展映射文件缺 mappings 数组。",)
    return validate_gene_complex_mapping_payloads(rows)


DRAFT_NOTE = "草稿：基因来自同源比对，角色与化学计量待人工判断"


def build_draft_mappings_from_candidates(
    candidate_rows: list[dict[str, Any]] | tuple[dict[str, Any], ...],
) -> tuple[GeneComplexMapping, ...]:
    """从策展候选行自动起草映射，供人工复核——把"从零查"变成"打勾/否决"。

    草稿一律取**最保守**的三档：`auxiliary` + `unknown` + `pending_review`，因此
    `may_enter_executable_single_gene_oe` 恒为 False——**草稿绝不会自己生效**，必须有人复核后
    才可能进入可执行路径。角色与化学计量是生物学判断，软件不猜。
    """
    drafts: list[GeneComplexMapping] = []
    seen: set[tuple[str, str]] = set()
    for row in candidate_rows:
        gene_id = str(row.get("gene_id") or "").strip()
        if not gene_id:
            continue  # 复合体/家族名（OST 复合体、KTR 等）没有单一基因，只能人工查
        reactions = [
            str(item).strip()
            for item in (
                row.get("review_reactions")
                or row.get("executable_oe_proxy_reactions")
                or row.get("executable_ko_reactions")
                or []
            )
            if str(item).strip()
        ]
        homology_status = str(row.get("homology_review_status") or "")
        common_name = str(row.get("source_common_name") or "").strip()
        for reaction_id in reactions:
            key = (gene_id, reaction_id)
            if key in seen:
                continue
            seen.add(key)
            drafts.append(
                GeneComplexMapping(
                    pichia_gene_id=gene_id,
                    complex_reaction_id=reaction_id,
                    subunit_role=SUBUNIT_ROLE_AUXILIARY,
                    stoichiometry_status=STOICHIOMETRY_UNKNOWN,
                    review_status=REVIEW_PENDING,
                    evidence_source="homology_rbh_draft",
                    evidence_citation=f"策展俗名 {common_name}" if common_name else "",
                    note=f"{DRAFT_NOTE}；同源状态={homology_status or '未知'}",
                )
            )
    return tuple(drafts)


def serialize_gene_complex_mappings(
    rows: tuple[GeneComplexMapping, ...] | list[GeneComplexMapping],
) -> dict[str, Any]:
    """按策展文件格式序列化（供导出 / 回填）。"""
    return {
        "schema_version": GENE_COMPLEX_MAPPING_SCHEMA_VERSION,
        "mappings": [row.to_dict() for row in rows],
    }


def genes_for_complex(
    rows: tuple[GeneComplexMapping, ...] | list[GeneComplexMapping],
    complex_reaction_id: str,
    *,
    reviewed_only: bool = True,
) -> tuple[GeneComplexMapping, ...]:
    """反向查询：模型说过表达这个复合体有效 → 实验室该动哪几个基因。"""
    reaction_id = str(complex_reaction_id or "").strip()
    if not reaction_id:
        return ()
    return tuple(
        row
        for row in rows
        if row.complex_reaction_id == reaction_id and (row.is_reviewed or not reviewed_only)
    )


def complexes_for_gene(
    rows: tuple[GeneComplexMapping, ...] | list[GeneComplexMapping],
    pichia_gene_id: str,
    *,
    reviewed_only: bool = True,
) -> tuple[GeneComplexMapping, ...]:
    """正向查询：这个基因参与了哪些分泌机器复合体（让基于基因的筛查够得到分泌层）。"""
    gene_id = str(pichia_gene_id or "").strip()
    if not gene_id:
        return ()
    return tuple(
        row for row in rows if row.pichia_gene_id == gene_id and (row.is_reviewed or not reviewed_only)
    )


def summarize_gene_complex_mapping_rows(
    rows: tuple[GeneComplexMapping, ...] | list[GeneComplexMapping],
) -> dict[str, Any]:
    all_rows = tuple(rows)
    return {
        "schema_version": GENE_COMPLEX_MAPPING_SCHEMA_VERSION,
        "mapping_count": len(all_rows),
        "reviewed_count": sum(1 for row in all_rows if row.is_reviewed),
        "pending_review_count": sum(1 for row in all_rows if row.review_status == REVIEW_PENDING),
        "distinct_genes": len({row.pichia_gene_id for row in all_rows}),
        "distinct_complexes": len({row.complex_reaction_id for row in all_rows}),
        "single_gene_oe_eligible_count": sum(1 for row in all_rows if row.may_enter_executable_single_gene_oe),
        "stoichiometry_unknown_count": sum(
            1 for row in all_rows if row.stoichiometry_status == STOICHIOMETRY_UNKNOWN
        ),
    }


__all__ = [
    "DRAFT_NOTE",
    "GENE_COMPLEX_MAPPING_SCHEMA_VERSION",
    "build_draft_mappings_from_candidates",
    "serialize_gene_complex_mappings",
    "REVIEW_PENDING",
    "REVIEW_REJECTED",
    "REVIEW_REVIEWED",
    "REVIEW_STATUSES",
    "STOICHIOMETRY_KNOWN",
    "STOICHIOMETRY_STATUSES",
    "STOICHIOMETRY_UNKNOWN",
    "SUBUNIT_ROLES",
    "SUBUNIT_ROLE_AUXILIARY",
    "SUBUNIT_ROLE_REPLACEABLE",
    "SUBUNIT_ROLE_REQUIRED",
    "GeneComplexMapping",
    "complexes_for_gene",
    "genes_for_complex",
    "load_gene_complex_mapping_file",
    "summarize_gene_complex_mapping_rows",
    "validate_gene_complex_mapping_payloads",
]
