"""筛查结果的基因名显示层：把模型 locus tag（如 `PAS_chr3_0199`）补上研究员看得懂的
数据库注释名，并**按注释档位门控 + 显式标注可信度**，不把推断当事实。

诚实边界（2026-07-28 对 catalog 命名准确性的评估结论，决定了本模块的门控设计）：

- 名字来源是 locus tag **精确字符串匹配** UniProt / KEGG（`gene_evidence.py` 里
  `gene_exact:` 查询 + Pichia proteome 限定），**全程没有序列相似 / BLAST / 模糊名匹配**
  这类会张冠李戴的回退：查得到才贴名，查不到就留空标 `low_model_only`（不编）。
  故 1017/1025 命中 `high_exact_locus_tag`，8 个无注释。
- **但精确匹配只保证【名字 ↔ 位点】对得上，不保证【位点 ↔ 研究员心里的那个基因】**。
  项目文档已记录一类反例：策展俗名 → 模型位点的对应本身低置信待复核
  （见 docs/pichia_current_architecture_and_requirements.md 的 PEP4/PRB1/YPS 条目）。
  这类位点单独打 `identity_review`，UI 必须显式提示。
- **描述性注释 ≠ 正式基因符号**：只有 ~46/1025 有真正的 symbol（如 SEC11），其余是蛋白
  描述串；另有约 56 个是 `uncharacterized / domain-containing protein` 这类泛化自动注释
  （诚实但没功能信息量）。三者用 `annotation_tier` 区分，不混作一谈。
- 名字只用于**显示**，绝不参与可执行 id：KO/OE 求解一律仍按 `gene_id` 走。

分层：本模块属 app/services（门面），只读已有标准化缓存、不触发重建（重建要读模型 / 联网，
不能放在页面渲染路径上）；缓存缺失时优雅降级为"只显示 locus tag"。
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from app import ensure_python_pichia_on_path
from app.services.pichia_gene_catalog_service import (
    pichia_gene_id_standardization_cache_path,
)

ensure_python_pichia_on_path()

from pcsec_pichia.services.gene_id_standardization import (  # noqa: E402 - engine import after path bootstrap (services are allowlisted)
    PichiaGeneIdStandardName,
    build_standard_name_lookup,
    load_pichia_gene_id_standardization_cache,
    standard_name_fields_for_gene,
)

# 只有这一档是"locus tag 精确命中外部库"，才允许把名字显示给研究员当参考。
EXACT_LOCUS_TAG_CONFIDENCE = "high_exact_locus_tag"

# 泛化自动注释标志词：命中即降档为 `generic_annotation`（诚实但无功能信息量，别让研究员
# 以为这是已知功能）。小写比对。
GENERIC_ANNOTATION_MARKERS: tuple[str, ...] = (
    "uncharacterized",
    "hypothetical",
    "putative",
    "probable",
    "predicted",
    "domain-containing protein",
    "duf",
)

# 策展俗名 → 模型位点的对应**本身**低置信待复核（不是名字错，是"这个位点是不是那个基因"存疑）。
# 来源：docs/pichia_current_architecture_and_requirements.md（PEP4/PRB1 模型基因 ID 标注为
# 低置信度待复核；YPS1-3 尚未进入基因目录，故此处无位点可标）。这类位点上外部库给出的
# 名字是**该位点真实的**功能注释，可能与团队口头说的俗名完全不同 —— 必须显式提示，
# 否则研究员会以为工具在说"这就是 PEP4"。
UNVERIFIED_IDENTITY_LOCI: dict[str, str] = {
    "PAS_chr2-2_0107": "PEP4",
    "PAS_chr2-1_0785": "PRB1",
}

NAME_ANNOTATION_COLUMNS: tuple[str, ...] = (
    "gene_display_name",
    "standard_symbol",
    "protein_name",
    "annotation_confidence",
    "standard_name_status",
    "annotation_tier",
    "annotation_accession",
    "identity_review",
)

# 研究员向的一句话免责说明（面板 caption 复用，保证口径一致）。
GENE_NAME_CAVEAT = (
    "基因名是**外部数据库（UniProt/KEGG）按位点号精确匹配**给出的描述性注释，"
    "**不是经验证的正式基因符号**，也不代表该位点就是你口头说的那个基因；"
    "括号里的位点号（如 PAS_chr3_0199）才是模型实际算的对象。"
    "标「泛化注释」的只是自动注释（如“某结构域蛋白”）、没有功能结论；"
    "标「身份待复核」的位点其俗名对应关系尚未坐实，下实验前请按位点号 / 数据库 ID 自行核对。"
)


def load_standard_name_lookup(paths: Any | None = None) -> dict[str, PichiaGeneIdStandardName]:
    """读基因命名标准化缓存（只读，不重建）。缺失 / 损坏 → 返回空表（显示降级为纯 locus tag）。

    刻意不用 `load_pichia_gene_id_standardization()`：那个在缓存缺失时会重建（读模型 + 联网），
    不能放在页面渲染路径上。
    """
    try:
        cache_path = pichia_gene_id_standardization_cache_path(paths)
    except Exception:  # noqa: BLE001 - 路径发现失败也只应降级、不该打断结果页
        return {}
    if not cache_path.exists():
        return {}
    try:
        rows = load_pichia_gene_id_standardization_cache(cache_path)
    except (OSError, ValueError, UnicodeDecodeError):
        return {}
    if not rows:
        return {}
    return build_standard_name_lookup(list(rows))


def classify_annotation_tier(standard_symbol: str, protein_name: str) -> str:
    """名字的信息量分档：正式符号 > 描述性注释 > 泛化自动注释 > 无注释。"""
    if str(standard_symbol or "").strip():
        return "verified_symbol"
    name = str(protein_name or "").strip()
    if not name:
        return "no_annotation"
    lowered = name.lower()
    if any(marker in lowered for marker in GENERIC_ANNOTATION_MARKERS):
        return "generic_annotation"
    return "descriptive_annotation"


def _accession(external_ids: object) -> str:
    """给研究员一个能点开自查的库 ID（优先 UniProt，退 KEGG）。"""
    if not isinstance(external_ids, dict):
        return ""
    uniprot = str(external_ids.get("uniprot") or "").strip()
    if uniprot:
        return f"UniProt:{uniprot}"
    kegg = str(external_ids.get("kegg") or "").strip()
    return f"KEGG:{kegg}" if kegg else ""


def annotate_gene_name_fields(
    gene_id: str,
    lookup: dict[str, PichiaGeneIdStandardName],
    *,
    candidate_kind: str = "gene",
    intervention_type: str = "",
) -> dict[str, Any]:
    """单个候选的显示用命名字段（已按档位门控 + 打身份复核标记）。

    非基因候选（策展反应 / 复合体假设）由引擎侧判为 `not_gene_candidate`，本函数原样透传、
    不硬塞基因名。
    """
    fields = standard_name_fields_for_gene(
        gene_id,
        lookup,
        candidate_kind=candidate_kind,
        intervention_type=intervention_type,
    )
    confidence = str(fields.get("annotation_confidence") or "")
    resolved_gene_id = str(gene_id or "").strip()

    # 门控：只有精确 locus tag 命中档才把名字显示出去；其余（含 low_model_only）一律不给名字，
    # 让 UI 老实显示"仅模型位点、无注释"，而不是给一个来路不明的名字。
    if confidence != EXACT_LOCUS_TAG_CONFIDENCE:
        return {
            "gene_display_name": "",
            "standard_symbol": "",
            "protein_name": "",
            "annotation_confidence": confidence,
            "standard_name_status": str(fields.get("standard_name_status") or ""),
            "annotation_tier": "no_annotation",
            "annotation_accession": "",
            "identity_review": _identity_review(resolved_gene_id),
        }

    standard_symbol = str(fields.get("standard_symbol") or "").strip()
    protein_name = str(fields.get("protein_name") or "").strip()
    display_name = str(fields.get("gene_display_name") or "").strip()
    return {
        "gene_display_name": display_name,
        "standard_symbol": standard_symbol,
        "protein_name": protein_name,
        "annotation_confidence": confidence,
        "standard_name_status": str(fields.get("standard_name_status") or ""),
        "annotation_tier": classify_annotation_tier(standard_symbol, protein_name),
        "annotation_accession": _accession(fields.get("external_ids")),
        "identity_review": _identity_review(resolved_gene_id),
    }


def _identity_review(gene_id: str) -> str:
    return "curated_identity_unverified" if gene_id in UNVERIFIED_IDENTITY_LOCI else ""


def annotate_screen_frame(
    frame: pd.DataFrame,
    lookup: dict[str, PichiaGeneIdStandardName] | None = None,
    *,
    paths: Any | None = None,
) -> pd.DataFrame:
    """给筛查 tradeoff 表补显示用命名列（不改任何可执行 id、不改行序）。

    对**旧筛查 CSV**（命名功能上线前跑的，只有 `gene_id`）和新 CSV 都适用：新 CSV 已带的
    非空名字原样保留，只补空缺，避免覆盖当次运行记录下来的值。
    """
    if frame.empty or "gene_id" not in frame.columns:
        return frame
    resolved_lookup = load_standard_name_lookup(paths) if lookup is None else lookup
    if not resolved_lookup:
        # 缓存不可用：补齐列（保持下游列契约稳定），值为空 → UI 只显示 locus tag。
        annotated = frame.copy()
        for column in NAME_ANNOTATION_COLUMNS:
            if column not in annotated.columns:
                annotated[column] = ""
        return annotated

    annotated = frame.copy()
    for column in NAME_ANNOTATION_COLUMNS:
        if column not in annotated.columns:
            annotated[column] = ""

    has_existing_name = annotated["gene_display_name"].astype(str).str.strip() != ""
    kinds = (
        annotated["candidate_kind"].astype(str)
        if "candidate_kind" in annotated.columns
        else pd.Series("gene", index=annotated.index)
    )
    interventions = (
        annotated["intervention_type"].astype(str)
        if "intervention_type" in annotated.columns
        else pd.Series("", index=annotated.index)
    )

    computed: dict[str, list[Any]] = {column: [] for column in NAME_ANNOTATION_COLUMNS}
    for position, index in enumerate(annotated.index):
        if bool(has_existing_name.iloc[position]):
            # 已带名字的新 CSV：保留原值，只补本模块新增的派生列（档位 / accession / 身份复核）。
            row_symbol = str(annotated.at[index, "standard_symbol"] or "").strip()
            row_protein = str(annotated.at[index, "protein_name"] or "").strip()
            fields = {
                "gene_display_name": str(annotated.at[index, "gene_display_name"] or ""),
                "standard_symbol": row_symbol,
                "protein_name": row_protein,
                "annotation_confidence": str(annotated.at[index, "annotation_confidence"] or ""),
                "standard_name_status": str(annotated.at[index, "standard_name_status"] or ""),
                "annotation_tier": classify_annotation_tier(row_symbol, row_protein),
                "annotation_accession": str(annotated.at[index, "annotation_accession"] or ""),
                "identity_review": _identity_review(str(annotated.at[index, "gene_id"] or "").strip()),
            }
        else:
            fields = annotate_gene_name_fields(
                str(annotated.at[index, "gene_id"] or ""),
                resolved_lookup,
                candidate_kind=kinds.iloc[position],
                intervention_type=interventions.iloc[position],
            )
        for column in NAME_ANNOTATION_COLUMNS:
            computed[column].append(fields.get(column, ""))

    for column in NAME_ANNOTATION_COLUMNS:
        annotated[column] = computed[column]
    return annotated


def safe_gene_display_label(row: Any) -> str:
    """一个候选的显示标签。名字梯队：正式符号 → 标准显示名 → 蛋白描述名 → 策展常用名。

    位点号只在名字来自**外部数据库注释**时附在括号里——研究员需要它来核对"模型实际算的是哪个
    位点"，也因为数据库注释名可能与俗名不是一回事。**策展候选**（curated catalog 的复合体 /
    反应级条目，如 `PDI1/ERO1`）不加后缀：那是人工挑的可读标签、本身就是该候选的身份，
    后面缀一串模型反应名（`sec_..._complex_formation`）只会变噪声。无任何名字时只给 id。
    """

    def value(key: str) -> str:
        try:
            raw = row.get(key)
        except AttributeError:
            return ""
        return str(raw or "").strip()

    gene_id = value("gene_id")
    annotation_name = value("standard_symbol") or value("gene_display_name") or value("protein_name")
    if annotation_name and annotation_name != gene_id:
        return f"{annotation_name}（{gene_id}）"
    curated_name = value("common_name")
    if curated_name and curated_name != gene_id:
        return curated_name
    return gene_id


__all__ = [
    "EXACT_LOCUS_TAG_CONFIDENCE",
    "GENE_NAME_CAVEAT",
    "GENERIC_ANNOTATION_MARKERS",
    "NAME_ANNOTATION_COLUMNS",
    "UNVERIFIED_IDENTITY_LOCI",
    "annotate_gene_name_fields",
    "annotate_screen_frame",
    "classify_annotation_tier",
    "load_standard_name_lookup",
    "safe_gene_display_label",
]
