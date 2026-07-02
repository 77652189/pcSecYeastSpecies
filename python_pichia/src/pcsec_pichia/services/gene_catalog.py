"""毕赤酵母分泌通路基因目录。

按分泌通路步骤组织的已知 KO/OE 靶点，来源于 pcSecPichia 模型
和已报道的毕赤酵母分泌工程文献。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


InterventionType = Literal["KO", "OE", "both"]


@dataclass(frozen=True)
class SecretionGeneEntry:
    """一条分泌相关基因或反应的目录条目。

    Attributes:
        category:        分泌通路分类
        common_name:     常用名（如 Kar2、PDI1）
        description:     中文描述
        gene_id:         模型基因 ID（PAS_chr...），用于 KO
        oe_reaction_id:  过表达反应 ID（sec_...），用于 OE
        intervention:    建议的扰动类型
        evidence:        参考来源
        homolog_note:    毕赤酵母中的同源基因说明
    """

    category: str
    common_name: str
    description: str
    gene_id: str = ""
    ko_reaction_id: str = ""  # 直接 KO 的反应 ID（用于无 gene_id 的复合体）
    oe_reaction_id: str = ""
    intervention: InterventionType = "both"
    evidence: str = ""
    homolog_note: str = ""


# ---------------------------------------------------------------------------
# 分泌通路分类
# ---------------------------------------------------------------------------
CAT_ER_TRANSLOCATION = "ER 转运"
CAT_ER_FOLDING = "ER 折叠与分子伴侣"
CAT_DSB = "二硫键 (DSB)"
CAT_N_GLYCAN = "N-糖基化"
CAT_O_GLYCAN = "O-糖基化"
CAT_ERAD = "错误折叠与 ERAD"
CAT_COPII = "COPII 囊泡转运"
CAT_GOLGI = "Golgi 加工"
CAT_COPI = "COPI 逆向转运"
CAT_EXOCYTOSIS = "胞吐与分泌"
CAT_PROTEASOME = "蛋白酶体与降解"
CAT_GENERAL = "通用/其他"

# ---------------------------------------------------------------------------
# 完整基因目录
# ---------------------------------------------------------------------------

SECRETION_GENE_CATALOG: tuple[SecretionGeneEntry, ...] = (
    # ======================== ER 转运 ========================
    SecretionGeneEntry(
        category=CAT_ER_TRANSLOCATION,
        common_name="SEC61",
        description="ER 易位通道核心亚基，蛋白进入 ER 的门户",
        oe_reaction_id="sec_SEC61SEC63C_complex_formation",
        intervention="OE",
        evidence="模型 sec_SEC61SEC63C 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ER_TRANSLOCATION,
        common_name="SSH1",
        description="SEC61 旁系同源易位通道",
        oe_reaction_id="sec_SSH1C_complex_formation",
        intervention="OE",
        evidence="模型 sec_SSH1C 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ER_TRANSLOCATION,
        common_name="SRP/SRP受体",
        description="信号肽识别颗粒受体，介导共翻译转运",
        oe_reaction_id="sec_SRPC_complex_formation",
        intervention="OE",
        evidence="模型 sec_SRPC/SRC 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ER_TRANSLOCATION,
        common_name="SPC",
        description="信号肽酶复合体，切除信号肽",
        oe_reaction_id="sec_SPC_complex_formation",
        intervention="OE",
        evidence="模型 sec_SPC 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ER_TRANSLOCATION,
        common_name="RAC",
        description="核糖体对接复合体，连接翻译与易位",
        oe_reaction_id="sec_RAC_complex_formation",
        intervention="OE",
        evidence="模型 sec_RAC 复合体",
    ),

    # ======================== ER 折叠与分子伴侣 ========================
    SecretionGeneEntry(
        category=CAT_ER_FOLDING,
        common_name="KAR2 / BiP",
        description="ER 分子伴侣，帮助蛋白正确折叠，**最重要的 OE 靶点之一**",
        oe_reaction_id="sec_Kar2p_complex_formation",
        intervention="OE",
        evidence="已报道 Kar2 过表达可提升毕赤酵母外源蛋白分泌",
    ),
    SecretionGeneEntry(
        category=CAT_ER_FOLDING,
        common_name="BIP/NEFS",
        description="BiP 核苷酸交换因子，调控 Kar2 活性",
        oe_reaction_id="sec_BIP_NEFS_complex_formation",
        intervention="OE",
        evidence="模型 sec_BIP_NEFS 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ER_FOLDING,
        common_name="SSA1",
        description="Hsp70 家族胞质分子伴侣，协助新生链折叠",
        oe_reaction_id="sec_Ssa1_Ydj1_Snl1_complex_formation",
        intervention="OE",
        evidence="模型 sec_Ssa1_Ydj1_Snl1 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ER_FOLDING,
        common_name="YDJ1",
        description="Hsp40 共伴侣，协助 SSA1 识别未折叠蛋白",
        oe_reaction_id="sec_Ssa1_Ydj1_Snl1_complex_formation",
        intervention="OE",
        evidence="模型 sec_Ssa1_Ydj1_Snl1 复合体",
    ),

    # ======================== 二硫键 ========================
    SecretionGeneEntry(
        category=CAT_DSB,
        common_name="PDI1",
        description="蛋白二硫键异构酶，促进正确二硫键形成",
        oe_reaction_id="sec_PDI1_ERV2_Ero1p_complex_formation",
        intervention="OE",
        evidence="已报道 PDI1 过表达可提升含 DSB 蛋白的分泌",
    ),
    SecretionGeneEntry(
        category=CAT_DSB,
        common_name="ERO1",
        description="ER 氧化还原酶，为 PDI1 提供氧化力",
        oe_reaction_id="sec_PDI1_ERV2_Ero1p_complex_formation",
        intervention="OE",
        evidence="模型 sec_PDI1_ERV2_Ero1p 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_DSB,
        common_name="ERV2",
        description="ER 氧化还原酶，与 ERO1 协同",
        oe_reaction_id="sec_PDI1_ERV2_Ero1p_complex_formation",
        intervention="OE",
        evidence="模型 sec_PDI1_ERV2_Ero1p 复合体",
    ),

    # ======================== N-糖基化 ========================
    SecretionGeneEntry(
        category=CAT_N_GLYCAN,
        common_name="OST 复合体",
        description="寡糖转移酶复合体，催化 N-糖基化",
        oe_reaction_id="sec_OSTC_complex_formation",
        intervention="OE",
        evidence="模型 sec_OSTC 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_N_GLYCAN,
        common_name="CWH41",
        description="ER 糖苷酶 I，N-聚糖加工",
        oe_reaction_id="sec_Cwh41p_complex_formation",
        intervention="OE",
        evidence="模型 sec_Cwh41p 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_N_GLYCAN,
        common_name="ROT2",
        description="ER 糖苷酶 II，N-聚糖加工",
        oe_reaction_id="sec_Rot2p_complex_formation",
        intervention="OE",
        evidence="模型 sec_Rot2p 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_N_GLYCAN,
        common_name="MNS1",
        description="ER 甘露糖苷酶 I，N-聚糖修剪",
        oe_reaction_id="sec_Mns1p_complex_formation",
        intervention="OE",
        evidence="模型 sec_Mns1p 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_N_GLYCAN,
        common_name="OCH1",
        description="Golgi α-1,6-甘露糖转移酶，毕赤酵母糖基化工程关键靶点",
        ko_reaction_id="sec_Och1p_complex_formation",
        oe_reaction_id="sec_Och1p_complex_formation",
        intervention="OE",
        evidence="毕赤酵母糖基化工程（人源化）；模型 sec_Och1p 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_N_GLYCAN,
        common_name="MPOLI/MPoLII",
        description="Golgi 甘露糖基转移酶复合体",
        oe_reaction_id="sec_MPOLI_complex_formation",
        intervention="OE",
        evidence="模型 sec_MPOLI/MPoLII 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_N_GLYCAN,
        common_name="MNN2",
        description="Golgi 甘露糖基转移酶",
        oe_reaction_id="sec_Mnn2pA_complex_formation",
        intervention="OE",
        evidence="模型 sec_Mnn2pA/Mnn2pB/Mnn2pC",
    ),

    # ======================== O-糖基化 ========================
    SecretionGeneEntry(
        category=CAT_O_GLYCAN,
        common_name="PMT1/PMT2/PMT4-6",
        description="O-甘露糖转移酶复合体，催化 O-糖基化起始",
        oe_reaction_id="sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_complex_formation",
        intervention="OE",
        evidence="模型 sec_Pmt 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_O_GLYCAN,
        common_name="KTR",
        description="Golgi 甘露糖转移酶，O-聚糖延伸",
        oe_reaction_id="sec_KTR_complex_formation",
        intervention="OE",
        evidence="模型 sec_KTR 复合体",
    ),

    # ======================== COPII 囊泡 ========================
    SecretionGeneEntry(
        category=CAT_COPII,
        common_name="SEC12",
        description="COPII 衣壳装配因子，囊泡出芽起始",
        oe_reaction_id="sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex_formation",
        intervention="OE",
        evidence="模型 COPII 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_COPII,
        common_name="SAR1",
        description="COPII 小 GTP 酶，调控囊泡出芽",
        oe_reaction_id="sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex_formation",
        intervention="OE",
        evidence="已报道 SAR1 过表达可提升分泌",
    ),
    SecretionGeneEntry(
        category=CAT_COPII,
        common_name="SEC23/SEC24",
        description="COPII 衣壳内层，货物选择",
        oe_reaction_id="sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_Bet1p_Bos1p_complex_formation",
        intervention="OE",
        evidence="模型 COPII 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_COPII,
        common_name="SEC13/SEC31",
        description="COPII 衣壳外层，囊泡形成",
        oe_reaction_id="sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex_formation",
        intervention="OE",
        evidence="模型 COPII 外层复合体",
    ),
    SecretionGeneEntry(
        category=CAT_COPII,
        common_name="YPT1/USO1",
        description="小 GTP 酶 + 拴系因子，囊泡与 Golgi 对接",
        oe_reaction_id="sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex_formation",
        intervention="OE",
        evidence="模型 COPII/Golgi 拴系复合体",
    ),
    SecretionGeneEntry(
        category=CAT_COPII,
        common_name="EMP24/ERP",
        description="COPII 货物受体，选择分泌蛋白进入囊泡",
        oe_reaction_id="sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex_formation",
        intervention="OE",
        evidence="模型 GPI-COPII 复合体",
    ),

    # ======================== ERAD ========================
    SecretionGeneEntry(
        category=CAT_ERAD,
        common_name="HRD1/HRD3/DER1",
        description="ERAD E3 连接酶复合体核心，错误折叠蛋白逆向转运",
        ko_reaction_id="sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex_formation",
        oe_reaction_id="sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex_formation",
        intervention="OE",
        evidence="敲除 HRD1 可减少 ERAD，提升外源蛋白积累",
    ),
    SecretionGeneEntry(
        category=CAT_ERAD,
        common_name="UBC6/UBC7",
        description="ERAD 泛素结合酶（模型仅支持复合体级 OE）",
        intervention="OE",
        evidence="模型 ERAD 复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ERAD,
        common_name="CDC48",
        description="AAA-ATP 酶，从 ER 膜提取错误折叠蛋白（模型仅支持复合体级 OE）",
        ko_reaction_id="sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex_formation",
        oe_reaction_id="sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex_formation",
        intervention="OE",
        evidence="模型 ERAD 逆向转运复合体",
    ),
    SecretionGeneEntry(
        category=CAT_ERAD,
        common_name="DOA10",
        description="ERAD E3 连接酶（M 分支），降解 ER 膜蛋白（模型仅支持复合体级 OE）",
        ko_reaction_id="sec_Ubc6p_Ubc7p_Doa10p_complex_formation",
        intervention="OE",
        evidence="模型 ERAD-M 复合体",
    ),

    # ======================== COPI 逆向转运 ========================
    SecretionGeneEntry(
        category=CAT_COPI,
        common_name="ARF1",
        description="COPI 小 GTP 酶，介导 Golgi→ER 逆向转运",
        oe_reaction_id="sec_Arf1p_Gea2p_Rer1p_Erd2p_Cop1p_Sec26p_Sec27p_Sec21p_Ret2p_Sec28p_Ret3p_complex_formation",
        intervention="OE",
        evidence="模型 COPI 复合体",
    ),

    # ======================== 胞吐 ========================
    SecretionGeneEntry(
        category=CAT_EXOCYTOSIS,
        common_name="SEC3/SEC5/SEC6/SEC8/SEC10/SEC15",
        description="胞吐复合体 (exocyst)，囊泡与质膜融合",
        oe_reaction_id="sec_Arf1p_Sec3p_Sec5p_Sec6p_Sec8p_Sec10p_Sec15p_Exo70p_Exo84p_Sec4p_Chc1p_Clc1p_complex_formation",
        intervention="OE",
        evidence="模型 exocyst 复合体",
    ),

    # ======================== 蛋白酶体/降解 ========================
    SecretionGeneEntry(
        category=CAT_PROTEASOME,
        common_name="PEP4",
        description="液泡蛋白酶 A，敲除可减少目标蛋白降解",
        gene_id="PAS_chr2-2_0107",
        intervention="KO",
        evidence="毕赤酵母蛋白表达常用 KO 靶点",
    ),
    SecretionGeneEntry(
        category=CAT_PROTEASOME,
        common_name="PRB1",
        description="液泡蛋白酶 B，与 PEP4 协同降解",
        gene_id="PAS_chr2-1_0785",
        intervention="KO",
        evidence="为默认 KO 候选之一",
    ),
    SecretionGeneEntry(
        category=CAT_PROTEASOME,
        common_name="蛋白酶体",
        description="26S 蛋白酶体复合体，调控 ERAD 后的降解",
        oe_reaction_id="Mach_proteasome_complex_formation",
        intervention="OE",
        evidence="模型蛋白酶体复合体",
    ),

    # ======================== GPI 锚定 ========================
    SecretionGeneEntry(
        category=CAT_GENERAL,
        common_name="GPI 锚定复合体",
        description="GPI 锚定蛋白修饰（对跨膜/GPI 蛋白重要）",
        oe_reaction_id="sec_GPIR_complex_formation",
        intervention="OE",
        evidence="模型 GPI 锚定复合体",
    ),
)


def get_catalog_by_category() -> dict[str, list[SecretionGeneEntry]]:
    """按分泌通路分类返回基因目录。"""
    result: dict[str, list[SecretionGeneEntry]] = {}
    for entry in SECRETION_GENE_CATALOG:
        result.setdefault(entry.category, []).append(entry)
    return result


def search_catalog(query: str = "") -> list[SecretionGeneEntry]:
    """搜索基因目录，返回匹配的条目。"""
    q = query.lower().strip()
    if not q:
        return list(SECRETION_GENE_CATALOG)
    results: list[SecretionGeneEntry] = []
    for entry in SECRETION_GENE_CATALOG:
        if (q in entry.common_name.lower()
                or q in entry.description.lower()
                or q in entry.category.lower()
                or q in entry.gene_id.lower()):
            results.append(entry)
    return results


def get_oe_reactions_for_selection(selected_names: list[str]) -> list[str]:
    """根据选中的常用名返回对应的 OE reaction IDs。"""
    name_set = set(selected_names)
    reactions: list[str] = []
    for entry in SECRETION_GENE_CATALOG:
        if entry.common_name in name_set and entry.oe_reaction_id:
            if entry.oe_reaction_id not in reactions:
                reactions.append(entry.oe_reaction_id)
    return reactions


def get_ko_reactions_for_selection(selected_names: list[str]) -> list[str]:
    """根据选中的常用名返回可直接 KO 的反应 ID（复合体级）。"""
    name_set = set(selected_names)
    reactions: list[str] = []
    for entry in SECRETION_GENE_CATALOG:
        if entry.common_name in name_set and entry.ko_reaction_id:
            if entry.ko_reaction_id not in reactions:
                reactions.append(entry.ko_reaction_id)
    return reactions


def get_ko_genes_for_selection(selected_names: list[str]) -> list[str]:
    """根据选中的常用名返回对应的 KO gene IDs。"""
    name_set = set(selected_names)
    genes: list[str] = []
    for entry in SECRETION_GENE_CATALOG:
        if entry.common_name in name_set and entry.gene_id:
            if entry.gene_id not in genes:
                genes.append(entry.gene_id)
    return genes


_SECRETION_GENE_EVIDENCE_CACHE: list[dict[str, object]] | None = None


def build_secretion_gene_evidence_map(model=None) -> list[dict[str, object]]:
    """Map curated secretion-engineering names to model GPR genes and reaction proxies.

    This is an evidence layer, not a new perturbation algorithm. A curated name such as
    PDI1 can be present as a MATLAB/secretory-pathway reaction proxy while still missing
    from ``model.genes`` and therefore not being executable as a gene-level GPR KO/OE.
    """
    global _SECRETION_GENE_EVIDENCE_CACHE
    use_cache = model is None
    if use_cache and _SECRETION_GENE_EVIDENCE_CACHE is not None:
        return _SECRETION_GENE_EVIDENCE_CACHE
    if model is None:
        from pcsec_pichia.loading import load_pcsec_pichia_inputs, repo_root

        inputs = load_pcsec_pichia_inputs(repo_root())
        model = inputs.prepared_model
    full_genes = load_full_model_genes(model)
    full_genes_by_id = {str(row.get("gene_id") or ""): row for row in full_genes}
    rxns = [str(rxn_id) for rxn_id in getattr(model, "rxns", ())]
    reaction_index = {rxn_id: idx for idx, rxn_id in enumerate(rxns)}
    rules = list(getattr(model, "rules", ()) or ())
    gr_rules = list(getattr(model, "gr_rules", ()) or ())
    rows: list[dict[str, object]] = []
    for entry in SECRETION_GENE_CATALOG:
        model_gene = _model_gene_row_for_curated_entry(entry, full_genes, full_genes_by_id)
        proxy_reactions = [value for value in (entry.ko_reaction_id, entry.oe_reaction_id) if value]
        reaction_evidence = [
            _reaction_proxy_evidence(reaction_id, reaction_index, rules=rules, gr_rules=gr_rules)
            for reaction_id in dict.fromkeys(proxy_reactions)
        ]
        has_gpr_gene = bool(model_gene)
        has_proxy = bool(reaction_evidence)
        proxy_exists = any(bool(item["exists_in_model"]) for item in reaction_evidence)
        proxy_has_gpr_rule = any(bool(item["has_gpr_rule"]) for item in reaction_evidence)
        if has_gpr_gene:
            mapping_status = "model_gpr_gene_available"
            recommended_use = "gene_level_gpr_perturbation"
        elif proxy_exists:
            mapping_status = "reaction_proxy_only"
            recommended_use = "reaction_level_proxy_requires_locus_review"
        elif has_proxy:
            mapping_status = "declared_proxy_missing_in_model"
            recommended_use = "manual_review_required"
        else:
            mapping_status = "literature_name_only"
            recommended_use = "manual_review_required"
        rows.append(
            {
                "common_name": entry.common_name,
                "category": entry.category,
                "description": entry.description,
                "curated_evidence": entry.evidence,
                "homolog_note": entry.homolog_note,
                "declared_model_gene_id": entry.gene_id,
                "mapped_model_gene_id": str(model_gene.get("gene_id") or "") if model_gene else "",
                "mapped_display_name": str(model_gene.get("display_name") or model_gene.get("protein_name") or "") if model_gene else "",
                "mapped_aliases": list(model_gene.get("aliases") or []) if model_gene else [],
                "mapping_status": mapping_status,
                "recommended_use": recommended_use,
                "ko_reaction_id": entry.ko_reaction_id,
                "oe_reaction_id": entry.oe_reaction_id,
                "reaction_evidence": reaction_evidence,
                "proxy_exists_in_model": proxy_exists,
                "proxy_has_gpr_rule": proxy_has_gpr_rule,
                "gene_level_ready": has_gpr_gene,
                "reaction_proxy_ready": proxy_exists,
            }
        )
    if use_cache:
        _SECRETION_GENE_EVIDENCE_CACHE = rows
    return rows


def search_secretion_gene_evidence(query: str = "", model=None) -> list[dict[str, object]]:
    """Search curated secretion gene evidence mapping rows."""
    rows = build_secretion_gene_evidence_map(model)
    q = str(query or "").strip().lower()
    if not q:
        return rows
    return [row for row in rows if _secretion_gene_evidence_matches(row, q)]


def build_lightweight_secretion_gene_evidence(model=None) -> list[dict[str, object]]:
    """Build curated evidence rows without the full 1025-gene capability scan."""
    gene_ids = {str(gene_id) for gene_id in getattr(model, "genes", ())} if model is not None else set()
    reaction_ids = {str(rxn_id) for rxn_id in getattr(model, "rxns", ())} if model is not None else set()
    rows: list[dict[str, object]] = []
    for entry in SECRETION_GENE_CATALOG:
        proxy_reactions = [value for value in (entry.ko_reaction_id, entry.oe_reaction_id) if value]
        proxy_exists = any(reaction_id in reaction_ids for reaction_id in proxy_reactions) if reaction_ids else bool(proxy_reactions)
        has_gpr_gene = bool(entry.gene_id and (not gene_ids or entry.gene_id in gene_ids))
        if has_gpr_gene:
            mapping_status = "model_gpr_gene_available"
            recommended_use = "gene_level_gpr_perturbation"
        elif proxy_exists:
            mapping_status = "reaction_proxy_only"
            recommended_use = "reaction_level_proxy_requires_locus_review"
        elif proxy_reactions:
            mapping_status = "declared_proxy_missing_in_model"
            recommended_use = "manual_review_required"
        else:
            mapping_status = "literature_name_only"
            recommended_use = "manual_review_required"
        rows.append(
            {
                "common_name": entry.common_name,
                "category": entry.category,
                "description": entry.description,
                "curated_evidence": entry.evidence,
                "homolog_note": entry.homolog_note,
                "declared_model_gene_id": entry.gene_id,
                "mapped_model_gene_id": entry.gene_id if has_gpr_gene else "",
                "mapped_display_name": entry.common_name if has_gpr_gene else "",
                "mapped_aliases": (),
                "mapping_status": mapping_status,
                "recommended_use": recommended_use,
                "ko_reaction_id": entry.ko_reaction_id,
                "oe_reaction_id": entry.oe_reaction_id,
                "reaction_evidence": tuple(
                    {
                        "reaction_id": reaction_id,
                        "exists_in_model": reaction_id in reaction_ids if reaction_ids else True,
                        "reaction_index_1based": None,
                        "has_gpr_rule": False,
                        "rule": "",
                        "gr_rule": "",
                    }
                    for reaction_id in dict.fromkeys(proxy_reactions)
                ),
                "proxy_exists_in_model": proxy_exists,
                "proxy_has_gpr_rule": False,
                "gene_level_ready": has_gpr_gene,
                "reaction_proxy_ready": proxy_exists,
            }
        )
    return rows


def _model_gene_row_for_curated_entry(
    entry: SecretionGeneEntry,
    full_genes: list[dict[str, object]],
    full_genes_by_id: dict[str, dict[str, object]],
) -> dict[str, object] | None:
    if entry.gene_id and entry.gene_id in full_genes_by_id:
        return full_genes_by_id[entry.gene_id]
    tokens = _curated_name_tokens(entry.common_name)
    for row in full_genes:
        searchable = {
            str(row.get("gene_id") or ""),
            str(row.get("canonical_gene_id") or ""),
            str(row.get("standard_gene_symbol") or ""),
            str(row.get("display_name") or ""),
            str(row.get("protein_name") or ""),
            str(row.get("ortholog_symbol") or ""),
            *[str(alias) for alias in row.get("aliases") or []],
        }
        lowered = {item.lower() for item in searchable if item}
        if any(token.lower() in lowered for token in tokens):
            return row
    return None


def _curated_name_tokens(common_name: str) -> list[str]:
    raw = str(common_name or "").replace("(", " ").replace(")", " ")
    tokens: list[str] = []
    for chunk in raw.replace("/", " ").replace(",", " ").split():
        token = chunk.strip()
        if token and token not in tokens:
            tokens.append(token)
    return tokens


def _reaction_proxy_evidence(
    reaction_id: str,
    reaction_index: dict[str, int],
    *,
    rules: list[object],
    gr_rules: list[object],
) -> dict[str, object]:
    index = reaction_index.get(reaction_id)
    rule = str(rules[index]) if index is not None and index < len(rules) else ""
    gr_rule = str(gr_rules[index]) if index is not None and index < len(gr_rules) else ""
    return {
        "reaction_id": reaction_id,
        "exists_in_model": index is not None,
        "reaction_index_1based": index + 1 if index is not None else None,
        "has_gpr_rule": bool(rule.strip() or gr_rule.strip()),
        "rule": rule,
        "gr_rule": gr_rule,
    }


def _secretion_gene_evidence_matches(row: dict[str, object], query: str) -> bool:
    reaction_text = " ".join(
        str(item.get("reaction_id") or "")
        for item in row.get("reaction_evidence") or []
        if isinstance(item, dict)
    )
    fields = (
        row.get("common_name"),
        row.get("category"),
        row.get("description"),
        row.get("curated_evidence"),
        row.get("declared_model_gene_id"),
        row.get("mapped_model_gene_id"),
        row.get("mapped_display_name"),
        " ".join(str(alias) for alias in row.get("mapped_aliases") or []),
        row.get("mapping_status"),
        row.get("recommended_use"),
        row.get("ko_reaction_id"),
        row.get("oe_reaction_id"),
        reaction_text,
    )
    return any(query in str(value or "").lower() for value in fields)


_FULL_GENE_CACHE: list[dict[str, object]] | None = None


def clear_full_model_gene_cache() -> None:
    global _FULL_GENE_CACHE
    global _SECRETION_GENE_EVIDENCE_CACHE
    _FULL_GENE_CACHE = None
    _SECRETION_GENE_EVIDENCE_CACHE = None


def load_full_model_genes(
    model=None,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    evidence_cache_path: Path | str | None = None,
) -> list[dict[str, object]]:
    global _FULL_GENE_CACHE
    use_cache = model is None and complex_subunits is None and evidence_cache_path is None
    if use_cache and _FULL_GENE_CACHE is not None:
        return _FULL_GENE_CACHE

    if model is None:
        from pcsec_pichia.loading import load_pcsec_pichia_inputs, repo_root
        inputs = load_pcsec_pichia_inputs(repo_root())
        model = inputs.prepared_model
        complex_subunits = inputs.secretory.complex_subunits
    from pcsec_pichia.screens import classify_secretory_process
    from pcsec_pichia.screens import build_gene_capability_profile
    from pcsec_pichia.screens import reactions_for_gene
    from pcsec_pichia.services.gene_evidence import evidence_for_gene, load_gene_evidence_cache

    evidence_by_gene = load_gene_evidence_cache(evidence_cache_path)

    process_by_reaction = {str(rxn_id): classify_secretory_process(str(rxn_id)) for rxn_id in model.rxns}
    gene_map: dict[str, dict[str, object]] = {}
    for gene_id in model.genes:
        matched = reactions_for_gene(model, str(gene_id))
        gene_map[str(gene_id)] = {
            "reactions": matched,
            "processes": {process_by_reaction.get(str(reaction_id), "unknown") for reaction_id in matched},
        }

    results = []
    for gene_idx, gene_id in enumerate(model.genes):
        data = gene_map[gene_id]
        external_evidence = evidence_for_gene(gene_id, evidence_by_gene)
        capability = build_gene_capability_profile(
            model,
            gene_id,
            complex_subunits=complex_subunits,
            aliases=external_evidence.aliases if external_evidence else (),
        ).to_dict()
        procs: set[str] = data["processes"]
        if procs - {"unknown", "metabolic_or_other"}:
            primary = "分泌相关"
        elif "metabolic_or_other" in procs:
            primary = "代谢"
        else:
            primary = "未分类"
        results.append({
            "gene_id": gene_id,
            "canonical_gene_id": capability["canonical_gene_id"],
            "aliases": capability["aliases"],
            "gene_index": gene_idx + 1,
            "n_reactions": len(data["reactions"]),
            "processes": ", ".join(sorted(procs)) if procs else "unknown",
            "primary_category": primary,
            "sample_reactions": data["reactions"][:5],
            "affected_reactions": capability["affected_reactions"],
            "inactive_reactions_if_ko": capability["inactive_reactions_if_ko"],
            "oe_executable_reactions": capability["oe_executable_reactions"],
            "oe_explain_only_reactions": capability["oe_explain_only_reactions"],
            "gpr_rules": capability["gpr_rules"],
            "gpr_role": capability["gpr_role"],
            "ko_support_status": capability["ko_support_status"],
            "oe_support_status": capability["oe_support_status"],
            "support_reason": capability["support_reason"],
            "missing_information": capability["missing_information"],
            "confidence": capability["confidence"],
            "external_ids": (external_evidence.external_ids or {}) if external_evidence else {},
            "standard_gene_symbol": external_evidence.standard_gene_symbol if external_evidence else "",
            "display_name": external_evidence.display_name if external_evidence else "",
            "protein_name": external_evidence.protein_name if external_evidence else "",
            "function_annotation": external_evidence.function_annotation if external_evidence else "",
            "subcellular_location": external_evidence.subcellular_location if external_evidence else "",
            "ec_numbers": list(external_evidence.ec_numbers) if external_evidence else [],
            "go_terms": list(external_evidence.go_terms) if external_evidence else [],
            "ortholog_symbol": external_evidence.ortholog_symbol if external_evidence else "",
            "wet_lab_readiness": external_evidence.wet_lab_readiness if external_evidence else "model_only_not_experiment_ready",
            "evidence_sources": list(external_evidence.evidence_sources) if external_evidence else [],
            "evidence_confidence": external_evidence.evidence_confidence if external_evidence else "",
            "last_refreshed": external_evidence.last_refreshed if external_evidence else "",
        })

    if use_cache:
        _FULL_GENE_CACHE = results
    return results


def search_full_catalog(query: str, model=None) -> list[dict[str, object]]:
    """搜索全部模型基因，返回匹配的条目。"""
    q = query.lower().strip()
    all_genes = load_full_model_genes(model)
    if not q:
        return all_genes
    return [
        g
        for g in all_genes
        if q in str(g["gene_id"]).lower()
        or q in str(g.get("primary_category") or "").lower()
        or q in " ".join(str(item) for item in g.get("aliases") or ()).lower()
        or q in str(g.get("protein_name") or "").lower()
        or q in str(g.get("function_annotation") or "").lower()
        or q in str(g.get("ko_support_status") or "").lower()
        or q in str(g.get("oe_support_status") or "").lower()
    ]


__all__ = [
    "CAT_ER_TRANSLOCATION",
    "CAT_ER_FOLDING",
    "CAT_DSB",
    "CAT_N_GLYCAN",
    "CAT_O_GLYCAN",
    "CAT_ERAD",
    "CAT_COPII",
    "CAT_GOLGI",
    "CAT_COPI",
    "CAT_EXOCYTOSIS",
    "CAT_PROTEASOME",
    "CAT_GENERAL",
    "SECRETION_GENE_CATALOG",
    "SecretionGeneEntry",
    "get_catalog_by_category",
    "get_ko_genes_for_selection",
    "get_ko_reactions_for_selection",
    "get_oe_reactions_for_selection",
    "build_lightweight_secretion_gene_evidence",
    "build_secretion_gene_evidence_map",
    "clear_full_model_gene_cache",
    "load_full_model_genes",
    "search_catalog",
    "search_secretion_gene_evidence",
    "search_full_catalog",
]
