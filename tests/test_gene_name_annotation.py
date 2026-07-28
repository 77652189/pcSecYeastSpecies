"""基因名显示层的诚实边界测试。

重点不是"名字能显示出来"，而是**不该显示的时候确实不显示**：非精确匹配档不给名字、
非基因候选不硬塞基因名、身份待复核位点必须打标、位点号永远在场。
"""

from __future__ import annotations

import pandas as pd

from app.services import gene_name_annotation as gene_names
from pcsec_pichia.services.gene_id_standardization import (
    PichiaGeneIdStandardName,
    build_standard_name_lookup,
)


def _lookup() -> dict[str, PichiaGeneIdStandardName]:
    return build_standard_name_lookup(
        [
            # 精确命中 + 有正式符号
            PichiaGeneIdStandardName(
                gene_id="PAS_chr1-4_0187",
                display_name="SEC11",
                standard_symbol="SEC11",
                protein_name="Signal peptidase complex catalytic subunit SEC11",
                external_ids={"uniprot": "C4QZZ0", "kegg": "ppa:PAS_chr1-4_0187"},
                annotation_sources=("UniProt", "KEGG"),
                annotation_confidence="high_exact_locus_tag",
            ),
            # 精确命中、无符号、泛化注释
            PichiaGeneIdStandardName(
                gene_id="PAS_chr3_0199",
                display_name="Glutamine amidotransferase type-2 domain-containing protein",
                protein_name="Glutamine amidotransferase type-2 domain-containing protein",
                external_ids={"uniprot": "C4R3U2"},
                annotation_sources=("UniProt",),
                annotation_confidence="high_exact_locus_tag",
            ),
            # 仅模型内部，无外部注释 —— 不允许显示任何名字
            PichiaGeneIdStandardName(
                gene_id="PAS_c034_0014",
                display_name="PAS_c034_0014",
                annotation_sources=("model_only",),
                annotation_confidence="low_model_only",
            ),
            # 身份待复核位点（策展俗名 PEP4，注释却是别的酶）
            PichiaGeneIdStandardName(
                gene_id="PAS_chr2-2_0107",
                display_name="Palmitoyl-protein thioesterase 1",
                protein_name="Palmitoyl-protein thioesterase 1",
                external_ids={"uniprot": "C4R2U6"},
                annotation_sources=("UniProt",),
                annotation_confidence="high_exact_locus_tag",
            ),
        ]
    )


def test_exact_locus_tag_hit_surfaces_name_with_accession() -> None:
    fields = gene_names.annotate_gene_name_fields("PAS_chr1-4_0187", _lookup())

    assert fields["standard_symbol"] == "SEC11"
    assert fields["annotation_tier"] == "verified_symbol"
    assert fields["annotation_accession"] == "UniProt:C4QZZ0"
    assert fields["identity_review"] == ""


def test_generic_auto_annotation_is_flagged_not_sold_as_function() -> None:
    fields = gene_names.annotate_gene_name_fields("PAS_chr3_0199", _lookup())

    assert fields["protein_name"].startswith("Glutamine amidotransferase")
    # "domain-containing protein" 属自动注释，必须降档，不能与已验证功能名混为一谈
    assert fields["annotation_tier"] == "generic_annotation"


def test_non_exact_tier_never_shows_a_name() -> None:
    """核心门控：查不到精确匹配就不给名字（宁可显示位点号，也不给来路不明的名字）。"""
    fields = gene_names.annotate_gene_name_fields("PAS_c034_0014", _lookup())

    assert fields["gene_display_name"] == ""
    assert fields["standard_symbol"] == ""
    assert fields["protein_name"] == ""
    assert fields["annotation_tier"] == "no_annotation"
    assert fields["annotation_accession"] == ""


def test_unknown_gene_id_degrades_without_inventing_a_name() -> None:
    fields = gene_names.annotate_gene_name_fields("PAS_chr9_9999", _lookup())

    assert fields["gene_display_name"] == ""
    assert fields["annotation_tier"] == "no_annotation"
    assert fields["standard_name_status"] == "missing_standard_name"


def test_curated_identity_unverified_locus_is_flagged() -> None:
    """PEP4/PRB1 类：名字忠实于位点，但"这个位点就是 PEP4"未坐实——必须打标让 UI 警告。"""
    fields = gene_names.annotate_gene_name_fields("PAS_chr2-2_0107", _lookup())

    assert fields["identity_review"] == "curated_identity_unverified"
    assert "PEP4" in gene_names.UNVERIFIED_IDENTITY_LOCI["PAS_chr2-2_0107"]
    # 显示的仍是该位点真实注释，而不是俗名 PEP4
    assert "PEP4" not in fields["protein_name"]


def test_reaction_level_candidate_gets_no_gene_name() -> None:
    fields = gene_names.annotate_gene_name_fields(
        "sec_Pdi1p_complex_formation",
        _lookup(),
        candidate_kind="catalog_reaction",
        intervention_type="OE_reaction",
    )

    assert fields["gene_display_name"] == ""
    assert fields["standard_name_status"] == "not_gene_candidate"


def test_annotate_screen_frame_fills_old_csv_and_keeps_locus_visible() -> None:
    """旧筛查 CSV 只有 gene_id；补名字后位点号仍须在标签里可见。"""
    frame = pd.DataFrame(
        [
            {"gene_id": "PAS_chr1-4_0187", "candidate_kind": "gene", "intervention_type": "OE"},
            {"gene_id": "PAS_c034_0014", "candidate_kind": "gene", "intervention_type": "KO"},
        ]
    )

    annotated = gene_names.annotate_screen_frame(frame, _lookup())

    assert annotated.loc[0, "standard_symbol"] == "SEC11"
    assert annotated.loc[1, "gene_display_name"] == ""
    label = gene_names.safe_gene_display_label(annotated.iloc[0])
    assert "SEC11" in label and "PAS_chr1-4_0187" in label
    # 无注释的候选只显示位点号，不编名字
    assert gene_names.safe_gene_display_label(annotated.iloc[1]) == "PAS_c034_0014"


def test_annotate_screen_frame_preserves_names_already_in_new_csv() -> None:
    """新 CSV 自带的名字不被覆盖（那是那次运行记录下来的值），但派生列要补上。"""
    frame = pd.DataFrame(
        [
            {
                "gene_id": "PAS_chr2-2_0107",
                "candidate_kind": "gene",
                "intervention_type": "KO",
                "gene_display_name": "运行时记录名",
                "standard_symbol": "",
                "protein_name": "Palmitoyl-protein thioesterase 1",
                "annotation_confidence": "high_exact_locus_tag",
                "standard_name_status": "annotated",
            }
        ]
    )

    annotated = gene_names.annotate_screen_frame(frame, _lookup())

    assert annotated.loc[0, "gene_display_name"] == "运行时记录名"
    assert annotated.loc[0, "annotation_tier"] == "descriptive_annotation"
    assert annotated.loc[0, "identity_review"] == "curated_identity_unverified"


def test_annotate_screen_frame_degrades_when_cache_missing() -> None:
    """命名缓存不可用时只补空列（UI 退回纯位点号显示），不得抛错打断结果页。"""
    frame = pd.DataFrame([{"gene_id": "PAS_chr1-4_0187", "candidate_kind": "gene"}])

    annotated = gene_names.annotate_screen_frame(frame, {})

    for column in gene_names.NAME_ANNOTATION_COLUMNS:
        assert column in annotated.columns
    assert annotated.loc[0, "gene_display_name"] == ""
    assert gene_names.safe_gene_display_label(annotated.iloc[0]) == "PAS_chr1-4_0187"


def test_curated_candidate_label_keeps_human_name_without_model_id_suffix() -> None:
    """策展候选（如 PDI1/ERO1）的可读名本身就是身份，后缀模型反应名只会变噪声；
    位点号后缀只用于**数据库注释**来的名字（研究员需要它核对算的是哪个位点）。"""
    curated = pd.Series(
        {"gene_id": "sec_Pdi1p_Ero1p_complex_formation", "common_name": "PDI1/ERO1", "standard_symbol": "", "gene_display_name": "", "protein_name": ""}
    )
    annotated_gene = pd.Series(
        {"gene_id": "PAS_chr1-4_0187", "common_name": "", "standard_symbol": "SEC11", "gene_display_name": "SEC11", "protein_name": ""}
    )

    assert gene_names.safe_gene_display_label(curated) == "PDI1/ERO1"
    assert gene_names.safe_gene_display_label(annotated_gene) == "SEC11（PAS_chr1-4_0187）"


def test_long_annotation_name_is_truncated_but_locus_stays_visible() -> None:
    """UniProt 会把整段功能描述当名字（库里最长 239 字符）；不截断会撑爆柱状图 Y 轴。
    截断后位点号必须仍然在场——研究员靠它认模型实际算的对象。"""
    long_name = "Arginine biosynthesis bifunctional protein ArgJ, mitochondrial [Cleaved into: fragment]"
    row = pd.Series(
        {"gene_id": "PAS_chr3_0176", "standard_symbol": "", "gene_display_name": long_name, "protein_name": "", "common_name": ""}
    )

    label = gene_names.safe_gene_display_label(row)

    assert label.endswith("（PAS_chr3_0176）")
    assert "…" in label
    assert len(label) < len(long_name)
    # 不截断时能拿到完整名字（表格列仍显示全名）
    assert gene_names.safe_gene_display_label(row, max_name_chars=0).startswith(long_name)


def test_annotation_display_version_is_wired_into_screen_caches() -> None:
    """嵌套 st.cache_data 的坑：只改富集逻辑不会让下游短名单缓存失效，页面会继续显示旧名字
    （2026-07-28 实测踩到）。四个缓存函数都必须收 name_version 参数。"""
    import inspect

    from app.ui.views import genome_wide_screen as view

    for function in (
        view._cached_tradeoff_frame,
        view._cached_single_target,
        view._cached_divergence,
        view._cached_shortlist_readout,
    ):
        target = getattr(function, "__wrapped__", function)
        assert "name_version" in inspect.signature(target).parameters, function


def test_classify_annotation_tier_boundaries() -> None:
    assert gene_names.classify_annotation_tier("SEC11", "anything") == "verified_symbol"
    assert gene_names.classify_annotation_tier("", "Phosphoenolpyruvate carboxykinase") == "descriptive_annotation"
    assert gene_names.classify_annotation_tier("", "Uncharacterized protein") == "generic_annotation"
    assert gene_names.classify_annotation_tier("", "STAS domain-containing protein") == "generic_annotation"
    assert gene_names.classify_annotation_tier("", "") == "no_annotation"
