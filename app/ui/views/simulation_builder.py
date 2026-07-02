from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import streamlit as st

from app.services.pichia_target_catalog_service import list_builtin_target_templates
from app.ui.common import PATHS
from app.ui.views.simulation_display import target_semantics_label


MEDIA_TYPE_LABELS: dict[int, str] = {
    2: "YNB 基础培养基（维生素，无氨基酸）",
    4: "YNB + 核心氨基酸（15 种，默认）",
    5: "YNB + 全氨基酸（20 种）",
}

MEDIA_TYPE_DESCRIPTIONS: dict[int, str] = {
    2: "碳源 + 无机盐 + YNB 维生素 7 种 @ -2.0 mmol·gDCW⁻¹·h⁻¹\n不开放氨基酸补料",
    4: "碳源 + 无机盐 + YNB 维生素\n开放 15 种核心氨基酸 @ -0.08\n不含: Ala, Asn, Cys, Gln, Pro, Ser",
    5: "碳源 + 无机盐 + YNB 维生素\n开放全 20 种氨基酸 @ -0.08",
}


def medium_type_label(media_type: int) -> str:
    return MEDIA_TYPE_LABELS.get(int(media_type), f"未知培养基配方（media_type={media_type}）")


@dataclass(frozen=True)
class TargetBuildFormState:
    build_mode: str
    target_id: str
    target_name: str
    custom_json_path: Path | None
    disulfide_sites: int
    n_glycosylation_sites: int
    o_glycosylation_sites: int
    signal_peptide_sequence: str
    leader_sequence: str
    mature_sequence: str
    enable_ribosome: bool
    enable_misfolding: bool
    enable_cost_slope_compatibility: bool
    cost_slope_medium_compatibility_mode: str
    mu: float
    media_type: int
    carbon_source_id: str


def render_target_build_form() -> TargetBuildFormState:
    build_mode = st.radio(
        "构建模式",
        ["快速选择（内置模板）", "三段式构建（自定义组合）", "自定义 JSON"],
        horizontal=True,
        key="pichia_draft_build_mode",
    )
    target_id, target_name, custom_json_path = "OPN", "", None
    disulfide_sites = n_glycosylation_sites = o_glycosylation_sites = 0
    signal_peptide_sequence = leader_sequence = mature_sequence = ""

    if build_mode == "快速选择（内置模板）":
        templates = {item.target_id: item for item in list_builtin_target_templates(PATHS)}
        choice = st.selectbox(
            "模板",
            list(templates.keys()),
            format_func=lambda key: templates[key].label,
            key="pichia_template",
        )
        target_id = choice
        target_name = choice
        selected_template = templates[choice]
        st.caption(
            f"引导肽 {selected_template.leader_length} aa；"
            f"信号肽 {selected_template.signal_peptide_length} aa；"
            f"成熟链 {selected_template.mature_sequence_length} aa；"
            f"全长 {selected_template.full_sequence_length} aa；"
            f"参数状态：{selected_template.parameter_status}"
        )
        st.caption(
            "目标语义："
            f"{target_semantics_label(selected_template.alignment_target_kind)}；"
            f"序列角色：{target_semantics_label(selected_template.sequence_role)}；"
            f"规范化：{target_semantics_label(selected_template.normalization_mode)}"
        )
        if selected_template.note:
            st.info(selected_template.note)
        if selected_template.target_warning:
            st.warning(selected_template.target_warning)

    elif build_mode == "三段式构建（自定义组合）":
        from app.services.pichia_target_catalog_service import (
            known_leaders,
            known_mature_proteins,
            known_signal_peptides,
        )

        signal_peptides = known_signal_peptides()
        leaders = known_leaders()
        mature_proteins = known_mature_proteins()
        signal_peptide_id = st.selectbox(
            "信号肽",
            list(signal_peptides.keys()),
            format_func=lambda key: signal_peptides[key].get("label", key),
            key="pichia_sp",
        )
        signal_peptide_sequence = str(signal_peptides.get(signal_peptide_id, {}).get("sequence", ""))
        leader_id = st.selectbox(
            "引导肽",
            list(leaders.keys()),
            format_func=lambda key: leaders[key].get("label", key),
            key="pichia_ld",
        )
        leader_sequence = str(leaders.get(leader_id, {}).get("sequence", ""))
        mature_id = st.selectbox(
            "成熟蛋白",
            list(mature_proteins.keys()),
            format_func=lambda key: mature_proteins[key].get("label", key),
            key="pichia_mt",
        )
        mature_info = mature_proteins.get(mature_id, {})
        mature_sequence = str(mature_info.get("sequence", ""))
        disulfide_sites = int(mature_info.get("disulfide_sites", 0))
        n_glycosylation_sites = int(mature_info.get("n_glycosylation_sites", 0))
        o_glycosylation_sites = int(mature_info.get("o_glycosylation_sites", 0))
        st.info(f"全长: {len(signal_peptide_sequence) + len(leader_sequence) + len(mature_sequence)} aa")
        st.caption(
            "自定义组合不会智能推断 PTM；当前 DSB/NG/OG 来自所选成熟蛋白模板："
            f"{disulfide_sites}/{n_glycosylation_sites}/{o_glycosylation_sites}。"
        )
        target_id = f"{signal_peptide_id}_{leader_id}_{mature_id}"
        target_name = target_id

    else:
        custom_json_path = Path(
            st.text_input(
                "自定义 JSON 文件路径",
                value=str(PATHS.repo_root / "local_runs" / "pichia_hlf_opn_probe" / "targets.example.json"),
                key="pichia_json",
            )
        )
        target_id = st.text_input("目标蛋白 ID", value="OPN_CUSTOM", key="pichia_json_target")
        target_name = target_id
        st.warning("自定义 JSON 需要显式提供成熟序列、leader/signal peptide 边界和 DSB/NG/OG 计数；当前不会自动推断。")

    enable_ribosome = st.checkbox("启用核糖体约束", value=True)
    enable_misfolding = st.checkbox("启用错误折叠约束", value=True)
    enable_cost_slope_compatibility = st.checkbox(
        "启用蛋白成本斜率对比（MATLAB 历史路线，可选，较慢）",
        value=False,
        help=(
            "当前默认路线：固定生长率 μ，在 corrected 培养基下最大化目标蛋白分泌通量，"
            "用于估计当前条件下的分泌能力。"
            "历史 MATLAB 成本路线：固定生长率 μ，再固定一组目标蛋白分泌比例，"
            "然后优化葡萄糖摄取反应 Ex_glc_D；通过葡萄糖摄取变化和核糖体通量变化估算 protein cost slope。"
            "打开此项只会额外运行历史成本路线用于对比/解释，不会替换或改变当前默认 corrected pipeline 的数值结果。"
        ),
    )
    cost_slope_medium_compatibility_mode = "corrected"
    if enable_cost_slope_compatibility:
        cost_slope_medium_compatibility_mode = st.selectbox(
            "蛋白成本对比使用的培养基边界",
            ["corrected", "matlab_legacy_cost"],
            format_func=lambda value: {
                "corrected": "Python corrected：使用当前修正后的培养基边界",
                "matlab_legacy_cost": "MATLAB 历史 artifact：仅为旧 Protein_cost_TP 对齐关闭 9 个 exchange",
            }[value],
            help=(
                "只影响上面的蛋白成本斜率对比，不影响默认分泌仿真。"
                "Python corrected 更适合当前模型解释；MATLAB 历史 artifact 用于复现旧 MATLAB 成本分析的培养基边界，"
                "不代表更推荐的默认科学设置。"
            ),
            key="pichia_cost_slope_medium_mode",
        )
    col_mu, col_media, col_carbon = st.columns(3)
    with col_mu:
        mu = st.number_input("μ (h⁻¹)", 0.01, 0.44, 0.10, 0.01, format="%.2f", key="pichia_mu")
    with col_media:
        media_type = int(
            st.selectbox(
                "培养基配方",
                list(MEDIA_TYPE_LABELS),
                index=1,
                format_func=medium_type_label,
                key="pichia_media",
                help="这里显示的是成分名称；内部仍映射到 MATLAB/Python 使用的 media_type 编号，数值行为不变。",
            )
        )
        with st.expander("成分", expanded=False):
            st.code(MEDIA_TYPE_DESCRIPTIONS.get(media_type, ""), language="text")
    with col_carbon:
        carbon_source_id = st.selectbox(
            "碳源",
            ["glucose", "methanol", "glycerol", "glucose_glycerol", "glycerol_methanol"],
            format_func=lambda value: {
                "glucose": "葡萄糖 glucose",
                "methanol": "甲醇 methanol",
                "glycerol": "甘油 glycerol",
                "glucose_glycerol": "葡萄糖 + 甘油",
                "glycerol_methanol": "甘油 + 甲醇",
            }[value],
            key="pichia_carbon_source",
            help="切换模型中允许摄取的主要碳源。葡萄糖是当前默认 corrected 条件；甲醇/甘油为 Python draft 边界配置，仍需按目标场景验证。",
        )

    return TargetBuildFormState(
        build_mode=build_mode,
        target_id=target_id,
        target_name=target_name,
        custom_json_path=custom_json_path,
        disulfide_sites=disulfide_sites,
        n_glycosylation_sites=n_glycosylation_sites,
        o_glycosylation_sites=o_glycosylation_sites,
        signal_peptide_sequence=signal_peptide_sequence,
        leader_sequence=leader_sequence,
        mature_sequence=mature_sequence,
        enable_ribosome=enable_ribosome,
        enable_misfolding=enable_misfolding,
        enable_cost_slope_compatibility=enable_cost_slope_compatibility,
        cost_slope_medium_compatibility_mode=cost_slope_medium_compatibility_mode,
        mu=float(mu),
        media_type=media_type,
        carbon_source_id=str(carbon_source_id),
    )


__all__ = ["TargetBuildFormState", "medium_type_label", "render_target_build_form"]
