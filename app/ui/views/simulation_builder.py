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
    target_is_custom_library: bool
    enable_ribosome: bool
    enable_misfolding: bool
    enable_cost_slope_compatibility: bool
    cost_slope_medium_compatibility_mode: str
    enable_solver_robustness_check: bool
    enable_oe_dose_response: bool
    mu: float
    media_type: int
    carbon_source_id: str


def _render_template_sequences(template) -> None:
    """显示模板各段的实际氨基酸序列（此前只给长度，看不出自己选了什么）。"""
    segments = [
        ("信号肽", template.signal_peptide_sequence),
        ("引导肽", template.leader_sequence),
        ("成熟蛋白", template.mature_sequence),
    ]
    present = [(name, sequence) for name, sequence in segments if sequence]
    if not present:
        st.caption(
            f"引导肽 {template.leader_length} aa；信号肽 {template.signal_peptide_length} aa；"
            f"成熟链 {template.mature_sequence_length} aa；全长 {template.full_sequence_length} aa"
        )
        return
    st.caption(
        f"全长 {template.full_sequence_length} aa　·　信号肽 {template.signal_peptide_length} aa"
        f"　·　引导肽 {template.leader_length} aa　·　成熟链 {template.mature_sequence_length} aa"
        f"　·　参数状态：{template.parameter_status}"
    )
    with st.expander("查看序列", expanded=False):
        for name, sequence in present:
            st.markdown(f"**{name}**（{len(sequence)} aa）")
            st.code(sequence, language="text")


def _render_custom_template_summary(entry: dict) -> None:
    signal_sequence = str(entry.get("signal_peptide_sequence", ""))
    leader_sequence = str(entry.get("leader_sequence", ""))
    mature_sequence = str(entry.get("mature_sequence", ""))
    summary = st.columns(4)
    summary[0].metric("全长", f"{len(signal_sequence) + len(leader_sequence) + len(mature_sequence)} aa")
    summary[1].metric("二硫键", int(entry.get("disulfide_sites", 0) or 0))
    summary[2].metric("N-糖位点", int(entry.get("n_glycosylation_sites", 0) or 0))
    summary[3].metric("O-糖位点", int(entry.get("o_glycosylation_sites", 0) or 0))
    st.caption("自建模板：按你保存的序列直接构建，不走内置目标的参数解析。")
    with st.expander("查看序列", expanded=False):
        for name, sequence in (("信号肽", signal_sequence), ("引导肽", leader_sequence), ("成熟蛋白", mature_sequence)):
            if sequence:
                st.markdown(f"**{name}**（{len(sequence)} aa）")
                st.code(sequence, language="text")


def render_target_selection() -> dict:
    """只负责"要做哪个目标蛋白"。条件与约束见 render_conditions_and_constraints——
    两者分开，好让上层把它们放进各自的步骤标签页，而不是堆在一屏。"""
    # 三段式排在最前：它把"信号肽 + 引导肽 + 成熟蛋白"三段拆开显式选择，研究员看得见自己在组装什么；
    # 快速模板是打包好的整体，适合复现既有目标，故退居其次。
    build_mode = st.radio(
        "构建模式",
        ["三段式构建（自定义组合）", "快速选择（内置模板）", "自定义 JSON"],
        horizontal=True,
        key="pichia_draft_build_mode",
    )
    target_id, target_name, custom_json_path = "OPN", "", None
    disulfide_sites = n_glycosylation_sites = o_glycosylation_sites = 0
    signal_peptide_sequence = leader_sequence = mature_sequence = ""
    target_is_custom_library = False

    if build_mode == "快速选择（内置模板）":
        from app.ui.views.target_library_manager import merged_entries

        builtin_templates = {item.target_id: item for item in list_builtin_target_templates(PATHS)}
        custom_templates = {
            key: value
            for key, value in merged_entries("templates").items()
            if value.get("source") == "custom"
        }
        options = [*builtin_templates.keys(), *custom_templates.keys()]

        def _template_label(key: str) -> str:
            if key in builtin_templates:
                return builtin_templates[key].label
            return f"{custom_templates[key].get('label', key)}（自建）"

        choice = st.selectbox("模板", options, format_func=_template_label, key="pichia_template")
        target_id = choice
        target_name = choice
        if choice in builtin_templates:
            selected_template = builtin_templates[choice]
            # 快速模板此前只给各段长度，看不到实际序列——三段式反而看得见，两种模式信息不对等。
            _render_template_sequences(selected_template)
            if selected_template.note:
                st.info(selected_template.note)
            # 目标语义与历史对照是溯源信息，研究员日常用不到，收进折叠区、别在主流程上制造噪声。
            with st.expander("数据来源与历史对照（溯源信息）", expanded=False):
                st.caption(
                    "目标语义："
                    f"{target_semantics_label(selected_template.alignment_target_kind)}；"
                    f"序列角色：{target_semantics_label(selected_template.sequence_role)}；"
                    f"规范化：{target_semantics_label(selected_template.normalization_mode)}"
                )
                if selected_template.target_warning:
                    st.caption(selected_template.target_warning)
        else:
            # 自建模板没有内置 spec，不能按 id 解析——必须把序列显式带上走 custom_sequence 路径。
            entry = custom_templates[choice]
            target_is_custom_library = True
            signal_peptide_sequence = str(entry.get("signal_peptide_sequence", ""))
            leader_sequence = str(entry.get("leader_sequence", ""))
            mature_sequence = str(entry.get("mature_sequence", ""))
            disulfide_sites = int(entry.get("disulfide_sites", 0) or 0)
            n_glycosylation_sites = int(entry.get("n_glycosylation_sites", 0) or 0)
            o_glycosylation_sites = int(entry.get("o_glycosylation_sites", 0) or 0)
            _render_custom_template_summary(entry)

    elif build_mode == "三段式构建（自定义组合）":
        from app.ui.views.target_library_manager import merged_entries

        signal_peptides = merged_entries("signal_peptides")
        leaders = merged_entries("leaders")
        mature_proteins = merged_entries("mature_proteins")
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
        summary = st.columns(4)
        summary[0].metric("全长", f"{len(signal_peptide_sequence) + len(leader_sequence) + len(mature_sequence)} aa")
        summary[1].metric("二硫键", disulfide_sites)
        summary[2].metric("N-糖位点", n_glycosylation_sites)
        summary[3].metric("O-糖位点", o_glycosylation_sites)
        st.caption("修饰位点数取自所选成熟蛋白条目，不会自动推断——换蛋白请确认这三个数。")
        with st.expander("查看已组装的序列", expanded=False):
            for name, sequence in (
                ("信号肽", signal_peptide_sequence),
                ("引导肽", leader_sequence),
                ("成熟蛋白", mature_sequence),
            ):
                if sequence:
                    st.markdown(f"**{name}**（{len(sequence)} aa）")
                    st.code(sequence, language="text")
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

    st.caption("要新增/修改信号肽、引导肽、蛋白或模板，去侧边栏「序列库与映射管理」。")

    return {
        "build_mode": build_mode,
        "target_id": target_id,
        "target_name": target_name,
        "custom_json_path": custom_json_path,
        "target_is_custom_library": target_is_custom_library,
        "disulfide_sites": disulfide_sites,
        "n_glycosylation_sites": n_glycosylation_sites,
        "o_glycosylation_sites": o_glycosylation_sites,
        "signal_peptide_sequence": signal_peptide_sequence,
        "leader_sequence": leader_sequence,
        "mature_sequence": mature_sequence,
    }


def render_conditions_and_constraints() -> dict:
    """培养条件 + 模型约束与可选分析。与目标选择分开，供上层放进各自的步骤标签页。"""
    st.markdown("**培养条件**")
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
    with st.expander("培养基成分", expanded=False):
        st.code(MEDIA_TYPE_DESCRIPTIONS.get(media_type, ""), language="text")

    st.divider()
    st.markdown("**模型约束与可选分析**")
    col_basic, col_optional = st.columns(2)
    enable_ribosome = col_basic.checkbox("启用核糖体约束", value=True)
    enable_misfolding = col_basic.checkbox("启用错误折叠约束", value=True)
    enable_cost_slope_compatibility = st.checkbox(
        "启用蛋白成本分析（固定生长率+分泌比例网格测算成本斜率，较慢）",
        value=False,
        help=(
            "这是目标蛋白成本分析功能本身：固定生长率 μ，再固定一组目标蛋白分泌比例，"
            "然后优化葡萄糖摄取反应 Ex_glc_D；通过葡萄糖摄取变化和核糖体通量变化的斜率，"
            "估算增加单位分泌量需要多少额外代谢成本（沿用 MATLAB 版的成本估算做法）。"
            "不勾选时不会展示任何蛋白成本分析——没有不用模型求解结果的简化替代版本，"
            "因为那类替代版本不代表真实成本，容易造成误导。"
            "勾选后会多跑一组模型求解（默认2个生长率×5个分泌比例=10次），不影响默认分泌仿真本身的数值结果。"
        ),
    )
    enable_solver_robustness_check = st.checkbox(
        "启用求解器稳健性检查（换 highs-ds/highs-ipm 重解，判断瓶颈归因是否为数值假象，较慢）",
        value=False,
        help=(
            "限制强度这类对偶解在退化最优解处并不唯一，同一个模型换个求解算法可能把瓶颈归到不同资源层。"
            "勾选后会用另外两种求解算法（highs-ds、highs-ipm）各重解一次，比对最靠前的“可过表达缓解”瓶颈是否跨算法一致："
            "一致=结论稳健；翻转=这个瓶颈多半是数值巧合、不是生物学结论。"
            "每多一个算法就多跑一次完整模型求解，不影响默认分泌仿真本身的数值结果，也不改默认求解器。"
        ),
    )
    enable_oe_dose_response = st.checkbox(
        "启用 OE 剂量响应形状（扫描多个过表达倍数，看分泌提升会很快到顶还是持续上升，较慢）",
        value=False,
        help=(
            "默认的过表达筛查只在固定 2× 一个点上测提升，看不出再加大表达量还有没有用。"
            "勾选后会对候选反应扫描一组过表达倍数（默认 1.25/1.5/2/3/5/8×），把分泌响应形状分成四类："
            "饱和型（适度过表达就够，再加收益递减）、线性型（还在涨，值得进一步加大）、"
            "阈值型（要超过某个最小倍数才起效）、无响应（任何倍数都几乎没提升，别过表达这个基因）。"
            "这是相对形状信号，不产出绝对产量或最优倍数；每个倍数都要多跑一次模型求解，不影响默认分泌仿真的数值结果。"
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
    return {
        "enable_ribosome": enable_ribosome,
        "enable_misfolding": enable_misfolding,
        "enable_cost_slope_compatibility": enable_cost_slope_compatibility,
        "cost_slope_medium_compatibility_mode": cost_slope_medium_compatibility_mode,
        "enable_solver_robustness_check": enable_solver_robustness_check,
        "enable_oe_dose_response": enable_oe_dose_response,
        "mu": float(mu),
        "media_type": media_type,
        "carbon_source_id": str(carbon_source_id),
    }


def render_target_build_form() -> TargetBuildFormState:
    """兼容入口：把目标选择与条件设置渲染在一起（不分步时使用）。"""
    return TargetBuildFormState(**render_target_selection(), **render_conditions_and_constraints())


__all__ = [
    "TargetBuildFormState",
    "medium_type_label",
    "render_conditions_and_constraints",
    "render_target_build_form",
    "render_target_selection",
]
