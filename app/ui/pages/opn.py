from __future__ import annotations

import pandas as pd
import streamlit as st

from app.adapters.matlab import MatlabAdapter
from app.core.i18n import status_label
from app.services.opn import (
    DEFAULT_OPN_CANDIDATE,
    DEFAULT_OPN_PRODUCTION_RATIO,
    OpnCandidateCatalog,
    OpnSimulationService,
    opn_category_label,
)
from app.ui.common import (
    PATHS,
    cached_opn_candidates,
    cached_opn_construct_designs,
    cached_opn_rankings,
)
from app.ui.pages.opn_cds import render_opn_cds_design_panel, render_opn_construct_design_table
from app.ui.pages.opn_library import render_opn_workflow_overview, render_signal_peptide_library_manager
from app.ui.pages.opn_results import (
    render_opn_latest_result_explanation,
    render_opn_recommendation_board,
    render_opn_result_explanation,
)


def render_opn_signal_peptides() -> None:
    st.markdown(
        """
        <div class="concept-box">
        这里分成两个清楚的步骤：先用 pcSecPichia 在蛋白分泌层面比较候选信号肽，再把已经选中的蛋白构建交给
        PichiaCLM 做下游 DNA/CDS 设计。PichiaCLM 不参与分泌模型评分，它只是筛选之后的密码子优化工具。
        </div>
        """,
        unsafe_allow_html=True,
    )
    candidates = pd.DataFrame(cached_opn_candidates())
    if candidates.empty:
        st.error("没有找到 OPN 候选表。请先生成 Data/pcSecPichia/TargetProtein_OPN_candidates_meta.csv。")
        return

    candidates["分类"] = candidates["category"].map(opn_category_label)
    candidates["leader长度aa"] = candidates["leader_length"]
    candidates["完整构建长度aa"] = candidates["construct_length"]

    metric1, metric2, metric3, metric4 = st.columns(4)
    metric1.metric("候选数量", len(candidates))
    metric2.metric("成熟 OPN 长度", "298 aa")
    metric3.metric("默认候选", DEFAULT_OPN_CANDIDATE)
    metric4.metric("默认产量约束", f"{DEFAULT_OPN_PRODUCTION_RATIO:.0e}")

    render_opn_workflow_overview()
    render_opn_recommendation_board(pd.DataFrame(cached_opn_rankings()))

    with st.expander("这个页面怎么看", expanded=True):
        st.markdown(
            """
            - **第一步：分泌模型筛选。** 比较候选 leader 接到成熟 OPN 后，在 pcSecPichia 中的模型可行性和资源成本。
            - **第二步：候选解释与验证。** 查看每个 leader 的来源、风险，并可运行一个小规模 LP/SoPlex 验证。
            - **第三步：下游 CDS 设计。** 只对已经决定进入首轮实验的候选调用 PichiaCLM，生成毕赤酵母 DNA/CDS 序列。
            """
        )

    display = candidates[
        [
            "candidate_id",
            "分类",
            "leader长度aa",
            "完整构建长度aa",
            "leader_sequence",
            "processing_route",
        ]
    ].rename(
        columns={
            "candidate_id": "候选 ID",
            "leader_sequence": "leader 序列",
            "processing_route": "加工路线",
        }
    )
    st.subheader("候选信号肽表")
    st.dataframe(display, use_container_width=True, hide_index=True)
    render_signal_peptide_library_manager()

    default_index = candidates.index[candidates["candidate_id"] == DEFAULT_OPN_CANDIDATE]
    selected_index = int(default_index[0]) if len(default_index) else 0
    selected_id = st.selectbox(
        "选择一个候选查看详情",
        candidates["candidate_id"].tolist(),
        index=selected_index,
    )
    selected = candidates.loc[candidates["candidate_id"] == selected_id].iloc[0]

    left, right = st.columns([1, 1])
    with left:
        st.subheader("候选详情")
        st.markdown(f"**分类：** {opn_category_label(selected['category'])}")
        st.markdown(f"**leader 长度：** {selected['leader_length']} aa")
        st.markdown(f"**信号肽长度：** {len(selected['signal_peptide_sequence'])} aa")
        st.markdown(f"**完整构建长度：** {selected['construct_length']} aa")
        st.markdown(f"**选择理由：** {selected['rationale']}")
        st.markdown(f"**注意事项：** {selected['caution']}")
        with st.expander("查看序列"):
            st.code(selected["leader_sequence"], language="text")
    with right:
        st.subheader("模型会看什么")
        st.markdown(
            """
            模型主要会把 leader/OPN 的序列长度、氨基酸组成、ER 分泌路径、翻译和分泌机器资源约束写入 LP。
            它不会直接预测真实滴度，也不会模拟宿主蛋白酶是否切割 OPN 内部位点。
            """
        )
        st.warning("成熟 OPN 内部检测到 `RR` 和 `KR` 二碱性位点；使用 alpha-factor pro 路线时需要额外关注 Kex2 类异常切割风险。")
        latest, summary = OpnSimulationService(PATHS, MatlabAdapter()).latest_candidate_result(selected_id)
        if latest and summary:
            st.info(
                f"最近一次验证：{status_label(summary.optimal)}。"
                f"目标函数值 {summary.objective_value or '未读取到'}，输出文件 {latest.name}。"
            )

    st.subheader("运行 OPN 小规模验证")
    with st.form("opn_candidate_smoke_form"):
        form_left, form_right = st.columns([1, 1])
        with form_left:
            run_candidate_id = st.selectbox(
                "要运行的候选",
                candidates["candidate_id"].tolist(),
                index=candidates["candidate_id"].tolist().index(selected_id),
                key="opn_run_candidate_id",
            )
            mu = st.number_input("生长速率 mu（h^-1）", min_value=0.01, max_value=0.44, value=0.10, step=0.01, format="%.2f")
        with form_right:
            production_ratio = st.selectbox(
                "OPN 生产通量约束",
                [1e-8, 1e-6],
                index=0,
                format_func=lambda value: f"{value:.0e}",
            )
            timeout = st.number_input("SoPlex 超时秒数", min_value=120, max_value=3600, value=300, step=60)
        st.caption("运行会调用 MATLAB R2020b+ 生成 LP，再调用 Docker 镜像 `pcsec-soplex:24.04` 求解，通常需要 1-3 分钟。")
        run_clicked = st.form_submit_button("运行 OPN 候选验证", type="primary")

    if run_clicked:
        service = OpnSimulationService(PATHS, MatlabAdapter())
        with st.spinner("正在生成 OPN/Pichia LP 并调用 SoPlex，请稍等..."):
            result = service.run_candidate_smoke(
                candidate_id=run_candidate_id,
                mu=float(mu),
                production_ratio=float(production_ratio),
                timeout_seconds=int(timeout),
            )
        st.session_state["last_opn_result"] = result.model_dump()

    result_data = st.session_state.get("last_opn_result")
    if result_data:
        render_opn_result_explanation(result_data)
    else:
        latest, summary = OpnSimulationService(PATHS, MatlabAdapter()).latest_candidate_result(selected_id)
        if latest and summary:
            render_opn_latest_result_explanation(selected_id, latest, summary)

    if PATHS.opn_candidate_csv.exists():
        st.download_button(
            "下载模型候选 CSV",
            PATHS.opn_candidate_csv.read_bytes(),
            file_name=PATHS.opn_candidate_csv.name,
            mime="text/csv",
        )

    render_opn_construct_design_table(pd.DataFrame(cached_opn_construct_designs()))
    render_opn_cds_design_panel()


