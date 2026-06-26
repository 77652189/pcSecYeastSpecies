from __future__ import annotations

import pandas as pd
import streamlit as st

from app.adapters.matlab import MatlabAdapter
from app.services.opn import (
    DEFAULT_OPN_CANDIDATE,
    DEFAULT_OPN_PRODUCTION_RATIO,
    OpnCandidateCatalog,
    OpnSimulationService,
)
from app.ui.common import PATHS


def render_matlab_reference() -> None:
    if OpnCandidateCatalog is None or OpnSimulationService is None:
        st.warning("历史 MATLAB OPN 参考验证需要本地 MATLAB + Docker SoPlex。")
        return

    st.markdown(
        """<div class="concept-box">历史 OPN/MATLAB 兼容验证入口。它只用于复查旧 OPN MATLAB LP/SoPlex 参考结果，不属于当前 Python corrected 分泌仿真主流程。</div>""",
        unsafe_allow_html=True,
    )
    catalog = OpnCandidateCatalog(PATHS)
    candidates = [candidate.candidate_id for candidate in catalog.list_candidates()] or [DEFAULT_OPN_CANDIDATE]
    left, right = st.columns([1, 2])
    with left:
        candidate_id = st.selectbox("OPN 候选", candidates, key="matlab_cid")
        mu = st.number_input("μ", 0.01, 0.44, 0.10, 0.01, format="%.2f", key="matlab_mu")
        production_ratio = st.number_input(
            "通量固定值",
            1e-10,
            1e-4,
            DEFAULT_OPN_PRODUCTION_RATIO,
            1e-8,
            format="%.0e",
            key="matlab_pr",
        )
        timeout_seconds = st.number_input("超时", 120, 3600, 300, 60, key="matlab_to")
        st.caption("会调用旧 MATLAB OPN smoke 工作流；不生成新的 hLF artifact，也不改变 Python corrected 结果。")
        run = st.button("运行历史 MATLAB OPN 参考验证", type="primary")
    with right:
        if run:
            service = OpnSimulationService(PATHS, MatlabAdapter())
            with st.spinner("运行中…"):
                result = service.run_candidate_smoke(
                    candidate_id=candidate_id,
                    mu=mu,
                    production_ratio=production_ratio,
                    timeout_seconds=timeout_seconds,
                    engine_mode="matlab",
                )
            st.session_state["last_matlab_ref"] = result.model_dump()
        result_data = st.session_state.get("last_matlab_ref")
        if result_data:
            st.success(result_data["message"]) if result_data["success"] else st.error(result_data["message"])
            st.dataframe(pd.DataFrame([result_data]), use_container_width=True, hide_index=True)
        else:
            latest, summary = OpnSimulationService(PATHS, MatlabAdapter()).latest_candidate_result()
            if latest and summary:
                st.info("最近一次结果：")
                st.dataframe(
                    pd.DataFrame(
                        [
                            {
                                "文件": str(latest),
                                "optimal": summary.optimal,
                                "objective": summary.objective_value,
                            }
                        ]
                    ),
                    use_container_width=True,
                    hide_index=True,
                )


__all__ = ["render_matlab_reference"]
