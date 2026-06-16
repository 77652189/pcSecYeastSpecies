from __future__ import annotations

import pandas as pd
import streamlit as st

from app.adapters.matlab import MatlabAdapter
from app.adapters.powershell import PowerShellAdapter
from app.core.i18n import SOPLEX_COLUMN_LABELS, status_label
from app.services.logs import RunLogService
from app.services.simulation import SimulationService
from app.ui.common import PATHS, rename_columns


def render_simulation() -> None:
    st.markdown(
        """
        <div class="concept-box">
        这里运行一个已验证的小规模示例计算：酿酒酵母（SCE）在葡萄糖（glucose）条件下，给定生长速率 mu，
        生成线性规划 LP 文件并调用 Docker 中的 SoPlex 求解器。
        </div>
        """,
        unsafe_allow_html=True,
    )
    left, right = st.columns([1, 2])
    with left:
        mu = st.number_input("生长速率 mu（h^-1）", min_value=0.01, max_value=0.44, value=0.10, step=0.01, format="%.2f")
        timeout = st.number_input("SoPlex 超时秒数", min_value=120, max_value=3600, value=300, step=60)
        st.caption("点击运行后会调用 MATLAB 生成 LP，再调用 Docker SoPlex 求解，可能需要数分钟。")
        run_clicked = st.button("运行示例仿真", type="primary")

    with st.expander("这些输入和输出是什么意思", expanded=True):
        st.markdown(
            """
            - **生长速率 mu：** 模型中固定的细胞生长速度，单位是 h^-1。
            - **线性规划 LP：** MATLAB 根据模型和条件生成的优化问题文件。
            - **SoPlex：** 用来求解 LP 的优化求解器。
            - **objective value：** 求解器返回的目标函数值，用来判断这次优化结果的数值结果；它不是直接的产量单位。
            """
        )

    if run_clicked:
        service = SimulationService(PATHS, MatlabAdapter(), PowerShellAdapter())
        with st.spinner("正在生成 LP 并调用 Docker SoPlex，请稍等..."):
            result = service.run_sce_glucose_smoke(mu=float(mu), timeout_seconds=int(timeout))
        st.session_state["last_simulation_result"] = result.model_dump()

    result_data = st.session_state.get("last_simulation_result")
    with right:
        if result_data:
            st.success(result_data["message"]) if result_data["success"] else st.error(result_data["message"])
            summary = {
                "生长速率 mu": result_data["mu"],
                "目标函数值（objective value）": result_data["objective_value"],
                "LP 文件": result_data["lp_file"],
                "SoPlex 输出文件": result_data["output_file"],
            }
            st.dataframe(pd.DataFrame([summary]), use_container_width=True, hide_index=True)
            if result_data["command_output"]:
                with st.expander("查看命令输出"):
                    st.code(result_data["command_output"][-12000:], language="text")
        else:
            latest, summary = RunLogService(PATHS).latest_soplex_summary()
            if latest and summary:
                st.info("下面显示最近一次 SoPlex 输出摘要。")
                display = pd.DataFrame(
                    [
                        {
                            "输出文件": str(latest.relative_to(PATHS.repo_root)),
                            "optimal_label": status_label(summary.optimal),
                            "objective": summary.objective_value,
                            "status": summary.status_line,
                        }
                    ]
                )
                st.dataframe(rename_columns(display, SOPLEX_COLUMN_LABELS), use_container_width=True, hide_index=True)


