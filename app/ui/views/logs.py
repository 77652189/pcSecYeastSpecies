from __future__ import annotations

import pandas as pd
import streamlit as st

from app.core.i18n import RUN_FILE_COLUMN_LABELS, SOPLEX_COLUMN_LABELS, status_label
from app.services.logs import RunLogService
from app.ui.common import PATHS, rename_columns


def render_logs() -> None:
    st.markdown(
        """
        <div class="concept-box">
        这里汇总本地运行产生的文件。`.lp` 是优化问题输入，`.lp.out` 是 SoPlex 求解输出，`.png` 是图形验证结果。
        </div>
        """,
        unsafe_allow_html=True,
    )
    service = RunLogService(PATHS)
    recent = pd.DataFrame(service.recent_files())
    st.subheader("最近运行文件")
    with st.expander("哪些文件最值得看", expanded=True):
        st.markdown(
            """
            - **`.lp.out`：** 最重要，包含 SoPlex 是否求解成功以及 objective value。
            - **`.lp`：** MATLAB 生成的线性规划输入文件，通常很大，不建议手工阅读。
            - **`.png`：** 第一阶段图形验证输出，可用于确认绘图链路正常。
            """
        )
    if recent.empty:
        st.info("local_runs 目录暂无文件。")
    else:
        recent["修改时间"] = pd.to_datetime(recent["修改时间"], unit="s").dt.strftime("%Y-%m-%d %H:%M:%S")
        st.dataframe(rename_columns(recent, RUN_FILE_COLUMN_LABELS), use_container_width=True, hide_index=True)

    latest, summary = service.latest_soplex_summary()
    if latest and summary:
        st.subheader("最近 SoPlex 求解摘要")
        display = pd.DataFrame(
            [
                {
                    "文件": str(latest.relative_to(PATHS.repo_root)),
                    "optimal_label": status_label(summary.optimal),
                    "objective": summary.objective_value,
                    "status": summary.status_line,
                }
            ]
        )
        st.dataframe(rename_columns(display, SOPLEX_COLUMN_LABELS), use_container_width=True, hide_index=True)
        with st.expander("查看输出文件末尾日志"):
            text = latest.read_text(encoding="utf-8", errors="replace")
            st.code(text[-12000:], language="text")


