from __future__ import annotations

import pandas as pd
import streamlit as st

from app.core.i18n import status_label
from app.services.logs import RunLogService
from app.ui.common import HEALTH_COLUMN_LABELS, PATHS, cached_health, dataset_frame, rename_columns, request_navigation


def render_overview() -> None:
    datasets = dataset_frame()
    health = cached_health()
    items = pd.DataFrame(health["items"])

    st.header("项目总览")
    st.markdown(
        """
        <div class="concept-box">
        这个工具用于查看和验证酵母蛋白分泌模型的计算结果。你可以把它理解成一个“细胞工厂模拟器”：
        选择物种、条件或结果文件后，观察模型预测的生长、代谢通量、蛋白成本和求解状态。
        </div>
        """,
        unsafe_allow_html=True,
    )

    start_col, note_col = st.columns([1, 3])
    if start_col.button("从这里开始 →", type="primary", key="overview_start_here_btn"):
        request_navigation("全基因组KO/OE筛查")
        st.rerun()
    note_col.caption(
        "想知道某个目标蛋白（如 hLF）该**敲除 / 过表达哪些基因**来提升分泌？从「全基因组KO/OE筛查」开始——"
        "那里出一份按相对提升排好序的候选短名单、可导出 CSV；再走到候选核实、证据复核确认。"
    )

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("结果文件", len(datasets))
    col2.metric("结果主题", datasets["category_label"].nunique() if not datasets.empty else 0)
    col3.metric("模型物种", 3)
    col4.metric("最近输出", len(RunLogService(PATHS).recent_files()))

    with st.expander("这个网站能做什么、输入是什么、输出是什么", expanded=True):
        st.markdown(
            """
            **能做什么：** 比较三种酵母在蛋白分泌、温度变化、碳源变化、目标蛋白生产等场景下的模型预测结果。

            **输入是什么：** 当前版本主要读取项目中已经计算好的 `Results/` 文件；仿真验证页额外允许输入一个生长速率 `mu`。

            **输出是什么：** 表格、趋势图、求解器 SoPlex 输出、线性规划 LP 文件，以及最近一次仿真的目标函数值（objective value）。
            """
        )

    with st.expander("部署状态（依赖 / 环境自检 · 研究员一般不用看）", expanded=False):
        if not items.empty:
            items["status_label"] = items["status"].map(status_label)
            status_order = {"ok": 0, "warning": 1, "missing": 2, "error": 3}
            items["sort"] = items["status"].map(status_order).fillna(9)
            display = items.sort_values(["sort", "name"])[["name", "status_label", "detail"]]
            st.dataframe(rename_columns(display, HEALTH_COLUMN_LABELS), width='stretch', hide_index=True)
        st.markdown(
            """
            - **正常：** 依赖、模型文件或结果目录已找到。
            - **缺失/错误：** 相关依赖没有安装、路径不对，或预检脚本运行失败。
            - Windows 原生 `soplex` 可以缺失；本项目当前通过 Docker 里的 SoPlex 求解。
            """
        )

    st.subheader("模型覆盖范围")
    st.dataframe(
        pd.DataFrame(
            [
                {"物种": "酿酒酵母（S. cerevisiae）", "代码": "SCE", "模型": "pcSecYeast", "当前功能": "结果浏览 + 已验证小规模仿真"},
                {"物种": "毕赤酵母（K. phaffii）", "代码": "PPA", "模型": "pcSecPichia", "当前功能": "结果浏览"},
                {"物种": "马克斯克鲁维酵母（K. marxianus）", "代码": "KMX", "模型": "pcSecKmarx", "当前功能": "结果浏览"},
            ]
        ),
        width='stretch',
        hide_index=True,
    )

    if PATHS.phase1_png.exists():
        st.subheader("最近一次图形验证")
        st.image(str(PATHS.phase1_png), width=620)


