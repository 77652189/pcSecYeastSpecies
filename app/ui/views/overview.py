from __future__ import annotations

import pandas as pd
import streamlit as st

from app.core.i18n import status_label
from app.services.logs import RunLogService
from app.ui.common import HEALTH_COLUMN_LABELS, PATHS, cached_health, dataset_frame, rename_columns


def render_overview() -> None:
    datasets = dataset_frame()
    health = cached_health()
    items = pd.DataFrame(health["items"])

    st.markdown(
        """
        <div class="concept-box">
        这个工具用于查看和验证酵母蛋白分泌模型的计算结果。你可以把它理解成一个“细胞工厂模拟器”：
        选择物种、条件或结果文件后，观察模型预测的生长、代谢通量、蛋白成本和求解状态。
        </div>
        """,
        unsafe_allow_html=True,
    )

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("可浏览结果文件", len(datasets))
    col2.metric("结果主题数量", datasets["category_label"].nunique() if not datasets.empty else 0)
    col3.metric("模型物种数量", 3)
    col4.metric("最近输出文件", len(RunLogService(PATHS).recent_files()))

    with st.expander("这个网站能做什么、输入是什么、输出是什么", expanded=True):
        st.markdown(
            """
            **能做什么：** 比较三种酵母在蛋白分泌、温度变化、碳源变化、目标蛋白生产等场景下的模型预测结果。

            **输入是什么：** 当前版本主要读取项目中已经计算好的 `Results/` 文件；仿真验证页额外允许输入一个生长速率 `mu`。

            **输出是什么：** 表格、趋势图、求解器 SoPlex 输出、线性规划 LP 文件，以及最近一次仿真的目标函数值（objective value）。
            """
        )

    st.subheader("部署状态")
    if not items.empty:
        items["status_label"] = items["status"].map(status_label)
        status_order = {"ok": 0, "warning": 1, "missing": 2, "error": 3}
        items["sort"] = items["status"].map(status_order).fillna(9)
        display = items.sort_values(["sort", "name"])[["name", "status_label", "detail"]]
        st.dataframe(rename_columns(display, HEALTH_COLUMN_LABELS), use_container_width=True, hide_index=True)
    with st.expander("怎么看部署状态"):
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
        use_container_width=True,
        hide_index=True,
    )

    if PATHS.phase1_png.exists():
        st.subheader("最近一次图形验证")
        st.image(str(PATHS.phase1_png), width=620)


