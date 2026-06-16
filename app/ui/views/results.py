from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from app.core.i18n import DATASET_COLUMN_LABELS, category_label, short_path, species_label
from app.services.results import PlotBuilder, ResultCatalog, ResultLoader
from app.ui.common import PATHS, cached_loaded_dataset, dataset_frame, rename_columns


def render_results_browser() -> None:
    datasets = dataset_frame()
    st.markdown(
        """
        <div class="concept-box">
        在这里浏览项目已经计算好的结果文件。建议先按“结果主题”和“物种”筛选，再打开一个数据集查看表格和图。
        </div>
        """,
        unsafe_allow_html=True,
    )
    if datasets.empty:
        st.warning("没有找到 Results 目录中的 Excel 或 MATLAB 结果文件。")
        return

    with st.sidebar:
        st.subheader("结果筛选")
        species_options = ["全部", *sorted(datasets["species_label"].unique())]
        category_options = ["全部", *sorted(datasets["category_label"].unique())]
        species = st.selectbox("物种", species_options, index=0)
        category = st.selectbox("结果主题", category_options, index=0)
        keyword = st.text_input("关键词", placeholder="例如 Insulin、CSource、SCE")

    filtered = datasets.copy()
    if species != "全部":
        filtered = filtered[filtered["species_label"] == species]
    if category != "全部":
        filtered = filtered[filtered["category_label"] == category]
    if keyword:
        mask = filtered["name"].str.contains(keyword, case=False, na=False) | filtered["id"].str.contains(keyword, case=False, na=False)
        filtered = filtered[mask]

    with st.expander("如何选择数据集", expanded=True):
        st.markdown(
            """
            - **结果主题** 对应原项目 `Results/` 下的文件夹，例如“碳源分析”“目标蛋白生长分析”。
            - **物种** 是从文件名推断的；如果显示“未知物种”，通常说明文件名没有明确物种缩写。
            - **Excel 结果表** 更适合直接阅读；**MATLAB 结果文件** 通常包含矩阵或结构体，页面会先显示变量摘要。
            """
        )

    st.subheader("结果数据集")
    display_columns = ["name", "category_label", "species_label", "suffix", "size_kb", "modified_at", "id"]
    st.dataframe(
        rename_columns(filtered[display_columns], DATASET_COLUMN_LABELS),
        use_container_width=True,
        hide_index=True,
    )
    if filtered.empty:
        st.info("当前筛选条件下没有结果文件。")
        return

    selected_id = st.selectbox(
        "打开数据集",
        filtered["id"].tolist(),
        format_func=lambda value: filtered.loc[filtered["id"] == value, "name"].iloc[0],
    )
    loaded_data = cached_loaded_dataset(selected_id)
    loaded = ResultLoader().load_dataset(ResultCatalog(PATHS).get_dataset(selected_id))

    st.subheader(loaded.info.name)
    st.caption(f"{category_label(loaded.info.category)} / {species_label(loaded.info.species)} / {short_path(selected_id)}")

    if loaded.variable_summary:
        st.markdown("**MATLAB 变量摘要**")
        st.dataframe(pd.DataFrame(loaded.variable_summary), use_container_width=True, hide_index=True)

    table_names = list(loaded.tables)
    if table_names:
        table_name = st.selectbox("表格或变量", table_names)
        table = pd.DataFrame(loaded_data["tables"][table_name])
        st.dataframe(table, use_container_width=True, hide_index=True)

        with st.expander("图怎么看"):
            st.markdown(
                """
                图表是对当前表格中数值列的快速预览。对于碳源分析，散点图会比较实验测得生长速率和模型预测生长速率；
                对于一般矩阵结果，页面会展示前几个数值变量的趋势，帮助你判断数据是否可读、是否需要下载原始文件进一步分析。
                """
            )
        figure = PlotBuilder().build_chart(loaded, table_name)
        if figure is not None:
            st.plotly_chart(figure, use_container_width=True)
        else:
            st.info("这个数据集没有适合自动绘图的数值列，建议下载原始文件进一步查看。")

    file_path = Path(loaded_data["info"]["path"])
    if file_path.exists():
        st.caption("下载的是原始结果文件，不是重新计算生成的报告。")
        st.download_button(
            "下载原始结果文件",
            file_path.read_bytes(),
            file_name=file_path.name,
            mime="application/octet-stream",
        )


