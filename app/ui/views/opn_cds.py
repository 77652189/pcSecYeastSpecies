from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from app.services.opn import OPN_SHORTLIST
from app.ui.common import cached_opn_candidates, cached_opn_cds_designs


def render_opn_construct_design_table(designs: pd.DataFrame) -> None:
    st.subheader("模型筛选后的蛋白构建清单")
    if designs.empty:
        st.info("暂时没有可导出的 OPN 构建设计。")
        return

    st.markdown(
        """
        这张表仍然是蛋白层面的构建清单：它把“推荐角色、信号肽序列、成熟 OPN、完整蛋白序列、Kex2 风险”放在一起。
        到这一步还没有做 DNA 密码子优化；PichiaCLM 是下面的下游步骤。
        """
    )
    preview = designs[
        [
            "candidate_id",
            "experimental_role",
            "leader_sequence",
            "signal_peptide_sequence",
            "full_protein_length",
            "contains_alpha_pro_region",
            "kex2_risk",
            "codon_optimization_next",
        ]
    ].rename(
        columns={
            "candidate_id": "候选 ID",
            "experimental_role": "实验角色",
            "leader_sequence": "leader 序列",
            "signal_peptide_sequence": "信号肽序列",
            "full_protein_length": "完整蛋白长度",
            "contains_alpha_pro_region": "含 alpha pro 区",
            "kex2_risk": "Kex2/切割风险",
            "codon_optimization_next": "下一步",
        }
    )
    st.dataframe(preview, use_container_width=True, hide_index=True)

    export = designs.rename(
        columns={
            "candidate_id": "候选ID",
            "experimental_role": "实验角色",
            "recommendation": "推荐级别",
            "leader_sequence": "leader序列",
            "signal_peptide_sequence": "信号肽序列",
            "mature_opn_sequence": "成熟OPN序列",
            "full_protein_sequence": "完整蛋白序列",
            "leader_length": "leader长度",
            "signal_peptide_length": "信号肽长度",
            "mature_opn_length": "成熟OPN长度",
            "full_protein_length": "完整蛋白长度",
            "contains_alpha_pro_region": "是否含alpha_pro区",
            "processing_route": "加工路线",
            "kex2_risk": "Kex2风险",
            "codon_optimization_next": "下一步密码子优化",
            "note": "备注",
        }
    )
    st.download_button(
        "下载实验构建设计 CSV",
        export.to_csv(index=False).encode("utf-8-sig"),
        file_name="OPN_Pichia_signal_peptide_construct_design.csv",
        mime="text/csv",
    )


def render_opn_cds_design_panel() -> None:
    st.subheader("下游步骤：PichiaCLM 生成毕赤酵母 CDS")
    st.markdown(
        """
        只有在上面的 pcSecPichia 分泌模型筛选完成后，才进入这一步。这里调用 PichiaCLM 为首轮候选生成
        毕赤酵母 DNA/CDS 序列，并导出 CSV、XLSX 和 FASTA，供后续载体设计和实验讨论。
        """
    )
    st.info("注意：PichiaCLM 做的是 DNA/CDS 层面的密码子优化，不参与前面的 pcSec 分泌模型评分，也不代表真实发酵产量。")
    with st.expander("这两个项目现在怎么通信？", expanded=True):
        st.markdown(
            """
            - **当前做法：下游函数级调用。** pcSecYeastSpecies 的 Streamlit 页面调用 `CdsDesignService`，服务层再调用 `PichiaClmAdapter`，适配器直接 import 本机的 PichiaCLM 核心模型。
            - **不推荐：Streamlit 调 Streamlit。** 两个网页都偏展示层，互相 HTTP 调用会让错误处理、状态管理和未来上线都变复杂。
            - **后期上线：服务级调用。** PichiaCLM 可以单独提供 FastAPI；pcSecYeastSpecies 只把这个下游适配器从“本地 import”换成“HTTP 请求”，前端和分泌模型逻辑不需要重写。
            """
        )

    candidates = pd.DataFrame(cached_opn_candidates())
    default_ids = [candidate_id for candidate_id in OPN_SHORTLIST if candidate_id in set(candidates["candidate_id"])]
    selected_ids = st.multiselect(
        "选择要生成 CDS 的 OPN 构建",
        candidates["candidate_id"].tolist(),
        default=default_ids,
    )
    left, right = st.columns([1, 1])
    with left:
        per_construct = st.number_input("每个构建生成几个 CDS 候选", min_value=1, max_value=10, value=3, step=1)
    with right:
        seed = st.number_input("随机种子", min_value=1, max_value=9999, value=42, step=1)

    if not selected_ids:
        st.info("请至少选择一个 OPN 构建。")
        return
    if not st.button("生成下游 CDS 候选", type="primary"):
        st.caption("为了避免页面打开时自动加载深度学习模型，点击按钮后才会调用 PichiaCLM。")
        return

    with st.spinner("正在调用 PichiaCLM 生成 CDS，首次加载模型可能需要几十秒..."):
        result = cached_opn_cds_designs(tuple(selected_ids), int(per_construct), int(seed))

    if not result["available"]:
        st.error(result["message"])
        return
    st.success(result["message"])
    records = pd.DataFrame(result["records"])
    if records.empty:
        return

    matched_count = int(records["translation_matches_input"].sum())
    stop_count = int(records["internal_stop_codons"].map(len).sum())
    summary_col1, summary_col2, summary_col3, summary_col4 = st.columns(4)
    summary_col1.metric("CDS 候选记录", len(records))
    summary_col2.metric("回译匹配蛋白", f"{matched_count}/{len(records)}")
    summary_col3.metric("内部终止密码子", stop_count)
    summary_col4.metric("导出文件", "CSV / XLSX / FASTA")

    preview = records[
        [
            "construct_id",
            "experimental_role",
            "recommendation",
            "cds_candidate_rank",
            "recommended_subset",
            "aa_length",
            "cds_length",
            "gc_percent",
            "gc_status",
            "cai_public",
            "quality_status",
            "warnings",
            "restriction_sites",
            "motif_hits",
            "length_multiple_of_three",
            "translation_matches_input",
            "internal_stop_codons",
        ]
    ].rename(
        columns={
            "construct_id": "OPN 构建 ID",
            "experimental_role": "实验角色",
            "recommendation": "推荐级别",
            "cds_candidate_rank": "CDS 候选序号",
            "recommended_subset": "推荐保留",
            "aa_length": "氨基酸长度",
            "cds_length": "CDS 长度 bp",
            "gc_percent": "GC%",
            "gc_status": "GC 状态",
            "cai_public": "CAI（公开毕赤参考）",
            "quality_status": "质控状态",
            "warnings": "提醒数",
            "restriction_sites": "默认/自定义酶切位点数",
            "motif_hits": "不期望 motif 数",
            "length_multiple_of_three": "长度为3倍数",
            "translation_matches_input": "回译匹配蛋白",
            "internal_stop_codons": "内部终止密码子位置",
        }
    )
    st.dataframe(preview, use_container_width=True, hide_index=True)

    selected_record = st.selectbox(
        "查看一条 CDS 序列",
        records.index.tolist(),
        format_func=lambda idx: f"{records.loc[idx, 'construct_id']} / CDS {records.loc[idx, 'cds_candidate_rank']}",
    )
    st.code(records.loc[selected_record, "cds"], language="text")

    st.markdown("**下载下游 CDS 设计文件**")
    col_csv, col_xlsx, col_fasta = st.columns(3)
    download_specs = [
        (col_csv, result.get("csv_file"), "下载 CSV 构建表", "text/csv"),
        (col_xlsx, result.get("xlsx_file"), "下载 XLSX 构建表", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
        (col_fasta, result.get("fasta_file"), "下载 FASTA 序列", "text/plain"),
    ]
    for column, path_value, label, mime in download_specs:
        if not path_value:
            column.warning("文件未生成")
            continue
        path = Path(path_value)
        if path.exists():
            column.download_button(label, path.read_bytes(), file_name=path.name, mime=mime)
        else:
            column.warning(f"未找到 {path.name}")


