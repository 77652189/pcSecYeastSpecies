from __future__ import annotations

import pandas as pd
import streamlit as st

from app.services.signal_peptide_screening import SignalPeptideScreeningService
from app.ui.common import (
    PATHS,
    _download_file_button,
    cached_signal_peptide_library,
    cached_uniprot_signal_peptides,
    signal_peptide_library_service,
)


def render_opn_workflow_overview() -> None:
    st.subheader("推荐工作流")
    step1, step2, step3 = st.columns(3)
    step1.markdown(
        """
        **1. pcSecPichia 分泌模型筛选**

        输入是 `leader + 成熟 OPN` 的氨基酸序列。这里比较的是蛋白分泌负担、ER 路径和 LP 求解可行性。
        """
    )
    step2.markdown(
        """
        **2. 确定首轮实验候选**

        保留 PAS_chr3_0030、DDDK18 和 alpha-factor 对照。这个决策来自模型参考、文献证据和加工风险。
        """
    )
    step3.markdown(
        """
        **3. PichiaCLM 下游 CDS 设计**

        输入仍是已选蛋白构建的氨基酸序列，输出 DNA/CDS。它不改变 pcSec 的分泌模型结果。
        """
    )


def render_signal_peptide_library_manager(*, include_external: bool = True) -> None:
    st.subheader("信号肽候选库管理")
    st.markdown(
        """
        这里管理的是“候选库”，不是最终实验清单。新增候选会先作为草案通过序列和字段检查，
        再由项目维护者审核后写入模型输入表；这样可以避免网页误改 Git 中的正式数据。
        """
    )
    service = signal_peptide_library_service()
    library = pd.DataFrame(cached_signal_peptide_library())
    if library.empty:
        st.info("候选库为空。")
        return

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("正式候选数", len(library))
    col2.metric("首轮推荐", int((library["library_stage"] == "首轮推荐").sum()))
    col3.metric("Pichia 来源", int((library["category"] == "pichia_native_signal").sum()))
    col4.metric("未入首轮候选", int((library["library_stage"] == "候选库").sum()))

    with st.expander("筛选和查看候选库", expanded=True):
        filter_left, filter_right = st.columns([1, 1])
        with filter_left:
            stage_options = sorted(library["library_stage"].unique())
            selected_stages = st.multiselect("候选阶段", stage_options, default=stage_options)
        with filter_right:
            category_options = sorted(library["category_label"].unique())
            selected_categories = st.multiselect("候选来源/类别", category_options, default=category_options)
        keyword = st.text_input("按 ID、来源说明或理由搜索", placeholder="例如 PASCHR3、DDDK、UniProt")

        filtered = library[
            library["library_stage"].isin(selected_stages)
            & library["category_label"].isin(selected_categories)
        ].copy()
        if keyword:
            keyword_mask = (
                filtered["candidate_id"].str.contains(keyword, case=False, na=False)
                | filtered["source_note"].str.contains(keyword, case=False, na=False)
                | filtered["rationale"].str.contains(keyword, case=False, na=False)
            )
            filtered = filtered[keyword_mask]
        st.dataframe(
            filtered[
                [
                    "candidate_id",
                    "library_stage",
                    "category_label",
                    "source_type",
                    "leader_length",
                    "signal_peptide_length",
                    "processing_route",
                    "source_note",
                    "rationale",
                ]
            ].rename(
                columns={
                    "candidate_id": "候选 ID",
                    "library_stage": "候选阶段",
                    "category_label": "类别",
                    "source_type": "来源类型",
                    "leader_length": "leader 长度",
                    "signal_peptide_length": "信号肽长度",
                    "processing_route": "加工路线",
                    "source_note": "来源说明",
                    "rationale": "入库理由",
                }
            ),
            use_container_width=True,
            hide_index=True,
        )

    with st.expander("新增候选草案", expanded=False):
        st.markdown(
            """
            下载模板后填入新候选，再上传 CSV 做预检。预检通过后，页面会生成一个“合并草案 CSV”下载；
            这个草案还不会自动进入正式模型，需要人工审核来源、切割位点和文献证据。
            """
        )
        st.download_button(
            "下载新增候选模板 CSV",
            service.template_csv(),
            file_name="signal_peptide_candidate_import_template.csv",
            mime="text/csv",
        )
        upload = st.file_uploader("上传新增候选 CSV 草案", type=["csv"], key="signal_peptide_candidate_upload")
        if upload is not None:
            validation = service.validate_import_csv(upload.getvalue())
            if validation.valid:
                st.success(f"预检通过：{len(validation.rows)} 条新候选可以进入人工审核。")
                preview = pd.DataFrame(validation.rows)
                st.dataframe(
                    preview[
                        [
                            "candidate_id",
                            "category_label",
                            "leader_length",
                            "signal_peptide_length",
                            "processing_route",
                            "source_note",
                        ]
                    ].rename(
                        columns={
                            "candidate_id": "候选 ID",
                            "category_label": "类别",
                            "leader_length": "leader 长度",
                            "signal_peptide_length": "信号肽长度",
                            "processing_route": "加工路线",
                            "source_note": "来源说明",
                        }
                    ),
                    use_container_width=True,
                    hide_index=True,
                )
                st.download_button(
                    "下载合并后的候选库草案 CSV",
                    service.merged_draft_csv(validation.rows),
                    file_name="signal_peptide_candidate_library_draft.csv",
                    mime="text/csv",
                )
            else:
                st.error("预检未通过，请修正后重新上传。")
                st.write(validation.errors)

    if include_external:
        render_external_signal_peptide_screening()


def render_external_signal_peptide_screening() -> None:
    render_uniprot_discovery_panel()
    render_uniprot_method_screening()


def render_uniprot_discovery_panel() -> None:
    with st.expander("从 UniProt API 发现更多候选", expanded=True):
        st.markdown(
            """
            这一步会联网查询 UniProt 中带 `signal peptide` 注释的 Komagataella/Pichia 蛋白，
            自动提取 N 端 signal peptide 作为候选草案，并进行重复检测。它只负责“发现候选”，不会自动加入正式模型。
            """
        )
        api_left, api_mid, api_right = st.columns([1, 1, 1])
        with api_left:
            taxon_id = st.number_input("UniProt organism/taxon ID", min_value=1, value=4922, step=1)
        with api_mid:
            size = st.number_input("最多拉取条目数", min_value=25, max_value=500, value=300, step=25)
        with api_right:
            reviewed_only = st.checkbox("只查 reviewed", value=False)
        if st.button("从 UniProt 获取候选", type="secondary"):
            with st.spinner("正在查询 UniProt 并提取 signal peptide..."):
                discovery = cached_uniprot_signal_peptides(int(taxon_id), int(size), bool(reviewed_only))
            if discovery["errors"] and not discovery["rows"]:
                st.error("没有获取到可用候选。")
                st.write(discovery["errors"])
                st.caption(discovery["source_url"])
            else:
                if discovery["errors"]:
                    st.warning("部分结果被跳过。")
                    st.write(discovery["errors"])
                discovered = pd.DataFrame(discovery["rows"])
                st.success(
                    f"发现 {len(discovered)} 条去重后的外部候选草案；"
                    f"检测到 {int(discovery.get('duplicate_count', 0))} 条重复记录。"
                )
                st.caption(f"已持久保存到：{PATHS.opn_signal_peptides_dir}")
                st.caption(f"来源：{discovery['source_url']}")
                st.dataframe(
                    discovered[
                        [
                            "candidate_id",
                            "leader_length",
                            "signal_peptide_sequence",
                            "source_note",
                            "caution",
                        ]
                    ].rename(
                        columns={
                            "candidate_id": "候选 ID",
                            "leader_length": "信号肽长度",
                            "signal_peptide_sequence": "信号肽序列",
                            "source_note": "UniProt 来源",
                            "caution": "注意事项",
                        }
                    ),
                    use_container_width=True,
                    hide_index=True,
                )
                import_columns = [
                    "candidate_id",
                    "leader_sequence",
                    "signal_peptide_sequence",
                    "category",
                    "processing_route",
                    "source_note",
                    "rationale",
                    "caution",
                ]
                st.download_button(
                    "下载 UniProt 发现候选 CSV",
                    discovered[import_columns].to_csv(index=False).encode("utf-8-sig"),
                    file_name="uniprot_signal_peptide_candidates.csv",
                    mime="text/csv",
                )
                duplicate_rows = pd.DataFrame(discovery.get("duplicate_rows", []))
                if not duplicate_rows.empty:
                    with st.expander("查看重复检测明细", expanded=False):
                        st.dataframe(
                            duplicate_rows[
                                [
                                    "candidate_id",
                                    "accession",
                                    "protein_name",
                                    "signal_peptide_sequence",
                                    "duplicate_reason",
                                    "duplicate_of",
                                ]
                            ].rename(
                                columns={
                                    "candidate_id": "候选 ID",
                                    "accession": "UniProt accession",
                                    "protein_name": "来源蛋白",
                                    "signal_peptide_sequence": "重复信号肽序列",
                                    "duplicate_reason": "重复原因",
                                    "duplicate_of": "重复于",
                                }
                            ),
                            use_container_width=True,
                            hide_index=True,
                        )
                        st.download_button(
                            "下载重复记录 CSV",
                            duplicate_rows.to_csv(index=False).encode("utf-8-sig"),
                            file_name="uniprot_signal_peptide_duplicates.csv",
                            mime="text/csv",
                        )

def render_uniprot_method_screening() -> None:
    with st.expander("从 UniProt 建库并比较筛选方法", expanded=True):
        st.markdown(
            """
            这一步用于扩大候选库并比较证据：先从 UniProt 拉取毕赤酵母相关、带信号肽注释的蛋白，
            再用自研透明规则做预筛；如果本机部署了 USPNet-fast，还会做机器学习复核。
            规则分数只用于解释和排序，不等同于真实表达量；通过候选仍需要进入 pcSecPichia 模型和小试验证。
            """
        )
        service = SignalPeptideScreeningService(PATHS)
        status = service.uspnet_adapter.status()
        if status.available:
            st.success(status.message)
        else:
            st.warning(status.message)
            st.markdown("[打开 USPNet 官方仓库](https://github.com/ml4bio/USPNet)")
        with st.expander("这些方法分别能说明什么", expanded=False):
            st.markdown(
                """
                - **UniProt 注释**：说明该来源蛋白在数据库中带有 signal peptide 注释，是候选发现证据。
                - **自研规则**：只检查长度、N 端电荷、疏水核心和切割位点附近特征，用于可解释预筛，不是训练模型。
                - **USPNet-fast**：MIT License 的外部机器学习模型；部署后可作为商业友好的复核方法。
                - **Razor**：可作为后续第三方参考，但许可不如 MIT 直接，因此不作为默认依赖。
                """
            )

        col_taxon, col_limit, col_reviewed = st.columns([1, 1, 1])
        with col_taxon:
            taxon_id = st.number_input(
                "UniProt taxon ID",
                min_value=1,
                value=4922,
                step=1,
                key="method_screen_taxon_id",
            )
        with col_limit:
            max_records = st.number_input(
                "最多拉取记录数",
                min_value=25,
                max_value=500,
                value=300,
                step=25,
                key="method_screen_max_records",
            )
        with col_reviewed:
            reviewed_only = st.checkbox("只查 reviewed", value=False, key="method_screen_reviewed")

        st.caption(
            "默认范围是 Pichia/毕赤相关 taxon 4922；当前 UniProt 预检约 265 条带 signal peptide 注释的条目。"
            "筛选阶段会优先复用上方 UniProt 查询已保存的候选文件，避免重复联网拉取。"
        )
        if st.button("建立候选库并比较方法", type="primary"):
            with st.spinner("正在读取已保存的 UniProt 候选，并运行规则/USPNet 对比..."):
                st.session_state["method_screening_result"] = service.screen_uniprot_candidates(
                    taxon_id=int(taxon_id),
                    max_records=int(max_records),
                    reviewed_only=bool(reviewed_only),
                )

        result = st.session_state.get("method_screening_result")
        loaded_from_disk = False
        if result is None:
            result = service.load_persisted_screening_result()
            if result is not None:
                st.session_state["method_screening_result"] = result
                loaded_from_disk = True
        if result is None:
            st.info(
                "点击按钮后会生成 UniProt 初始候选、规则评分表和 FASTA；如果 USPNet 已安装，还会给出多方法一致通过结果。"
                "生成后的 CSV/FASTA/JSON 会保存在本地，下次打开页面会自动加载。"
            )
            return

        if loaded_from_disk:
            st.info("已自动加载上次保存的方法比较结果；如果需要更新 UniProt 或 USPNet 结果，再点击按钮重新运行。")
        if result.success:
            st.success(result.message)
        elif result.available:
            st.warning(result.message)
        else:
            st.warning(result.message)

        summary = result.summary
        metric1, metric2, metric3, metric4, metric5, metric6 = st.columns(6)
        metric1.metric("UniProt 初始命中", int(summary.get("uniprot_initial_hits", 0)))
        metric2.metric("去重后候选", int(summary.get("deduplicated_candidates", 0)))
        metric3.metric("重复检测", int(summary.get("uniprot_duplicate_count", 0)))
        metric4.metric("规则高优先", int(summary.get("rules_high_priority", 0)))
        metric5.metric("USPNet 通过", int(summary.get("uspnet_passed", 0)))
        metric6.metric("一致通过", int(summary.get("consensus_passed", 0)))
        st.caption(f"候选来源：{summary.get('uniprot_candidate_source', '本地或本次运行结果')}")

        with st.expander("规则分数和 USPNet 预测怎么读", expanded=True):
            st.markdown(
                """
                - **规则分数**是一个可解释的质控分数，检查长度、N 端电荷、疏水核心、切割位点附近小残基等典型 signal peptide 特征。
                - 很多 UniProt 已注释 signal peptide 会接近满分，这只说明“像标准信号肽”，**不代表分泌产量更高**。
                - **USPNet=SP** 表示机器学习模型也判断它是信号肽；**NO_SP** 表示 USPNet 不支持，需要降级或人工复核。
                - 真正用于 OPN 生产的优先级，还要结合 pcSecPichia 模型、小试表达和加工风险。
                """
            )
            dist1, dist2, dist3, dist4 = st.columns(4)
            dist1.metric("规则分数 ≥95", int(summary.get("rules_score_95_plus", 0)))
            dist2.metric("80-94", int(summary.get("rules_score_80_to_94", 0)))
            dist3.metric("65-79", int(summary.get("rules_score_65_to_79", 0)))
            dist4.metric("<65", int(summary.get("rules_score_below_65", 0)))

        rows = pd.DataFrame(result.rows)
        if not rows.empty:
            passed = rows[rows["recommended_for_draft_library"] == True].copy()
            if passed.empty:
                st.info("当前没有候选通过规则预筛；可以检查 UniProt 来源、长度和疏水核心等规则列。")
            else:
                st.dataframe(
                    passed[
                        [
                            "candidate_id",
                            "accession",
                            "protein_name",
                            "signal_peptide_sequence",
                            "rules_score",
                            "rules_priority",
                            "rules_score_note",
                            "uspnet_prediction_label",
                            "screening_status",
                            "recommended_for_draft_library",
                        ]
                    ].rename(
                        columns={
                            "candidate_id": "候选 ID",
                            "accession": "UniProt accession",
                            "protein_name": "来源蛋白",
                            "signal_peptide_sequence": "UniProt 注释信号肽",
                            "rules_score": "规则分数",
                            "rules_priority": "规则优先级",
                            "rules_score_note": "规则分数说明",
                            "uspnet_prediction_label": "USPNet 预测说明",
                            "screening_status": "综合状态",
                            "recommended_for_draft_library": "建议进入草案库",
                        }
                    ),
                    use_container_width=True,
                    hide_index=True,
                )
                with st.expander("查看规则解释列", expanded=False):
                    st.dataframe(
                        passed[
                            [
                                "candidate_id",
                                "rules_reasons",
                                "rules_risks",
                                "rules_score_note",
                                "uspnet_interpretation",
                                "rules_h_region_max_hydrophobicity",
                                "rules_n_region_positive_count",
                                "rules_c_region_small_neutral",
                            ]
                        ].rename(
                            columns={
                                "candidate_id": "候选 ID",
                                "rules_reasons": "规则支持理由",
                                "rules_risks": "规则风险提示",
                                "rules_score_note": "规则分数说明",
                                "uspnet_interpretation": "USPNet 解释",
                                "rules_h_region_max_hydrophobicity": "最大疏水窗口",
                                "rules_n_region_positive_count": "N 端正电残基数",
                                "rules_c_region_small_neutral": "C 端小残基规则",
                            }
                        ),
                        use_container_width=True,
                        hide_index=True,
                    )
        if result.errors:
            with st.expander("运行提示和错误详情", expanded=False):
                st.write(result.errors)

        st.markdown("**下载候选库和方法对比结果**")
        download_cols = st.columns(4)
        _download_file_button(download_cols[0], result.uniprot_csv, "下载 UniProt 初始候选 CSV", "text/csv")
        _download_file_button(download_cols[1], result.duplicate_csv, "下载重复记录 CSV", "text/csv")
        _download_file_button(download_cols[2], result.comparison_csv, "下载方法对比 CSV", "text/csv")
        _download_file_button(download_cols[3], result.recommended_fasta, "下载推荐候选 FASTA", "text/plain")
        st.caption(f"输出目录：{result.output_dir}")


