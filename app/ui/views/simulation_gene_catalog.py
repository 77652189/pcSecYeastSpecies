from __future__ import annotations

import hashlib
import math

import pandas as pd
import streamlit as st

from app.services.pichia_gene_catalog_service import (
    build_pichia_gene_evidence_cache,
    get_pichia_ko_genes_for_selection,
    get_pichia_ko_reactions_for_selection,
    get_pichia_oe_reactions_for_selection,
    list_pichia_gene_rule_evidence,
    list_pichia_secretion_gene_evidence,
    list_verified_secretion_gene_library,
    load_pichia_full_model_gene_catalog,
    pichia_full_model_gene_catalog_cache_path,
)
from app.ui.views.simulation_gene_text import merge_candidate_text

KO_RUNNABLE_STATUS = "ko_runnable_gpr_gene_deletion"
OE_RUNNABLE_STATUS = "oe_runnable_reaction_proxy"


def render_gene_lookup_panel() -> None:
    st.markdown("**已验证分泌工程候选库**")
    st.caption(
        "默认展示 37 条分泌工程候选，并合并本地数据库/模型证据。"
        "这里的“已验证”指数据库证据支持与模型可解释，不等于湿实验已完成验证。"
    )
    query = st.text_input(
        "搜索候选基因 / locus / 功能 / 反应",
        value=st.session_state.get("pichia_gene_lookup_query", ""),
        placeholder="例如：PDI、ERO1、Kar2、ERAD、PAS_chr",
        key="pichia_gene_lookup_query",
    )
    _render_verified_secretion_gene_library(query)

    with st.expander("高级：全模型 GPR 基因库（1025 个模型基因）", expanded=False):
        st.caption("这里是原始模型 GPR 视图，用于严格 gene-level KO 和 GPR-aware OE proxy 检索；默认不加载，避免首屏卡顿。")
        show_full = st.checkbox("加载全模型 GPR 基因库", value=False, key="pichia_gene_show_full")
        if show_full:
            _render_full_model_gene_lookup(query)
        else:
            st.info("打开上方开关后加载全模型基因库。")

    with st.expander("高级：反应级代理（sec_* 复合体 / 路径反应）", expanded=False):
        st.caption("反应级代理可用于模型解释或 OE proxy，但不是湿实验基因名，也不是 gene-level GPR。")
        show_reaction_proxies = st.checkbox(
            "加载反应级代理",
            value=False,
            key="pichia_gene_show_reaction_proxies",
        )
        if show_reaction_proxies:
            _render_reaction_proxy_lookup(query)
        else:
            st.info("打开上方开关后查看 sec_* 反应代理。")

    with st.expander("高级：外部证据 GPR overlay / 证据维护", expanded=False):
        st.caption("外部证据 overlay 是实验性证据层，默认不进入仿真；当前无可执行补充规则时只作人工复核参考。")
        maintenance_col, overlay_col = st.columns([1.0, 1.0])
        with maintenance_col:
            refresh_lightweight_cache = st.button(
                "刷新常用基因证据缓存",
                key="pichia_gene_refresh_lightweight_cache",
                help="刷新 37 条分泌工程候选及反应代理缓存；不会重建全模型湿实验注释。",
            )
            if refresh_lightweight_cache:
                refreshed = list_pichia_secretion_gene_evidence("", force_refresh=True)
                st.success(f"常用基因证据缓存已刷新：{len(refreshed)} 条。")
        with overlay_col:
            show_overlay = st.checkbox("显示外部证据 GPR overlay", value=False, key="pichia_gene_show_overlay")
        if show_overlay:
            _render_gene_rule_overlay_lookup(query)
        else:
            st.info("打开上方开关后查看 PDI1/ERO1/KAR2/OCH1/PEP4/PRB1 等 overlay 证据。")


def _render_verified_secretion_gene_library(query: str) -> None:
    rows = list_verified_secretion_gene_library(query)
    if not rows:
        st.info("已验证分泌工程候选库中未找到匹配。可在高级区加载全模型 GPR 基因库继续搜索。")
        return
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "标准名": row["display_name"],
                    "locus tag / 模型基因 ID": row["locus_tag"] or row["model_gene_id"] or "需人工确认",
                    "功能": row["function_annotation"],
                    "可执行操作": row["operation_status"],
                    "证据等级": row["evidence_tier"],
                    "推荐用途": row["recommended_use"],
                    "分类": row["category"],
                }
                for row in rows
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )
    selected = st.multiselect(
        "选择候选基因",
        [str(row["display_name"]) for row in rows if row.get("display_name")],
        key="pichia_verified_gene_sel",
    )
    action_col_ko, action_col_oe, detail_col = st.columns([1.0, 1.0, 1.0])
    with action_col_ko:
        if st.button("添加到敲除输入", key="pichia_verified_gene_add_ko") and selected:
            _add_curated_knockout_selection(selected)
    with action_col_oe:
        if st.button("添加到过表达反应代理", key="pichia_verified_gene_add_oe") and selected:
            _add_curated_oe_reaction_selection(selected)
    with detail_col:
        st.caption("PDI1/ERO1 等无可执行 GPR overlay 时不会被写入 gene-level KO。")
    message = st.session_state.pop("pichia_gene_catalog_message", "")
    if message:
        st.info(message)
    with st.expander("查看候选证据详情", expanded=False):
        st.dataframe(
            pd.DataFrame(
                [
                    {
                        "标准名": row["display_name"],
                        "模型 gene": row["model_gene_id"] or "",
                        "候选 locus": row["locus_tag"] or "",
                        "映射状态": _curated_mapping_status_label(row.get("mapping_status")),
                        "GPR overlay": _rule_overlay_status_label(row.get("rule_status")),
                        "外部置信度": row.get("rule_confidence") or "",
                        "KO 代理反应": row.get("ko_reaction_id") or "",
                        "OE 代理反应": row.get("oe_reaction_id") or "",
                        "证据摘要": row["source_summary"],
                    }
                    for row in rows
                ]
            ),
            use_container_width=True,
            hide_index=True,
        )


def _render_full_model_gene_lookup(query: str) -> None:
    filter_col, oe_col, wet_lab_col, refresh_col, page_size_col = st.columns([1.1, 1.2, 1.4, 1.1, 1.0])
    with filter_col:
        only_ko = st.checkbox("只显示可敲除基因", value=False, key="pichia_gene_only_ko")
    with oe_col:
        only_oe = st.checkbox("只显示可过表达代理", value=False, key="pichia_gene_only_oe")
    with wet_lab_col:
        wet_lab_filter = st.selectbox(
            "湿实验注释",
            ["全部", "可直接推进湿实验", "需人工确认", "仅模型级候选"],
            index=0,
            key="pichia_gene_wet_lab_filter",
        )
    with refresh_col:
        force_refresh = st.button("刷新基因目录缓存", key="pichia_gene_refresh_cache")
    with page_size_col:
        page_size = int(
            st.number_input(
                "每页最大行数",
                min_value=25,
                max_value=500,
                value=100,
                step=25,
                key="pichia_gene_page_size",
            )
        )

    if st.button("在线刷新湿实验注释缓存", key="pichia_gene_refresh_evidence"):
        with st.spinner("正在从 UniProt/KEGG 构建湿实验注释缓存..."):
            summary = build_pichia_gene_evidence_cache()
        st.success(
            "湿实验注释缓存已更新："
            f"{summary.get('database_supported_count', 0)} / {summary.get('total_genes', 0)} 个基因有数据库支持。"
        )
        force_refresh = False

    full_rows = load_pichia_full_model_gene_catalog(force_refresh=force_refresh)
    filtered_rows = _filter_full_model_gene_rows(
        full_rows,
        query=query,
        only_ko=only_ko,
        only_oe=only_oe,
        wet_lab_filter=wet_lab_filter,
    )
    if not filtered_rows:
        st.info("模型 GPR 基因目录中未找到匹配。可以打开分泌工程基因名或反应级代理继续检索。")
        return
    page_signature = (query.strip(), only_ko, only_oe, wet_lab_filter, page_size, len(filtered_rows))
    if st.session_state.get("pichia_gene_page_signature") != page_signature:
        st.session_state["pichia_gene_page"] = 1
        st.session_state["pichia_gene_page_signature"] = page_signature
    page_number = int(st.session_state.get("pichia_gene_page", 1))
    display_rows, page_number, total_pages = _paginate_full_model_gene_rows(
        filtered_rows,
        page_number=page_number,
        page_size=page_size,
    )
    st.session_state["pichia_gene_page"] = page_number
    cache_path = pichia_full_model_gene_catalog_cache_path()
    st.caption(
        f"模型 GPR 基因缓存：`{cache_path.parent}`。"
        f"当前第 {page_number} / {total_pages} 页，显示 {len(display_rows)} 个；"
        f"共 {len(filtered_rows)} 个匹配基因。"
    )
    previous_col, page_col, next_col = st.columns([1.0, 1.2, 1.0])
    with previous_col:
        if st.button("上一页", disabled=page_number <= 1, key="pichia_gene_previous_page"):
            st.session_state["pichia_gene_page"] = max(1, page_number - 1)
            st.rerun()
    with page_col:
        requested_page = int(
            st.number_input(
                "页码",
                min_value=1,
                max_value=total_pages,
                value=page_number,
                step=1,
                key=_page_input_widget_key(page_signature, page_number, total_pages),
            )
        )
        if requested_page != page_number:
            st.session_state["pichia_gene_page"] = requested_page
            st.rerun()
    with next_col:
        if st.button("下一页", disabled=page_number >= total_pages, key="pichia_gene_next_page"):
            st.session_state["pichia_gene_page"] = min(total_pages, page_number + 1)
            st.rerun()
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "显示名称": _full_model_gene_display_name(gene),
                    "模型基因 ID": gene["gene_id"],
                    "可用操作": _gene_action_label(gene),
                    "湿实验状态": _wet_lab_readiness_label(gene.get("wet_lab_readiness")),
                    "分类": gene["primary_category"],
                    "功能 / 依据": _full_model_gene_function_summary(gene),
                    "数据库 ID": _external_id_summary(gene),
                    "通路": _process_label(gene.get("processes")),
                    "反应数": gene["n_reactions"],
                    "KO 状态": _ko_status_label(gene.get("ko_support_status")),
                    "OE 状态": _oe_status_label(gene.get("oe_support_status")),
                    "GPR 角色": _gpr_role_label(gene.get("gpr_role")),
                    "置信度": gene.get("confidence", ""),
                    "证据等级": gene.get("evidence_confidence", ""),
                    "证据来源": ", ".join(str(item) for item in gene.get("evidence_sources") or []),
                    "别名": ", ".join(str(item) for item in gene.get("aliases") or []),
                    "缺失信息": ", ".join(str(item) for item in gene.get("missing_information") or []),
                }
                for gene in display_rows
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )
    gene_name_by_id = {str(gene["gene_id"]): _full_model_gene_display_name(gene) for gene in display_rows}
    selected = st.multiselect(
        "选择模型基因",
        [str(gene["gene_id"]) for gene in display_rows],
        key="pichia_full_sel",
        format_func=lambda gene_id: f"{gene_name_by_id.get(str(gene_id), str(gene_id))}（{gene_id}）",
    )
    add_ko_col, add_oe_col = st.columns(2)
    with add_ko_col:
        if st.button("添加到敲除基因") and selected:
            current = str(st.session_state.get("pichia_draft_ko_genes", ""))
            st.session_state["pichia_draft_ko_genes"] = merge_candidate_text(current, selected)
            st.rerun()
    with add_oe_col:
        if st.button("添加到过表达基因代理") and selected:
            current = str(st.session_state.get("pichia_draft_oe_genes", ""))
            st.session_state["pichia_draft_oe_genes"] = merge_candidate_text(current, selected)
            st.rerun()


def _paginate_full_model_gene_rows(
    rows: list[dict[str, object]],
    *,
    page_number: int,
    page_size: int,
) -> tuple[list[dict[str, object]], int, int]:
    safe_page_size = max(1, int(page_size))
    total_pages = max(1, math.ceil(len(rows) / safe_page_size))
    current_page = min(max(1, int(page_number)), total_pages)
    start = (current_page - 1) * safe_page_size
    end = start + safe_page_size
    return rows[start:end], current_page, total_pages


def _page_input_widget_key(page_signature: tuple[object, ...], page_number: int, total_pages: int) -> str:
    digest = hashlib.sha1(repr((page_signature, page_number, total_pages)).encode("utf-8")).hexdigest()[:12]
    return f"pichia_gene_page_input_{digest}"


def _filter_full_model_gene_rows(
    rows: list[dict[str, object]],
    *,
    query: str = "",
    only_ko: bool = False,
    only_oe: bool = False,
    wet_lab_filter: str = "全部",
) -> list[dict[str, object]]:
    query_text = query.strip().lower()
    filtered_rows = [
        gene
        for gene in rows
        if (not only_ko or gene.get("ko_support_status") == KO_RUNNABLE_STATUS)
        and (not only_oe or gene.get("oe_support_status") == OE_RUNNABLE_STATUS)
        and _matches_wet_lab_filter(gene, wet_lab_filter)
    ]
    if not query_text:
        return filtered_rows
    return [gene for gene in filtered_rows if _full_model_gene_matches_query(gene, query_text)]


def _full_model_gene_matches_query(gene: dict[str, object], query_text: str) -> bool:
    return (
        query_text in str(gene.get("gene_id") or "").lower()
        or query_text in _full_model_gene_display_name(gene).lower()
        or query_text in _full_model_gene_function_summary(gene).lower()
        or query_text in str(gene.get("primary_category") or "").lower()
        or query_text in str(gene.get("processes") or "").lower()
        or query_text in " ".join(str(item) for item in gene.get("aliases") or []).lower()
        or query_text in str(gene.get("protein_name") or "").lower()
        or query_text in str(gene.get("function_annotation") or "").lower()
        or query_text in str(gene.get("ko_support_status") or "").lower()
        or query_text in str(gene.get("oe_support_status") or "").lower()
    )


def _full_model_gene_display_name(gene: dict[str, object]) -> str:
    display_name = str(gene.get("display_name") or "").strip()
    if display_name and display_name != str(gene.get("gene_id") or "").strip():
        return display_name
    standard_symbol = str(gene.get("standard_gene_symbol") or "").strip()
    if standard_symbol:
        return standard_symbol
    protein_name = str(gene.get("protein_name") or "").strip()
    if protein_name:
        return protein_name
    aliases = [str(item).strip() for item in gene.get("aliases") or [] if str(item).strip()]
    if aliases:
        return aliases[0]
    reactions = _reaction_tokens(gene)
    if reactions:
        return f"{'/'.join(reactions[:3])} 相关酶（未注释）"
    return "未注释模型基因"


def _full_model_gene_function_summary(gene: dict[str, object]) -> str:
    annotation = str(gene.get("function_annotation") or "").strip()
    if annotation:
        return annotation
    reactions = _reaction_tokens(gene)
    if reactions:
        joined = ", ".join(reactions[:5])
        suffix = " 等" if len(reactions) > 5 else ""
        return f"按模型 GPR 关联到反应：{joined}{suffix}。尚无外部基因名/功能注释。"
    return "尚无外部基因名/功能注释；仅可按模型 locus ID 和 GPR 关系使用。"


def _matches_wet_lab_filter(gene: dict[str, object], wet_lab_filter: str) -> bool:
    readiness = str(gene.get("wet_lab_readiness") or "model_only_not_experiment_ready")
    if wet_lab_filter == "可直接推进湿实验":
        return readiness == "database_supported_experiment_candidate"
    if wet_lab_filter == "需人工确认":
        return readiness == "manual_review_required"
    if wet_lab_filter == "仅模型级候选":
        return readiness == "model_only_not_experiment_ready"
    return True


def _wet_lab_readiness_label(readiness: object) -> str:
    labels = {
        "database_supported_experiment_candidate": "可直接推进：数据库精确支持",
        "manual_review_required": "需人工确认：有部分数据库证据",
        "model_only_not_experiment_ready": "仅模型级候选：不建议直接实验",
    }
    return labels.get(str(readiness or ""), str(readiness or ""))


def _external_id_summary(gene: dict[str, object]) -> str:
    external_ids = gene.get("external_ids") if isinstance(gene.get("external_ids"), dict) else {}
    if not external_ids:
        return ""
    ordered_keys = ("uniprot", "ncbi_gene", "kegg", "refseq")
    return "; ".join(f"{key}: {external_ids[key]}" for key in ordered_keys if external_ids.get(key))


def _reaction_tokens(gene: dict[str, object]) -> list[str]:
    raw_reactions = gene.get("sample_reactions") or gene.get("affected_reactions") or []
    tokens: list[str] = []
    seen: set[str] = set()
    for item in raw_reactions:
        token = str(item).strip()
        if not token:
            continue
        token = token.split("_no_", 1)[0]
        token = token.removesuffix("_fwd").removesuffix("_rvs")
        if token and token not in seen:
            seen.add(token)
            tokens.append(token)
    return tokens


def _gene_action_label(gene: dict[str, object]) -> str:
    actions: list[str] = []
    if gene.get("ko_support_status") == KO_RUNNABLE_STATUS:
        actions.append("可敲除")
    if gene.get("oe_support_status") == OE_RUNNABLE_STATUS:
        actions.append("可过表达代理")
    if not actions:
        return "仅解释 / 暂不可运行"
    return " / ".join(actions)


def _ko_status_label(status: object) -> str:
    labels = {
        "ko_runnable_gpr_gene_deletion": "可运行：基因级 KO",
        "ko_no_gpr_effect": "不可运行：无 GPR 影响",
        "ko_explain_only_complex_subunit": "仅解释：复合体亚基",
    }
    return labels.get(str(status or ""), str(status or ""))


def _oe_status_label(status: object) -> str:
    labels = {
        "oe_runnable_reaction_proxy": "可运行：反应级 OE 代理",
        "oe_no_gpr_effect": "不可运行：无 GPR 影响",
        "oe_explain_only_complex_subunit": "仅解释：复合体亚基",
    }
    return labels.get(str(status or ""), str(status or ""))


def _gpr_role_label(role: object) -> str:
    labels = {
        "single_gene": "单基因",
        "complex_subunit": "复合体亚基",
        "isozyme": "同工酶",
        "no_gpr_effect": "无 GPR 影响",
    }
    return labels.get(str(role or ""), str(role or ""))


def _process_label(process: object) -> str:
    labels = {
        "metabolic_or_other": "代谢 / 其他",
        "translation": "翻译",
        "er_translocation": "ER 转运",
        "folding_dsb": "折叠 / DSB",
        "glycosylation": "糖基化",
        "misfolding_erad": "错误折叠 / ERAD",
        "transport_secretion": "运输 / 分泌",
    }
    text = str(process or "")
    return labels.get(text, text)


def _render_matlab_gene_target_lookup(query: str, *, force_refresh: bool = False) -> None:
    evidence_rows = list_pichia_secretion_gene_evidence(query, force_refresh=force_refresh)
    if not evidence_rows:
        st.info("分泌工程基因名中未找到匹配。")
        return
    st.markdown("**分泌工程基因名（带证据映射）**")
    st.caption(
        "这里显示分泌工程常用名、MATLAB 路径名称与模型证据的对应关系。"
        "若没有模型 GPR gene ID，则不能作为严格 gene-level 输入。"
    )
    rule_evidence_by_name = {
        str(row.get("common_name") or ""): row
        for row in list_pichia_gene_rule_evidence(query)
    }
    rows, current_category = [], None
    for row in evidence_rows:
        if row["category"] != current_category:
            current_category = str(row["category"])
            rows.append(
                {
                    "常用基因名": f"▸ {current_category}",
                    "说明": "",
                    "模型 GPR gene ID": "",
                    "证据状态": "",
                    "推荐用途": "",
                    "候选 locus tag": "",
                    "外部证据": "",
                    "GPR 补充状态": "",
                    "补充建议": "",
                    "可用代理反应": "",
                    "代理反应证据": "",
                }
            )
        rule_evidence = rule_evidence_by_name.get(str(row["common_name"]), {})
        rows.append(
            {
                "常用基因名": row["common_name"],
                "说明": row["description"],
                "模型 GPR gene ID": row["mapped_model_gene_id"] or row["declared_model_gene_id"] or "无模型 GPR gene ID",
                "证据状态": _curated_mapping_status_label(row.get("mapping_status")),
                "推荐用途": _curated_recommended_use_label(row.get("recommended_use")),
                "候选 locus tag": rule_evidence.get("candidate_locus_tag") or "",
                "外部证据": _rule_evidence_source_label(rule_evidence),
                "GPR 补充状态": _rule_overlay_status_label(rule_evidence.get("rule_status")),
                "补充建议": _rule_overlay_action_label(rule_evidence.get("recommended_action")),
                "可用代理反应": row.get("oe_reaction_id") or row.get("ko_reaction_id") or "",
                "代理反应证据": _reaction_proxy_evidence_label(row),
            }
        )
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    selected = st.multiselect(
        "选择分泌工程基因名",
        [str(row["common_name"]) for row in evidence_rows if row.get("common_name")],
        key="pichia_matlab_gene_sel",
    )
    action_col_ko, action_col_oe = st.columns(2)
    with action_col_ko:
        if st.button("添加到敲除输入", key="pichia_matlab_gene_add_ko") and selected:
            _add_curated_knockout_selection(selected)
    with action_col_oe:
        if st.button("添加到过表达反应代理", key="pichia_matlab_gene_add_oe") and selected:
            _add_curated_oe_reaction_selection(selected)
    message = st.session_state.pop("pichia_gene_catalog_message", "")
    if message:
        st.info(message)
    st.caption(
        "如果模型 GPR gene ID 为空，这个名称不能直接作为 gene-level KO/OE 输入；"
        "请使用其对应反应级代理，或先人工确认 K. phaffii locus ID。"
    )


def _render_reaction_proxy_lookup(query: str, *, force_refresh: bool = False) -> None:
    proxies = [
        row
        for row in list_pichia_secretion_gene_evidence(query, force_refresh=force_refresh)
        if row.get("oe_reaction_id") or row.get("ko_reaction_id")
    ]
    if not proxies:
        st.info("反应级代理中未找到匹配。")
        return
    st.markdown("**反应级代理**")
    st.caption(
        "这些条目直接写入 KO/OE 反应代理输入，代表模型中的复合体形成或分泌路径反应；"
        "它们不是 gene-level 扰动，也不能直接等同于湿实验基因名。"
    )
    rows, current_category = [], None
    for row in proxies:
        if row["category"] != current_category:
            current_category = str(row["category"])
            rows.append({"代理名称": f"▸ {current_category}", "代理反应 ID": "", "来源基因名": "", "用途": "", "说明": ""})
        rows.append(
            {
                "代理名称": row["common_name"],
                "代理反应 ID": row.get("oe_reaction_id") or row.get("ko_reaction_id") or "",
                "来源基因名": row["common_name"],
                "用途": _curated_recommended_use_label(row.get("recommended_use")),
                "说明": f"{row['description']}；{_reaction_proxy_evidence_label(row)}",
            }
        )
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    selected = st.multiselect(
        "选择反应级代理",
        [str(row["common_name"]) for row in proxies if row.get("common_name")],
        key="pichia_reaction_proxy_sel",
    )
    action_col_ko, action_col_oe = st.columns(2)
    with action_col_ko:
        if st.button("添加到敲除反应代理", key="pichia_reaction_proxy_add_ko") and selected:
            _add_curated_knockout_selection(selected)
    with action_col_oe:
        if st.button("添加到过表达反应代理", key="pichia_reaction_proxy_add_oe") and selected:
            _add_curated_oe_reaction_selection(selected)
    message = st.session_state.pop("pichia_gene_catalog_message", "")
    if message:
        st.info(message)


def _render_gene_rule_overlay_lookup(query: str) -> None:
    rows = list_pichia_gene_rule_evidence(query)
    if not rows:
        st.info("外部证据 GPR overlay 中未找到匹配。")
        return
    executable_rows = [
        row
        for row in rows
        if str(row.get("rule_status") or "").startswith("overlay_executable")
    ]
    if not executable_rows:
        st.warning("暂无可执行补充规则；这些证据只能用于人工复核或反应级代理解释。")
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "常用名": row.get("common_name") or "",
                    "候选 locus tag": row.get("candidate_locus_tag") or "未解析",
                    "蛋白 / 功能": row.get("protein_name") or "",
                    "证据置信度": row.get("confidence") or "",
                    "GPR 补充状态": _rule_overlay_status_label(row.get("rule_status")),
                    "推荐动作": _rule_overlay_action_label(row.get("recommended_action")),
                    "目标反应": ", ".join(str(item) for item in row.get("target_reaction_ids") or []),
                    "证据来源": _rule_evidence_source_label(row),
                }
                for row in rows
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )


def _curated_mapping_status_label(status: object) -> str:
    labels = {
        "model_gpr_gene_available": "已映射到模型 GPR gene",
        "reaction_proxy_only": "无 GPR gene；仅反应级代理",
        "declared_proxy_missing_in_model": "声明了代理但模型未找到反应",
        "literature_name_only": "仅文献/路径名称；需人工确认",
    }
    return labels.get(str(status or ""), str(status or ""))


def _curated_recommended_use_label(use: object) -> str:
    labels = {
        "gene_level_gpr_perturbation": "可用于 gene-level GPR 扰动",
        "reaction_level_proxy_requires_locus_review": "可用于反应级代理；湿实验需确认 locus ID",
        "manual_review_required": "需人工确认后使用",
    }
    return labels.get(str(use or ""), str(use or ""))


def _rule_evidence_source_label(row: dict[str, object]) -> str:
    if not row:
        return ""
    sources = ", ".join(str(item) for item in row.get("evidence_sources") or [])
    confidence = str(row.get("confidence") or "")
    external_ids = row.get("external_ids") if isinstance(row.get("external_ids"), dict) else {}
    ids = "; ".join(f"{key}: {value}" for key, value in external_ids.items() if value)
    parts = [part for part in (confidence, sources, ids) if part]
    return " | ".join(parts)


def _rule_overlay_status_label(status: object) -> str:
    labels = {
        "overlay_executable": "可执行 overlay（实验性）",
        "overlay_executable_complex_rule": "可执行复合体 overlay（实验性）",
        "display_only_requires_multi_source_confirmation": "仅展示：需要多源确认",
        "display_only_name_context_not_exact_kar2_locus": "仅展示：名称上下文不是精确 locus",
        "display_only_conflicts_with_existing_model_gene_annotation": "仅展示：与当前模型注释冲突",
        "display_only_multiple_or_indirect_candidates": "仅展示：多候选或间接证据",
        "not_executable": "不可执行：证据不足",
    }
    return labels.get(str(status or ""), str(status or ""))


def _rule_overlay_action_label(action: object) -> str:
    labels = {
        "enable_only_for_explicit_analysis": "仅在显式实验模式中使用",
        "keep_reaction_level_proxy_until_locus_is_confirmed": "确认 locus 前保留反应级代理",
        "manual_locus_review_required": "需要人工确认 K. phaffii locus",
        "do_not_replace_existing_model_gene_without_review": "不要在未复核前替换模型基因",
        "manual_review_required": "需要人工复核",
    }
    return labels.get(str(action or ""), str(action or ""))


def _reaction_proxy_evidence_label(row: dict[str, object]) -> str:
    evidence_rows = [item for item in row.get("reaction_evidence") or [] if isinstance(item, dict)]
    if not evidence_rows:
        return "无代理反应"
    existing = [str(item.get("reaction_id")) for item in evidence_rows if item.get("exists_in_model")]
    missing = [str(item.get("reaction_id")) for item in evidence_rows if not item.get("exists_in_model")]
    has_gpr = any(bool(item.get("has_gpr_rule")) for item in evidence_rows)
    parts: list[str] = []
    if existing:
        parts.append(f"模型中存在：{', '.join(existing)}")
    if missing:
        parts.append(f"模型中未找到：{', '.join(missing)}")
    parts.append("代理反应含 GPR 规则" if has_gpr else "代理反应无 GPR 规则")
    return "；".join(parts)


def _add_curated_knockout_selection(selected: list[str]) -> None:
    genes = get_pichia_ko_genes_for_selection(selected)
    if genes:
        current = str(st.session_state.get("pichia_draft_ko_genes", ""))
        st.session_state["pichia_draft_ko_genes"] = merge_candidate_text(current, genes)
        st.session_state["pichia_gene_catalog_message"] = f"已加入敲除基因：{', '.join(genes)}"
        st.rerun()
    reactions = get_pichia_ko_reactions_for_selection(selected)
    if reactions:
        current = str(st.session_state.get("pichia_draft_ko_reactions", ""))
        st.session_state["pichia_draft_ko_reactions"] = merge_candidate_text(current, reactions)
        st.session_state["pichia_gene_catalog_message"] = f"已加入复合体级 KO 反应：{', '.join(reactions)}"
        st.rerun()
    st.session_state["pichia_gene_catalog_message"] = "所选策略热点没有可靠的敲除模型基因或 KO 反应 ID；可尝试添加到过表达反应代理。"


def _add_curated_oe_reaction_selection(selected: list[str]) -> None:
    reactions = get_pichia_oe_reactions_for_selection(selected)
    if reactions:
        current = str(st.session_state.get("pichia_draft_oe_reactions", ""))
        st.session_state["pichia_draft_oe_reactions"] = merge_candidate_text(current, reactions)
        st.session_state["pichia_gene_catalog_message"] = f"已加入过表达反应代理：{', '.join(reactions)}"
        st.rerun()
    st.session_state["pichia_gene_catalog_message"] = "所选条目没有可用的过表达反应代理 ID；请勾选显示全部模型基因后选择模型基因 ID。"


__all__ = ["render_gene_lookup_panel"]
