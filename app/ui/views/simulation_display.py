from __future__ import annotations

import pandas as pd


CANDIDATE_DISPLAY_COLUMNS = {
    "input_gene_id": "输入基因",
    "canonical_gene_id": "模型基因",
    "standard_gene_symbol": "标准基因名",
    "display_name": "显示名称",
    "protein_name": "蛋白名称",
    "function_annotation": "功能注释",
    "external_ids": "数据库 ID",
    "evidence_sources": "证据来源",
    "evidence_confidence": "证据等级",
    "wet_lab_readiness": "湿实验状态",
    "resolved_reaction_id": "解析反应",
    "reaction_id": "反应",
    "intervention_type": "扰动类型",
    "effect_label": "效果",
    "objective_value": "目标值",
    "baseline_objective_value": "基线",
    "candidate_id": "候选ID",
    "delta_objective": "Δ目标值",
    "secretory_process": "分泌环节",
    "solver_status_label": "求解状态",
    "failure_reason": "失败原因",
    "mapping_level": "映射层级",
    "mapping_confidence": "解释置信度",
    "mapping_interpretation": "解释",
    "complex_id": "复合体",
    "gpr_role": "GPR 角色",
    "capacity_effect": "容量影响",
    "simulation_basis": "模拟依据",
    "ko_support_status": "KO 支持状态",
    "oe_support_status": "OE 支持状态",
    "support_reason": "支持依据",
    "missing_information": "缺失信息",
    "database_annotation_sources": "数据库注释来源",
    "database_annotation_confidence": "数据库注释置信度",
    "model_gpr_executable": "模型 GPR 可执行",
    "oe_reaction_proxy": "OE 反应代理",
    "phenotype_evidence": "表型证据",
    "recommendation_tier": "证据分级",
    "recommendation_tier_reason": "分级理由",
    "inactive_reaction_count": "失活反应数",
    "inactive_reactions_preview": "失活反应预览",
    "affected_reactions": "关联反应",
    "warnings": "提示",
}


TARGET_SEMANTICS_LABELS = {
    "project_defined_hLF": "项目定义 hLF（用户提供序列）",
    "product_target": "产品草案目标",
    "opn_leader_candidate": "OPN leader 候选",
    "custom": "自定义目标",
    "native_signal_plus_mature_hLF": "人源天然信号肽 + mature hLF",
    "mature_secreted_with_leader_candidate": "成熟蛋白 + leader 候选",
    "mature_secreted": "成熟分泌蛋白",
    "custom_user_sequence": "用户自定义序列",
    "user_provided_as_provided": "用户提供序列，按原样使用",
    "as_provided": "按提供序列记录",
}


MAPPING_LEVEL_LABELS = {
    "direct_gpr": "GPR 直接关联",
    "complex_subunit": "复合体亚基",
    "reaction_proxy": "反应代理",
    "metabolic_or_other": "代谢/其它",
    "unresolved": "未解析",
}


MAPPING_CONFIDENCE_LABELS = {
    "high": "高",
    "medium": "中",
    "low": "低",
    "unresolved": "未解析",
}


GPR_ROLE_LABELS = {
    "single_gene": "单基因酶",
    "isoenzyme": "同工酶",
    "complex_subunit": "复合体亚基",
    "mixed": "混合 GPR",
    "reaction_level": "反应级",
    "unresolved": "未解析",
}


CAPACITY_EFFECT_LABELS = {
    "disables_reactions": "失活反应",
    "no_reaction_disabled": "未失活反应",
    "reaction_disabled": "反应敲除",
    "reaction_capacity_proxy": "反应容量代理",
    "partial_reaction_capacity_proxy": "部分反应容量代理",
    "complex_subunit_limited": "复合体亚基受限",
    "manual_review_required": "需要人工检查",
    "no_gpr_effect": "模型无 GPR 影响",
    "unresolved": "未解析",
}


SIMULATION_BASIS_LABELS = {
    "gpr_gene_deletion": "GPR 基因敲除",
    "reaction_deletion": "反应敲除",
    "reaction_level_capacity_proxy": "反应级容量代理",
    "explain_only": "仅解释",
    "unresolved": "未解析",
}


KO_SUPPORT_STATUS_LABELS = {
    "ko_runnable_gpr_gene_deletion": "可运行：GPR 基因敲除",
    "ko_no_reaction_disabled": "无明显 GPR 失活",
    "ko_no_gpr_effect": "模型无 GPR 影响",
    "unresolved_gene": "基因未解析",
    "reaction_level_diagnostic": "反应级诊断入口",
}


OE_SUPPORT_STATUS_LABELS = {
    "oe_runnable_reaction_proxy": "可运行：反应容量代理",
    "oe_explain_only_complex_subunit": "仅解释：复合体亚基受限",
    "oe_explain_only_no_capacity_model": "仅解释：缺少容量模型",
    "oe_no_gpr_effect": "模型无 GPR 影响",
    "unresolved_gene": "基因未解析",
    "reaction_level_diagnostic": "反应级诊断入口",
}


WET_LAB_READINESS_LABELS = {
    "database_supported_experiment_candidate": "可直接推进：数据库精确支持",
    "manual_review_required": "需人工确认：有部分数据库证据",
    "model_only_not_experiment_ready": "仅模型级候选：不建议直接实验",
}


def display_value(value: object, fallback: str = "") -> str:
    if value is None:
        return fallback
    try:
        if pd.isna(value):
            return fallback
    except TypeError:
        pass
    text = str(value)
    if text.lower() in {"nan", "none", "nat"}:
        return fallback
    return text


def target_semantics_label(value: object) -> str:
    text = display_value(value)
    return TARGET_SEMANTICS_LABELS.get(text, text or "未声明")


def screen_status_label(raw_status: object, success: object) -> str:
    status = display_value(raw_status)
    if status == "no_reaction_disabled":
        return "未运行：GPR 未失活任何反应"
    if status == "not_run_no_gpr_effect":
        return "仅解释，模型无 GPR 影响"
    success_text = str(success).strip().lower()
    if success is True or success_text == "true":
        return "求解成功"
    return {
        "2": "约束不可行",
        "3": "目标无界",
        "4": "求解器数值错误",
        "missing_reaction": "反应未找到",
        "unresolved_gene": "基因未解析",
        "unresolved_reaction": "反应未解析",
        "not_run_complex_subunit_limited": "仅解释，未求解",
        "not_run_gene_oe_proxy": "仅解释，未求解",
        "not_run_no_gpr_effect": "仅解释，模型无 GPR 影响",
        "no_reaction_disabled": "未运行：GPR 未失活任何反应",
        "missing_objective": "目标反应未找到",
    }.get(status, "求解失败")


def normalise_candidate_frame_for_display(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    if "solver_status_label" not in frame.columns:
        frame["solver_status_label"] = ""
    if "failure_reason" not in frame.columns:
        frame["failure_reason"] = ""
    if "mapping_level" in frame.columns:
        frame["mapping_level"] = frame["mapping_level"].map(
            lambda value: MAPPING_LEVEL_LABELS.get(display_value(value), display_value(value))
        )
    if "mapping_confidence" in frame.columns:
        frame["mapping_confidence"] = frame["mapping_confidence"].map(
            lambda value: MAPPING_CONFIDENCE_LABELS.get(display_value(value), display_value(value))
        )
    if "gpr_role" in frame.columns:
        frame["gpr_role"] = frame["gpr_role"].map(
            lambda value: GPR_ROLE_LABELS.get(display_value(value), display_value(value))
        )
    if "capacity_effect" in frame.columns:
        frame["capacity_effect"] = frame["capacity_effect"].map(
            lambda value: CAPACITY_EFFECT_LABELS.get(display_value(value), display_value(value))
        )
    if "simulation_basis" in frame.columns:
        frame["simulation_basis"] = frame["simulation_basis"].map(
            lambda value: SIMULATION_BASIS_LABELS.get(display_value(value), display_value(value))
        )
    if "ko_support_status" in frame.columns:
        frame["ko_support_status"] = frame["ko_support_status"].map(
            lambda value: KO_SUPPORT_STATUS_LABELS.get(display_value(value), display_value(value))
        )
    if "oe_support_status" in frame.columns:
        frame["oe_support_status"] = frame["oe_support_status"].map(
            lambda value: OE_SUPPORT_STATUS_LABELS.get(display_value(value), display_value(value))
        )
    if "wet_lab_readiness" in frame.columns:
        frame["wet_lab_readiness"] = frame["wet_lab_readiness"].map(
            lambda value: WET_LAB_READINESS_LABELS.get(display_value(value), display_value(value))
        )
    for idx, row in frame.iterrows():
        label = screen_status_label(row.get("status"), row.get("success"))
        if not display_value(row.get("solver_status_label")):
            frame.at[idx, "solver_status_label"] = label
        if label != "求解成功" and not display_value(row.get("failure_reason")):
            frame.at[idx, "failure_reason"] = label
        if (
            display_value(row.get("status")) == "2"
            and display_value(row.get("effect_label")) in {"", "求解失败"}
        ):
            frame.at[idx, "effect_label"] = "约束不可行"
        if (
            display_value(row.get("status")) in {"not_run_complex_subunit_limited", "not_run_gene_oe_proxy"}
            and display_value(row.get("effect_label")) in {"", "求解失败", "未解析"}
        ):
            frame.at[idx, "effect_label"] = "未运行"
        if (
            display_value(row.get("status")) == "not_run_no_gpr_effect"
            and display_value(row.get("effect_label")) in {"", "求解失败", "未解析"}
        ):
            frame.at[idx, "effect_label"] = "未运行"
    return frame


def candidate_effect_counts(frame: pd.DataFrame) -> dict[str, int]:
    counts = {
        "提升分泌": 0,
        "降低分泌": 0,
        "约束不可行": 0,
        "求解失败": 0,
        "未运行": 0,
        "未解析": 0,
        "无明显变化": 0,
    }
    for _, row in frame.iterrows():
        effect = display_value(row.get("effect_label"), "未解析")
        status = display_value(row.get("status"))
        solver_status = display_value(row.get("solver_status_label") or row.get("failure_reason"))
        if status == "2" or "不可行" in effect or "不可行" in solver_status:
            counts["约束不可行"] += 1
        elif "提升" in effect:
            counts["提升分泌"] += 1
        elif "降低" in effect:
            counts["降低分泌"] += 1
        elif "无明显" in effect:
            counts["无明显变化"] += 1
        elif "失败" in effect or "失败" in solver_status:
            counts["求解失败"] += 1
        elif "未运行" in effect or "仅解释" in solver_status:
            counts["未运行"] += 1
        else:
            counts["未解析"] += 1
    return {key: value for key, value in counts.items() if value}


def candidate_row_label(index: int, row: pd.Series) -> str:
    gene = display_value(
        row.get("standard_gene_symbol")
        or row.get("display_name")
        or row.get("input_gene_id")
        or row.get("gene_id")
        or row.get("candidate_id")
    )
    effect = display_value(row.get("effect_label"), "未解析")
    delta = display_value(row.get("delta_objective"), "无可行目标值")
    return f"{index + 1}. {gene} | {effect} | Δ={delta}"


__all__ = [
    "CANDIDATE_DISPLAY_COLUMNS",
    "MAPPING_CONFIDENCE_LABELS",
    "MAPPING_LEVEL_LABELS",
    "GPR_ROLE_LABELS",
    "CAPACITY_EFFECT_LABELS",
    "SIMULATION_BASIS_LABELS",
    "KO_SUPPORT_STATUS_LABELS",
    "OE_SUPPORT_STATUS_LABELS",
    "WET_LAB_READINESS_LABELS",
    "candidate_effect_counts",
    "candidate_row_label",
    "display_value",
    "normalise_candidate_frame_for_display",
    "screen_status_label",
    "target_semantics_label",
]
