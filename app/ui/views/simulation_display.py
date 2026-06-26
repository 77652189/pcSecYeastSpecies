from __future__ import annotations

import pandas as pd


CANDIDATE_DISPLAY_COLUMNS = {
    "input_gene_id": "输入基因",
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
    success_text = str(success).strip().lower()
    if success is True or success_text == "true":
        return "求解成功"
    status = display_value(raw_status)
    return {
        "2": "约束不可行",
        "3": "目标无界",
        "4": "求解器数值错误",
        "missing_reaction": "反应未找到",
        "unresolved_gene": "基因未解析",
        "unresolved_reaction": "反应未解析",
        "missing_objective": "目标反应未找到",
    }.get(status, "求解失败")


def normalise_candidate_frame_for_display(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    if "solver_status_label" not in frame.columns:
        frame["solver_status_label"] = ""
    if "failure_reason" not in frame.columns:
        frame["failure_reason"] = ""
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
    return frame


def candidate_effect_counts(frame: pd.DataFrame) -> dict[str, int]:
    counts = {
        "提升分泌": 0,
        "降低分泌": 0,
        "约束不可行": 0,
        "求解失败": 0,
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
        else:
            counts["未解析"] += 1
    return {key: value for key, value in counts.items() if value}


def candidate_row_label(index: int, row: pd.Series) -> str:
    gene = display_value(row.get("gene_id") or row.get("input_gene_id") or row.get("candidate_id"))
    effect = display_value(row.get("effect_label"), "未解析")
    delta = display_value(row.get("delta_objective"), "无可行目标值")
    return f"{index + 1}. {gene} | {effect} | Δ={delta}"


__all__ = [
    "CANDIDATE_DISPLAY_COLUMNS",
    "candidate_effect_counts",
    "candidate_row_label",
    "display_value",
    "normalise_candidate_frame_for_display",
    "screen_status_label",
    "target_semantics_label",
]
