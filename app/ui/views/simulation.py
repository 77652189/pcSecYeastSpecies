from __future__ import annotations

from pathlib import Path

import streamlit as st

from app.services.pichia_background_tasks import clear_last_result, submit_background_simulation
from app.services.pichia_secretion_schema import SecretionRunRequest
from app.ui.common import PATHS, request_navigation
from app.ui.views.simulation_builder import (
    TargetBuildFormState,
    render_conditions_and_constraints,
    render_target_selection,
)
from app.ui.views.simulation_gene_inputs import render_gene_perturbation_form
from app.ui.views.simulation_matlab_reference import render_matlab_reference
from app.ui.views.simulation_results import render_pichia_results


def _prefill_field_values(candidate_id: str, intervention_type: str, candidate_kind: str) -> dict[str, str]:
    """Pure routing decision behind apply_simulation_prefill - which draft field gets candidate_id.

    "gene" (full-genome screen rows, OE capacity gene ids) goes to the gene-ID inputs, which
    resolve through GPR; everything else ("catalog_reaction" from the curated catalog screen,
    "complex_oe_hypothesis" from the hypothetical whole-complex OE test, and any future
    non-gene candidate_kind) goes to the reaction-ID inputs, since those rows' id field
    actually holds a direct reaction id with no gene to resolve. Allowlisting "gene" rather
    than denylisting "catalog_reaction" so a new non-gene candidate_kind fails safe (reaction
    routing) instead of silently landing in the wrong box.
    """
    is_gene = candidate_kind == "gene"
    return {
        "pichia_draft_ko_genes": candidate_id if (is_gene and intervention_type == "KO") else "",
        "pichia_draft_ko_reactions": candidate_id if (not is_gene and intervention_type == "KO") else "",
        "pichia_draft_oe_genes": candidate_id if (is_gene and intervention_type == "OE") else "",
        "pichia_draft_oe_reactions": candidate_id if (not is_gene and intervention_type == "OE") else "",
    }


def apply_simulation_prefill(
    target_id: str, candidate_id: str, intervention_type: str, candidate_kind: str = "gene"
) -> None:
    """Pre-fill this page's session_state with one candidate and jump here.

    Shared by every page that wants a "verify this candidate in simulation"
    button (genome-wide screen, OE capacity comparison, ...), so there is one
    owner for how the draft build form gets pre-filled instead of each caller
    re-deriving the same routing rules.

    Replaces (not appends to) the KO/OE draft fields so the simulation isolates
    exactly this one candidate instead of mixing in whatever was left over from a
    previous exploration.
    """
    st.session_state["pichia_tab_selector"] = "仿真构建"
    st.session_state["pichia_draft_build_mode"] = "快速选择（内置模板）"
    st.session_state["pichia_template"] = target_id
    st.session_state.update(_prefill_field_values(candidate_id, intervention_type, candidate_kind))
    request_navigation("仿真验证")
    st.rerun()


def render_simulation() -> None:
    st.header("仿真验证")
    tab_key = "pichia_tab_selector"
    options = ["仿真构建", "仿真结果", "历史 MATLAB OPN 参考"]

    if st.session_state.get("pichia_draft_task_status_path"):
        st.radio("切换页面", options, index=1, horizontal=True, disabled=True, key="pichia_running_tab_display")
        st.caption("当前有仿真任务正在运行，完成前保持在结果页面。")
        render_pichia_results()
        return

    # Keep an active background task in the results view so users see progress.
    if st.session_state.pop("pichia_switch_to_results", False):
        st.session_state[tab_key] = "仿真结果"

    tab = st.radio("切换页面", options, horizontal=True, key=tab_key)
    if tab == "仿真构建":
        _render_pichia_builder()
    elif tab == "仿真结果":
        render_pichia_results()
    else:
        render_matlab_reference()


_BUILDER_STEP_KEY = "pichia_builder_step"
_BUILDER_STEPS = ("目标蛋白", "改造候选（KO/OE）", "培养条件与分析")


def _render_pichia_builder() -> None:
    # 分步向导：标签页地位平等、看不出先后，用户反馈"太不明显"。改成一次只展开当前步 +
    # 显式的上一步/下一步，运行按钮只在最后一步出现。
    #
    # 为什么三步全部渲染、只是折叠起来：Streamlit 在某个控件本轮未渲染时会清掉它的 session_state，
    # 只渲染当前步会导致"回到上一步时填过的内容全没了"，也会打断从筛查页跳来的预填。
    # 全渲染 + 只展开当前步，既有向导的引导感，又不丢状态。
    step = int(st.session_state.get(_BUILDER_STEP_KEY, 1))
    step = min(max(step, 1), len(_BUILDER_STEPS))
    st.caption(f"第 {step} / {len(_BUILDER_STEPS)} 步 — {_BUILDER_STEPS[step - 1]}")
    st.progress(step / len(_BUILDER_STEPS))

    with st.expander(f"① {_BUILDER_STEPS[0]}", expanded=step == 1):
        target_fields = render_target_selection()
    with st.expander(f"② {_BUILDER_STEPS[1]}", expanded=step == 2):
        gene_state = render_gene_perturbation_form(str(target_fields["target_id"]))
    with st.expander(f"③ {_BUILDER_STEPS[2]}", expanded=step == 3):
        condition_fields = render_conditions_and_constraints()
    build_state = TargetBuildFormState(**target_fields, **condition_fields)

    col_back, col_next, _ = st.columns([1, 1, 3])
    if step > 1 and col_back.button("← 上一步", key="pichia_builder_prev_step"):
        st.session_state[_BUILDER_STEP_KEY] = step - 1
        st.rerun()
    if step < len(_BUILDER_STEPS) and col_next.button("下一步 →", key="pichia_builder_next_step", type="primary"):
        st.session_state[_BUILDER_STEP_KEY] = step + 1
        st.rerun()

    if step < len(_BUILDER_STEPS):
        st.caption("走到第 3 步后会出现运行按钮。")
        return

    st.divider()
    out_dir = st.text_input("输出目录", value=str(PATHS.local_runs_dir/"streamlit_pichia_runs"), key="pichia_out")

    task_sp = st.session_state.get("pichia_draft_task_status_path")
    run_clicked = st.button(
        "运行 Python 分泌仿真",
        type="primary",
        disabled=task_sp is not None,
        key="pichia_run_simulation_button",
    )
    if st.button("清除上次结果", key="pichia_clear_last_result_button"):
        for k in ["last_pichia_secretion_draft_response","pichia_draft_task_status_path","pichia_draft_task_id"]:
            st.session_state.pop(k, None)
        clear_last_result(PATHS); st.rerun()

    if run_clicked and not task_sp:
        common = dict(enable_ribosome_translation_constraint=build_state.enable_ribosome, enable_misfolding_constraint=build_state.enable_misfolding,
                      enable_cost_slope_compatibility=build_state.enable_cost_slope_compatibility,
                      cost_slope_medium_compatibility_mode=build_state.cost_slope_medium_compatibility_mode,
                      enable_solver_robustness_check=build_state.enable_solver_robustness_check,
                      enable_oe_dose_response=build_state.enable_oe_dose_response,
                      mu=build_state.mu, media_type=build_state.media_type, carbon_source_id=build_state.carbon_source_id,
                      ko_gene_ids=gene_state.ko_gene_ids,
                      ko_reaction_ids=gene_state.ko_reaction_ids,
                      oe_gene_ids=gene_state.oe_gene_ids, oe_reaction_ids=gene_state.oe_reaction_ids,
                      screen_candidate_limit=gene_state.candidate_limit,
                      enable_gene_rule_overlay=gene_state.enable_gene_rule_overlay,
                      output_dir=Path(out_dir) if out_dir.strip() else None)
        # 自建模板没有内置 spec，不能按 id 解析——它保存的是显式序列，必须走 custom_sequence，
        # 否则运行时会按 id 找不到目标而失败。
        if build_state.build_mode == "快速选择（内置模板）" and not build_state.target_is_custom_library:
            req = SecretionRunRequest(target_source="builtin", target_id=build_state.target_id, target_name=build_state.target_name, **common)
        elif build_state.build_mode == "自定义 JSON":
            req = SecretionRunRequest(target_source="custom_json", target_id=build_state.target_id, target_name=build_state.target_name,
                custom_json_path=build_state.custom_json_path, **common)
        else:
            req = SecretionRunRequest(target_source="custom_sequence", target_id=build_state.target_id, target_name=build_state.target_name,
                sequence=build_state.mature_sequence, leader_sequence=build_state.leader_sequence, signal_peptide_sequence=build_state.signal_peptide_sequence,
                sequence_role="mature_secreted", normalization_mode="as_provided",
                disulfide_sites=build_state.disulfide_sites, n_glycosylation_sites=build_state.n_glycosylation_sites,
                o_glycosylation_sites=build_state.o_glycosylation_sites, **common)
        # Stash the strain-defining modifications so the results page can offer "下一步 OE 候选"
        # (a modified-strain re-solve). Reaction/complex-level KO/OE are what the re-solve applies;
        # gene-level entries are recorded so the results page can flag they are not yet re-solved.
        st.session_state["pichia_last_run_modifications"] = {
            "target_id": build_state.target_id,
            "target_is_builtin": build_state.build_mode == "快速选择（内置模板）",
            "carbon_source_id": build_state.carbon_source_id,
            "media_type": build_state.media_type,
            "mu": build_state.mu,
            "enable_ribosome": bool(build_state.enable_ribosome),
            "enable_misfolding": bool(build_state.enable_misfolding),
            "ko_reaction_ids": list(gene_state.ko_reaction_ids),
            "oe_reaction_ids": list(gene_state.oe_reaction_ids),
            "ko_gene_ids": list(gene_state.ko_gene_ids),
            "oe_gene_ids": list(gene_state.oe_gene_ids),
        }
        st.session_state.pop("pichia_next_oe_candidates_result", None)  # invalidate stale readout
        tid, tsp = submit_background_simulation(req, PATHS)
        st.session_state["pichia_draft_task_id"] = tid
        st.session_state["pichia_draft_task_status_path"] = tsp
        st.session_state.pop("last_pichia_secretion_draft_response", None)
        clear_last_result(PATHS)
        st.toast("任务已提交，跳转到结果页面…")
        st.session_state["pichia_switch_to_results"] = True
        st.rerun()
