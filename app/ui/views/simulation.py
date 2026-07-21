from __future__ import annotations

from pathlib import Path

import streamlit as st

from app.services.pichia_background_tasks import clear_last_result, submit_background_simulation
from app.services.pichia_secretion_schema import SecretionRunRequest
from app.ui.common import PATHS, request_navigation
from app.ui.views.simulation_builder import render_target_build_form
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


def _render_pichia_builder() -> None:
    st.markdown("""<div class="concept-box">毕赤酵母分泌仿真工作台。三段式构建目标蛋白，可选培养基、基因扰动。</div>""", unsafe_allow_html=True)

    build_state = render_target_build_form()
    gene_state = render_gene_perturbation_form(build_state.target_id)

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
        if build_state.build_mode == "快速选择（内置模板）":
            req = SecretionRunRequest(target_source="builtin", target_id=build_state.target_id, target_name=build_state.target_name, **common)
        elif build_state.build_mode == "三段式构建（自定义组合）":
            req = SecretionRunRequest(target_source="custom_sequence", target_id=build_state.target_id, target_name=build_state.target_name,
                sequence=build_state.mature_sequence, leader_sequence=build_state.leader_sequence, signal_peptide_sequence=build_state.signal_peptide_sequence,
                sequence_role="mature_secreted", normalization_mode="as_provided",
                disulfide_sites=build_state.disulfide_sites, n_glycosylation_sites=build_state.n_glycosylation_sites,
                o_glycosylation_sites=build_state.o_glycosylation_sites, **common)
        else:
            req = SecretionRunRequest(target_source="custom_json", target_id=build_state.target_id, target_name=build_state.target_name,
                custom_json_path=build_state.custom_json_path, **common)
        tid, tsp = submit_background_simulation(req, PATHS)
        st.session_state["pichia_draft_task_id"] = tid
        st.session_state["pichia_draft_task_status_path"] = tsp
        st.session_state.pop("last_pichia_secretion_draft_response", None)
        clear_last_result(PATHS)
        st.toast("任务已提交，跳转到结果页面…")
        st.session_state["pichia_switch_to_results"] = True
        st.rerun()
