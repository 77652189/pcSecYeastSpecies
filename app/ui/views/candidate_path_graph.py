from __future__ import annotations

import streamlit as st

from app.ui.views.simulation_display import display_value


def _target_ptm_counts(summary: dict[str, object]) -> dict[str, object]:
    metadata = summary.get("target_metadata") if isinstance(summary, dict) else {}
    if isinstance(metadata, dict):
        return {
            "DSB": metadata.get("disulfide_sites", 0),
            "NG": metadata.get("n_glycosylation_sites", 0),
            "OG": metadata.get("o_glycosylation_sites", 0),
        }
    plan = summary.get("secretion_plan") if isinstance(summary, dict) else {}
    ptm = (plan or {}).get("ptm_counts") if isinstance(plan, dict) else {}
    if isinstance(ptm, dict):
        return {
            "DSB": ptm.get("disulfide_sites", 0),
            "NG": ptm.get("n_glycosylation_sites", 0),
            "OG": ptm.get("o_glycosylation_sites", 0),
        }
    return {"DSB": 0, "NG": 0, "OG": 0}


def render_secretion_path_graph(row: dict[str, object], summary: dict[str, object]) -> None:
    """用 graphviz 展示选定候选的分泌路径影响图。"""
    gene = display_value(row.get("input_gene_id") or row.get("gene_id") or row.get("candidate_id"))
    intervention = display_value(row.get("intervention_type"))
    effect = display_value(row.get("effect_label"), "未解析")
    delta = display_value(row.get("delta_objective"), "无可行目标值")
    process = display_value(row.get("secretory_process"))
    reaction = display_value(row.get("resolved_reaction_id") or row.get("reaction_id"))
    ptm = _target_ptm_counts(summary)
    ptm_label = f"目标 PTM\nDSB={ptm.get('DSB')} / NG={ptm.get('NG')} / OG={ptm.get('OG')}"

    nodes = {
        "Gene": {"label": f"基因/反应\n{gene[:30]}", "fill": "#e5e7eb"},
        "Intervention": {"label": f"扰动\n{intervention}", "fill": "#e5e7eb"},
        "Reaction": {"label": f"反应/复合体\n{reaction[:30]}", "fill": "#e5e7eb"},
        "Process": {"label": f"分泌环节\n{process[:30]}", "fill": "#dbeafe"},
        "PTM": {"label": ptm_label, "fill": "#fef3c7"},
        "Secretion": {
            "label": f"分泌通量变化\n{effect}",
            "fill": "#bbf7d0"
            if "提升" in effect
            else ("#fecaca" if "降低" in effect else ("#fde68a" if "不可行" in effect or "失败" in effect else "#e5e7eb")),
        },
        "Result": {"label": f"Δ目标值\n{delta}", "fill": "#e5e7eb"},
    }
    edges = [
        ("Gene", "Intervention"),
        ("Intervention", "Reaction"),
        ("Reaction", "Process"),
        ("Process", "PTM"),
        ("PTM", "Secretion"),
        ("Secretion", "Result"),
    ]

    dot = "digraph G {\n  rankdir=LR;\n  node [style=filled, shape=box, fontsize=11];\n"
    for node_id, info in nodes.items():
        dot += f'  {node_id} [label="{info["label"]}", fillcolor="{info["fill"]}"];\n'
    for src, dst in edges:
        dot += f"  {src} -> {dst};\n"
    dot += "}"

    try:
        import graphviz

        graph = graphviz.Source(dot)
        st.graphviz_chart(graph, use_container_width=True)
    except ImportError:
        st.graphviz_chart(dot, use_container_width=True)

    st.caption("路径图展示扰动从基因/反应到分泌通量的影响链路，非实际发酵网络。")


__all__ = ["render_secretion_path_graph"]
