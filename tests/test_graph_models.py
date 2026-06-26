from __future__ import annotations

from pcsec_pichia.core.graph_models import SecretionGraphEdge, minimal_target_secretion_graph


def test_minimal_target_secretion_graph_is_traceable() -> None:
    graph = minimal_target_secretion_graph("OPN", "骨桥蛋白 OPN")

    assert graph.target_id == "OPN"
    assert "不是完整细胞代谢网络" in graph.note
    assert graph.validation_errors() == []
    assert any(node.kind == "target" and node.evidence_id == "OPN" for node in graph.nodes)
    assert any("ER translocation" in node.label for node in graph.nodes)


def test_graph_validation_reports_missing_edge_nodes() -> None:
    graph = minimal_target_secretion_graph("hLF")
    broken = graph.__class__(
        graph_id=graph.graph_id,
        target_id=graph.target_id,
        title=graph.title,
        nodes=graph.nodes,
        edges=graph.edges + [SecretionGraphEdge(source="missing", target="target:hLF", kind="affects")],
    )

    assert broken.validation_errors() == ["edge source not found: missing"]
