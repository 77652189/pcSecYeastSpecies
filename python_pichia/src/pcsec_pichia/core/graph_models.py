from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


NodeKind = Literal["target", "process", "module", "reaction", "enzyme", "gene", "constraint"]
EdgeKind = Literal["flows_to", "uses", "limited_by", "encoded_by", "affects"]
BottleneckLevel = Literal["unknown", "low", "medium", "high"]


@dataclass(frozen=True)
class SecretionGraphNode:
    node_id: str
    label: str
    kind: NodeKind
    evidence_type: str = ""
    evidence_id: str = ""
    bottleneck_level: BottleneckLevel = "unknown"
    note: str = ""


@dataclass(frozen=True)
class SecretionGraphEdge:
    source: str
    target: str
    kind: EdgeKind
    label: str = ""
    evidence_id: str = ""


@dataclass(frozen=True)
class SecretionGraph:
    graph_id: str
    target_id: str
    title: str
    nodes: list[SecretionGraphNode] = field(default_factory=list)
    edges: list[SecretionGraphEdge] = field(default_factory=list)
    note: str = "这是模型解释图，不是完整细胞代谢网络。"

    def node_ids(self) -> set[str]:
        return {node.node_id for node in self.nodes}

    def validation_errors(self) -> list[str]:
        node_ids = self.node_ids()
        errors: list[str] = []
        for edge in self.edges:
            if edge.source not in node_ids:
                errors.append(f"edge source not found: {edge.source}")
            if edge.target not in node_ids:
                errors.append(f"edge target not found: {edge.target}")
        return errors


def minimal_target_secretion_graph(target_id: str, target_label: str | None = None) -> SecretionGraph:
    label = target_label or target_id
    steps = [
        ("translation", "翻译"),
        ("er_translocation", "ER translocation"),
        ("folding", "folding / disulfide"),
        ("glycosylation", "N/O-glycosylation"),
        ("golgi_processing", "Golgi processing"),
        ("vesicle_transport", "vesicle transport"),
        ("extracellular_secretion", "extracellular secretion"),
    ]
    nodes = [
        SecretionGraphNode(
            node_id=f"target:{target_id}",
            label=label,
            kind="target",
            evidence_type="target_input",
            evidence_id=target_id,
        )
    ]
    nodes.extend(
        SecretionGraphNode(
            node_id=f"process:{step_id}",
            label=step_label,
            kind="process",
            evidence_type="pcSecPichia module",
            evidence_id=step_id,
        )
        for step_id, step_label in steps
    )
    edges: list[SecretionGraphEdge] = []
    previous = f"target:{target_id}"
    for step_id, step_label in steps:
        current = f"process:{step_id}"
        edges.append(SecretionGraphEdge(source=previous, target=current, kind="flows_to", label=step_label))
        previous = current
    return SecretionGraph(
        graph_id=f"{target_id}_minimal_secretion_path",
        target_id=target_id,
        title=f"{label} 分泌路径总览",
        nodes=nodes,
        edges=edges,
    )
