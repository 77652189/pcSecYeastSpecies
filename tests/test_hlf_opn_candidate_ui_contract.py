from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_hlf_opn_candidate_panel_is_embedded_before_gene_textareas() -> None:
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_inputs.py").read_text(encoding="utf-8")

    assert "render_hlf_opn_candidate_panel" in source
    assert source.index("render_hlf_opn_candidate_panel(target_id)") < source.index(
        'key="pichia_draft_ko_genes"'
    )


def test_hlf_opn_candidate_panel_uses_app_facade_not_runtime_blast() -> None:
    view_path = REPO_ROOT / "app" / "ui" / "views" / "hlf_opn_candidate_panel.py"
    module_ast = ast.parse(view_path.read_text(encoding="utf-8"))
    imported_modules: set[str] = set()
    called_names: set[str] = set()
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                called_names.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                called_names.add(node.func.attr)

    assert "app.services.pichia_gene_catalog_service" in imported_modules
    assert not any(module.startswith(("pcsec_pichia", "python_pichia")) for module in imported_modules)
    assert "pcsec_pichia.homology.blast_runner" not in imported_modules
    assert "pcsec_pichia.homology.external_fetch" not in imported_modules
    assert not {"run_blastp", "make_blast_db", "fetch_external_name_references", "urlopen"} & called_names


def test_hlf_opn_target_context_normalizes_builtin_target_ids() -> None:
    from app.ui.views.hlf_opn_candidate_panel import target_context_for_hlf_opn_candidates

    assert target_context_for_hlf_opn_candidates("hLF") == "hLF"
    assert target_context_for_hlf_opn_candidates("HLF_PROJECT_710") == "hLF"
    assert target_context_for_hlf_opn_candidates("OPN_ALPHA_FULL_PROJECT") == "OPN"
    assert target_context_for_hlf_opn_candidates("custom_secreted_target") is None


def test_hlf_opn_executable_inputs_update_only_ko_and_oe_gene_boxes(monkeypatch) -> None:
    import app.ui.views.hlf_opn_candidate_panel as view

    fake_st = _FakeStreamlit()
    fake_st.session_state["pichia_draft_ko_genes"] = "PAS_EXISTING_KO"
    fake_st.session_state["pichia_draft_oe_genes"] = "PAS_EXISTING_OE"
    monkeypatch.setattr(view, "st", fake_st)

    added = view._apply_executable_candidate_inputs(
        {
            "ko_gene_ids": ["PAS_EXISTING_KO", "PAS_NEW_KO"],
            "oe_gene_ids": ["PAS_EXISTING_OE", "PAS_NEW_OE"],
            "excluded_count": 99,
            "warnings": ["not_in_model candidates are excluded"],
        }
    )

    assert added == {"ko": ["PAS_EXISTING_KO", "PAS_NEW_KO"], "oe": ["PAS_EXISTING_OE", "PAS_NEW_OE"]}
    assert fake_st.session_state["pichia_draft_ko_genes"] == "PAS_EXISTING_KO\nPAS_NEW_KO"
    assert fake_st.session_state["pichia_draft_oe_genes"] == "PAS_EXISTING_OE\nPAS_NEW_OE"
    assert "pichia_draft_ko_reactions" not in fake_st.session_state
    assert "pichia_draft_oe_reactions" not in fake_st.session_state


def test_hlf_opn_candidate_frame_shows_required_evidence_boundaries() -> None:
    from app.ui.views.hlf_opn_candidate_panel import _candidate_frame

    frame = _candidate_frame(
        [
            {
                "target_context": "hLF",
                "gene_id": "PAS_KAR2",
                "display_name": "KAR2",
                "standard_symbol": "KAR2",
                "protein_name": "ER chaperone",
                "external_ids": {"uniprot": "U-KAR2"},
                "operability_status": "model_oe_proxy_executable",
                "recommended_intervention": "OE",
                "evidence_type": "homology_auxiliary",
                "evidence_confidence": "model_executable_annotation_supported",
                "homology_review_status": "model_ready_rbh_high_confidence",
                "rule_transfer_status": "rule_transfer_ready",
                "warnings": ["OE candidates are reaction-level proxies"],
            }
        ]
    )

    assert {
        "目标",
        "gene_id",
        "标准命名",
        "标准符号",
        "蛋白注释",
        "外部ID",
        "模型可操作性",
        "扰动",
        "证据类型",
        "证据置信度",
        "同源状态",
        "规则迁移",
        "限制",
    }.issubset(frame.columns)
    rendered = frame.to_csv(index=False)
    assert "PAS_KAR2" in rendered
    assert "reaction-level proxies" in rendered
    assert "mg/L" not in rendered


def test_hlf_opn_panel_filters_context_and_adds_only_executable_inputs(monkeypatch) -> None:
    import app.ui.views.hlf_opn_candidate_panel as view

    calls: dict[str, Any] = {}
    fake_st = _FakeStreamlit(clicked_keys={"hlf_opn_add_executable_hLF"})

    def fake_load_candidates(**kwargs):
        calls["candidate_kwargs"] = kwargs
        return [
            _candidate_row("PAS_KO", "model_ko_executable", "KO", model_operable=True),
            _candidate_row("PAS_OE", "model_oe_proxy_executable", "OE", model_operable=True),
            _candidate_row("PAS_NOT_MODEL", "not_in_model", "OE", model_operable=False),
        ]

    def fake_executable(**kwargs):
        calls["executable_kwargs"] = kwargs
        return {
            "target_context": "hLF",
            "ko_gene_ids": ["PAS_KO"],
            "oe_gene_ids": ["PAS_OE"],
            "excluded_count": 1,
            "warnings": ["not_in_model candidates are excluded from executable KO/OE inputs."],
        }

    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(view, "load_hlf_opn_candidate_genes", fake_load_candidates)
    monkeypatch.setattr(view, "hlf_opn_candidate_gene_summary", lambda: {"total_candidates": 3})
    monkeypatch.setattr(view, "hlf_opn_executable_candidate_inputs", fake_executable)
    monkeypatch.setattr(view, "load_hlf_opn_gpr_overlay_review", lambda **kwargs: [])

    view.render_hlf_opn_candidate_panel("hLF")

    assert calls["candidate_kwargs"]["target_context"] == "hLF"
    assert calls["candidate_kwargs"]["include_shared"] is True
    assert calls["executable_kwargs"]["target_context"] == "hLF"
    assert fake_st.session_state["pichia_draft_ko_genes"] == "PAS_KO"
    assert fake_st.session_state["pichia_draft_oe_genes"] == "PAS_OE"
    rendered = "\n".join(frame.to_csv(index=False) for frame in fake_st.dataframes)
    assert "PAS_NOT_MODEL" in rendered
    assert "not_in_model" in rendered


def _candidate_row(
    gene_id: str,
    operability_status: str,
    intervention: str,
    *,
    model_operable: bool,
) -> dict[str, Any]:
    return {
        "target_context": "hLF",
        "gene_id": gene_id,
        "display_name": gene_id,
        "standard_symbol": gene_id,
        "protein_name": "",
        "external_ids": {},
        "operability_status": operability_status,
        "recommended_intervention": intervention,
        "evidence_type": "curated_review",
        "evidence_confidence": "manual_review_required",
        "model_operable": model_operable,
        "homology_review_status": "",
        "rule_transfer_status": "",
        "warnings": [],
    }


class _FakeStreamlit:
    def __init__(self, *, clicked_keys: set[str] | None = None) -> None:
        self.clicked_keys = clicked_keys or set()
        self.session_state: dict[str, Any] = {}
        self.dataframes: list[pd.DataFrame] = []
        self.downloads: list[dict[str, Any]] = []
        self.metrics: list[tuple[str, Any]] = []
        self.info_calls: list[str] = []
        self.warning_calls: list[str] = []
        self.toast_calls: list[str] = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False

    def expander(self, *args, **kwargs):
        return self

    def columns(self, spec):
        count = spec if isinstance(spec, int) else len(spec)
        return [self for _ in range(count)]

    def caption(self, *args, **kwargs) -> None:
        pass

    def markdown(self, *args, **kwargs) -> None:
        pass

    def subheader(self, *args, **kwargs) -> None:
        pass

    def metric(self, label: str, value: Any, *args, **kwargs) -> None:
        self.metrics.append((label, value))

    def selectbox(self, label: str, options, index: int = 0, *args, **kwargs):
        return list(options)[index]

    def dataframe(self, frame: pd.DataFrame, *args, **kwargs) -> None:
        self.dataframes.append(frame)

    def button(self, label: str, *args, **kwargs) -> bool:
        return bool(kwargs.get("key") in self.clicked_keys and not kwargs.get("disabled", False))

    def download_button(self, label: str, data: bytes, *args, **kwargs) -> None:
        self.downloads.append({"label": label, "data": data, **kwargs})

    def info(self, value: str, *args, **kwargs) -> None:
        self.info_calls.append(value)

    def warning(self, value: str, *args, **kwargs) -> None:
        self.warning_calls.append(value)

    def toast(self, value: str, *args, **kwargs) -> None:
        self.toast_calls.append(value)
