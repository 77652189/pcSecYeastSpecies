from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import pandas as pd

from app.ui.common import HOMOLOGY_AUDIT_PAGE


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_homology_audit_page_is_registered_in_navigation_and_entrypoint() -> None:
    common_text = (REPO_ROOT / "app" / "ui" / "common.py").read_text(encoding="utf-8")
    entrypoint_text = (REPO_ROOT / "app" / "ui" / "streamlit_app.py").read_text(encoding="utf-8")

    assert HOMOLOGY_AUDIT_PAGE in common_text
    assert "render_homology_audit" in entrypoint_text
    assert "elif page == HOMOLOGY_AUDIT_PAGE" in entrypoint_text


def test_homology_audit_view_uses_service_facade_not_engine_runtime() -> None:
    view_path = REPO_ROOT / "app" / "ui" / "views" / "homology_audit.py"
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

    assert "app.services.pichia_homology_audit_service" in imported_modules
    assert not any(module.startswith(("pcsec_pichia", "python_pichia")) for module in imported_modules)
    assert not {"run_blastp", "make_blast_db", "makeblastdb"} & called_names
    assert "pcsec_pichia.homology.external_fetch" not in imported_modules
    assert not {"fetch_external_name_references", "fetch_uniprot_name_reference"} & called_names


def test_homology_audit_name_table_includes_external_crosscheck_columns() -> None:
    import app.ui.views.homology_audit as view

    assert "external_crosscheck_status" in view.NAME_AUDIT_COLUMNS
    assert "external_crosscheck_sources" in view.NAME_AUDIT_COLUMNS
    assert "external_crosscheck_warnings" in view.NAME_AUDIT_COLUMNS


def test_homology_audit_external_reference_columns_include_review_fields() -> None:
    import app.ui.views.homology_audit as view

    assert "evidence_kind" in view.EXTERNAL_REFERENCE_COLUMNS
    assert "function_description" in view.EXTERNAL_REFERENCE_COLUMNS
    assert "gpr_transfer_status" in view.EXTERNAL_REFERENCE_COLUMNS
    assert "external_model_sources" in view.EXTERNAL_REFERENCE_COLUMNS
    assert "gpr_source_priority" in view.EXTERNAL_REFERENCE_COLUMNS
    assert "external_gpr_mapping_status" in view.EXTERNAL_REFERENCE_COLUMNS
    assert "external_gpr_conflict_warnings" in view.EXTERNAL_REFERENCE_COLUMNS
    assert "manual_review_reasons" in view.EXTERNAL_REFERENCE_COLUMNS


def test_ready_rule_transfer_rows_only_include_model_operable_ready_candidates() -> None:
    import app.ui.views.homology_audit as view

    rows = [
        _rule_row("PAS_chr2-1_0140", "rule_transfer_ready", True),
        _rule_row("PAS_chr4_0844", "rule_transfer_supported_not_model_operable", True),
        _rule_row("PAS_chr4_0156", "rule_transfer_ready", False),
        _rule_row("", "rule_transfer_ready", True),
        _rule_row("PAS_chr2-1_0140", "rule_transfer_ready", True),
    ]

    ready = view._ready_rule_transfer_rows(rows)

    assert [row["pichia_model_gene_id"] for row in ready] == ["PAS_chr2-1_0140"]


def test_homology_rule_transfer_controls_add_ko_candidate_and_navigate(monkeypatch) -> None:
    import app.ui.views.homology_audit as view

    fake_st = _FakeStreamlit(
        multiselect_values=["PAS_chr2-1_0140 — KAR2 / BiP"],
        clicked_keys={"homology_add_ready_ko"},
    )
    fake_st.session_state["pichia_draft_ko_genes"] = "PAS_chr1-1_0001\nPAS_chr2-1_0140"
    navigation: list[str] = []
    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(view, "request_navigation", navigation.append)

    view._render_add_to_simulation_controls(
        [
            _rule_row("PAS_chr2-1_0140", "rule_transfer_ready", True, common_name="KAR2 / BiP"),
            _rule_row("PAS_chr4_0844", "rule_transfer_supported_not_model_operable", True),
        ]
    )

    assert fake_st.multiselect_calls[0]["options"] == ["PAS_chr2-1_0140 — KAR2 / BiP"]
    assert fake_st.session_state["pichia_draft_ko_genes"] == "PAS_chr1-1_0001\nPAS_chr2-1_0140"
    assert "pichia_draft_oe_genes" not in fake_st.session_state
    assert navigation == ["仿真验证"]
    assert fake_st.rerun_called is True
    assert "已加入 KO 输入" in fake_st.toasts[0]


def test_homology_rule_transfer_selection_adds_oe_candidate_with_dedupe(monkeypatch) -> None:
    import app.ui.views.homology_audit as view

    fake_st = _FakeStreamlit()
    fake_st.session_state["pichia_draft_oe_genes"] = "PAS_chr3_0001"
    navigation: list[str] = []
    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(view, "request_navigation", navigation.append)

    view._apply_rule_transfer_selection(
        [
            _rule_row("PAS_chr3_0001", "rule_transfer_ready", True),
            _rule_row("PAS_chr3_0001", "rule_transfer_ready", True),
            _rule_row("PAS_chr4_0844", "rule_transfer_low_confidence", True),
            _rule_row("PAS_chr5_0005", "rule_transfer_ready", True),
        ],
        action="oe",
    )

    assert fake_st.session_state["pichia_draft_oe_genes"] == "PAS_chr3_0001\nPAS_chr5_0005"
    assert "pichia_draft_ko_genes" not in fake_st.session_state
    assert navigation == ["仿真验证"]
    assert fake_st.rerun_called is True
    assert "已加入 OE 输入" in fake_st.toasts[0]


def test_homology_audit_page_handles_missing_cache_without_crashing(monkeypatch) -> None:
    import app.ui.views.homology_audit as view

    fake_st = _FakeStreamlit()
    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(
        view,
        "load_homology_audit_browser_data",
        lambda **filters: {
            "cache_status": {
                "cache_available": False,
                "recommended_build_command": "python scripts\\build_pichia_homology_cache.py --catalog-only",
                "missing_files": ["sce_to_pichia_name_audit.jsonl"],
            },
            "summary": {},
            "name_audit_rows": [],
            "rule_transfer_audit_rows": [],
        },
    )

    view.render_homology_audit()

    assert any("尚未找到可用的同源审计 cache" in item for item in fake_st.warning_calls)
    assert fake_st.code_calls == ["python scripts\\build_pichia_homology_cache.py --catalog-only"]


def test_homology_audit_page_gets_rows_from_service_and_renders_tables(monkeypatch) -> None:
    import app.ui.views.homology_audit as view

    calls: list[dict[str, Any]] = []
    fake_st = _FakeStreamlit(text_input_value="kar")

    def fake_load(**filters):
        calls.append(filters)
        return _payload()

    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(view, "load_homology_audit_browser_data", fake_load)

    view.render_homology_audit()

    assert calls[0]["query"] == "kar"
    assert calls[0]["min_identity"] is None
    assert len(fake_st.dataframes) == 3
    rendered_text = "\n".join(frame.to_csv(index=False) for frame in fake_st.dataframes)
    assert "KAR2 / BiP" in rendered_text
    assert "规则迁移·就绪" in rendered_text  # rule_transfer_status 已汉化（原 rule_transfer_ready）
    assert "gpr_candidate" in rendered_text  # evidence_kind 值未汉化，保持原样


def test_homology_audit_export_button_uses_service_export(monkeypatch) -> None:
    import app.ui.views.homology_audit as view

    exported: dict[str, Any] = {}
    fake_st = _FakeStreamlit()

    def fake_export(rows: list[dict[str, Any]], *, file_format: str = "tsv") -> bytes:
        exported["rows"] = rows
        exported["file_format"] = file_format
        return b"exported-by-service"

    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(view, "load_homology_audit_browser_data", lambda **filters: _payload())
    monkeypatch.setattr(view, "export_homology_audit_rows", fake_export)

    view.render_homology_audit()

    assert exported["rows"][0]["internal_common_name"] == "KAR2 / BiP"
    assert exported["file_format"] == "tsv"
    assert fake_st.downloads[0]["data"] == b"exported-by-service"


def test_homology_audit_external_reference_export_uses_external_service(monkeypatch) -> None:
    import app.ui.views.homology_audit as view

    exported: dict[str, Any] = {}
    fake_st = _FakeStreamlit(selectbox_values={"导出表": "外部数据库证据"})

    def fake_export(rows: list[dict[str, Any]], *, file_format: str = "tsv") -> bytes:
        exported["rows"] = rows
        exported["file_format"] = file_format
        return b"external-exported-by-service"

    monkeypatch.setattr(view, "st", fake_st)
    monkeypatch.setattr(view, "load_homology_audit_browser_data", lambda **filters: _payload())
    monkeypatch.setattr(view, "export_external_reference_rows", fake_export)

    view.render_homology_audit()

    assert exported["rows"][0]["evidence_kind"] == "gpr_candidate"
    assert exported["file_format"] == "tsv"
    assert fake_st.downloads[0]["data"] == b"external-exported-by-service"
    assert fake_st.downloads[0]["file_name"].endswith("external_reference_evidence.tsv")


def test_homology_audit_cache_tab_shows_external_cache_status(monkeypatch) -> None:
    import app.ui.views.homology_audit as view

    payload = _payload()
    fake_st = _FakeStreamlit()
    monkeypatch.setattr(view, "st", fake_st)

    view._render_cache_and_export(
        payload["cache_status"],
        payload["name_audit_rows"],
        payload["rule_transfer_audit_rows"],
    )

    rendered_markdown = "\n".join(fake_st.markdown_calls)
    assert "外部命名参考缓存" in rendered_markdown  # 原 External name reference cache
    assert "外部参考缓存" in rendered_markdown  # 原 External reference cache
    assert "UniProt:1" in rendered_markdown
    assert "gene_function:1" in rendered_markdown
    assert "yeast-GEM" in rendered_markdown
    assert "homology_supported_yeast_gpr" in rendered_markdown
    assert "external_match_confirmed: 1" in rendered_markdown
    assert any("build_pichia_external_name_reference_cache.py" in call for call in fake_st.code_calls)
    assert any("build_pichia_external_reference_cache.py" in call for call in fake_st.code_calls)


def _payload() -> dict[str, Any]:
    return {
        "cache_status": {
            "cache_available": True,
            "cache_root": "local_runs/pichia_homology_cache/smoke",
            "generated_at": "2026-07-07T12:00:00",
            "row_count": 1,
            "missing_files": [],
            "recommended_build_command": "python scripts\\build_pichia_homology_cache.py --catalog-only",
            "external_cache_available": True,
            "external_cache_path": "local_runs/pichia_homology_cache/smoke/external_name_references.jsonl",
            "external_cache_generated_at": "2026-07-08T12:00:00",
            "external_reference_count": 1,
            "external_sources": ["UniProt"],
            "external_source_counts": {"UniProt": 1},
            "external_cache_warnings": [],
            "recommended_external_build_command": (
                "python scripts\\build_pichia_external_name_reference_cache.py "
                "--name-audit-jsonl local_runs\\pichia_homology_cache\\smoke\\sce_to_pichia_name_audit.jsonl "
                "--output-path local_runs\\pichia_homology_cache\\smoke\\external_name_references.jsonl"
            ),
            "external_reference_cache": {
                "cache_available": True,
                "records_path": "local_runs/pichia_homology_cache/smoke/external_reference_records.jsonl",
                "record_count": 1,
                "source_counts": {"yeast-gem": 1},
                "record_type_counts": {"gene_function": 1, "gpr_candidate": 1},
                "external_model_sources": ["yeast-GEM"],
                "gpr_source_priority": {"best_priority_tier": "homology_supported_yeast_gpr"},
                "external_gpr_candidate_count": 1,
                "best_external_gpr_source": "yeast-gem:yeast-GEM",
                "external_gpr_mapping_status": {"gene_mapping_required": 1},
                "external_gpr_conflict_warnings": [],
                "retrieved_at_range": {
                    "first": "2026-07-09T00:00:00Z",
                    "last": "2026-07-09T00:00:00Z",
                },
                "recommended_refresh_command": (
                    "python scripts\\build_pichia_external_reference_cache.py "
                    "--sources uniprot,sgd --output-dir local_runs\\pichia_external_reference_cache\\smoke"
                ),
                "warnings": [],
            },
        },
        "summary": {
            "rule_transfer_row_count": 1,
            "rule_transfer_status_counts": {
                "rule_transfer_ready": 1,
                "rule_transfer_supported_not_model_operable": 0,
                "rule_transfer_low_confidence": 0,
                "rule_transfer_unresolved": 0,
                "rule_transfer_paralog_risk": 0,
                "rule_transfer_not_supported": 0,
            },
        },
        "name_audit_rows": [
            {
                "internal_common_name": "KAR2 / BiP",
                "internal_gene_id": "PAS_chr2-1_0140",
                "internal_sequence_id": "YJL034W",
                "external_gene_name": "KAR2",
                "external_locus_tag": "PAS_chr2-1_0140",
                "external_accession": "C4R",
                "identity_pct": 75.0,
                "query_coverage": 95.0,
                "subject_coverage": 95.0,
                "evalue": 1e-100,
                "is_rbh": True,
                "in_model_gene_index": True,
                "name_consistency_status": "name_confirmed_by_rbh",
                "review_status": "model_ready_rbh_high_confidence",
                "external_crosscheck_status": "external_match_confirmed",
                "external_crosscheck_sources": ["UniProt:2026_01:C4R"],
                "external_crosscheck_warnings": [],
                "warnings": [],
            }
        ],
        "rule_transfer_audit_rows": [
            {
                "internal_common_name": "KAR2 / BiP",
                "query_symbol": "KAR2",
                "sce_orf": "YJL034W",
                "pichia_gene_id": "PAS_chr2-1_0140",
                "pichia_model_gene_id": "PAS_chr2-1_0140",
                "identity_pct": 75.0,
                "query_coverage": 95.0,
                "subject_coverage": 95.0,
                "evalue": 1e-100,
                "is_rbh": True,
                "in_model_gene_index": True,
                "homology_review_status": "model_ready_rbh_high_confidence",
                "rule_transfer_status": "rule_transfer_ready",
                "warnings": [],
            }
        ],
        "external_reference_rows": [
            {
                "evidence_kind": "gpr_candidate",
                "source_database": "yeast-gem",
                "query_gene_id": "YJL034W",
                "external_model_sources": ["yeast-GEM"],
                "gpr_source_priority": "homology_supported_yeast_gpr",
                "external_gpr_mapping_status": "gene_mapping_required",
                "external_gpr_conflict_warnings": [],
                "gpr_transfer_status": "gene_mapping_required",
                "manual_review_reasons": ["external gene rule is not mapped"],
            }
        ],
    }


def _rule_row(
    model_gene_id: str,
    rule_transfer_status: str,
    in_model_gene_index: bool,
    *,
    common_name: str = "KAR2",
) -> dict[str, Any]:
    return {
        "internal_common_name": common_name,
        "query_symbol": common_name.split()[0],
        "sce_orf": "YJL034W",
        "pichia_gene_id": model_gene_id or "PAS_chr_unmodelled",
        "pichia_model_gene_id": model_gene_id,
        "identity_pct": 75.0,
        "query_coverage": 95.0,
        "subject_coverage": 95.0,
        "evalue": 1e-100,
        "is_rbh": True,
        "in_model_gene_index": in_model_gene_index,
        "homology_review_status": "model_ready_rbh_high_confidence",
        "rule_transfer_status": rule_transfer_status,
        "warnings": [],
    }


class _FakeStreamlit:
    def __init__(
        self,
        *,
        text_input_value: str = "",
        multiselect_values: list[str] | None = None,
        clicked_keys: set[str] | None = None,
        selectbox_values: dict[str, Any] | None = None,
    ) -> None:
        self.text_input_value = text_input_value
        self.multiselect_values = multiselect_values or []
        self.clicked_keys = clicked_keys or set()
        self.selectbox_values = selectbox_values or {}
        self.session_state: dict[str, Any] = {}
        self.warning_calls: list[str] = []
        self.info_calls: list[str] = []
        self.code_calls: list[str] = []
        self.markdown_calls: list[str] = []
        self.dataframes: list[pd.DataFrame] = []
        self.downloads: list[dict[str, Any]] = []
        self.metrics: list[tuple[str, Any]] = []
        self.tab_labels: list[str] = []
        self.multiselect_calls: list[dict[str, Any]] = []
        self.button_calls: list[dict[str, Any]] = []
        self.toasts: list[str] = []
        self.rerun_called = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False

    def header(self, *args, **kwargs) -> None:
        pass

    def caption(self, *args, **kwargs) -> None:
        pass

    def markdown(self, value: str, *args, **kwargs) -> None:
        self.markdown_calls.append(value)

    def subheader(self, *args, **kwargs) -> None:
        pass

    def warning(self, value: str, *args, **kwargs) -> None:
        self.warning_calls.append(value)

    def info(self, value: str, *args, **kwargs) -> None:
        self.info_calls.append(value)

    def code(self, value: str, *args, **kwargs) -> None:
        self.code_calls.append(value)

    def text_input(self, *args, **kwargs) -> str:
        return self.text_input_value

    def selectbox(self, label: str, options, index: int = 0, *args, **kwargs):
        if label in self.selectbox_values:
            return self.selectbox_values[label]
        return list(options)[index]

    def multiselect(self, label: str, options, *args, **kwargs):
        self.multiselect_calls.append({"label": label, "options": list(options), **kwargs})
        return self.multiselect_values

    def slider(self, *args, **kwargs) -> float:
        return float(kwargs.get("value", 0.0))

    def columns(self, spec):
        count = spec if isinstance(spec, int) else len(spec)
        return [self for _ in range(count)]

    def metric(self, label: str, value: Any, *args, **kwargs) -> None:
        self.metrics.append((label, value))

    def tabs(self, labels: list[str]):
        self.tab_labels = labels
        return [self for _ in labels]

    def dataframe(self, frame: pd.DataFrame, *args, **kwargs) -> None:
        self.dataframes.append(frame)

    def button(self, label: str, *args, **kwargs) -> bool:
        self.button_calls.append({"label": label, **kwargs})
        return bool(kwargs.get("key") in self.clicked_keys and not kwargs.get("disabled", False))

    def download_button(self, label: str, data: bytes, *args, **kwargs) -> None:
        self.downloads.append({"label": label, "data": data, **kwargs})

    def toast(self, value: str, *args, **kwargs) -> None:
        self.toasts.append(value)

    def rerun(self) -> None:
        self.rerun_called = True
