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
    assert len(fake_st.dataframes) == 2
    rendered_text = "\n".join(frame.to_csv(index=False) for frame in fake_st.dataframes)
    assert "KAR2 / BiP" in rendered_text
    assert "rule_transfer_ready" in rendered_text


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


def _payload() -> dict[str, Any]:
    return {
        "cache_status": {
            "cache_available": True,
            "cache_root": "local_runs/pichia_homology_cache/smoke",
            "generated_at": "2026-07-07T12:00:00",
            "row_count": 1,
            "missing_files": [],
            "recommended_build_command": "python scripts\\build_pichia_homology_cache.py --catalog-only",
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
    }


class _FakeStreamlit:
    def __init__(self, *, text_input_value: str = "") -> None:
        self.text_input_value = text_input_value
        self.warning_calls: list[str] = []
        self.info_calls: list[str] = []
        self.code_calls: list[str] = []
        self.dataframes: list[pd.DataFrame] = []
        self.downloads: list[dict[str, Any]] = []
        self.metrics: list[tuple[str, Any]] = []
        self.tab_labels: list[str] = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False

    def header(self, *args, **kwargs) -> None:
        pass

    def caption(self, *args, **kwargs) -> None:
        pass

    def markdown(self, *args, **kwargs) -> None:
        pass

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
        return list(options)[index]

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

    def download_button(self, label: str, data: bytes, *args, **kwargs) -> None:
        self.downloads.append({"label": label, "data": data, **kwargs})
