from __future__ import annotations

import ast
import json
import os
import time
from pathlib import Path

import pandas as pd

from app.services.pichia_request_mapping_service import (
    request_warnings,
    sequence_contract_for_engine,
    target_input_payload,
)
from app.services.pichia_background_tasks import (
    BACKGROUND_TASK_STALE_SECONDS,
    load_latest_completed_background_result,
    load_last_result,
    poll_background_simulation,
    response_to_summary,
    status_path_for_background_task,
)
from app.services.pichia_screen_preview_service import _preview_screen_inputs_for_model
from app.services.pichia_secretion_schema import SecretionRunRequest, SecretionRunResponse
from app.services.pichia_secretion_runner import _ensure_pcsec_pichia_analysis_api
from app.services.pichia_target_catalog_service import (
    _builtin_target_semantics,
)
from app.services.pichia_target_catalog_service import (
    known_mature_proteins,
    known_signal_peptides,
)
from app.ui.views.simulation_display import (
    candidate_effect_counts,
    candidate_row_label,
    normalise_candidate_frame_for_display,
    target_semantics_label,
)
from app.ui.views.simulation_builder import medium_type_label
from app.ui.views.simulation_gene_inputs import gene_mapping_rows_for_display
from app.ui.views.simulation_gene_text import merge_candidate_text, parse_candidate_text
from app.services.pichia_secretion_service import (
    discover_project_paths,
    run_pichia_secretion_draft,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_public_pichia_secretion_facade_exposes_run_entrypoint_and_paths() -> None:
    assert callable(run_pichia_secretion_draft)
    paths = discover_project_paths()
    assert paths.repo_root == REPO_ROOT


def test_runner_refreshes_stale_pcsec_analysis_module(monkeypatch) -> None:
    import sys
    import types

    import app.services.pichia_secretion_runner as runner

    stale_module = types.ModuleType("pcsec_pichia.analysis")
    actual_module = sys.modules.get("pcsec_pichia.analysis")

    def fake_reload(module):
        module.analyze_target_growth_impact = lambda *args, **kwargs: None
        module.analyze_yield_improvement_candidates = lambda *args, **kwargs: None
        module.summarize_protein_cost_slope_compatibility = lambda *args, **kwargs: {}
        module.summarize_yield_improvement_recommendations = lambda *args, **kwargs: {}
        return module

    monkeypatch.setitem(sys.modules, "pcsec_pichia.analysis", stale_module)
    monkeypatch.setattr(runner.importlib, "reload", fake_reload)
    try:
        _ensure_pcsec_pichia_analysis_api()
        assert hasattr(stale_module, "analyze_target_growth_impact")
    finally:
        if actual_module is not None:
            monkeypatch.setitem(sys.modules, "pcsec_pichia.analysis", actual_module)


def test_pichia_secretion_facade_exports_only_reviewed_public_symbols() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_secretion_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    exported: list[str] | None = None
    for node in module_ast.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__all__":
                    exported = ast.literal_eval(node.value)

    assert exported == [
        "BuiltinTargetTemplate",
        "NormalizationMode",
        "SecretionRunRequest",
        "SecretionRunResponse",
        "SequenceRole",
        "TargetSource",
        "TerminalStopPolicy",
        "discover_project_paths",
        "poll_background_simulation",
        "run_pichia_secretion_draft",
        "status_path_for_background_task",
        "submit_background_simulation",
    ]


def test_pichia_secretion_facade_stays_thin_and_imports_owner_modules_only() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_secretion_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    function_names = [
        node.name for node in module_ast.body if isinstance(node, ast.FunctionDef)
    ]
    class_names = [node.name for node in module_ast.body if isinstance(node, ast.ClassDef)]
    imported_modules: set[str] = set()
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)

    assert function_names == [
        "discover_project_paths",
        "run_pichia_secretion_draft",
    ]
    assert class_names == []
    assert imported_modules == {
        "__future__",
        "pathlib",
        "app",
        "pcsec_pichia.core.paths",
        "app.services.pichia_background_tasks",
        "app.services.pichia_secretion_schema",
        "app.services.pichia_request_mapping_service",
        "app.services.pichia_secretion_runner",
    }
    assert not any(module_name.startswith("app.ui") for module_name in imported_modules)
    assert not any(module_name.startswith("app.api") for module_name in imported_modules)
    assert not any(module_name.startswith("app.engines") for module_name in imported_modules)


def test_pichia_app_services_use_central_python_pichia_bootstrap() -> None:
    offenders: list[str] = []
    for path in (REPO_ROOT / "app" / "services").glob("pichia_*.py"):
        source = path.read_text(encoding="utf-8")
        if "sys.path" in source or "import sys" in source:
            offenders.append(str(path.relative_to(REPO_ROOT)))

    assert offenders == []


def test_target_catalog_service_uses_formal_targets_not_probe_private_module() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_target_catalog_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    imported_modules: set[str] = set()
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)

    signal_peptides = known_signal_peptides()
    mature_proteins = known_mature_proteins()

    assert "pcsec_pichia.probe" not in imported_modules
    assert "pcsec_pichia.probe._prototype" not in imported_modules
    assert "pcsec_pichia.targets" in imported_modules
    assert signal_peptides["native_hLF"]["length"] == 19
    assert mature_proteins["hLF"]["length"] == 691
    assert "用户提供" in str(mature_proteins["hLF"]["source"])
    assert mature_proteins["OPN_ALPHA_FULL_PROJECT"]["length"] == 298


def test_pichia_app_services_do_not_import_probe_private_modules() -> None:
    offenders: list[str] = []
    for path in (REPO_ROOT / "app" / "services").glob("pichia_*.py"):
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(module_ast):
            imported: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
            elif isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            for module_name in imported:
                if module_name.startswith("pcsec_pichia.probe"):
                    offenders.append(f"{path.relative_to(REPO_ROOT)}: {module_name}")

    assert offenders == []


def test_service_contract_uses_facade_for_public_entrypoints_only() -> None:
    test_path = Path(__file__)
    module_ast = ast.parse(test_path.read_text(encoding="utf-8"))
    facade_imports: list[str] = []
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module == "app.services.pichia_secretion_service":
            facade_imports.extend(alias.name for alias in node.names)

    assert sorted(facade_imports) == [
        "discover_project_paths",
        "run_pichia_secretion_draft",
    ]


def test_streamlit_ui_does_not_import_engine_directly() -> None:
    direct_engine_imports: list[str] = []
    for path in (REPO_ROOT / "app" / "ui").rglob("*.py"):
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(module_ast):
            imported: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
            elif isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            for module_name in imported:
                if module_name.startswith(("pcsec_pichia", "python_pichia")):
                    direct_engine_imports.append(f"{path.relative_to(REPO_ROOT)}: {module_name}")

    assert direct_engine_imports == []


def test_python_draft_streamlit_views_do_not_import_legacy_opn_service() -> None:
    draft_view_paths = [
        REPO_ROOT / "app" / "ui" / "views" / "simulation.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_builder.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_display.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_catalog.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_inputs.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_text.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py",
        REPO_ROOT / "app" / "ui" / "views" / "candidate_path_graph.py",
    ]
    legacy_imports: list[str] = []
    for path in draft_view_paths:
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(module_ast):
            imported: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
            elif isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            for module_name in imported:
                if module_name in {"app.services.opn", "app.adapters.matlab"}:
                    legacy_imports.append(f"{path.relative_to(REPO_ROOT)}: {module_name}")

    assert legacy_imports == []


def test_python_draft_streamlit_views_use_owner_services_not_fat_facade() -> None:
    draft_view_paths = [
        REPO_ROOT / "app" / "ui" / "views" / "simulation.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_builder.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_display.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_catalog.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_inputs.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_text.py",
        REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py",
        REPO_ROOT / "app" / "ui" / "views" / "candidate_path_graph.py",
    ]
    fat_facade_imports: list[str] = []
    for path in draft_view_paths:
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(module_ast):
            if (
                isinstance(node, ast.ImportFrom)
                and node.module == "app.services.pichia_secretion_service"
            ):
                fat_facade_imports.append(str(path.relative_to(REPO_ROOT)))

    assert fat_facade_imports == []


def test_background_task_cache_ignores_corrupt_utf8_json(tmp_path: Path) -> None:
    from pcsec_pichia.core.paths import ProjectPaths

    paths = ProjectPaths(repo_root=tmp_path)
    cache_path = tmp_path / "local_runs" / "streamlit_pichia_runs" / ".last_result" / "result.json"
    cache_path.parent.mkdir(parents=True)
    cache_path.write_bytes(b"\xff\xfe\x00bad")

    status_path = tmp_path / "status.json"
    status_path.write_bytes(b"\xff\xfe\x00bad")

    assert load_last_result(paths) is None
    status, message, result = poll_background_simulation(status_path)
    assert status == "lost"
    assert message
    assert result is None


def test_background_task_poll_marks_old_running_status_stale(tmp_path: Path) -> None:
    status_path = tmp_path / "status.json"
    status_path.write_text(
        json.dumps({"status": "running", "message": "running"}, ensure_ascii=False),
        encoding="utf-8",
    )
    stale_time = time.time() - BACKGROUND_TASK_STALE_SECONDS - 5
    os.utime(status_path, (stale_time, stale_time))

    status, message, result = poll_background_simulation(status_path)

    assert status == "stale"
    assert "长时间未更新" in message
    assert result is None


def test_background_task_loader_recovers_latest_completed_result(tmp_path: Path) -> None:
    from pcsec_pichia.core.paths import ProjectPaths

    paths = ProjectPaths(repo_root=tmp_path)
    task_root = tmp_path / "local_runs" / "streamlit_pichia_runs" / ".background_tasks"
    old_dir = task_root / "old"
    new_dir = task_root / "new"
    old_dir.mkdir(parents=True)
    new_dir.mkdir(parents=True)
    (old_dir / "status.json").write_text(
        json.dumps({"status": "done", "result": {"target_id": "old"}}, ensure_ascii=False),
        encoding="utf-8",
    )
    (new_dir / "status.json").write_text(
        json.dumps({"status": "done", "result": {"target_id": "new"}}, ensure_ascii=False),
        encoding="utf-8",
    )
    old_time = time.time() - 100
    new_time = time.time()
    os.utime(old_dir / "status.json", (old_time, old_time))
    os.utime(new_dir / "status.json", (new_time, new_time))

    latest = load_latest_completed_background_result(paths)

    assert latest == {"target_id": "new"}


def test_background_task_status_path_is_scoped_to_project_paths(tmp_path: Path) -> None:
    from pcsec_pichia.core.paths import ProjectPaths

    first = ProjectPaths(repo_root=tmp_path / "first")
    second = ProjectPaths(repo_root=tmp_path / "second")

    first_status = status_path_for_background_task("task-a", first)
    second_status = status_path_for_background_task("task-b", second)

    assert first_status.parent.parent == first.local_runs_dir / "streamlit_pichia_runs" / ".background_tasks"
    assert second_status.parent.parent == second.local_runs_dir / "streamlit_pichia_runs" / ".background_tasks"
    assert first_status != second_status


def test_curated_gene_catalog_supports_advanced_oe_reaction_proxy_inputs() -> None:
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_catalog.py").read_text(encoding="utf-8")

    assert "get_pichia_oe_reactions_for_selection" in source
    assert "pichia_draft_oe_reactions" in source
    assert "添加到过表达反应代理" in source
    assert "不是 gene-level 扰动" in source
    assert "没有可靠的敲除模型基因或 KO 反应 ID" in source


def test_simulation_view_reaches_legacy_matlab_only_through_reference_tab() -> None:
    view_path = REPO_ROOT / "app" / "ui" / "views" / "simulation.py"
    module_ast = ast.parse(view_path.read_text(encoding="utf-8"))
    imported_modules: list[str] = []
    imported_names_by_module: dict[str, list[str]] = {}
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.append(node.module)
            imported_names_by_module.setdefault(node.module, []).extend(
                alias.name for alias in node.names
            )
        elif isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)

    assert "app.ui.views.simulation_matlab_reference" in imported_modules
    assert imported_names_by_module["app.ui.views.simulation_matlab_reference"] == [
        "render_matlab_reference"
    ]
    assert "app.adapters.matlab" not in imported_modules
    assert "app.services.opn" not in imported_modules
    assert not any(module_name.startswith("app.engines") for module_name in imported_modules)


def test_simulation_run_button_switches_to_results_page() -> None:
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation.py").read_text(encoding="utf-8")

    assert 'key="pichia_run_simulation_button"' in source
    assert 'key="pichia_clear_last_result_button"' in source
    assert 'st.session_state.get("pichia_draft_task_status_path")' in source
    assert 'st.session_state.pop("pichia_switch_to_results", False)' in source
    assert 'st.session_state[tab_key] = "仿真结果"' in source
    assert 'st.session_state["pichia_switch_to_results"] = True' in source
    assert "st.rerun()" in source
    results_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py").read_text(encoding="utf-8")
    assert "刷新仿真状态" in results_source
    assert "load_latest_completed_background_result(PATHS)" in results_source
    assert "time.sleep" not in results_source


def test_streamlit_display_helpers_localize_candidate_status_without_engine_logic() -> None:
    frame = pd.DataFrame(
        [
            {"status": "2", "success": False, "effect_label": "求解失败"},
            {
                "status": "optimal",
                "success": True,
                "effect_label": "提升分泌",
                "mapping_level": "complex_subunit",
                "mapping_confidence": "medium",
                "gpr_role": "complex_subunit",
                "capacity_effect": "complex_subunit_limited",
                "simulation_basis": "explain_only",
            },
        ]
    )

    display_frame = normalise_candidate_frame_for_display(frame)
    counts = candidate_effect_counts(display_frame)

    assert display_frame.loc[0, "solver_status_label"] == "约束不可行"
    assert display_frame.loc[0, "effect_label"] == "约束不可行"
    assert display_frame.loc[1, "mapping_level"] == "复合体亚基"
    assert display_frame.loc[1, "mapping_confidence"] == "中"
    assert display_frame.loc[1, "gpr_role"] == "复合体亚基"
    assert display_frame.loc[1, "capacity_effect"] == "复合体亚基受限"
    assert display_frame.loc[1, "simulation_basis"] == "仅解释"
    assert counts == {"提升分泌": 1, "约束不可行": 1}
    assert target_semantics_label("project_defined_hLF") == "项目定义 hLF（用户提供序列）"


def test_streamlit_gene_input_text_helpers_dedupe_multiline_candidates() -> None:
    parsed = parse_candidate_text("G1, G2\nG1\uFF1BG3\uFF0CG2")

    assert parsed == ("G1", "G2", "G3")
    assert merge_candidate_text("G1\nG2", ["G2", "G4"]) == "G1\nG2\nG4"


def test_streamlit_gene_perturbation_help_marks_ko_as_gene_level() -> None:
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_inputs.py").read_text(encoding="utf-8")

    assert "正式基因级 KO" in source
    assert "按 GPR 规则关闭会失活的反应" in source
    assert "reaction-level OE proxy" in source


def test_gene_rule_overlay_is_explicit_experimental_request_option() -> None:
    schema_request = SecretionRunRequest(target_source="builtin", target_id="OPN_ALPHA_FULL_PROJECT")
    ui_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_inputs.py").read_text(encoding="utf-8")
    preview_source = (REPO_ROOT / "app" / "services" / "pichia_screen_preview_service.py").read_text(encoding="utf-8")
    catalog_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_catalog.py").read_text(encoding="utf-8")

    assert schema_request.enable_gene_rule_overlay is False
    assert "使用外部证据补充 GPR（实验性，默认关闭）" in ui_source
    assert "不会写回原始模型，也不是 MATLAB 原始 GPR" in ui_source
    assert "build_gpr_overlay" in preview_source
    assert "apply_gpr_overlay_for_analysis" in preview_source
    assert "proposed_rule" not in preview_source
    assert "proposed_gr_rule" not in preview_source
    assert "候选 locus tag" in catalog_source
    assert "GPR 补充状态" in catalog_source
    assert "补充建议" in catalog_source


def test_streamlit_medium_type_labels_use_composition_names_not_internal_numbers() -> None:
    assert medium_type_label(2) == "YNB 基础培养基（维生素，无氨基酸）"
    assert medium_type_label(4) == "YNB + 核心氨基酸（15 种，默认）"
    assert medium_type_label(5) == "YNB + 全氨基酸（20 种）"
    assert "media_type=99" in medium_type_label(99)


def test_streamlit_cost_slope_option_explains_matlab_compatibility_route() -> None:
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_builder.py").read_text(encoding="utf-8")
    results_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py").read_text(encoding="utf-8")

    assert "启用蛋白成本斜率对比（MATLAB 历史路线，可选，较慢）" in source
    assert "当前默认路线" in source
    assert "最大化目标蛋白分泌通量" in source
    assert "历史 MATLAB 成本路线" in source
    assert "固定生长率 μ" in source
    assert "固定一组目标蛋白分泌比例" in source
    assert "优化葡萄糖摄取反应 Ex_glc_D" in source
    assert "不会替换或改变当前默认 corrected pipeline 的数值结果" in source
    assert "capacity_fraction_ratios" in results_source
    assert "按当前 corrected 分泌 capacity" in results_source


def test_python_draft_service_does_not_depend_on_legacy_app_engines() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_secretion_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    imported_modules: list[str] = []
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.append(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)

    assert not any(module_name.startswith("app.engines") for module_name in imported_modules)


def test_legacy_matlab_runtime_imports_stay_in_reference_boundaries() -> None:
    matlab_adapter_imports: list[str] = []
    legacy_engine_imports: list[str] = []
    for path in (REPO_ROOT / "app").rglob("*.py"):
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(module_ast):
            imported: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
            elif isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            for module_name in imported:
                relative_path = str(path.relative_to(REPO_ROOT))
                if module_name == "app.adapters.matlab":
                    matlab_adapter_imports.append(relative_path)
                if module_name.startswith("app.engines"):
                    legacy_engine_imports.append(relative_path)

    assert sorted(set(matlab_adapter_imports)) == [
        "app\\engines\\matlab_pichia_engine.py",
        "app\\services\\opn.py",
        "app\\services\\simulation.py",
        "app\\ui\\views\\simulation_matlab_reference.py",
    ]
    assert sorted(set(legacy_engine_imports)) == [
        "app\\engines\\__init__.py",
        "app\\services\\opn.py",
    ]


def test_custom_sequence_payload_strips_terminal_stop_and_records_contract() -> None:
    request = SecretionRunRequest(
        target_source="custom_sequence",
        target_id="CUSTOM_STOP",
        target_name="custom stop",
        sequence="ACD*",
        leader_sequence="MM AA",
        signal_peptide_sequence="MM",
        sequence_role="mature_secreted",
        normalization_mode="remove_terminal_stop",
        contains_signal_peptide=True,
        contains_leader=True,
        terminal_stop_policy="strip",
        disulfide_sites=1,
    )

    payload = target_input_payload(request)
    contract = sequence_contract_for_engine(request)

    assert payload["mature_sequence"] == "ACD"
    assert payload["leader_sequence"] == "MMAA"
    assert contract["sequence_role"] == "mature_secreted"
    assert contract["normalization_mode"] == "remove_terminal_stop"
    assert contract["terminal_stop_policy"] == "strip"
    assert contract["contains_signal_peptide"] is True
    assert contract["contains_leader"] is True
    assert contract["original_sequence_length"] == 4
    assert contract["normalized_sequence_length"] == 3
    assert contract["original_full_sequence_length"] == 8
    assert contract["normalized_full_sequence_length"] == 7
    assert contract["terminal_stop_present"] is True
    assert contract["terminal_stop_removed"] is True


def test_custom_sequence_request_warnings_explain_ambiguous_input() -> None:
    request = SecretionRunRequest(
        target_source="custom_sequence",
        target_id="CUSTOM_AMBIGUOUS",
        sequence="AC D?*",
        sequence_role="unknown",
        normalization_mode="as_provided",
        terminal_stop_policy="allow_for_record_only",
        disulfide_sites=0,
        n_glycosylation_sites=0,
        o_glycosylation_sites=0,
    )

    warnings = request_warnings(request)

    assert any("DSB/NG/OG" in item and "不做智能推断" in item for item in warnings)
    assert any("序列角色为「未知」" in item for item in warnings)
    assert any("包含空白字符" in item for item in warnings)
    assert any("非标准氨基酸字符" in item and "?" in item for item in warnings)
    assert any("DSB/NG/OG 均为 0" in item for item in warnings)
    assert any("序列末尾包含终止符 *" in item for item in warnings)


def test_screen_input_preview_resolves_manual_ko_oe_candidates() -> None:
    class TinyModel:
        rxns = ["R1", "R2"]
        rules = ["x(1)", "x(1) | x(2)"]
        gr_rules = ["G1", "G1 or G2"]
        gene_index = {"G1": 0, "G2": 1}
        reaction_index = {"R1": 0, "R2": 1}

    request = SecretionRunRequest(
        target_source="builtin",
        target_id="OPN_ALPHA_FULL_PROJECT",
        ko_gene_ids=("G1", "NO_SUCH_GENE"),
        ko_reaction_ids=("R1", "NO_SUCH_KO_RXN"),
        oe_gene_ids=("G1", "NO_SUCH_OE_GENE"),
        oe_reaction_ids=("R2", "NO_SUCH_OE_RXN"),
        screen_candidate_limit=2,
    )

    preview = _preview_screen_inputs_for_model(TinyModel(), request)

    assert preview["candidate_limit"] == 2
    assert preview["ko_genes"][0]["status"] == "resolved"
    assert preview["ko_genes"][1]["status"] == "unresolved_gene"
    assert preview["ko_reactions"][0]["status"] == "resolved"
    assert preview["ko_reactions"][1]["status"] == "unresolved_reaction"
    assert preview["oe_genes"][0]["intervention_type"] == "OE_gene_proxy"
    assert preview["oe_genes"][0]["resolved_reactions_preview"] == ["R1", "R2"]
    assert preview["oe_genes"][0]["simulation_basis"] == "reaction_level_capacity_proxy"
    assert preview["oe_genes"][0]["capacity_effect"] == "reaction_capacity_proxy"
    assert preview["oe_genes"][1]["status"] == "unresolved_gene"
    assert preview["oe_reactions"][0]["status"] == "resolved"
    assert preview["oe_reactions"][1]["status"] == "unresolved_reaction"
    assert preview["gene_mapping"]["genes"][0]["gene_id"] == "G1"
    assert preview["gene_mapping"]["genes"][0]["reaction_count"] == 2
    assert preview["gene_capabilities"][0]["gene_id"] == "G1"
    assert preview["gene_capabilities"][0]["ko_support_status"] == "ko_runnable_gpr_gene_deletion"
    assert preview["gene_capabilities"][0]["oe_support_status"] == "oe_runnable_reaction_proxy"
    assert preview["ko_genes"][0]["ko_support_status"] == "ko_runnable_gpr_gene_deletion"
    assert preview["oe_genes"][0]["oe_support_status"] == "oe_runnable_reaction_proxy"
    assert preview["ko_genes"][0]["recommendation_tier"] == "model_executable"
    assert preview["oe_genes"][0]["recommendation_tier"] == "model_executable"
    assert preview["oe_genes"][0]["oe_reaction_proxy"] is True
    assert "phenotype_evidence" in preview["ko_genes"][0]
    assert "database_annotation_sources" in preview["gene_capabilities"][0]
    assert preview["oe_genes"][0]["support_reason"]
    assert any(
        row["gene_id"] == "NO_SUCH_GENE"
        and row["mapping_level"] == "unresolved"
        and row["mapping_confidence"] == "unresolved"
        for row in preview["gene_mapping_rows"]
    )
    assert "GPR-aware" in preview["semantics"]["OE_gene_proxy"]
    assert any("敲除基因未在模型中找到" in item for item in preview["warnings"])
    assert any("GPR-aware planning + reaction-level proxy" in item for item in preview["warnings"])

    display_frame = gene_mapping_rows_for_display(preview["gene_mapping_rows"])
    assert {"基因", "反应", "分泌环节", "映射层级", "置信度", "解释"}.issubset(display_frame.columns)
    assert "未解析" in set(display_frame["置信度"])


def test_screen_input_preview_resolves_gene_aliases_with_offline_evidence(monkeypatch) -> None:
    from pcsec_pichia.services import gene_evidence
    from pcsec_pichia.services.gene_evidence import GeneExternalEvidence

    class TinyModel:
        rxns = ["R1"]
        rules = ["x(1)"]
        gr_rules = ["G1"]
        gene_index = {"G1": 0}
        reaction_index = {"R1": 0}

    monkeypatch.setattr(
        gene_evidence,
        "load_gene_evidence_cache",
        lambda *args, **kwargs: {
            "G1": GeneExternalEvidence(
                gene_id="G1",
                canonical_gene_id="G1",
                aliases=("ALIAS1",),
                evidence_sources=("offline_cache",),
            )
        },
    )
    request = SecretionRunRequest(
        target_source="builtin",
        target_id="OPN_ALPHA_FULL_PROJECT",
        ko_gene_ids=("ALIAS1",),
        oe_gene_ids=("ALIAS1",),
        screen_candidate_limit=2,
    )

    preview = _preview_screen_inputs_for_model(TinyModel(), request)

    assert preview["ko_genes"][0]["input_id"] == "ALIAS1"
    assert preview["ko_genes"][0]["canonical_gene_id"] == "G1"
    assert preview["ko_genes"][0]["status"] == "resolved"
    assert preview["oe_genes"][0]["input_id"] == "ALIAS1"
    assert preview["oe_genes"][0]["canonical_gene_id"] == "G1"
    assert preview["oe_genes"][0]["resolved_reactions_preview"] == ["R1"]
    assert preview["gene_capabilities"][0]["gene_id"] == "ALIAS1"
    assert preview["gene_capabilities"][0]["canonical_gene_id"] == "G1"
    assert preview["gene_mapping"]["genes"][0]["gene_id"] == "G1"
    assert preview["gene_mapping"]["genes"][0]["input_gene_ids"] == ["ALIAS1"]
    assert preview["gene_mapping_rows"][0]["input_gene_id"] == "ALIAS1"
    assert preview["gene_mapping_rows"][0]["canonical_gene_id"] == "G1"
    assert preview["gene_mapping_rows"][0]["mapping_level"] != "unresolved"
    assert any("基因别名 `ALIAS1` 已解析为模型基因 ID `G1`" in item for item in preview["warnings"])
    display_frame = gene_mapping_rows_for_display(preview["gene_mapping_rows"])
    assert display_frame.loc[0, "基因"] == "ALIAS1"
    assert display_frame.loc[0, "模型基因"] == "G1"


def test_screen_input_preview_loads_gene_evidence_from_repo_root_when_cwd_differs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    class TinyModel:
        rxns = ["R1"]
        rules = ["x(1)"]
        gr_rules = ["G1"]
        gene_index = {"G1": 0}
        reaction_index = {"R1": 0}

    cache_path = tmp_path / "local_runs" / "gene_evidence_cache" / "gene_evidence.json"
    cache_path.parent.mkdir(parents=True)
    cache_path.write_text(
        json.dumps(
            {
                "genes": [
                    {
                        "gene_id": "G1",
                        "canonical_gene_id": "G1",
                        "aliases": ["ALIAS1"],
                        "evidence_sources": ["offline_cache"],
                        "evidence_confidence": "high_exact_locus_tag",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    other_cwd = tmp_path / "not_repo_root"
    other_cwd.mkdir()
    monkeypatch.chdir(other_cwd)
    request = SecretionRunRequest(
        target_source="builtin",
        target_id="OPN_ALPHA_FULL_PROJECT",
        ko_gene_ids=("ALIAS1",),
        screen_candidate_limit=1,
    )

    preview = _preview_screen_inputs_for_model(TinyModel(), request, repo_root=tmp_path)

    assert preview["ko_genes"][0]["input_id"] == "ALIAS1"
    assert preview["ko_genes"][0]["canonical_gene_id"] == "G1"
    assert preview["ko_genes"][0]["status"] == "resolved"
    assert preview["ko_genes"][0]["database_annotation_sources"] == ["offline_cache"]
    assert preview["gene_capabilities"][0]["canonical_gene_id"] == "G1"


def test_screen_input_preview_uses_canonical_overlay_locus_for_capability(monkeypatch) -> None:
    from pcsec_pichia.services import gene_rule_overlay
    from pcsec_pichia.services.gene_rule_overlay import HIGH_CONFIDENCE, GeneRuleEvidence

    class TinyModel:
        rxns = ["sec_PDI1_ERV2_Ero1p_complex_formation"]
        rules = [""]
        gr_rules = [""]
        genes = ["G_ORIGINAL"]
        gene_index = {"G_ORIGINAL": 0}
        reaction_index = {"sec_PDI1_ERV2_Ero1p_complex_formation": 0}

    evidence = {
        "PDI1": GeneRuleEvidence(
            common_name="PDI1",
            candidate_locus_tag="PAS_PDI1",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        ),
        "ERO1": GeneRuleEvidence(
            common_name="ERO1",
            candidate_locus_tag="PAS_ERO1",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        ),
        "ERV2": GeneRuleEvidence(
            common_name="ERV2",
            candidate_locus_tag="PAS_ERV2",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        ),
    }
    monkeypatch.setattr(gene_rule_overlay, "load_gene_rule_evidence_cache", lambda *args, **kwargs: evidence)
    request = SecretionRunRequest(
        target_source="builtin",
        target_id="OPN_ALPHA_FULL_PROJECT",
        ko_gene_ids=("PDI1",),
        oe_gene_ids=("ERO1",),
        screen_candidate_limit=2,
        enable_gene_rule_overlay=True,
    )

    preview = _preview_screen_inputs_for_model(TinyModel(), request)

    assert preview["gene_rule_overlay"]["entry_count"] == 1
    assert preview["ko_genes"][0]["input_id"] == "PDI1"
    assert preview["ko_genes"][0]["canonical_gene_id"] == "PAS_PDI1"
    assert preview["ko_genes"][0]["status"] == "resolved"
    assert preview["ko_genes"][0]["ko_support_status"] == "ko_runnable_gpr_gene_deletion"
    assert preview["oe_genes"][0]["input_id"] == "ERO1"
    assert preview["oe_genes"][0]["canonical_gene_id"] == "PAS_ERO1"
    assert preview["oe_genes"][0]["status"] == "not_run_complex_subunit_limited"
    assert preview["oe_genes"][0]["oe_support_status"] == "oe_explain_only_complex_subunit"
    assert preview["gene_mapping_rows"][0]["input_gene_id"] == "PDI1"
    assert preview["gene_mapping_rows"][0]["canonical_gene_id"] == "PAS_PDI1"


def test_screen_preview_and_pipeline_share_engine_candidate_resolution_helpers() -> None:
    preview_path = REPO_ROOT / "app" / "services" / "pichia_screen_preview_service.py"
    pipeline_path = REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "pipeline.py"

    def module_ast_for(path: Path) -> ast.Module:
        return ast.parse(path.read_text(encoding="utf-8"))

    def imported_names(path: Path, module_name: str) -> set[str]:
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        names: set[str] = set()
        for node in ast.walk(module_ast):
            if isinstance(node, ast.ImportFrom) and node.module == module_name:
                names.update(alias.name for alias in node.names)
        return names

    preview_source = preview_path.read_text(encoding="utf-8")
    pipeline_source = pipeline_path.read_text(encoding="utf-8")
    planning_source = (REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "screens" / "planning.py").read_text(
        encoding="utf-8"
    )
    preview_ast = module_ast_for(preview_path)
    preview_names = {node.id for node in ast.walk(preview_ast) if isinstance(node, ast.Name)}
    preview_attributes = {node.attr for node in ast.walk(preview_ast) if isinstance(node, ast.Attribute)}

    assert "build_screen_plan" in imported_names(preview_path, "pcsec_pichia.screens.planning")
    assert "build_screen_plan" in imported_names(pipeline_path, "pcsec_pichia.screens.planning")
    assert {"plan_gene_overexpression", "split_existing_genes", "split_existing_reactions"}.isdisjoint(
        imported_names(preview_path, "pcsec_pichia.screens")
    )
    assert "import re" not in preview_source
    assert "gene_index" not in preview_names | preview_attributes
    assert "gr_rules" not in preview_names | preview_attributes
    assert "x\\(" not in preview_source
    assert "_build_screen_plan" not in pipeline_source
    assert "过表达基因先进行 GPR-aware 规划" in planning_source


def test_app_gene_catalog_facade_reuses_formal_engine_catalog() -> None:
    catalog_path = REPO_ROOT / "app" / "services" / "pichia_gene_catalog_service.py"
    module_ast = ast.parse(catalog_path.read_text(encoding="utf-8"))
    imported_names: dict[str, set[str]] = {}
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_names.setdefault(node.module, set()).update(alias.name for alias in node.names)

    source = catalog_path.read_text(encoding="utf-8")
    identifiers = {node.id for node in ast.walk(module_ast) if isinstance(node, ast.Name)}
    attributes = {node.attr for node in ast.walk(module_ast) if isinstance(node, ast.Attribute)}

    assert "search_full_catalog" in imported_names["pcsec_pichia.services.gene_catalog"]
    assert "load_full_model_genes" in imported_names["pcsec_pichia.services.gene_catalog"]
    assert "get_catalog_by_category" in imported_names["pcsec_pichia.services.gene_catalog"]
    assert "search_secretion_gene_evidence" in imported_names["pcsec_pichia.services.gene_catalog"]
    assert "import re" not in source
    assert "rules" not in identifiers | attributes
    assert "gr_rules" not in identifiers | attributes
    assert "x\\(" not in source


def test_secretion_gene_evidence_map_separates_gpr_genes_from_reaction_proxies() -> None:
    from pcsec_pichia.services.gene_catalog import search_secretion_gene_evidence

    class TinyModel:
        genes = ["PAS_chr2-2_0107"]
        rxns = ["sec_PDI1_ERV2_Ero1p_complex_formation"]
        rules = [""]
        gr_rules = [""]
        gene_index = {"PAS_chr2-2_0107": 0}
        reaction_index = {"sec_PDI1_ERV2_Ero1p_complex_formation": 0}

    pdi_rows = search_secretion_gene_evidence("PDI", TinyModel())
    pdi1 = next(row for row in pdi_rows if row["common_name"] == "PDI1")

    assert pdi1["mapped_model_gene_id"] == ""
    assert pdi1["mapping_status"] == "reaction_proxy_only"
    assert pdi1["recommended_use"] == "reaction_level_proxy_requires_locus_review"
    assert pdi1["proxy_exists_in_model"] is True
    assert pdi1["proxy_has_gpr_rule"] is False
    assert pdi1["gene_level_ready"] is False

    pep4 = next(row for row in search_secretion_gene_evidence("PEP4", TinyModel()) if row["common_name"] == "PEP4")

    assert pep4["mapped_model_gene_id"] == "PAS_chr2-2_0107"
    assert pep4["mapping_status"] == "model_gpr_gene_available"
    assert pep4["recommended_use"] == "gene_level_gpr_perturbation"
    assert pep4["gene_level_ready"] is True


def test_app_full_model_gene_catalog_uses_persistent_cache(tmp_path: Path, monkeypatch) -> None:
    from pcsec_pichia.core.paths import ProjectPaths
    from pcsec_pichia.services import gene_catalog

    from app.services import pichia_gene_catalog_service as service

    calls = {"count": 0}

    def fake_load_full_model_genes() -> list[dict[str, object]]:
        calls["count"] += 1
        return [
            {
                "gene_id": "G1",
                "primary_category": "分泌相关",
                "processes": "ER",
                "n_reactions": 1,
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "oe_support_status": "oe_runnable_reaction_proxy",
            }
        ]

    monkeypatch.setattr(gene_catalog, "load_full_model_genes", fake_load_full_model_genes)
    paths = ProjectPaths(repo_root=tmp_path)

    first_rows = service.load_pichia_full_model_gene_catalog(paths=paths)
    second_rows = service.load_pichia_full_model_gene_catalog(paths=paths)

    assert first_rows == second_rows
    assert calls["count"] == 1
    assert service.pichia_full_model_gene_catalog_cache_path(paths).exists()
    assert "local_runs" in str(service.pichia_full_model_gene_catalog_cache_path(paths))
    assert "gene_catalog_cache" in str(service.pichia_full_model_gene_catalog_cache_path(paths))

    service.load_pichia_full_model_gene_catalog(force_refresh=True, paths=paths)

    assert calls["count"] == 2


def test_app_secretion_gene_evidence_uses_persistent_cache(tmp_path: Path, monkeypatch) -> None:
    from pcsec_pichia.core.paths import ProjectPaths
    from pcsec_pichia.services import gene_catalog

    from app.services import pichia_gene_catalog_service as service

    calls = {"count": 0}

    def fake_search_secretion_gene_evidence(query: str = "") -> list[dict[str, object]]:
        calls["count"] += 1
        assert query == ""
        return [
            {
                "common_name": "PDI1",
                "category": "DSB",
                "description": "Protein disulfide isomerase",
                "mapped_model_gene_id": "",
                "declared_model_gene_id": "",
                "oe_reaction_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
                "ko_reaction_id": "",
                "reaction_evidence": [
                    {
                        "reaction_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
                        "exists_in_model": True,
                        "has_gpr_rule": False,
                    }
                ],
            }
        ]

    monkeypatch.setattr(gene_catalog, "search_secretion_gene_evidence", fake_search_secretion_gene_evidence)
    paths = ProjectPaths(repo_root=tmp_path)

    static_rows = service.list_pichia_secretion_gene_evidence("PDI", paths=paths)
    assert any(row["common_name"] == "PDI1" for row in static_rows)
    assert calls["count"] == 0
    assert not service.pichia_secretion_gene_evidence_cache_path(paths).exists()

    first_rows = service.list_pichia_secretion_gene_evidence("PDI", force_refresh=True, paths=paths)
    second_rows = service.list_pichia_secretion_gene_evidence("PDI", paths=paths)

    assert first_rows == second_rows
    assert first_rows[0]["common_name"] == "PDI1"
    assert calls["count"] == 1
    assert service.pichia_secretion_gene_evidence_cache_path(paths).exists()
    assert "local_runs" in str(service.pichia_secretion_gene_evidence_cache_path(paths))
    assert "gene_catalog_cache" in str(service.pichia_secretion_gene_evidence_cache_path(paths))

    service.list_pichia_secretion_gene_evidence(force_refresh=True, paths=paths)

    assert calls["count"] == 2


def test_verified_secretion_gene_library_classifies_execution_status() -> None:
    from app.services.pichia_gene_catalog_service import list_verified_secretion_gene_library

    rows = list_verified_secretion_gene_library()
    by_name = {str(row["display_name"]): row for row in rows}

    for name in ("PDI1", "ERO1", "KAR2 / BiP", "OCH1", "PEP4", "PRB1"):
        assert name in by_name

    assert "基因级 KO" in by_name["PEP4"]["operation_status"]
    assert "基因级 KO" in by_name["PRB1"]["operation_status"]
    assert "基因级 KO" not in by_name["PDI1"]["operation_status"]
    assert "基因级 KO" not in by_name["ERO1"]["operation_status"]
    assert "反应级 OE proxy" in by_name["PDI1"]["operation_status"]
    assert by_name["PDI1"]["model_gene_id"] == ""
    assert by_name["PEP4"]["model_gene_id"]
    assert by_name["PDI1"]["detail_payload"]["curated"]["mapping_status"] == "reaction_proxy_only"


def test_full_model_gene_catalog_ui_exposes_cache_and_capability_filters() -> None:
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_catalog.py").read_text(encoding="utf-8")
    input_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_inputs.py").read_text(encoding="utf-8")

    assert "已验证分泌工程候选库" in source
    assert "list_verified_secretion_gene_library" in source
    assert "高级：全模型 GPR 基因库" in source
    assert "高级：反应级代理" in source
    assert "高级：外部证据 GPR overlay / 证据维护" in source
    assert "默认不加载" in source
    assert "暂无可执行补充规则" in source
    assert "加载全模型 GPR 基因库" in source
    assert "加载反应级代理" in source
    assert "显示外部证据 GPR overlay" in source
    assert "显示名称" in source
    assert "模型基因 ID" in source
    assert "功能 / 依据" in source
    assert "可用操作" in source
    assert "湿实验状态" in source
    assert "pichia_gene_show_full" in source
    assert 'key="pichia_gene_show_reaction_proxies"' in source
    assert 'key="pichia_gene_show_matlab_genes"' not in source
    assert 'value=True' not in source[source.find("def render_gene_lookup_panel"):source.find("def _render_full_model_gene_lookup")]
    assert "分泌工程基因名（带证据映射）" in source
    assert "反应级代理" in source
    assert "模型 GPR gene ID" in source
    assert "刷新常用基因证据缓存" in source
    assert "pichia_secretion_gene_evidence_cache_path" in (
        REPO_ROOT / "app" / "services" / "pichia_gene_catalog_service.py"
    ).read_text(encoding="utf-8")
    assert "无模型 GPR gene ID" in source
    assert "不能直接作为 gene-level KO/OE 输入" in source
    assert "可用于反应级代理；湿实验需确认 locus ID" in source
    assert "代理反应无 GPR 规则" in source
    assert "它们不是 gene-level 扰动" in source
    assert "在线刷新湿实验注释缓存" in source
    assert "只显示可敲除基因" in source
    assert "只显示可过表达代理" in source
    assert "刷新基因目录缓存" in source
    assert "基因目录缓存" in source
    assert "每页最大行数" in source
    assert "上一页" in source
    assert "下一页" in source
    assert "页码" in source
    assert "_page_input_widget_key" in source
    assert "pichia_gene_page_input_" in source
    assert 'key="pichia_gene_page_input"' not in source
    assert 'st.session_state["pichia_gene_page_input"]' not in source
    assert "最多显示行数" not in source
    assert "ko_runnable_gpr_gene_deletion" in source
    assert "oe_runnable_reaction_proxy" in source
    assert "未注释模型基因" in source
    assert "database_supported_experiment_candidate" in source
    assert "render_gene_lookup_panel()" in input_source
    assert "pichia_gene_lookup_enabled" not in input_source


def test_full_model_gene_catalog_filter_helper_supports_ko_and_oe_modes() -> None:
    from app.ui.views.simulation_gene_catalog import _filter_full_model_gene_rows

    rows = [
        {
            "gene_id": "G_KO",
            "processes": "ER folding",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_no_gpr_effect",
            "wet_lab_readiness": "database_supported_experiment_candidate",
        },
        {
            "gene_id": "G_OE",
            "processes": "translation",
            "ko_support_status": "ko_no_gpr_effect",
            "oe_support_status": "oe_runnable_reaction_proxy",
            "wet_lab_readiness": "manual_review_required",
        },
        {
            "gene_id": "G_BOTH",
            "processes": "Golgi secretion",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_runnable_reaction_proxy",
            "wet_lab_readiness": "model_only_not_experiment_ready",
        },
    ]

    assert [row["gene_id"] for row in _filter_full_model_gene_rows(rows, only_ko=True)] == ["G_KO", "G_BOTH"]
    assert [row["gene_id"] for row in _filter_full_model_gene_rows(rows, only_oe=True)] == ["G_OE", "G_BOTH"]
    assert [row["gene_id"] for row in _filter_full_model_gene_rows(rows, only_ko=True, only_oe=True)] == ["G_BOTH"]
    assert [row["gene_id"] for row in _filter_full_model_gene_rows(rows, query="golgi")] == ["G_BOTH"]
    assert [row["gene_id"] for row in _filter_full_model_gene_rows(rows, wet_lab_filter="可直接推进湿实验")] == ["G_KO"]
    assert [row["gene_id"] for row in _filter_full_model_gene_rows(rows, wet_lab_filter="需人工确认")] == ["G_OE"]
    assert [row["gene_id"] for row in _filter_full_model_gene_rows(rows, wet_lab_filter="仅模型级候选")] == ["G_BOTH"]


def test_full_model_gene_catalog_pagination_helper_clamps_pages() -> None:
    from app.ui.views.simulation_gene_catalog import _paginate_full_model_gene_rows

    rows = [{"gene_id": f"G{i}"} for i in range(1, 251)]

    page_rows, page_number, total_pages = _paginate_full_model_gene_rows(rows, page_number=1, page_size=100)
    assert page_number == 1
    assert total_pages == 3
    assert len(page_rows) == 100
    assert page_rows[0]["gene_id"] == "G1"

    page_rows, page_number, total_pages = _paginate_full_model_gene_rows(rows, page_number=3, page_size=100)
    assert page_number == 3
    assert total_pages == 3
    assert len(page_rows) == 50
    assert page_rows[0]["gene_id"] == "G201"

    page_rows, page_number, total_pages = _paginate_full_model_gene_rows(rows, page_number=99, page_size=100)
    assert page_number == 3
    assert total_pages == 3
    assert len(page_rows) == 50


def test_full_model_gene_catalog_display_helpers_prefer_names_and_reaction_evidence() -> None:
    from app.ui.views.simulation_gene_catalog import (
        _external_id_summary,
        _full_model_gene_display_name,
        _full_model_gene_function_summary,
        _gene_action_label,
        _gpr_role_label,
        _ko_status_label,
        _oe_status_label,
        _process_label,
        _wet_lab_readiness_label,
    )

    annotated = {
        "gene_id": "G1",
        "protein_name": "Protein disulfide-isomerase",
        "display_name": "Protein disulfide-isomerase",
        "function_annotation": "Catalyzes disulfide bond formation.",
        "external_ids": {"uniprot": "P12345", "kegg": "ppa:G1"},
        "ko_support_status": "ko_runnable_gpr_gene_deletion",
        "oe_support_status": "oe_runnable_reaction_proxy",
    }
    reaction_only = {
        "gene_id": "PAS_chr1-4_0141",
        "sample_reactions": ["HMPK1_no_2_fwd", "PMPK_no_1_fwd"],
        "ko_support_status": "ko_runnable_gpr_gene_deletion",
        "oe_support_status": "oe_no_gpr_effect",
    }

    assert _full_model_gene_display_name(annotated) == "Protein disulfide-isomerase"
    assert _full_model_gene_function_summary(annotated) == "Catalyzes disulfide bond formation."
    assert _full_model_gene_display_name(reaction_only) == "HMPK1/PMPK 相关酶（未注释）"
    assert "按模型 GPR 关联到反应：HMPK1, PMPK" in _full_model_gene_function_summary(reaction_only)
    assert _gene_action_label(annotated) == "可敲除 / 可过表达代理"
    assert _gene_action_label(reaction_only) == "可敲除"
    assert _ko_status_label("ko_runnable_gpr_gene_deletion") == "可运行：基因级 KO"
    assert _oe_status_label("oe_runnable_reaction_proxy") == "可运行：反应级 OE 代理"
    assert _gpr_role_label("single_gene") == "单基因"
    assert _process_label("metabolic_or_other") == "代谢 / 其他"
    assert _external_id_summary(annotated) == "uniprot: P12345; kegg: ppa:G1"
    assert _wet_lab_readiness_label("database_supported_experiment_candidate") == "可直接推进：数据库精确支持"


def test_app_gene_catalog_option_row_preserves_capability_fields() -> None:
    from app.services.pichia_gene_catalog_service import _gene_option_row

    row = _gene_option_row(
        {
            "gene_id": "G1",
            "canonical_gene_id": "G1",
            "aliases": ["ALIAS1"],
            "n_reactions": 2,
            "sample_reactions": ["R1"],
            "processes": "ER",
            "primary_category": "分泌相关",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_runnable_reaction_proxy",
            "gpr_role": "single_gene",
            "support_reason": "model evidence",
            "missing_information": ["gene_expression_to_capacity_model"],
            "confidence": "high",
            "protein_name": "Example protein",
            "function_annotation": "Example function",
            "evidence_sources": ["offline_cache"],
        }
    )

    assert row["aliases"] == ["ALIAS1"]
    assert row["ko_support_status"] == "ko_runnable_gpr_gene_deletion"
    assert row["oe_support_status"] == "oe_runnable_reaction_proxy"
    assert row["gpr_role"] == "single_gene"
    assert row["missing_information"] == ["gene_expression_to_capacity_model"]
    assert row["function_annotation"] == "Example function"


def test_full_gene_catalog_loads_offline_external_evidence_without_network(tmp_path: Path) -> None:
    from pcsec_pichia.services.gene_catalog import load_full_model_genes

    class TinyModel:
        genes = ["G1", "G2"]
        rxns = ["R1"]
        rules = ["x(1)"]
        gr_rules = ["G1"]
        gene_index = {"G1": 0, "G2": 1}
        reaction_index = {"R1": 0}

    cache_path = tmp_path / "gene_evidence.json"
    cache_path.write_text(
        json.dumps(
            {
                "genes": [
                    {
                        "gene_id": "G1",
                        "canonical_gene_id": "G1",
                        "aliases": ["ALIAS1"],
                        "protein_name": "Example protein",
                        "function_annotation": "Example function",
                        "evidence_sources": ["UniProt", "NCBI"],
                        "evidence_confidence": "reviewed",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    rows = load_full_model_genes(TinyModel(), evidence_cache_path=cache_path)
    g1 = next(row for row in rows if row["gene_id"] == "G1")
    g2 = next(row for row in rows if row["gene_id"] == "G2")

    assert g1["aliases"] == ["ALIAS1"]
    assert g1["protein_name"] == "Example protein"
    assert g1["function_annotation"] == "Example function"
    assert g1["evidence_sources"] == ["UniProt", "NCBI"]
    assert g1["ko_support_status"] == "ko_runnable_gpr_gene_deletion"
    assert g2["ko_support_status"] == "ko_no_gpr_effect"
    assert g2["oe_support_status"] == "oe_no_gpr_effect"


def test_full_gene_catalog_uses_exact_gene_matching_for_text_gr_rules() -> None:
    from pcsec_pichia.services.gene_catalog import load_full_model_genes

    class TinyModel:
        genes = ["G1", "G10"]
        rxns = ["R10"]
        rules = ["[]"]
        gr_rules = ["G10"]
        gene_index = {"G1": 0, "G10": 1}
        reaction_index = {"R10": 0}

    rows = load_full_model_genes(TinyModel())
    g1 = next(row for row in rows if row["gene_id"] == "G1")
    g10 = next(row for row in rows if row["gene_id"] == "G10")

    assert g1["n_reactions"] == 0
    assert g1["sample_reactions"] == []
    assert g1["affected_reactions"] == []
    assert g10["n_reactions"] == 1
    assert g10["sample_reactions"] == ["R10"]
    assert g10["affected_reactions"] == ["R10"]


def test_ui_display_helpers_prefer_input_gene_alias_over_canonical_id() -> None:
    frame = normalise_candidate_frame_for_display(
        pd.DataFrame(
            [
                {
                    "gene_id": "G1",
                    "canonical_gene_id": "G1",
                    "input_gene_id": "ALIAS1",
                    "candidate_id": "G1",
                    "effect_label": "提升分泌",
                    "delta_objective": 0.1,
                    "success": True,
                    "status": "0",
                }
            ]
        )
    )

    assert candidate_row_label(0, frame.iloc[0]).startswith("1. ALIAS1 |")


def test_hlf_builtin_target_semantics_use_project_defined_sequence() -> None:
    hlf = _builtin_target_semantics("hLF")

    assert hlf["alignment_target_kind"] == "project_defined_hLF"
    assert hlf["sequence_role"] == "native_signal_plus_mature_hLF"
    assert hlf["normalization_mode"] == "user_provided_as_provided"
    assert "用户提供" in hlf["target_warning"]
    assert "hLF_PROJECT_710" in hlf["target_warning"]
    assert "aligned_except_known_matlab_compatibility_differences" in hlf["target_warning"]
    assert "matlab_failed" in hlf["target_warning"]
    assert "fully aligned" in hlf["target_warning"]


def test_hlf_request_warnings_separate_project_artifact_from_historical_matlab_failure() -> None:
    warnings = request_warnings(
        SecretionRunRequest(target_source="builtin", target_id="hLF")
    )

    assert any(
        "hLF_PROJECT_710" in item
        and "aligned_except_known_matlab_compatibility_differences" in item
        and "fully aligned" in item
        for item in warnings
    )
    assert any("historical matlab_failed" in item for item in warnings)
    assert all("失败" not in item or "不代表当前项目 hLF 710aa 失败" in item for item in warnings)


def test_response_summary_exposes_target_metadata_and_warnings() -> None:
    response = SecretionRunResponse(
        success=True,
        target_id="hLF",
        result_status="corrected_condition",
        matlab_alignment_status="aligned_except_known_matlab_compatibility_differences",
        alignment_summary={
            "target_id": "hLF_PROJECT_710",
            "python_target_id": "hLF",
            "alignment_artifact_target_id": "hLF_PROJECT_710",
            "matlab_alignment_status": "aligned_except_known_matlab_compatibility_differences",
            "is_fully_aligned": False,
        },
        target_metadata={
            "alignment_target_kind": "project_defined_hLF",
            "sequence_role": "native_signal_plus_mature_hLF",
            "normalization_mode": "user_provided_as_provided",
        },
        target_warnings=["hLF 使用用户提供的 710aa 目标序列。"],
        protein_cost_analysis={
            "result_status": "draft_explanatory",
            "total_relative_score": 100.0,
            "dominant_cost_categories": ["translation"],
            "lp_attribution": {
                "result_status": "draft_lp_sensitivity",
                "top_constraint_marginals": [{"block": "protein_mass", "marginal": 1.0}],
            },
        },
        target_growth_analysis={
            "result_status": "draft_explanatory",
            "growth_sensitivity_label": "increasing",
            "growth_sensitivity_reason": "monotonic_increasing_successful_grid",
            "valid_point_count": 1,
        },
        yield_improvement_recommendations={
            "result_status": "draft_model_recommendation",
            "summary_counts": {"recommended": 1, "not_recommended": 0, "unresolved": 0},
            "recommended_candidates": [{"display_name": "PEP4"}],
        },
        medium_condition={
            "condition_id": "glucose_glycerol_ynb_core_aa_corrected",
            "carbon_source_id": "glucose_glycerol",
            "scientific_status": "draft_co_carbon_boundary_requires_promoter_context",
        },
    )

    summary = response_to_summary(response)

    assert summary["target_metadata"]["alignment_target_kind"] == "project_defined_hLF"
    assert summary["target_metadata"]["sequence_role"] == "native_signal_plus_mature_hLF"
    assert summary["target_metadata"]["normalization_mode"] == "user_provided_as_provided"
    assert summary["alignment_summary"]["python_target_id"] == "hLF"
    assert summary["alignment_summary"]["alignment_artifact_target_id"] == "hLF_PROJECT_710"
    assert summary["alignment_summary"]["matlab_alignment_status"] == "aligned_except_known_matlab_compatibility_differences"
    assert summary["alignment_summary"]["is_fully_aligned"] is False
    assert summary["target_warnings"] == ["hLF 使用用户提供的 710aa 目标序列。"]
    assert summary["protein_cost_analysis"]["result_status"] == "draft_explanatory"
    assert summary["protein_cost_analysis"]["total_relative_score"] == 100.0
    assert summary["protein_cost_analysis"]["lp_attribution"]["result_status"] == "draft_lp_sensitivity"
    assert summary["target_growth_analysis"]["result_status"] == "draft_explanatory"
    assert summary["target_growth_analysis"]["growth_sensitivity_label"] == "increasing"
    assert summary["target_growth_analysis"]["growth_sensitivity_reason"] == "monotonic_increasing_successful_grid"
    assert summary["yield_improvement_recommendations"]["result_status"] == "draft_model_recommendation"
    assert summary["yield_improvement_recommendations"]["recommended_candidates"][0]["display_name"] == "PEP4"
    assert summary["medium_condition"]["condition_id"] == "glucose_glycerol_ynb_core_aa_corrected"
    assert summary["medium_condition"]["scientific_status"] == "draft_co_carbon_boundary_requires_promoter_context"
