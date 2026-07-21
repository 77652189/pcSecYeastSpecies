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
        module.classify_oe_dose_response_sweep = lambda *args, **kwargs: ()
        module.compare_solver_robustness = lambda *args, **kwargs: None
        module.summarize_oe_dose_response_shape = lambda *args, **kwargs: {}
        module.summarize_protein_cost_slope_compatibility = lambda *args, **kwargs: {}
        module.summarize_solver_robustness = lambda *args, **kwargs: {}
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

    # Selection routing is built dynamically (f"pichia_draft_{action}_reactions"/"_genes"),
    # not as separate literal keys per action - check for the template pieces instead.
    assert 'f"pichia_draft_{action}_genes"' in source
    assert 'f"pichia_draft_{action}_reactions"' in source
    assert "添加到过表达输入" in source
    assert "添加到敲除输入" in source
    # A curated entry's oe_reaction_id/ko_reaction_id are surfaced as "反应" kind rows so
    # they route to the *_reactions input, distinct from gene-level *_genes input.
    assert '"类型": "基因" if kind == "gene" else "反应"' in source


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
    assert "推荐动作" in catalog_source


def test_streamlit_medium_type_labels_use_composition_names_not_internal_numbers() -> None:
    assert medium_type_label(2) == "YNB 基础培养基（维生素，无氨基酸）"
    assert medium_type_label(4) == "YNB + 核心氨基酸（15 种，默认）"
    assert medium_type_label(5) == "YNB + 全氨基酸（20 种）"
    assert "media_type=99" in medium_type_label(99)


def test_streamlit_cost_slope_option_explains_it_is_the_protein_cost_analysis_feature() -> None:
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_builder.py").read_text(encoding="utf-8")
    results_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py").read_text(encoding="utf-8")

    assert "启用蛋白成本分析（固定生长率+分泌比例网格测算成本斜率，较慢）" in source
    assert "这是目标蛋白成本分析功能本身" in source
    assert "固定生长率 μ" in source
    assert "固定一组目标蛋白分泌比例" in source
    assert "优化葡萄糖摄取反应 Ex_glc_D" in source
    assert "不勾选时不会展示任何蛋白成本分析" in source
    assert "capacity_fraction_ratios" in results_source
    assert "按当前 corrected 分泌 capacity" in results_source


def test_streamlit_solver_robustness_option_and_result_panel_are_localized() -> None:
    builder_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_builder.py").read_text(encoding="utf-8")
    results_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py").read_text(encoding="utf-8")
    runner_source = (REPO_ROOT / "app" / "services" / "pichia_secretion_runner.py").read_text(encoding="utf-8")

    # builder checkbox is present and explains what solver-robustness means
    assert "启用求解器稳健性检查（换 highs-ds/highs-ipm 重解，判断瓶颈归因是否为数值假象，较慢）" in builder_source
    assert "对偶解在退化最优解处并不唯一" in builder_source
    # result panel surfaces the OE-actionable vs floor split and the solver-robustness verdict
    assert "OE 可缓解瓶颈（binding 上限，按复合体）" in results_source
    # the floor block was reframed as the 'why is it limited' answer (largest shadow prices, OE cannot relax)
    assert "为什么受限：最强约束层（下界/最低要求，OE 动不了）" in results_source
    assert "求解器稳健性（瓶颈归因是否跨求解器稳定）" in results_source
    assert "ranking-sensitive-to-solver" in results_source
    # service facade threads the flag through to the engine request
    assert "enable_solver_robustness_check=bool(request.enable_solver_robustness_check)" in runner_source


def test_streamlit_oe_dose_response_option_and_result_panel_are_localized() -> None:
    builder_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_builder.py").read_text(encoding="utf-8")
    results_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py").read_text(encoding="utf-8")
    runner_source = (REPO_ROOT / "app" / "services" / "pichia_secretion_runner.py").read_text(encoding="utf-8")

    # builder checkbox is present and explains that it replaces the single fixed 2x point
    assert "启用 OE 剂量响应形状（扫描多个过表达倍数，看分泌提升会很快到顶还是持续上升，较慢）" in builder_source
    assert "默认的过表达筛查只在固定 2× 一个点上测提升" in builder_source
    # result panel surfaces the shape section in plain Chinese for the research staff
    assert "OE 剂量响应形状（过表达越多，分泌是持续上升还是很快到顶）" in results_source
    # relative-signal framing must be explicit; no absolute capacity/optimal-dose claims
    assert "不产出绝对产量或最优倍数" in results_source
    # the shape-taxonomy labels now live in the single central i18n dict (not inline)
    from app.core.i18n import sim_result_value_label

    assert sim_result_value_label("saturating") == "饱和型（适度过表达就够，再加收益递减）"
    assert sim_result_value_label("flat_no_response") == "无响应（任何倍数都几乎没提升，别过表达）"
    # service facade threads the flag through to the engine request
    assert "enable_oe_dose_response=bool(request.enable_oe_dose_response)" in runner_source


def test_simulation_result_localization_goes_through_one_central_dictionary() -> None:
    # 用户要求：结果页的英文字段名/枚举值统一走一个集中字典翻译（app.core.i18n），避免散落、不一致、漏改。
    from app.core.i18n import sim_result_column_label, sim_result_value_label, sim_result_warning_label
    from app.ui.views.simulation_results import _localized_frame

    # 列名 -> 中文；未知列回退原文
    assert sim_result_column_label("secretory_process") == "分泌资源层"
    assert sim_result_column_label("abs_marginal") == "影子价格绝对值"
    assert sim_result_column_label("totally_unknown_field") == "totally_unknown_field"
    # 枚举/编码值 -> 中文（边界类型 / 分类 / 形状 / 布尔）；未知回退原文；None -> —
    assert sim_result_value_label("lower").startswith("下限")
    assert "跨求解器翻转" in sim_result_value_label("ranking-sensitive-to-solver")
    assert sim_result_value_label("saturating") == "饱和型（适度过表达就够，再加收益递减）"
    assert sim_result_value_label(True) == "是" and sim_result_value_label(False) == "否"
    assert sim_result_value_label("some_new_unmapped_code") == "some_new_unmapped_code"
    assert sim_result_value_label(None) == "—"
    # _localized_frame 是唯一漏斗：重命名列 + 映射 value_columns 里的枚举值，payload 英文不外泄
    frame = _localized_frame(
        [{"reaction_id": "sec_X", "bound_type": "lower", "abs_marginal": 5.0}],
        value_columns=("bound_type",),
    )
    assert list(frame.columns) == ["反应", "边界类型", "影子价格绝对值"]
    assert frame["边界类型"].iloc[0].startswith("下限")
    # 产量提升推荐表的单元格枚举也走同一字典
    assert sim_result_value_label("model_executable") == "模型可执行"
    assert sim_result_value_label("OE_reaction") == "反应级过表达（OE）"
    # 引擎英文警告按标志性子串翻译；未命中回退原文
    lp_warning = (
        "LP sensitivity is a Python draft based on SciPy HiGHS marginals; "
        "it is not MATLAB/SoPlex fully aligned shadow pricing."
    )
    assert sim_result_warning_label(lp_warning).startswith("LP 灵敏度")
    assert sim_result_warning_label("a brand new unmapped warning") == "a brand new unmapped warning"


def test_streamlit_relative_signal_charts_render_biologist_facing_figures() -> None:
    # End users are biologists, not engineers: R1 bottlenecks and R2 dose-response are charted
    # (Plotly), not just tabulated. These pure frame helpers back those charts.
    from app.ui.views.simulation_results import (
        _lp_floor_bottleneck_frame,
        _lp_oe_bottleneck_frame,
        _oe_dose_response_curve_frame,
        _resource_layer_label,
        _short_reaction_label,
    )

    # biologist-facing labels: strip boilerplate suffix, middle-truncate huge ids, and localize the
    # engine-provided secretory_process code via the central i18n dict (the engine now classifies
    # sec_* complexes itself, so this is a straight lookup with no name-based inference).
    assert _short_reaction_label("sec_Pdi1p_complex_formation") == "sec_Pdi1p"
    long_id = "PAS_chr2-2_0475_COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex"
    assert len(_short_reaction_label(long_id)) <= 34 and "…" in _short_reaction_label(long_id)
    assert _resource_layer_label("disulfide_folding") == "二硫键折叠 / DSB"  # PDI floor, classified by the engine
    assert _resource_layer_label("ribosome") == "翻译（核糖体）"
    assert _resource_layer_label("unknown") == "未解析"  # unmapped code degrades gracefully, no guessing

    # R2 dose-response -> factor on x, relative gain (%) on y, baseline factor 1.0 anchors at 0%
    oe_payload = {
        "reaction_shapes": [
            {
                "reaction_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
                "shape": "saturating",
                "point_deltas": [
                    {"factor": 1.0, "relative_gain": 0.0},
                    {"factor": 2.0, "relative_gain": 0.08},
                    {"factor": 4.0, "relative_gain": 0.11},
                ],
            }
        ]
    }
    curve = _oe_dose_response_curve_frame(oe_payload)
    assert list(curve.columns) == ["过表达倍数", "分泌相对提升(%)", "反应｜形状"]
    assert len(curve) == 3
    assert curve.loc[curve["过表达倍数"] == 2.0, "分泌相对提升(%)"].iloc[0] == 8.0
    assert curve.loc[curve["过表达倍数"] == 1.0, "分泌相对提升(%)"].iloc[0] == 0.0
    assert "饱和型" in curve["反应｜形状"].iloc[0]

    # R1 OE-actionable bottlenecks -> horizontal bar with the resource layer localized to Chinese.
    # secretory_process codes are the engine's own (sec_Pdi1p -> disulfide_folding), not inferred.
    lp_payload = {
        "oe_actionable_bottlenecks": [
            {"reaction_id": "sec_Pdi1p_complex_formation", "abs_marginal": 0.92, "secretory_process": "disulfide_folding"},
            {"reaction_id": "Mach_Ribosome_complex_formation", "abs_marginal": 0.5, "secretory_process": "ribosome"},
        ]
    }
    bars = _lp_oe_bottleneck_frame(lp_payload)
    assert set(bars["分泌资源层"]) == {"二硫键折叠 / DSB", "翻译（核糖体）"}
    assert bars["影子价格(绝对值)"].max() == 0.92

    # the large lower-bound floors are the 'why is it limited' answer and are charted separately.
    # hLF's dominant floor is the PDI disulfide-folding complex: the engine tags it disulfide_folding
    # in the payload (it used to fall through to 'unknown'), so it charts as the folding layer here.
    floor_payload = {
        "floor_constraints_not_oe_addressable": [
            {"reaction_id": "sec_Pdi1p_complex_formation", "abs_marginal": 5073.9, "secretory_process": "disulfide_folding"},
            {"reaction_id": "Mach_Ribosome_complex_formation", "abs_marginal": 180.8, "secretory_process": "ribosome"},
        ]
    }
    floors = _lp_floor_bottleneck_frame(floor_payload)
    assert floors["影子价格(绝对值)"].max() == 5073.9
    assert set(floors["分泌资源层"]) == {"二硫键折叠 / DSB", "翻译（核糖体）"}


def test_streamlit_value_of_information_panel_is_localized_and_chart_backed() -> None:
    # R4 (ADR-004): the value-of-information panel is a ranking-confidence + what-to-measure product,
    # framed as relative (never absolute yield), and reads the pipeline's value_of_information payload.
    results_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_results.py").read_text(encoding="utf-8")
    assert "排序可信度 & 该测什么（价值-of-information）" in results_source
    assert "只排测量优先级，不预测结果、不自动认定谁更好" in results_source
    assert "候选排序：分数越接近越难区分" in results_source

    from app.ui.views.simulation_results import _value_of_information_payload

    assert _value_of_information_payload({}) == {}
    payload = _value_of_information_payload(
        {"value_of_information": {"has_actionable_ambiguity": False, "ranked_candidates": []}}
    )
    assert payload["has_actionable_ambiguity"] is False


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
    assert preview["ko_genes"][0]["external_model_sources"] == []
    assert preview["ko_genes"][0]["gpr_source_priority"] == {}
    assert preview["ko_genes"][0]["external_gpr_candidate_count"] == 0
    assert preview["ko_genes"][0]["external_gpr_mapping_status"] == {}
    assert preview["ko_genes"][0]["external_gpr_conflict_warnings"] == []
    assert preview["ko_genes"][0]["manual_review_reasons"] == []
    assert preview["ko_reactions"][0]["external_model_sources"] == []
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


def test_screen_input_preview_passes_homology_evidence_without_changing_tier(monkeypatch) -> None:
    from pcsec_pichia.services import homology_evidence
    from pcsec_pichia.services.homology_evidence import GeneHomologyEvidence

    class TinyModel:
        rxns = ["R1"]
        rules = ["x(1)"]
        gr_rules = ["G1"]
        genes = ["G1"]
        gene_index = {"G1": 0}
        reaction_index = {"R1": 0}

    monkeypatch.setattr(
        homology_evidence,
        "load_homology_evidence_cache",
        lambda cache_dir=None: {
            "g1": GeneHomologyEvidence(
                gene_id="G1",
                internal_common_name="KAR2 / BiP",
                query_symbol="KAR2",
                pichia_gene_id="G1",
                pichia_model_gene_id="G1",
                homology_review_status="model_ready_rbh_high_confidence",
                rule_transfer_status="rule_transfer_ready",
                name_consistency_status="name_confirmed_by_rbh",
                is_rbh=True,
                in_model_gene_index=True,
            )
        },
    )
    request = SecretionRunRequest(
        target_source="builtin",
        target_id="OPN_ALPHA_FULL_PROJECT",
        ko_gene_ids=("G1",),
        oe_gene_ids=("G1",),
    )

    preview = _preview_screen_inputs_for_model(TinyModel(), request)

    assert preview["gene_capabilities"][0]["homology_review_status"] == "model_ready_rbh_high_confidence"
    assert preview["ko_genes"][0]["rule_transfer_status"] == "rule_transfer_ready"
    assert preview["oe_genes"][0]["homology_evidence"]["query_symbol"] == "KAR2"
    assert preview["oe_genes"][0]["recommendation_tier"] == "model_executable"
    assert preview["oe_genes"][0]["recommendation_tier"] != "experiment_calibrated"


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


def test_gene_lookup_panel_is_a_single_search_not_three_separate_browsers() -> None:
    """候选库合并成一个搜索入口（策展库+全模型基因），只保留外部证据overlay作为独立高级区。

    genome-wide screen 已经系统覆盖了全部模型基因和策展反应级候选，候选库不再需要自己的
    发现型大浏览器（分页/多重筛选器/单独的"加载"开关）——那些已经被筛查结果+核实跳转取代。
    """
    source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_catalog.py").read_text(encoding="utf-8")
    input_source = (REPO_ROOT / "app" / "ui" / "views" / "simulation_gene_inputs.py").read_text(encoding="utf-8")

    # One merged search, not three separately-toggled browser panels.
    assert "list_verified_secretion_gene_library" in source
    assert "load_pichia_full_model_gene_catalog" in source
    assert "全基因组KO/OE筛查已经系统覆盖了全部1025个模型基因和策展库的反应级候选" in source
    assert "在仿真验证中核实" in source
    assert '"高级：全模型 GPR 基因库' not in source
    assert '"高级：反应级代理' not in source
    assert "加载全模型 GPR 基因库" not in source
    assert "加载反应级代理" not in source
    assert "pichia_gene_show_full" not in source
    assert "pichia_gene_show_reaction_proxies" not in source

    # No more pagination/filter-toolbar chrome for browsing all genes.
    assert "_paginate_full_model_gene_rows" not in source
    assert "_page_input_widget_key" not in source
    assert "只显示可敲除基因" not in source
    assert "只显示可过表达代理" not in source
    assert "每页最大行数" not in source
    assert "上一页" not in source
    assert "下一页" not in source

    # External evidence GPR overlay is kept as its own advanced section - it is the one
    # thing genome-wide screening does not cover at all.
    assert "高级：外部证据 GPR overlay / 候选库维护" in source
    assert "暂无可执行补充规则" in source
    assert "显示外部证据 GPR overlay" in source
    assert "刷新常用基因证据缓存" in source
    assert "在线重建全模型湿实验注释缓存" in source

    # Combination testing (multi-select -> add to KO/OE) is preserved.
    assert "选择候选（可多选，用于组合测试）" in source
    assert 'key="pichia_gene_search_add_ko"' in source
    assert 'key="pichia_gene_search_add_oe"' in source

    assert "render_gene_lookup_panel()" in input_source


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


def test_gene_search_selection_routes_genes_and_reactions_to_separate_inputs() -> None:
    from app.ui.views.simulation_gene_catalog import _partition_selection_by_kind

    rows = [
        {"来源": "策展库", "ID": "PAS_chr2-2_0107", "kind": "gene"},
        {"来源": "策展库", "ID": "sec_PDI1_ERV2_Ero1p_complex_formation", "kind": "reaction"},
        {"来源": "全模型", "ID": "PAS_chr1-1_0013", "kind": "gene"},
    ]

    genes, reactions = _partition_selection_by_kind(rows)

    assert genes == ["PAS_chr2-2_0107", "PAS_chr1-1_0013"]
    assert reactions == ["sec_PDI1_ERV2_Ero1p_complex_formation"]


def test_gene_search_merges_curated_and_full_model_sources_without_duplicates(monkeypatch) -> None:
    import app.ui.views.simulation_gene_catalog as catalog_ui

    monkeypatch.setattr(
        catalog_ui,
        "list_verified_secretion_gene_library",
        lambda query: [
            {
                "display_name": "PEP4",
                "model_gene_id": "PAS_chr2-2_0107",
                "oe_reaction_id": "",
                "ko_reaction_id": "",
                "function_annotation": "液泡蛋白酶 A",
            }
        ],
    )
    monkeypatch.setattr(
        catalog_ui,
        "load_pichia_full_model_gene_catalog",
        lambda: [
            # Same gene_id the curated row already returned, plus one full-model-only
            # match the curated library doesn't know about - both should still be
            # matched by _filter_full_model_gene_rows's own query check.
            {
                "gene_id": "PAS_chr2-2_0107",
                "display_name": "PEP4",
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "oe_support_status": "oe_no_gpr_effect",
            },
            {
                "gene_id": "PAS_chr1-1_0099",
                "display_name": "PEP4-like paralog",
                "ko_support_status": "ko_no_gpr_effect",
                "oe_support_status": "oe_runnable_reaction_proxy",
            },
        ],
    )

    rows, truncated_count = catalog_ui._collect_search_rows("pep4")

    assert truncated_count == 0
    ids = [(row["kind"], row["ID"]) for row in rows]
    # PAS_chr2-2_0107 appears in both sources but must be deduplicated to one row.
    assert ids.count(("gene", "PAS_chr2-2_0107")) == 1
    assert ("gene", "PAS_chr1-1_0099") in ids


def test_full_model_gene_catalog_display_helpers_prefer_names_and_reaction_evidence() -> None:
    from app.ui.views.simulation_gene_catalog import (
        _full_model_gene_display_name,
        _full_model_gene_function_summary,
    )

    annotated = {
        "gene_id": "G1",
        "protein_name": "Protein disulfide-isomerase",
        "display_name": "Protein disulfide-isomerase",
        "function_annotation": "Catalyzes disulfide bond formation.",
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
            "result_status": "draft_cost_slope_analysis",
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
        value_of_information={
            "result_status": "draft_value_of_information",
            "has_actionable_ambiguity": False,
            "ranked_candidates": [{"rank": 1, "candidate_id": "PEP4", "score": 0.5}],
            "information_items": [],
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
    assert summary["protein_cost_analysis"]["result_status"] == "draft_cost_slope_analysis"
    assert summary["protein_cost_analysis"]["lp_attribution"]["result_status"] == "draft_lp_sensitivity"
    assert summary["target_growth_analysis"]["result_status"] == "draft_explanatory"
    assert summary["target_growth_analysis"]["growth_sensitivity_label"] == "increasing"
    assert summary["target_growth_analysis"]["growth_sensitivity_reason"] == "monotonic_increasing_successful_grid"
    assert summary["yield_improvement_recommendations"]["result_status"] == "draft_model_recommendation"
    # R4 (value-of-information) must flow through the cached response dict the results page reads,
    # not only via the summary-file fallback, so the panel renders consistently with the others.
    assert summary["value_of_information"]["result_status"] == "draft_value_of_information"
    assert summary["value_of_information"]["ranked_candidates"][0]["candidate_id"] == "PEP4"
    assert summary["yield_improvement_recommendations"]["recommended_candidates"][0]["display_name"] == "PEP4"
    assert summary["medium_condition"]["condition_id"] == "glucose_glycerol_ynb_core_aa_corrected"
    assert summary["medium_condition"]["scientific_status"] == "draft_co_carbon_boundary_requires_promoter_context"
