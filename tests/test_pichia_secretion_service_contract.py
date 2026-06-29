from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd

from app.services.pichia_request_mapping_service import (
    request_warnings,
    sequence_contract_for_engine,
    target_input_payload,
)
from app.services.pichia_background_tasks import response_to_summary
from app.services.pichia_screen_preview_service import _preview_screen_inputs_for_model
from app.services.pichia_secretion_schema import SecretionRunRequest, SecretionRunResponse
from app.services.pichia_target_catalog_service import (
    _builtin_target_semantics,
)
from app.services.pichia_target_catalog_service import (
    known_mature_proteins,
    known_signal_peptides,
)
from app.ui.views.simulation_display import (
    candidate_effect_counts,
    normalise_candidate_frame_for_display,
    target_semantics_label,
)
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


def test_streamlit_display_helpers_localize_candidate_status_without_engine_logic() -> None:
    frame = pd.DataFrame(
        [
            {"status": "2", "success": False, "effect_label": "求解失败"},
            {"status": "optimal", "success": True, "effect_label": "提升分泌"},
        ]
    )

    display_frame = normalise_candidate_frame_for_display(frame)
    counts = candidate_effect_counts(display_frame)

    assert display_frame.loc[0, "solver_status_label"] == "约束不可行"
    assert display_frame.loc[0, "effect_label"] == "约束不可行"
    assert counts == {"提升分泌": 1, "约束不可行": 1}
    assert target_semantics_label("project_defined_hLF") == "项目定义 hLF（用户提供序列）"


def test_streamlit_gene_input_text_helpers_dedupe_multiline_candidates() -> None:
    parsed = parse_candidate_text("G1, G2\nG1\uFF1BG3\uFF0CG2")

    assert parsed == ("G1", "G2", "G3")
    assert merge_candidate_text("G1\nG2", ["G2", "G4"]) == "G1\nG2\nG4"


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
    assert preview["oe_genes"][1]["status"] == "unresolved_gene"
    assert preview["oe_reactions"][0]["status"] == "resolved"
    assert preview["oe_reactions"][1]["status"] == "unresolved_reaction"
    assert "reaction-level OE proxy" in preview["semantics"]["OE_gene_proxy"]
    assert any("敲除基因未在模型中找到" in item for item in preview["warnings"])
    assert any("过表达基因当前按 reaction-level OE proxy" in item for item in preview["warnings"])


def test_screen_preview_and_pipeline_share_engine_candidate_resolution_helpers() -> None:
    preview_path = REPO_ROOT / "app" / "services" / "pichia_screen_preview_service.py"
    pipeline_path = REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "pipeline.py"
    expected_helpers = {
        "resolve_oe_gene_reactions",
        "split_existing_genes",
        "split_existing_reactions",
    }

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
    preview_ast = module_ast_for(preview_path)
    preview_names = {node.id for node in ast.walk(preview_ast) if isinstance(node, ast.Name)}
    preview_attributes = {node.attr for node in ast.walk(preview_ast) if isinstance(node, ast.Attribute)}

    assert expected_helpers.issubset(imported_names(preview_path, "pcsec_pichia.screens"))
    assert expected_helpers.issubset(imported_names(pipeline_path, "pcsec_pichia.screens"))
    assert "import re" not in preview_source
    assert "gene_index" not in preview_names | preview_attributes
    assert "gr_rules" not in preview_names | preview_attributes
    assert "x\\(" not in preview_source
    assert "过表达基因当前按 reaction-level OE proxy" in pipeline_source


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
    assert "import re" not in source
    assert "rules" not in identifiers | attributes
    assert "gr_rules" not in identifiers | attributes
    assert "x\\(" not in source


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
