from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

ACTIVE_DOCS = {
    "EXECUTION_PLAN.md",
    "data_and_results_policy.md",
    "handoff.md",
    "pichia_current_architecture_and_requirements.md",
    "pichia_next_plan.md",
    "README.md",
}


def test_execution_plan_is_the_project_priority_control() -> None:
    index = (REPO_ROOT / "docs" / "README.md").read_text(encoding="utf-8")
    technical_plan = (REPO_ROOT / "docs" / "pichia_next_plan.md").read_text(
        encoding="utf-8"
    )
    handoff = (REPO_ROOT / "docs" / "handoff.md").read_text(encoding="utf-8")

    assert "项目级执行与预算计划" in index
    assert "技术计划不得绕过它扩大范围" in index
    assert "docs/EXECUTION_PLAN.md" in technical_plan
    assert "项目级执行与预算计划" in handoff

NON_ACTIVE_REFERENCE_DOCS = {
    "cobrapy_phase0_baseline_assessment_2026-07-06.md",
    "cobrapy_phase3_installed_shadow_validation_2026-07-06.md",
    "opn_pichia_signal_peptide_candidates.md",
    "pichia_cobrapy_import_qa_shadow_plan.md",
    "pichia_homology_crosswalk_architecture.md",
    "pichia_ko_oe_genome_screen_design_2026-07-02.md",
    "pichia_medium_mixed_carbon_objective_plan_2026-06-30.md",
    "pichia_online_external_reference_architecture.md",
    "pichia_python_hlf_design_decisions.md",
    "pichia_python_hlf_project_710_alignment_status_2026-06-26.md",
    "pichia_sce_homology_feasibility_20260708.md",
}

DELETED_OBSOLETE_MIGRATION_DOCS = {
    "pichia_python_architecture.md",
    "pichia_python_next_development_slices_2026-06-26.md",
    "pichia_python_release_validation_2026-06-25.md",
    "pichia_python_migration_strategy.md",
    "pichia_python_refactor_plan.md",
    "migration_progress.md",
}


def test_docs_root_contains_only_reviewed_active_pichia_docs() -> None:
    docs_root = REPO_ROOT / "docs"
    root_markdown_files = {
        path.name for path in docs_root.glob("*.md") if path.is_file()
    }

    assert root_markdown_files == ACTIVE_DOCS


def test_completed_reference_docs_are_not_active_entries() -> None:
    docs_root = REPO_ROOT / "docs"
    root_markdown_files = {path.name for path in docs_root.glob("*.md")}

    assert root_markdown_files.isdisjoint(NON_ACTIVE_REFERENCE_DOCS)


def test_obsolete_migration_plans_are_deleted_not_kept_as_active_debt() -> None:
    docs_root = REPO_ROOT / "docs"
    archive_root = docs_root / "archive"
    python_pichia_docs_root = REPO_ROOT / "python_pichia" / "docs"
    all_doc_names = {
        path.name
        for root in (docs_root, archive_root, python_pichia_docs_root)
        for path in root.glob("*.md")
    }

    assert all_doc_names.isdisjoint(DELETED_OBSOLETE_MIGRATION_DOCS)


def test_next_plan_routes_to_oe_product_tiering_closure() -> None:
    text = (REPO_ROOT / "docs" / "pichia_next_plan.md").read_text(encoding="utf-8")
    current_state = text.split("```yaml", 1)[1].split("```", 1)[0]
    current_next_step = text.split("## 当前下一步", 1)[1]

    assert "current_program: mvp_directions_1_to_3" in current_state
    assert "current_slice: direction_2_oe_product_tiering_closure" in current_state
    assert "slice_status: ready" in current_state
    assert "从 Phase 1 Round 0 开始" not in text
    assert "Round 6A：外部 baseline capacity 候选与审核提升" in text
    assert "Round 6B：hLF/OPN 绝对容量重验收（条件路径）" in text
    assert "direction_2_oe_product_tiering_closure" in current_next_step
    assert "docs/EXECUTION_PLAN.md" in current_next_step
    assert "不得进入方向 3" in current_next_step
    assert "执行 Round 6A 的 A0a" not in current_next_step


def test_docs_readme_routes_to_oe_product_tiering_closure() -> None:
    text = (REPO_ROOT / "docs" / "README.md").read_text(encoding="utf-8")

    assert "方向 1 的研发发酵模板回填已经验收" in text
    assert "当前授权切片按 ADR-002" in text
    assert "当前授权切片是研发组真实发酵模板回填收口" not in text


def test_handoff_points_to_oe_product_tiering_closure() -> None:
    text = (REPO_ROOT / "docs" / "handoff.md").read_text(encoding="utf-8")

    assert "current_slice: direction_2_oe_product_tiering_closure" in text
    assert "slice_status: ready" in text
    assert "absolute_capacity_status: unavailable_waiting_for_qualified_evidence" in text
    assert "只执行 `direction_2_oe_product_tiering_closure`" in text
    assert "不得新增外部容量来源" in text
    assert "不得自动进入方向 3" in text
    assert "不能只修改文档" in text


def test_active_architecture_indexes_layered_oe_decision() -> None:
    architecture = (
        REPO_ROOT / "docs" / "pichia_current_architecture_and_requirements.md"
    ).read_text(encoding="utf-8")
    adr_index = (REPO_ROOT / "docs" / "adr" / "README.md").read_text(
        encoding="utf-8"
    )
    adr = (
        REPO_ROOT
        / "docs"
        / "adr"
        / "002-relative-oe-and-absolute-capacity-layers.md"
    ).read_text(encoding="utf-8")

    assert "## 产品验收分层" in architecture
    assert "ADR-002" in architecture
    assert "ADR-002" in adr_index
    assert "相对、未校准的 OE 决策层" in adr
    assert "绝对 gene-capacity 研究层" in adr
    assert "补充 ADR-001" in adr
