from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

ACTIVE_DOCS = {
    "EXECUTION_PLAN.md",
    "handoff.md",
    "pichia_current_architecture_and_requirements.md",
    "README.md",
}


def test_execution_plan_is_the_project_priority_control() -> None:
    index = (REPO_ROOT / "docs" / "README.md").read_text(encoding="utf-8")
    handoff = (REPO_ROOT / "docs" / "handoff.md").read_text(encoding="utf-8")

    assert "项目级执行计划" in index
    assert "技术计划不得绕过它扩大范围" in index
    assert "项目级执行计划" in handoff

NON_ACTIVE_REFERENCE_DOCS = {
    "data_and_results_policy.md",
    "cobrapy_phase0_baseline_assessment_2026-07-06.md",
    "cobrapy_phase3_installed_shadow_validation_2026-07-06.md",
    "opn_pichia_signal_peptide_candidates.md",
    "pichia_cobrapy_import_qa_shadow_plan.md",
    "pichia_homology_crosswalk_architecture.md",
    "pichia_ko_oe_genome_screen_design_2026-07-02.md",
    "pichia_medium_mixed_carbon_objective_plan_2026-06-30.md",
    "pichia_next_plan.md",
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
    """只扫仓库里真实存在、且被版本控制跟踪的文档目录。

    此前还扫了 docs/archive/ 和 python_pichia/docs/，两个都是空转：
    前者被 .gitignore 排除（内容因机器而异，干净 clone 上不存在），
    后者压根没有这个目录。两个 glob 都静默返回空，于是这条断言在别人机器上
    “通过”的理由和在本机完全不同——依赖未跟踪本地状态的测试没有检测力，只有摩擦。
    """
    scanned_roots = (REPO_ROOT / "docs", REPO_ROOT / "docs" / "adr")

    # 防空转：扫描目标消失时必须红，而不是悄悄通过。上面那个 bug 就是这么藏住的。
    for root in scanned_roots:
        assert root.is_dir(), f"扫描目标不存在，测试会静默通过：{root}"

    doc_names = {path.name for root in scanned_roots for path in root.glob("*.md")}

    assert doc_names.isdisjoint(DELETED_OBSOLETE_MIGRATION_DOCS)


def test_private_archive_stays_out_of_version_control() -> None:
    """docs/archive/ 放的是不对外的设计文档（2026-07-07 保密清理 0df8f92 的产物）。

    保护它的只有 .gitignore 里的一行。那行没了，下一次 `git add -A` 就会把这批
    私有文档提交进公开仓库。把“私有文件保持 gitignore”这个决定焊死在这里。
    """
    ignored_lines = {
        line.strip()
        for line in (REPO_ROOT / ".gitignore").read_text(encoding="utf-8").splitlines()
    }

    assert "docs/archive/" in ignored_lines


def test_docs_readme_routes_to_current_slice() -> None:
    text = (REPO_ROOT / "docs" / "README.md").read_text(encoding="utf-8")

    assert "碳源条件标定 + 短名单跨条件稳健性" in text
    assert "方向 4 组合设计与目标蛋白降解通路建模明确不做" in text


def test_handoff_points_to_current_slice() -> None:
    # handoff 精简为"当前目标 + 下一步 + 必读 + 验证"；current_slice 随推进更新
    # （direction_3_erad → direction_5 碳源标定 → 改造后分层短名单 → 当前的可用性/可达性）。
    text = (REPO_ROOT / "docs" / "handoff.md").read_text(encoding="utf-8")

    # 2026-07-28：迭代2 D1–D6 已全 push，slice 推进到阶段④ 可用性 + 可达性（ADR-007）。
    # 此前"剩余全部数据门控"的判断已被证伪——分泌机器复合体可跑但界面上不可达，属非数据门控。
    assert "current_slice: usability_and_secretory_machinery_reachability" in text
    assert "slice_status: in_progress" in text
    assert "previous_slice: modified_strain_ko_oe_layered_shortlist" in text
    assert "absolute_capacity_status: unavailable_waiting_for_qualified_evidence" in text
    assert "碳源条件标定 + 跨条件稳健性" in text
    # 硬边界必须在场
    assert "glucose 的 corrected_reference 结果不得改动" in text
    assert "保密湿实验数据只存仓库外本地私有区" in text


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
