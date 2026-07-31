from __future__ import annotations

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

# handoff 的 yaml 状态块允许的取值。断言「取值在枚举内」而不是「等于当前值」：
# slice 从 in_progress 推进到 done 不该让测试变红，只有冒出一个没见过的状态词才该红。
SLICE_STATUS_VALUES = {"not_started", "in_progress", "blocked", "done"}
REQUIRED_HANDOFF_STATE_KEYS = (
    "current_slice",
    "slice_status",
    "previous_slice",
    "previous_slice_status",
    "absolute_capacity_status",
)


def _handoff_state_block() -> dict[str, str]:
    """取 handoff 里那个 ```yaml 状态块，解析成 key -> value。"""
    text = (REPO_ROOT / "docs" / "handoff.md").read_text(encoding="utf-8")
    match = re.search(r"```yaml\n(.*?)```", text, re.DOTALL)
    assert match is not None, "handoff 里找不到 yaml 状态块，后续断言会全部落空"

    fields = {}
    for line in match.group(1).splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
    assert fields, "yaml 状态块解析出 0 个字段，解析口径可能已失效"
    return fields

ACTIVE_DOCS = {
    "EXECUTION_PLAN.md",
    "handoff.md",
    "requirements.md",
    "architecture.md",
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
    assert docs_root.is_dir(), f"扫描目标不存在，测试会静默通过：{docs_root}"

    root_markdown_files = {
        path.name for path in docs_root.glob("*.md") if path.is_file()
    }

    assert root_markdown_files == ACTIVE_DOCS


def test_completed_reference_docs_are_not_active_entries() -> None:
    docs_root = REPO_ROOT / "docs"
    assert docs_root.is_dir(), f"扫描目标不存在，测试会静默通过：{docs_root}"

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


def test_docs_readme_routes_instead_of_restating_state() -> None:
    """索引只负责路由。

    当前 slice、能力清单、门控项分别属于 handoff / 架构 / 执行计划；索引各抄一份
    就是双权威，而那份副本没人会记得更新——这不是假设，上一轮审计当场抓到过。
    """
    text = (REPO_ROOT / "docs" / "README.md").read_text(encoding="utf-8")

    for doc in ACTIVE_DOCS - {"README.md"} | {"adr/README.md"}:
        assert f"({doc})" in text, f"文档索引里没有指向 {doc} 的链接"

    # 永久边界可以留（不衰减）；当前状态不行。
    assert "方向 4" in text and "明确不做" in text
    assert "当前 slice" not in text, "索引复制了 handoff 的当前 slice —— 双权威，改成指针"


def test_handoff_state_block_is_wellformed() -> None:
    """断言状态块的**形状**合法，不断言它当前是什么。

    早先这里写的是 `assert "slice_status: in_progress" in text` 这类当前值断言：
    slice 一推进就红，而修法永远是改测试去迎合文档——那是家务，检测力为零。
    现在只有「字段缺失」或「冒出没见过的状态词」才会红。
    """
    fields = _handoff_state_block()

    missing = [key for key in REQUIRED_HANDOFF_STATE_KEYS if key not in fields]
    assert missing == [], f"handoff 状态块缺字段：{missing}"

    for key in ("slice_status", "previous_slice_status"):
        assert fields[key] in SLICE_STATUS_VALUES, (
            f"{key}={fields[key]!r} 不在合法枚举内 {sorted(SLICE_STATUS_VALUES)}；"
            "若确实新增了状态词，改枚举而不是删断言"
        )

    for key in ("current_slice", "previous_slice"):
        assert re.fullmatch(r"[a-z0-9_]+", fields[key]), f"{key} 不是 snake_case 标识：{fields[key]!r}"

    # 绝对容量恒 unavailable 是 ADR-002 定下的**永久边界**，不是当前值——
    # 具体后缀（等哪种证据）可以变，"unavailable" 这个前缀不能变。
    assert fields["absolute_capacity_status"].startswith("unavailable"), (
        f"绝对容量声称成了 {fields['absolute_capacity_status']!r}，"
        "这违反 ADR-002 的永久边界；要改先改 ADR"
    )


def test_handoff_states_the_hard_boundaries_verbatim() -> None:
    """硬约束必须以**原文**在场——这些是不变量，不是状态。"""
    text = (REPO_ROOT / "docs" / "handoff.md").read_text(encoding="utf-8")

    assert "glucose 的 corrected_reference 结果不得改动" in text
    assert "保密湿实验数据只存仓库外本地私有区" in text
    assert "方向 4 组合搜索、目标蛋白降解通路建模、换默认 solver：明确不做" in text


def test_active_architecture_indexes_layered_oe_decision() -> None:
    architecture = (
        REPO_ROOT / "docs" / "requirements.md"
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


# 「随决策变」的小节：回答"我们要什么、不能碰什么"。它们消失通常意味着
# 有人把需求或边界当成过期状态一起删了。
# 2026-07-31 需求与架构按衰减率拆开后，这批小节分属两份文档。
SLOW_LAYER_SECTIONS = {
    "requirements.md": (
        "## 原始研发目标",
        "## 核心证据边界",
        "## 产品验收分层",
        "## 项目成功标准",
    ),
    "architecture.md": (
        "## 架构边界",
        "## 数据与产物治理",
    ),
}


def _read_doc(name: str) -> str:
    return (REPO_ROOT / "docs" / name).read_text(encoding="utf-8")


def test_architecture_doc_keeps_its_slow_layer_sections() -> None:
    """慢层小节必须在场。

    快层（当前状态 / 进度 / 待办）**故意不做**标题断言：那类内容本来就该随进度变，
    锁住它只会退化成家务——每次推进都红，而修法永远是改测试去迎合文档，检测力为零。
    """
    # 必须按**整行**比对，不能用子串：`"## 核心证据边界" in text` 会被
    # `## 核心证据边界_DELETED` 满足，于是"改名"这种最常见的失效方式漏网。
    # （这条是变异检验当场抓出来的——第一版就是子串匹配。）
    missing: list[str] = []
    for doc_name, sections in SLOW_LAYER_SECTIONS.items():
        headings = {line.strip() for line in _read_doc(doc_name).splitlines()}
        missing += [f"{doc_name}::{s}" for s in sections if s not in headings]

    assert missing == [], f"需求/边界类小节不见了或被改名：{missing}"


def test_current_state_has_exactly_one_authority() -> None:
    """状态只能有一个权威，否则必然分叉。

    2026-07-31 之前：架构文档写"此处为权威状态"、执行计划写"只列当前状态"——两处双权威。
    同一件事各记一份的后果已经在 handoff 里现形过（被证伪的结论没删干净，与更正并存）。
    现在状态权威归执行计划，架构文档显式交出。
    """
    for doc_name in ("requirements.md", "architecture.md"):
        text = _read_doc(doc_name)
        assert "不持有当前状态的权威" in text, f"{doc_name} 没有显式交出状态权威"
        assert "此处为权威状态" not in text, f"{doc_name} 双权威回流了"

    assert "状态的唯一权威" in _read_doc("EXECUTION_PLAN.md")
