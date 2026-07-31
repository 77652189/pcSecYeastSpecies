from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _read_repo_text(relative_path: str) -> str:
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


def test_runtime_outputs_and_solver_artifacts_are_ignored() -> None:
    gitignore = _read_repo_text(".gitignore")

    for pattern in (
        "local_runs/",
        "python_pichia/local_runs/",
        "*.lp",
        "*.lp.out",
        "*.soplex.out",
        "*.float.out",
    ):
        assert pattern in gitignore

    for protected_dir in ("Data/", "Results/", "Model/", "Enzymedata/"):
        assert protected_dir not in {
            line.strip()
            for line in gitignore.splitlines()
            if line.strip() and not line.strip().startswith("#")
        }


def test_data_results_policy_documents_current_ownership() -> None:
    policy = _read_repo_text("docs/architecture.md")

    assert "`Results/` 是 legacy MATLAB results" in policy
    assert "不是当前 Python 或 Streamlit 的默认输出目录" in policy
    assert "`local_runs/` 是当前 Python、Streamlit、MATLAB harness" in policy
    assert "新生成的 LP、solver output" in policy
    assert "Git LFS" in policy


def test_current_docs_keep_runtime_outputs_out_of_science_assets() -> None:
    architecture = _read_repo_text("docs/architecture.md")

    assert "`Results/` 保留为 legacy MATLAB results" in architecture
    assert "统一落地目录，默认 ignored" in architecture
    assert "历史 `Results/` 迁移、Git LFS 改造或仓库历史瘦身" in architecture


def test_runtime_write_lines_do_not_target_protected_science_asset_dirs() -> None:
    protected_tokens = tuple(
        token
        for directory in ("Data", "Results", "Model", "Enzymedata")
        for token in (
            f'"{directory}"',
            f"'{directory}'",
            f' / "{directory}"',
            f"/{directory}/",
            f"\\{directory}\\",
            f'Path("{directory}")',
            f"Path('{directory}')",
            f" / {directory!r}",
        )
    )
    write_tokens = (
        ".write_text(",
        ".open(\"w\"",
        ".open('w'",
        "to_csv(",
        "to_excel(",
        "savemat(",
        "Export-Csv",
    )
    scanned_roots = (REPO_ROOT / "app", REPO_ROOT / "python_pichia", REPO_ROOT / "scripts")

    # 防空转：任一根目录被改名/移走时 rglob 静默返回空，"不许写保护目录"这条安全约束
    # 就会在无人察觉的情况下失效。边界测试最不能静默变绿。
    for root in scanned_roots:
        assert root.is_dir(), f"扫描目标不存在，测试会静默通过：{root}"

    suspicious: list[str] = []

    for root in scanned_roots:
        for path in root.rglob("*.py"):
            relative = path.relative_to(REPO_ROOT).as_posix()
            if "/tests/" in relative or relative.startswith("scripts/archive/"):
                continue
            text = path.read_text(encoding="utf-8")
            for line_number, line in enumerate(text.splitlines(), start=1):
                if any(token in line for token in write_tokens) and any(
                    token in line for token in protected_tokens
                ):
                    suspicious.append(f"{relative}:{line_number}:{line.strip()}")

    assert suspicious == []
