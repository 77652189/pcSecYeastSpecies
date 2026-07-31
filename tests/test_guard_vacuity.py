"""元守卫：遍历仓库目录的测试，必须先确认扫描目标在场。

一个测试若 `glob`/`rglob` 一个仓库目录、再断言 `violations == []`，那么该目录被改名或
移走时 glob 静默返回空，断言就 vacuously 通过——它守的边界无声失效，而且没有任何信号。

2026-07-31 实测：15 处目录遍历里有 10 处没有防护，其中包括
`test_streamlit_ui_does_not_import_engine_directly`（UI 不得直接 import 引擎）这类分层守卫。
本文件把"补完之后不许回退"焊死，新写的同类测试也必须自带防护。

只检查引用了 `REPO_ROOT` 的函数：`tmp_path` 之类由测试自己创建的目录，存在性由构造保证。
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent

_DIRECTORY_ITERATION = re.compile(r"\.(?:r?glob|iterdir)\(")
_VACUITY_GUARD = re.compile(r"\.is_dir\(\)|\.exists\(\)")


def _functions_iterating_repo_dirs() -> list[tuple[str, str, bool]]:
    """返回 (文件名, 函数名, 是否已防护)。

    跳过本文件：扫描器不扫自己。它的函数体里必然出现 `REPO_ROOT`、`glob` 这些字面量
    （那是判断规则本身），扫自己只会命中自己——首次运行就踩了这个自指陷阱。
    """
    assert TESTS_DIR.is_dir(), f"扫描目标不存在，测试会静默通过：{TESTS_DIR}"

    self_name = Path(__file__).name
    found: list[tuple[str, str, bool]] = []
    for path in sorted(TESTS_DIR.glob("test_*.py")):
        if path.name == self_name:
            continue
        source = path.read_text(encoding="utf-8")
        lines = source.splitlines()
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node, ast.FunctionDef):
                continue
            body = "\n".join(lines[node.lineno - 1:(node.end_lineno or node.lineno)])
            if not _DIRECTORY_ITERATION.search(body):
                continue
            if "REPO_ROOT" not in body:
                continue          # tmp_path 等自建目录，存在性由构造保证
            found.append((path.name, node.name, bool(_VACUITY_GUARD.search(body))))
    return found


def test_the_scan_itself_is_not_vacuous() -> None:
    """先证明这条元守卫自己没在扫空气——它正是要防的那种错误。"""
    assert TESTS_DIR.is_dir(), f"扫描目标不存在：{TESTS_DIR}"

    found = _functions_iterating_repo_dirs()

    assert len(found) >= 10, f"只找到 {len(found)} 处目录遍历，扫描口径可能已失效"


def test_every_repo_directory_scan_asserts_its_target_exists() -> None:
    """遍历仓库目录的测试必须先 assert 扫描目标在场。

    修法：把根路径提出来命名，然后
        assert some_root.is_dir(), f"扫描目标不存在，测试会静默通过：{some_root}"
    """
    unguarded = [
        f"{module}::{function}"
        for module, function, guarded in _functions_iterating_repo_dirs()
        if not guarded
    ]

    assert unguarded == [], "这些测试遍历仓库目录却没有防空转断言：\n  " + "\n  ".join(unguarded)
