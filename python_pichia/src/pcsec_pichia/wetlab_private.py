"""在手湿实验私有数据的本地读取护栏（保密关键基础设施）。

原始湿实验数据（发酵 OD 曲线 / 菌株 / titer / 温度·pH 明细）为机密，只存**仓库外**本地
私有区、不入 git / 不上云（见 ADR-003 / ADR-006 + 保密约定）。本模块是读取这批数据的
唯一受护栏入口：

- **路径配置（不入 git）**：环境变量 `PCSEC_WETLAB_PRIVATE_DIR`；未设则回退到约定的仓库
  同级目录 `../pcSec_wetlab_private`（存在才用）；都没有 = 优雅降级（返回 None / 空）。
- **护栏**：拒绝任何落在仓库树内的私有路径（防私有数据误入 repo / 被提交）；拒绝越出私有
  区的路径穿越。危险配置宁可响亮失败（抛 `PrivateDataGuardError`），不静默读。
- **只读**：本模块不写任何东西。调用方读到私有数据后，**只准往 repo 写机制层抽象**
  （μ 量级、相对结论等），绝不写原始 OD / 菌株 / titer 明细——这条契约由调用方遵守。
"""

from __future__ import annotations

import os
from pathlib import Path

DEFAULT_PRIVATE_DIR_ENV = "PCSEC_WETLAB_PRIVATE_DIR"
CONVENTIONAL_SIBLING_DIRNAME = "pcSec_wetlab_private"


class PrivateDataGuardError(RuntimeError):
    """私有数据读取护栏被触发（如私有路径落在仓库树内、或路径越出私有区）。"""


def _repo_root() -> Path:
    # 本模块位于 <repo>/python_pichia/src/pcsec_pichia/wetlab_private.py
    return Path(__file__).resolve().parents[3]


def _is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def resolve_private_data_dir(
    env: dict[str, str] | None = None,
    *,
    repo_root: Path | None = None,
) -> Path | None:
    """解析仓库外私有数据目录；未配置 / 不存在 = None（优雅降级）。

    危险配置（路径落在仓库树内）抛 `PrivateDataGuardError`——宁可响亮失败，也不让私有数据误入 git。
    `env` / `repo_root` 仅供测试注入。
    """
    environ = os.environ if env is None else env
    repo = Path(repo_root).resolve() if repo_root is not None else _repo_root()

    raw = environ.get(DEFAULT_PRIVATE_DIR_ENV)
    if raw and raw.strip():
        path = Path(raw).expanduser().resolve()
        if _is_within(path, repo):
            raise PrivateDataGuardError(
                f"{DEFAULT_PRIVATE_DIR_ENV} 指向仓库树内（{path}）；私有数据必须在仓库外，"
                "拒绝读取以防误入 git。请把私有区移到仓库外并重设该环境变量。"
            )
        return path if path.is_dir() else None

    # 约定回退：仓库同级 ../pcSec_wetlab_private（存在且在仓库外才用）
    sibling = (repo.parent / CONVENTIONAL_SIBLING_DIRNAME).resolve()
    if sibling.is_dir() and not _is_within(sibling, repo):
        return sibling
    return None


def private_data_available(
    env: dict[str, str] | None = None,
    *,
    repo_root: Path | None = None,
) -> bool:
    return resolve_private_data_dir(env, repo_root=repo_root) is not None


def resolve_private_file(
    relative_path: str | Path,
    env: dict[str, str] | None = None,
    *,
    repo_root: Path | None = None,
) -> Path | None:
    """把私有区内的相对路径安全解析为绝对路径；未配置 / 不存在 = None；越出私有区抛护栏错。"""
    base = resolve_private_data_dir(env, repo_root=repo_root)
    if base is None:
        return None
    target = (base / relative_path).resolve()
    if not _is_within(target, base):
        raise PrivateDataGuardError(f"私有文件路径越出私有区（{target}）；拒绝。")
    return target if target.is_file() else None


def read_private_text(
    relative_path: str | Path,
    *,
    encoding: str = "utf-8",
    env: dict[str, str] | None = None,
    repo_root: Path | None = None,
) -> str | None:
    """读私有区内文本文件；未配置 / 不存在 = None（优雅降级）；越界抛护栏错。"""
    target = resolve_private_file(relative_path, env, repo_root=repo_root)
    return None if target is None else target.read_text(encoding=encoding)


def read_private_bytes(
    relative_path: str | Path,
    *,
    env: dict[str, str] | None = None,
    repo_root: Path | None = None,
) -> bytes | None:
    """读私有区内二进制文件（如发酵 Excel）；未配置 / 不存在 = None；越界抛护栏错。"""
    target = resolve_private_file(relative_path, env, repo_root=repo_root)
    return None if target is None else target.read_bytes()


__all__ = [
    "CONVENTIONAL_SIBLING_DIRNAME",
    "DEFAULT_PRIVATE_DIR_ENV",
    "PrivateDataGuardError",
    "private_data_available",
    "read_private_bytes",
    "read_private_text",
    "resolve_private_data_dir",
    "resolve_private_file",
]
