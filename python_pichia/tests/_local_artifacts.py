"""依赖 gitignored 本地探针产物的测试共用的守卫。

这些产物在 `local_runs/` 下，按架构文档的数据治理规定**不进版本控制**，
因此干净 clone 上不存在。依赖它们的测试必须**显式跳过**而不是让断言崩掉：

- 崩掉 → 干净 clone 上一串莫名其妙的失败，看不出是环境问题还是真回归；
- 静默通过 → 更糟，「通过」的理由在两台机器上完全不同（见 canon 的空转守卫反模式）；
- 显式跳过 → 「这台机器少测了什么」在测试输出里可见。

不把产物提交成 fixture：那会破坏「运行产物一律进 local_runs/」这条已记录的边界。
"""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PROBE_DIR = REPO_ROOT / "local_runs" / "pichia_hlf_opn_probe"

MATLAB_STAGE3_ARTIFACT = PROBE_DIR / "matlab_stage3_alignment" / "matlab_stage3_alignment_summary.json"
HLF_PROJECT_710_ARTIFACT = (
    PROBE_DIR
    / "hlf_project_sequence_matlab_harness_2026-06-26"
    / "hlf_project_sequence_matlab_harness_summary.json"
)
TARGETS_EXAMPLE = PROBE_DIR / "targets.example.json"


def require_local_probe_artifact(path: Path) -> None:
    """产物缺失时跳过并说明原因；存在但为空则视为真故障。"""
    if not path.exists():
        pytest.skip(f"本地探针产物不在（gitignored，干净 clone 上没有）：{path}")
    assert path.stat().st_size > 0, f"本地探针产物是空文件：{path}"
