# pcSecPichia hLF 目标定义决策

日期：2026-06-26
状态：active decision note

## 当前 active hLF

Python builtin `hLF` 使用用户提供的 710 aa 项目序列：

- signal peptide：`MKLVFLVLLFLGALGLCLA`，19 aa。
- mature hLF：691 aa。
- full sequence：710 aa。
- DSB / NG / OG：21 / 4 / 0。
- sequence role：`native_signal_plus_mature_hLF`。
- normalization mode：`user_provided_as_provided`。

## 不再作为 active target 的旧入口

- `hLF_CLEAN`
- `hLF_MATURE_SECRETED`
- `hLF_NATIVE_SIGNAL`
- old `hLF` alpha-factor/pro-leader prototype

这些名字只允许出现在历史 artifact、诊断报告或 regression test 中，不应出现在 UI builtin target 列表。

## Alignment 语义

- Python target 仍叫 `hLF`。
- MATLAB harness artifact target 使用唯一 id：`hLF_PROJECT_710`。
- `hLF_PROJECT_710` 可报告为 `aligned_except_known_matlab_compatibility_differences`。
- 这不是 old MATLAB `hLF` fully aligned。
- old MATLAB `hLF` baseline 仍是 historical `matlab_failed`，原因是 harness/input mapping failure。

## 实现要求

- UI/service warning 必须同时说明：
  - 当前项目 hLF 710aa 有 `hLF_PROJECT_710` artifact。
  - 旧 MATLAB `hLF` baseline failure 是历史诊断，不等于当前项目 hLF 失败。
- report/summary 必须保留：
  - `python_target_id = hLF`
  - `alignment_artifact_target_id = hLF_PROJECT_710`
  - `is_fully_aligned = false`
- 不建议修改通用 DSB/misfolding 公式；此前差异来自 target parameter / harness artifact 语义。

## 后续可做

- 统一 UI/service 文案。
- 增加自定义 target sequence-role 校验。
- 如果未来做人源化糖基化，再单独引入 glycosylation mode，不混入当前 native hLF 定义。
