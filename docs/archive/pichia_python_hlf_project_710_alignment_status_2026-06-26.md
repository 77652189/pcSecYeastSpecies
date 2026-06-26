# hLF_PROJECT_710 Alignment Status

日期：2026-06-26
状态：active evidence note

## 结论

当前项目 hLF 710aa 序列对应的 MATLAB harness artifact target 是 `hLF_PROJECT_710`。

`hLF_PROJECT_710` 可以报告为：

```text
aligned_except_known_matlab_compatibility_differences
```

但它不是 old MATLAB `hLF` fully aligned。

## 关键证据

- MATLAB harness summary：
  - `local_runs/pichia_hlf_opn_probe/hlf_project_sequence_matlab_harness_2026-06-26/hlf_project_sequence_matlab_harness_summary.json`
- MATLAB LP：
  - `local_runs/pichia_hlf_opn_probe/hlf_project_sequence_matlab_harness_2026-06-26/Simulation_dilutionhLF_PROJECT_710_stage3_mu0p10_optional_PP.lp`
- LP diff summary：
  - `local_runs/pichia_hlf_opn_probe/hlf_project_710_lp_diff_2026-06-26/hLF_project_710_lp_diff_summary.json`

## 已知 compatibility exceptions

| exception | count | meaning |
|---|---:|---|
| corrected medium exchange bounds | 9 | Python corrected medium 与 MATLAB artifact medium bounds 不同 |
| misfolding dilution bounds | 1418 | MATLAB artifact 将 dilution_misfolding 变量封为 0；Python corrected 保持开放 |
| ribosome optional row mapping | 2 | probe-only replacement 后 row coefficient diff 可降为 0 |

## 当前 schema/reporting 要求

- Python builtin target：`hLF`。
- MATLAB artifact target：`hLF_PROJECT_710`。
- summary/report 应同时暴露两者，避免把 artifact target 当成用户输入 target。
- `is_fully_aligned` 必须为 `False`。
- old MATLAB `hLF` baseline 仍保留为 historical `matlab_failed`。

## 不再保留在 active docs 的细节

详细 row-level diff、probe-only LP 变体、旧 `hLF_CLEAN` 6 行差异等内容只作为 `local_runs/` artifact 证据保留，不再在 active docs 里展开。
