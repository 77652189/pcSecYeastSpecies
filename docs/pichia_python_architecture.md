# pcSecPichia Python 架构边界

日期：2026-06-26
状态：active architecture note

## 目标

只构建毕赤酵母 `pcSecPichia` Python engine，不全量迁移三物种 MATLAB 项目。

核心能力：

- OPN / hLF / custom target 分泌仿真。
- optional constraints。
- 小批量 KO/OE smoke。
- growth tradeoff smoke。
- summary/report/candidates/tradeoff 输出。
- MATLAB alignment schema/reporting。

## 分层

```text
Streamlit UI
  -> app service facade
    -> python_pichia pipeline
      -> loading / targets / secretion_plan / constraints / simulation / screens / reports / alignment

Future FastAPI
  -> app service facade
    -> python_pichia pipeline
```

## 边界规则

- `python_pichia/`：正式 engine，只放模型加载、目标构建、约束、求解、screen、report、alignment、pipeline。
- `app/services/`：service facade，负责 request/response、output_dir、warnings/errors、后台任务；不写核心模型算法。
- `app/ui/`：Streamlit 表单和展示；不直接调用 engine 内部 helper。
- `app/api/`：当前只允许作为 experimental facade，不直接调用 engine 内部模块。
- `Code/`、`Model/`、`Enzymedata/`、`Results/`：reference-only。
- `local_runs/`：原型和证据产物，不作为正式包结构。

## 当前已知架构债

| 债务 | 处理方向 |
|---|---|
| `app/ui/views/simulation.py` 过胖 | 拆 builder / results / candidate graph / MATLAB reference |
| `sys.path` bootstrap 分散 | 保留 `app/__init__.py` 为短期方案；长期用 editable install |
| FastAPI 定位不清 | 明确 experimental 或删除 |

已完成的架构收口：

- `app/services/pichia_secretion_service.py` 已拆成 schema / runner / catalog / background 等 owner modules。
- gene/reaction 解析已收敛到 engine helpers，preview 与 pipeline 复用同一解析边界。
- 旧迁移大文档和旧 route-by-route OPN/hLF draft 测试不再作为 active 维护对象。

## Alignment 语义

- `corrected_condition` 不等于旧 MATLAB fully aligned。
- OPN 与 hLF 710aa 当前可报告为 `aligned_except_known_matlab_compatibility_differences`，但 `is_fully_aligned` 必须为 `False`。
- Python builtin `hLF` 使用用户提供的 710 aa 项目序列。
- 对应 MATLAB harness artifact target 为 `hLF_PROJECT_710`。
- 旧 MATLAB `hLF` baseline 仍是 historical `matlab_failed`。

## 下一步入口

后续执行以此文档为入口：

- `docs/pichia_python_next_development_slices_2026-06-26.md`
