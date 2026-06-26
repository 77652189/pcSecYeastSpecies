# pcSecPichia Python 下一阶段开发切片

日期：2026-06-26
状态：active working plan

本文是接下来开发的执行入口。目标是把当前剩余问题拆成每轮对话可以完成、可以验证、不会顺手扩大范围的小任务。

## Active 文档

当前只保留 3 个 active 文档：

- `docs/pichia_python_next_development_slices_2026-06-26.md`：开发入口、当前状态和小目标切片。
- `docs/pichia_python_architecture.md`：架构边界和模块职责。
- `docs/pichia_python_release_validation_2026-06-25.md`：日常、聚焦、slow、release 验证命令。

历史材料放入 `docs/archive/`，只作为证据库，不作为继续开发计划。

## 当前状态快照

| 能力 | 状态 |
|---|---|
| OPN / hLF / custom target 输入 | 可用 |
| Python corrected 分泌仿真 | 可用，仍是 draft/corrected condition |
| optional constraints | 可用 |
| KO/OE smoke | 可用，小批量，不是全模型筛选 |
| Streamlit LAN UI | 可用 |
| FastAPI 最小 run/status/result | experimental |
| OPN alignment | aligned except known MATLAB compatibility differences，不是 fully aligned |
| hLF 710aa alignment | Python target `hLF` 对应 artifact target `hLF_PROJECT_710`，不是 old MATLAB hLF fully aligned |

## 固定边界

- 不修改 `Code/`、`Model/`、`Enzymedata/`、`Results/`。
- 不启动 MATLAB，除非任务明确要求生成新 harness artifact。
- 不跑全量测试；日常只跑相关聚焦测试。
- 不把核心模型逻辑写进 `app/`。`app/` 只做 UI/service/API 编排。
- Python corrected 默认行为不回退到旧 MATLAB bug/legacy 行为。
- `hLF` 当前指用户提供的 710 aa 项目序列；对应 MATLAB artifact target 是 `hLF_PROJECT_710`。

## 已完成的收口切片

- 提交边界清理：feature branch / Draft PR 已建立，保护目录和 local artifacts 不进入源码提交。
- 历史文档冻结：active docs 只保留 architecture / next slices / release validation。
- service facade 拆分：schema、runner、target catalog、gene catalog、background tasks 已拆出 owner modules。
- gene/reaction 解析统一：pipeline 与 preview 复用 engine-side 解析 helper。
- hLF 710aa 对齐展示统一：`hLF` 与 `hLF_PROJECT_710` 的语义已在 service/report warning 中区分。
- 慢测与 release 验证分层：日常、focused、slow、release 命令已固化。
- Streamlit 结果展示 helper：候选表中文列名、状态归一化、路径图空值显示已收敛到 `simulation_display.py`。

## 一轮对话可完成的小目标

### 1. UI 页面拆分

目标：把 `app/ui/views/simulation.py` 从单一大文件拆成可维护的视图模块。

建议拆分：
- `simulation_builder.py`
- `simulation_results.py`
- `candidate_path_graph.py`
- `matlab_reference.py`

每轮只拆一个模块，不改 UI 行为。

验证：

```powershell
python -m py_compile app\ui\views\simulation.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
```

### 2. 文档和测试继续瘦身

目标：删除对下一阶段开发没有帮助的历史材料，保留 release/review 需要的最小边界测试。

可完成内容：
- 删除旧 route-by-route OPN/hLF draft 输入测试。
- 删除只描述已废弃迁移路线的大文档。
- 保留 active docs 边界测试和 release slow gate 测试。
- 不删除当前 service/UI/pipeline/alignment 的聚焦门禁。

验证：

```powershell
python -m pytest -q tests\test_docs_active_boundary.py tests\test_review_package_boundaries.py tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

### 3. Streamlit 结果展示瘦身

目标：减少 UI 展示层中重复的中文列名、warning 组装和 markdown 拼接。

可完成内容：
- 把候选表列名映射收敛到一个小 helper。
- 把 alignment/status badge 文案收敛到一个小 helper。
- 不改 result payload，不改 pipeline。

验证：

```powershell
python -m py_compile app\ui\views\simulation_results.py app\ui\views\simulation.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
```

### 4. FastAPI 取舍（暂缓）

目标：决定 `app/api/` 是 active experimental，还是暂时移除。

可完成内容：
- 若保留：文档标注 experimental，并确保只调用 app service facade。
- 若删除：删除 `app/api/` 与对应测试。
- 不新增认证、队列、下载等新功能。

验证：

```powershell
python -m pytest -q tests\test_pichia_fastapi_entrypoints.py
```

### 5. Release 前验证 rehearsal

目标：在不启动 MATLAB 的前提下，跑一次 release 文档里的非 slow 验证命令，确认 Draft PR 可 review。

可完成内容：
- 跑 compileall / focused pytest。
- 记录任何跳过项或环境阻塞。
- 确认保护目录 diff 为空。

验证：

```powershell
python -m compileall app python_pichia tests
python -m pytest -q tests\test_pichia_secretion_service_contract.py tests\test_review_package_boundaries.py tests\test_docs_active_boundary.py
git diff --name-only -- Code Model Enzymedata Results
```

## 推荐顺序

1. UI 页面拆分。
2. 文档和测试继续瘦身。
3. Streamlit 结果展示瘦身。
4. Release 前验证 rehearsal。
5. FastAPI 取舍（暂缓）。

## 暂不进入下一阶段的内容

- Humanized 糖基化。
- 温度敏感性。
- 三物种全量迁移。
- 论文 figure 复现。
- 全模型 KO/OE 批量筛选。
- 新 MATLAB baseline 生成。
