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

## 一轮对话可完成的小目标

### 1. 提交边界清理

目标：让工作区能被 review，而不是一团 untracked 文件。

可完成内容：
- 列出应进入正式源码的新增文件。
- 列出应保留为 local artifact 的目录。
- 确认 `python_pichia/`、`app/services/pichia_secretion_service.py`、关键 tests/docs 的状态。
- 不提交，只产出清单或最小 `.gitignore` 修正。

验证：

```powershell
git status --short
git diff --name-only -- Code Model Enzymedata Results
```

### 2. 历史文档冻结

目标：保留历史决策证据，但不再把旧迁移路线作为 active plan。

可完成内容：
- 对已删除/过时迁移文档做最终决定：删除、恢复为 frozen，或移入 archive。
- 新增一个短索引说明哪些文档是 active，哪些是 frozen。
- 不再扩写 stage/route-by-route 旧路线。

验证：

```powershell
git status --short docs
```

### 3. Service facade 拆分

目标：降低 `app/services/pichia_secretion_service.py` 的耦合。

建议拆分：
- `pichia_secretion_schema.py`：request/response dataclass。
- `pichia_secretion_runner.py`：调用 `python_pichia.pipeline`。
- `pichia_gene_catalog_service.py`：KO/OE catalog facade。
- `pichia_background_tasks.py`：后台任务与 last-result cache。

每轮只拆一块，保持 import 兼容。

验证：

```powershell
python -m py_compile app\services\pichia_secretion_service.py app\ui\views\simulation.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
```

### 4. UI 页面拆分

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

### 5. Gene/reaction 解析统一

目标：避免“预检能解析，正式运行不能解析”的分叉。

可完成内容：
- 在 `python_pichia` 中建立唯一解析 helper。
- pipeline 与 service preview 都调用同一 helper。
- 保留 OE gene proxy 的明确警告。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_pipeline_entrypoints.py tests\test_pichia_secretion_service_contract.py
```

### 6. hLF 710aa 对齐展示统一

目标：用户看到的状态不再混淆旧 hLF MATLAB failure 与当前项目 hLF 710 artifact。

可完成内容：
- UI/service warning 改成双层表述：
  - 旧 MATLAB `hLF` baseline：historical `matlab_failed`。
  - 当前项目 `hLF` 710aa：使用 `hLF_PROJECT_710` artifact，状态为 `aligned_except_known_matlab_compatibility_differences`，但不是 fully aligned。
- report/summary 中保留 `python_target_id=hLF` 与 `alignment_artifact_target_id=hLF_PROJECT_710`。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_alignment_entrypoints.py python_pichia\tests\test_pipeline_entrypoints.py tests\test_pichia_secretion_service_contract.py
```

### 7. FastAPI 取舍

目标：决定 `app/api/` 是 active experimental，还是暂时移除。

可完成内容：
- 若保留：文档标注 experimental，并确保只调用 app service facade。
- 若删除：删除 `app/api/` 与对应测试。
- 不新增认证、队列、下载等新功能。

验证：

```powershell
python -m pytest -q tests\test_pichia_fastapi_entrypoints.py
```

### 8. 慢测与 release 验证分层

目标：日常测试保持快，release 前仍可跑真实求解。

可完成内容：
- 把慢求解测试统一标记为环境变量开关。
- 文档中给出 daily / focused / slow / release 四级命令。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_pipeline_entrypoints.py
$env:PCSEC_RUN_SLOW_PIPELINE_TESTS="1"; python -m pytest -q python_pichia\tests\test_pipeline_entrypoints.py
```

## 推荐顺序

1. 提交边界清理。
2. 历史文档冻结。
3. hLF 710aa 对齐展示统一。
4. service facade 拆分。
5. UI 页面拆分。
6. gene/reaction 解析统一。
7. FastAPI 取舍。
8. 慢测/release 验证分层完善。

## 暂不进入下一阶段的内容

- Humanized 糖基化。
- 温度敏感性。
- 三物种全量迁移。
- 论文 figure 复现。
- 全模型 KO/OE 批量筛选。
- 新 MATLAB baseline 生成。
