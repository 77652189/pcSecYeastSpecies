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
- Streamlit 基因扰动 UI：基因库展示、候选文本解析、KO/OE 表单已拆成独立视图 helper。

## 下一阶段优先级

后续不再围绕 MATLAB 项目全量迁移推进，而是按产品分析能力扩展。优先级固定为：

1. 目标蛋白蛋白成本分析。
2. 目标蛋白生长分析。
3. 代谢工程靶点分析。
4. 人源化糖基化蛋白成本分析。
5. 人源化糖基化生长分析。

## 一轮对话可完成的小目标

### 1. 目标蛋白蛋白成本分析

目标：基于当前 target/secretion plan/constraint summary，给出目标蛋白造成的主要蛋白成本来源，不先改变求解数值行为。

可完成内容：
- 在 `python_pichia/` 下建立正式 cost analysis 入口，例如 `pcsec_pichia.analysis` 或 `pcsec_pichia.costs`。
- 汇总 target sequence length、MW、DSB/NG/OG/GPI、translation、ER folding、glycosylation、misfolding 等成本项。
- 输出结构化 payload，供 report/UI 消费。
- 第一轮只做 OPN/hLF/custom target 的 deterministic summary，不做新 LP 求解和参数拟合。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_target_entrypoints.py python_pichia\tests\test_secretion_plan_entrypoints.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

### 2. 目标蛋白生长分析

目标：把已有 small-grid growth tradeoff 提升为面向目标蛋白解释的生长影响分析。

可完成内容：
- 复用 `simulation.run_growth_tradeoff`，不新增慢批量任务。
- 对 OPN/hLF/custom target 输出 growth vs secretion 的小表和解释字段。
- 在 report/UI 中区分 draft tradeoff 与真实发酵增长预测。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_simulation_entrypoints.py python_pichia\tests\test_pipeline_entrypoints.py
git diff --name-only -- Code Model Enzymedata Results
```

### 3. 代谢工程靶点分析

目标：在小批量 KO/OE 基础上提升解释能力，而不是立即做全模型筛选。

可完成内容：
- 保留手动小候选集。
- 对 KO/OE rows 增加 pathway/process grouping、effect summary 和 unresolved diagnostics。
- 继续明确 OE gene 是 reaction-level proxy。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

### 4. 人源化糖基化蛋白成本分析

目标：在 native OPN/hLF 稳定基础上引入 humanized glycosylation 的成本解释，不先承诺完整 pathway engineering。

可完成内容：
- 明确 humanized glycosylation mode 的输入 schema。
- 先做 cost/report 层面的新增成本项，不改变 corrected 默认 pipeline。
- 输出 warning：humanized 参数仍需实验/文献校准。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_target_entrypoints.py python_pichia\tests\test_reports_entrypoints.py
git diff --name-only -- Code Model Enzymedata Results
```

### 5. 人源化糖基化生长分析

目标：在人源化糖基化成本项稳定后，再分析其 growth tradeoff。

可完成内容：
- 复用目标蛋白生长分析入口。
- 增加 humanized glycosylation mode 对 growth tradeoff summary 的影响解释。
- 不做三物种迁移或新 MATLAB baseline。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_simulation_entrypoints.py python_pichia\tests\test_reports_entrypoints.py
git diff --name-only -- Code Model Enzymedata Results
```

## 明确不做的内容

- 三物种 MATLAB 项目全量迁移。
- 论文 figure 复现。
- 新 MATLAB baseline 自动生成。
- 温度敏感性。
- 全模型 KO/OE 批量筛选，除非后续单独立项。
