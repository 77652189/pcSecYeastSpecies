# pcSecPichia 下一步规划

状态：active

## 工作原则

- 每一步都围绕生物研发同事的产量提升问题展开。
- 先做可解释、可验证的小切片，再扩大模拟范围。
- Python 侧以 `python_pichia` 为 engine，`app/services` 和 `app/ui` 只做编排与展示。
- 不修改 `Code/`、`Model/`、`Enzymedata/`、`Results/`。
- 不声称预测绝对产量或实验成功率。

## 验证边界

日常只跑与当前切片相关的 focused tests。真实求解类慢测必须显式打开环境变量：

```powershell
$env:PCSEC_RUN_SLOW_PIPELINE_TESTS="1"
$env:PCSEC_RUN_SLOW_SCREEN_TESTS="1"
$env:PCSEC_RUN_SLOW_PROBE_TESTS="1"
```

慢测和本地 artifact 不属于默认门禁，包括：

- MATLAB harness / baseline 生成。
- 新 LP artifact 生成。
- 全模型 KO/OE 批量筛选。
- 长网格 growth tradeoff。
- 论文 figure pipeline 重建。

## 已完成的基础

- hLF / OPN / custom target 的 Python corrected 分泌仿真入口。
- 密码子优化和信号肽筛选结果可作为研发输入背景。
- 培养基和碳源条件层，包括 mixed-carbon objective probe。
- 目标蛋白成本、growth tradeoff、LP attribution 的解释型输出。
- KO/OE 候选预览和 screen 结果解释。
- KO/OE 外部表型证据层和 recommendation tier。
- Streamlit 展示、后台任务恢复、API/service facade 和 gene evidence cache 工具。
- 数据与结果治理 checkpoint：`Results/` 明确为 legacy MATLAB reference，当前运行产物统一进入 ignored `local_runs/`。

## 近期优先级

### 1. KO/OE 模拟工作流收口

目标：让研发同事能够围绕 hLF/OPN 提出候选 KO/OE，并得到清晰的模型解释和实验优先级。

范围：

- 范围已扩大为全模型批量筛查（约1025个基因，同时评估生长速率与分泌产量两方面影响），设计细节和取舍见 `pichia_ko_oe_genome_screen_design.md`。
- 强化 KO/OE 输入、预览、结果表和报告中的证据分层。
- 继续明确 OE reaction proxy 的限制。
- 对 essential KO、表型证据冲突、annotation-only 候选给出人工复核提示。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py python_pichia\tests\test_pipeline_entrypoints.py python_pichia\tests\test_reports_entrypoints.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

### 2. MATLAB 原有功能迁移盘点

目标：把“必须迁移”和“暂不迁移”的 MATLAB 功能列清楚，避免盲目全量搬运。

范围：

- 对照当前研发工作流，只盘点与 hLF/OPN 产量提升相关的 MATLAB 功能。
- 标注每项功能的 Python 状态：已迁移、部分迁移、只保留 reference、暂不做。
- 对需要迁移的能力定义最小验收测试。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_pipeline_entrypoints.py python_pichia\tests\test_reports_entrypoints.py
git diff --name-only -- Code Model Enzymedata Results
```

### 3. 研发可读报告整理

目标：让输出更像研发讨论材料，而不是求解器日志。

范围：

- 对 hLF/OPN 输出同一套摘要结构。
- 把候选分为推荐、证据支持、需人工复核、增长风险和不推荐。
- 报告中保留模型边界和实验限制。
- 不加入 mg/L 绝对产量预测。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_reports_entrypoints.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

### 4. 证据数据维护

目标：把当前内置/fixture evidence 逐步变成可维护的小型 curated 数据层。

范围：

- 先维护小型人工 curated 数据，不做大规模联网抓取。
- 记录 evidence source、confidence、target context 和 recommended use。
- 对同一 gene 的 KO 与 OE 分开维护。
- 保留 cache build 脚本作为辅助工具，不让联网构建成为默认工作流。

验证：

```powershell
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

### 5. 数据资产后续治理

目标：在不误伤历史 MATLAB 数据的前提下，逐步决定哪些大文件需要 Git LFS、GitHub Release 或外部存储。

范围：

- 第一轮只保留盘点和边界测试，不迁移历史结果。
- 如需移动 `Results/` 或大文件瘦身，单独建立 checkpoint。
- 新生成运行产物继续写入 `local_runs/`，不得漂移到科学资产目录。

## 暂不做

- 全模型 KO/OE 自动筛选。
- 新 MATLAB baseline 自动生成。
- 三物种完整迁移。
- 历史 `Results/` 迁移、Git LFS 改造或仓库历史瘦身。
- 温度敏感性和人源化糖基化完整 pathway engineering。
- 面向论文复现的 figure pipeline 重建。

这些内容可以单独立项，但不进入当前 hLF/OPN 产量提升工作流的近期计划。
