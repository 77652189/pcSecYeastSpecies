# pcSecPichia 下一步计划

状态：active  
最后更新：2026-07-07

## 工作原则

- 一次只做一个任务分支，完成后 review、测试、commit、合并回 `main`，再开始下一项。
- Python 核心逻辑放在 `python_pichia`，`app/services` 只做 facade，`app/ui` 只做展示。
- 不修改 `Code/`、`Model/`、`Enzymedata/`、`Results/`，除非明确声明为科学资产变更。
- 新运行产物、cache、报告和验证证据默认写入 ignored `local_runs/`。
- 输出必须保留模型边界，不声明 mg/L 绝对产量或实验成功率。

## 当前推荐优先级

### 1. BLAST/RBH Streamlit 同源审计

目标：把酿酒酵母的分泌工程知识安全迁移到 Pichia 模型前，建立可审计的同源证据层，并让研发用户能在 Streamlit 中查看“基因命名标准化 + 同源规则迁移评估”结果。

范围：

- 从本地 pcSecYeast / pcSecPichia protein sequence 资产导出 FASTA。
- 用本地 BLAST+ 运行 SCE -> Pichia 和 Pichia -> SCE 双向 blastp。
- 计算 reciprocal best hit、identity、evalue、query coverage、subject coverage。
- 合并当前 Pichia GEM `gene_index`，明确候选是否是模型可操作 gene。
- 输出 homology cache、name audit、rule-transfer audit 和 summary。
- 通过 `app/services` 只读 cache，通过 `app/ui` 展示、筛选、导出和解释结果。
- 先覆盖 `SECRETION_GENE_CATALOG` 相关 query，但不自动写回 catalog，也不把同源命中直接当作 KO/OE 模型 gene。
- CLI / scripts 只作为离线 cache builder，不是最终产品入口；Streamlit runtime 默认不跑 BLAST。

设计文档：[BLAST/RBH 同源映射架构](pichia_homology_crosswalk_architecture.md)

固定执行轮次：

1. Round 0：目标契约和现状审计，只更新 active 文档，明确 Streamlit 产品目标。
2. Round 1：核心 schema 和 audit 输出落地，生成 homology / name / rule-transfer 三类 cache。
3. Round 2：新增 app service facade，只读 cache、过滤、summary、导出，不运行 BLAST。
4. Round 3：新增 Streamlit 页面“基因命名与同源规则审计”。
5. Round 4：把 homology evidence 作为并行解释层接入 KO/OE 预览和 evidence summary，不改变 recommendation tier 边界。
6. Round 5：预留离线外部数据库名称校对契约，不做 runtime 联网。
7. Round 6：最终文档收束，说明已实现内容、边界和使用方式。

每轮完成后都要 review 本轮 diff、运行 focused tests / compileall / 保护目录检查，并将本轮相关改动 commit + push。

当前 cache builder 验收：

```powershell
python -m pytest -q python_pichia\tests\test_homology_sequence_sources.py python_pichia\tests\test_homology_rbh.py python_pichia\tests\test_homology_crosswalk.py
python -m compileall -q python_pichia\src
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### 2. Shadow LP service/backend toggle

目标：把 `solve_shadow_secretion_capacity(...)` 接入 service / pipeline 的 opt-in 后端选择，让 reference 和 shadow 能双轨对比。

范围：

- 新增 `solver_mode="reference" | "shadow"` 或等价配置。
- 默认仍为 reference，避免立刻改变用户工作流。
- 仅在 opt-in 时调用 shadow path。
- 报告中显示 `solver_mode`、backend、validation status 和 non-mg/L warning。
- reference solver 继续保留为 validation boundary。

建议验收：

```powershell
python -m pytest -q python_pichia\tests\test_shadow_lp_secretion_capacity_wrapper.py python_pichia\tests\test_shadow_lp_compare_mode.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
python -m compileall -q app python_pichia\src
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### 3. KO/OE evidence integration

目标：把同源证据、模型 GPR 可执行性、OE reaction proxy 和 phenotype evidence 合成一张解释表，而不是让 UI 或 service 自行判断。

范围：

- 输入：homology cache、gene capability profile、phenotype evidence、screen rows。
- 输出：candidate evidence summary。
- 明确区分 `homology_supported`、`model_executable`、`phenotype_supported`、`experiment_calibrated`。
- 不能仅靠 RBH 或 annotation 升级为 experiment calibrated。
- 保留 proxy warning 和 manual review reason。

建议验收：

```powershell
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py python_pichia\tests\test_yield_improvement_recommendations.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

## 暂不做

- 不把 BLAST/RBH 命中自动写入 curated catalog。
- 不把同源命中直接当作 Pichia 模型 gene_id。
- 不联网做实时数据库查询。
- 不把 COBRApy / optlang 设为 full pcSec 默认后端。
- 不做三物种全量迁移或论文 figure pipeline 重建。
- 不迁移历史 `Results/`、不改 Git LFS 或仓库历史。
