# pcSecPichia 当前技术执行计划

状态：active  
最后更新：2026-07-16

## 当前执行位置

```yaml
current_program: mvp_directions_1_to_3
current_slice: direction_3_erad_constraint_activation
slice_status: complete_kept_optional
direction_3_round0_status: complete
```

项目优先级和范围以 `docs/EXECUTION_PLAN.md` 为准；实际执行范围以 `docs/handoff.md` 为准。本文件只描述当前技术顺序，不恢复已完成阶段的逐轮日志。本文件下方大部分内容记录的是已完成的 `direction_3_secretory_resource_round_0` 技术细节，作为历史参考保留；`direction_3_erad_constraint_activation` 切片（验证与激活决定）也已完成，具体技术顺序未在本文件展开，验收证据读 `docs/EXECUTION_PLAN.md`。当前运行产物统一进入 ignored `local_runs/`，具体规则见 `docs/data_and_results_policy.md`；历史 `Results/` 迁移、Git LFS 改造或仓库历史瘦身不属于当前技术计划范围。

## 已完成基线

- 方向 1：研发发酵宽表 CSV/XLSX/JSONL、质量门禁、prediction linkage、校准资格、报告和脱敏回放已验收；等待获批真实数据——一份周报级别的进度汇报不构成这里说的"真实数据"，量级不够支撑回填。
- 方向 2：reaction proxy、relative uncalibrated、absolute unavailable 和 not executable 已由 `python_pichia` 核心层统一判定。
- relative uncalibrated 使用独立 current-model enzyme-coupling 成对求解；缺审核 anchor 的 absolute 请求在求解前拒绝。
- hLF/OPN relative smoke、旧 proxy/feature-off 回归、service/UI 透传和产品状态验收已通过。
- 绝对容量研究层仍缺 `reviewed_baseline_capacity`，保持 unavailable。
- 方向 3 Round 0 已完成：`pcsec_pichia.secretory_resources` 冻结了全部七类资源的架构与可执行契约，并逐一核实了每类是否有真实 kcat 数据支撑（详见下方"Round 0 完成后的核实结论"）。
- 基因命名、BLAST/RBH、外部 GPR 候选、Shadow LP、COBRApy import/QA 和 LLM judged report 已存在，不重复立项。

## 当前目标

`direction_3_secretory_resource_round_0` 已完成并通过验收：Round 0 只建立可执行架构与验证契约，不实现完整分泌机制求解，不产生新的 KO/OE 排名，不接 Streamlit 产品页面——这些边界在 Round 0 期间没有被突破。`direction_3_erad_constraint_activation` 也已完成并通过验收：验证了 hLF/OPN 求解可行性，决定 ERAD/misfolding 约束保持可选、不改默认值（见 `docs/EXECUTION_PLAN.md` 第 6-7 节）；目标蛋白降解通路（PEP4/PRB1/YPS）建模明确不做。当前没有处于进行中状态的技术工作；下一步范围由用户决定何时推进。本文件不自行扩大范围、不生成后续实现提示词。

## Round 0 完成后的核实结论

冻结架构之外，额外核实了每类资源的 handle 是否已有真实（非占位）kcat 数据：

- 转运、二硫键、糖基化、囊泡运输、folding/chaperone：确认有真实 kcat；前四类的约束已经无条件参与现有模型的每次求解，与 Round 0 新架构层无关。
- ER quality control/ERAD/proteasome：确认有真实 kcat，但约束生成默认关闭（`enable_misfolding_constraint`/`write_misfolding_constraints`，default False）。用小范围候选（约10个 ERAD/蛋白酶体相关基因/反应）做过一次开关前后对照测试：打开后，蛋白酶体、CDC48、HRD1 复合体、DSK2/RAD23 穿梭复合体等敲除开始显示 5%-14% 的真实产量差异——说明打开这个开关不是空转，会改变候选排名。
- target-specific translation/degradation cost：唯一确认真正缺数据的一类——具体到 hLF/OPN 自身的降解速率（kdeg），全仓库没有真实值。进一步查证还发现一个独立于 kdeg 数据缺口之外的结构性缺口：`r_{protein_id}_subunit_degradation` 这类降解反应完全没有 GPR，任何基因（含 PEP4/PRB1 这类液泡蛋白酶）敲除都无法影响它；PEP4/PRB1 在现有基因目录里用的模型基因 ID 本身也已经被标注为低置信度、待人工复核（`services/gene_rule_overlay.py`），不建议在没有可靠基因身份和真实动力学数据前把这条通路强行接入模型。

## 层级边界

```text
metabolic model
  -> 代谢计量、GPR、培养条件和生长约束
protein resource
  -> 酶容量、翻译、蛋白质量和通用蛋白资源
secretory resource
  -> 转运、折叠、二硫键、糖基化、质量控制和囊泡运输资源
experimental calibration
  -> 实验方向、排序和风险校准，不修改前三层科学资产
```

- 现有 `secretion_plan`、secretory enzyme data、target plan 和 Shadow LP 是输入，不由新层复制实现。
- 文献基因名单、数据库注释、BLAST/RBH 或同源关系只能形成候选证据，不能直接生成可执行约束。
- 没有合格参数、当前模型 handle 或单位转换链时，状态必须为 unavailable/not executable。
- hLF 与 OPN 的目标特异成本必须分别表达，不能把一个 target 的参数复制成另一个 target 的证据。

## Round 0 数据契约

每个资源定义至少包含：

- 稳定 `resource_id`、资源类别和生物过程。
- canonical unit 与允许的转换链。
- 当前模型 variable/reaction/constraint handle。
- source ref、version、hash/license 和 evidence class。
- host、target、培养条件和 model fingerprint 适用范围。
- nominal/lower/upper 或明确 unavailable 状态。
- uncertainty、warnings、limitations 和人工复核状态。
- feature flag、baseline 行为和 feature-off 回归要求。

首批资源类别只冻结语义：

1. ER translocation。
2. folding/chaperone。
3. disulfide bond formation。
4. glycosylation。
5. ER quality control、UPR/ERAD/proteasome。
6. vesicle trafficking/exocytosis。
7. target-specific translation、modification 和 degradation cost。

资源类别存在不表示当前可执行；每项必须独立报告 `executable`、`evidence_only`、`unavailable` 或 `conflict`。

## 模块与公共边界

Round 0 优先在 `python_pichia/src/pcsec_pichia/` 内建立可移除的 `secretory_resources` 核心包；最终目录名可根据现有调用链微调，但不得把核心判断写入 service/UI。

建议职责：

- `schema.py`：frozen dataclass、枚举、validation 和序列化。
- `catalog.py`：汇总当前模型/target/secretory data 的资源候选，不联网。
- `validation.py`：单位、适用范围、来源、冲突和可执行性门禁。
- `planning.py`：生成 backend-neutral resource plan；Round 0 不调用 solver。

首批公共契约应覆盖以下能力，具体命名以代码审计后的最小一致接口为准：

```python
build_secretory_resource_catalog(prepared_model, target, evidence=None)
validate_secretory_resource_catalog(catalog)
plan_secretory_resource_constraints(catalog, config)
summarize_secretory_resource_coverage(catalog)
```

`app/services`、Streamlit、pipeline 和完整 solver wiring 不属于 Round 0，除非为了证明核心契约可导入而增加最小 facade contract test。

## 固定执行循环

1. 读取 `docs/EXECUTION_PLAN.md`、架构文档、handoff 和数据治理策略。
2. 使用 codebase-memory-mcp 审计现有 secretory/target/constraint 调用链，确认所有权和可复用 handle。
3. 冻结本轮输入、输出、状态、单位和验收标准。
4. 先实现 `python_pichia` 数据契约、catalog/validation/planning 最小闭环。
5. 增加 focused tests，覆盖 hLF、OPN、缺参数、冲突、evidence-only、feature-off 和错误单位。
6. 运行只读 smoke，产物只写 ignored `local_runs/secretory_resources/`。
7. review/fix/verify，最多 3 轮。
8. 只在实现状态实质变化时更新 handoff；不得用文档修改代替产品代码。
9. 检查保护目录、依赖、密钥和 ignore 边界。

可以把现有调用链审计、测试缺口分析和 hLF/OPN 只读证据盘点并行交给子 Agent；schema、单位、状态、公共 API 和最终整合由主 Agent 决定。

## 验收与停止条件

Round 0 完成必须满足：

- 四层所有权清楚，现有模型逻辑没有被复制或改名。
- 每类资源有单位、来源、适用条件、不确定性、开关和 baseline/feature-off 契约。
- evidence-only 不能被提升为 executable。
- hLF/OPN 能生成独立 coverage 和 unavailable/conflict 解释。
- backend-neutral plan 不调用 solver，不修改稳定模型资产。
- focused tests、compileall、保护目录、依赖、密钥和 ignore 检查通过。

出现以下情况立即停止：

- 必须用默认值、伪参数或未审核实验数据才能声明资源可执行。
- 需要新增默认 solver、重写稳定约束层或修改 `Code/Model/Enzymedata/Results`。
- 资源单位或当前模型 handle 无法定义且没有诚实 unavailable 路径。
- 工作开始转向完整机制求解、组合搜索或方向 4。

## 慢速测试网关

三类真实求解回归默认跳过，需要显式设置对应环境变量才会运行：

- `PCSEC_RUN_SLOW_PIPELINE_TESTS="1"`：单候选全流程求解（`run_pichia_secretion_simulation`），覆盖 `python_pichia/tests/test_pipeline_entrypoints.py`。
- `PCSEC_RUN_SLOW_SCREEN_TESTS="1"`：全模型 KO/OE 批量筛选，覆盖 `python_pichia/tests/test_screens_entrypoints.py`。
- `PCSEC_RUN_SLOW_PROBE_TESTS="1"`：probe 迁移回归，对照 MATLAB harness / baseline 生成的既有产物，覆盖 `python_pichia/tests/test_probe_migration.py`。

## 后续顺序

Round 0 验收后由项目级计划重新确认方向 3 的实现轮次——已完成，结论是继续做 `direction_3_erad_constraint_activation`，范围只到验证与激活决定，不是完整求解；这一轮验证与激活工作本身也已完成（决定：保持可选）。除非 handoff 明确记录范围扩大，不实现完整 UPR/ERAD/糖基化/囊泡资源求解，不进入组合改造。

## 当前下一步

`direction_3_secretory_resource_round_0` 已完成架构和可执行契约，已按要求停在验收点。`direction_3_erad_constraint_activation` 也已完成：验证了 hLF 求解可行性、理解了已知 MATLAB 兼容性差异、决定 ERAD/misfolding 约束保持可选（不改默认值），验收证据读 `docs/EXECUTION_PLAN.md`。目标蛋白降解通路（PEP4/PRB1/YPS）建模和方向 4 仍不做，不实现完整机制求解。当前列出的技术工作已全部交付完成，没有处于进行中状态的技术工作；下一步范围由用户决定何时推进，本文件不生成后续实现提示词。当前边界读取 `docs/handoff.md`。
