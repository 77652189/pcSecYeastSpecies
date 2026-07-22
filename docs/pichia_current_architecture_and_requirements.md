# pcSecPichia 当前架构与能力边界

状态：active  
最后更新：2026-07-20

## 原始研发目标

项目要解决的问题是：围绕 hLF、OPN 等目标蛋白，通过 KO、OE 和分泌路径改造，找到更可能提高 Pichia 分泌表现的工程方案，并降低研发同事阅读模型结果和选择实验候选的成本。

当前项目的正确定位是“KO/OE 候选生成、模型验证、证据排序和实验反馈系统”，不是实验产量预测器。

OE 产品能力分为两个独立层级：相对、未校准的决策层用于候选比较和风险解释；绝对 gene-capacity 研究层只有在存在经审核的 baseline capacity 时才可执行。两层不能通过默认值或改名相互替代。

## 当前结论

当前系统可以回答：

- 哪些单基因 KO 在模型 GPR 中可执行。
- 哪些 OE 方向在反应容量层面可能改善目标蛋白分泌。
- 候选对生长、分泌和模型可行性的相对影响。
- 候选是否存在 essentiality、growth risk、proxy、同源映射或外部证据风险。
- hLF 和 OPN 是否出现不同的候选优先级。
- 哪些候选值得进入人工复核和小规模实验。
- 哪些相对 OE 场景可用于候选比较，哪些只能保持 unavailable。
- 当前哪个约束（反应容量、蛋白资源等）对某个 target 的分泌影响最大（draft 级证据，见 LP 敏感度归因）。

当前系统不能可靠回答：

- 某个改造会提高多少 mg/L。
- 某个候选的真实实验成功概率。
- reaction-level OE proxy 等同于真实 gene-level OE 后的表达变化。
- 缺少审核 baseline capacity 时的绝对 gene-level OE 容量。
- 模型外基因、调控网络和全部分泌通路基因的定量作用。
- 多基因组合的上位性和长期培养稳定性。
- 多亚基复合体内部哪个亚基是真正限速步骤——模型把复合体形成建模成一个反应，无法拆分归因到具体某个基因。
- 同一参数（如 kcat）按基因区分的真实测量不确定性——现有实现对所有基因使用统一假设的不确定性宽度，不是逐基因的真实置信区间。

## 研发工作流

```text
目标蛋白与培养条件
  -> pcSecPichia / Shadow LP 模型准备
  -> 单基因 KO 与 OE proxy 筛查
  -> 分泌、生长、可行性和资源代价计算
  -> gene/GPR/同源/外部数据库/表型证据合并
  -> recommendation tier 与风险分层
  -> 相对 OE 决策层 / 绝对容量可用性门禁
  -> fact pack + 程序 validator + LLM writer + Judge
  -> hLF / OPN 实验候选清单
  -> 实验反馈与下一轮排序校准
```

实验反馈闭环和 OE 产品分层已经完成。reaction proxy、relative uncalibrated、absolute unavailable 和 not executable 由核心层统一判定；绝对 hLF/OPN gene-capacity 因缺审核 baseline capacity 保持不可执行。当前架构工作转向独立 secretory resource layer。

## 已有能力

### 模型与筛查

- 加载 pcSecPichia 参考模型和目标蛋白输入。
- 构建内置目标和 custom target 的分泌路径。
- 设置培养基、碳源、固定生长率和目标蛋白负担。
- 执行 genome-wide / catalog-level 单基因 KO 筛查。
- 生成 OE reaction-capacity proxy，并保留“不是 gene-level OE”的 warning。
- 输出 secretion ratio、growth retention、最大可行生长率、protein cost 和 solver 状态。
- 对 essential KO、growth risk、不可解析和模型不可执行候选进行降级。

### 基因、GPR 与证据

- 模型 gene_index、蛋白序列和 UniProt/NCBI/KEGG 等命名注释可交叉查询。
- BLAST/RBH cache 可提供 SCE 到 Pichia 的同源候选、name audit 和 rule-transfer audit。
- 外部数据库 fetcher 可受控生成 UniProt、NCBI、SGD 和外部 GEM/GPR 证据缓存。
- 外部 GPR 只作为候选证据；映射到当前 reaction/gene 后才可成为模型可执行规则。
- phenotype evidence、数据库注释、模型 GPR、OE proxy 和同源证据保持分层。
- **模型的 GPR（基因-反应关联）只覆盖代谢反应（约 2732 个：中心代谢、氨基酸合成等，基因为 `PAS_chr*` 位点）。分泌机器层（chaperone、translocon、糖基化、COPII、核糖体、蛋白酶体等复合体形成反应）在参考模型里完全没有 GPR，其酶本身也不作为基因节点存在（已直接加载 `Model/pcSecPichia.mat` 核实：2793 个复合体形成反应，0 个有基因关联）。** 因此存在两条不同的筛查路径：走 GPR 的全基因组基因筛查只能触达代谢基因；分泌机器的 KO/OE 干预是在**复合体/反应层面**进行（直接把复合体形成反应流量压到 0，或把其 kcat 乘以倍数），通过分泌耦合约束传导到目标蛋白分泌，不经过"基因→GPR→反应"链路。这套复合体层面干预是真实有效、且针对具体目标蛋白的——例：过表达 PDI/ERO1/ERV2 复合体使 hLF 分泌 +8.15%，而无二硫键的 OPN 几乎不变。使用时需要人工把"过表达某复合体"翻译成"过表达对应的那几个基因"，且这类分泌机器候选只能通过 curated catalog 进入筛查，全基因组基因筛查发现不了它们。

### 求解与质量验证

- 当前大模型主路径使用结构化约束和 SciPy HiGHS / reference solve。
- Shadow LP 已承载 pcSec 约束并与内置 hLF/OPN reference 结果高精度对齐。
- COBRApy 用于外部 GEM import、基础 QA 和小模型语义验证，不作为 full pcSec 默认后端。
- Streamlit 可手动运行和读取 Shadow LP cross-check 报告。

### 结果解释

- 单候选 pipeline 无条件计算 LP 级敏感度归因（`analyze_target_protein_lp_attribution`，基于 SciPy HiGHS marginals），可给出当前解在哪个约束块、哪个具体反应的 bound 上最吃紧；这是 draft 级证据，不是 MATLAB/SoPlex 对齐的 shadow price，符号需谨慎解读，且仅覆盖单候选路径——批量筛查（genome-wide/catalog screen）目前不透出这项信息。**下界（`bound_type=lower`）上的大 marginal 不能当作 OE 候选线索**：下界是最低要求类约束，OE 放宽的是上限产能，不会缓解下界卡点；这条规则已写进函数自身的 `warnings`，因为 hLF 的 PDI1 单独反应和 hLF/OPN 的核糖体装配反应都出现过"下界 marginal 很大、但实测 OE 效果为零"的假信号。
- KO/OE preview、screen rows、recommendation、历史结果和 Streamlit 均可透传标准命名与证据字段。
- LLM 只读取程序生成的 fact pack。
- 程序 validator 检查 schema、evidence_id、target、gene/reaction、数值、KO/OE 分组和最低信息覆盖。
- Judge 读取 fact pack、writer 输出和 validator 结果，不合格时反馈重写。
- 最终报告按 hLF / OPN 分区，引用表由程序生成。

### 实验反馈

- CSV/XLSX/JSONL 实验记录可导入统一 schema，并保留原始值、单位、重复、失败和排除状态。
- prediction linkage 显式区分 matched、ambiguous、missing 和 context mismatch。
- hLF/OPN 回放可生成校准记录、方向一致性、Top-K 描述指标和证据不足提示。
- 实验导入不会自动修改 recommendation tier、curated evidence 或模型约束。

## 核心证据边界

| 层级 | 能说明什么 | 不能说明什么 |
| --- | --- | --- |
| 数据库注释 | 名称、功能和来源一致性 | 真实 KO/OE 表型 |
| BLAST/RBH | 序列层同源候选 | 当前 GEM 可执行性或功能完全等价 |
| 模型 GPR | KO 可以在当前模型中执行 | 实验一定可行 |
| OE reaction proxy | 关联反应容量变化方向 | 启动子、拷贝数和 gene-level 表达量 |
| 相对 gene-capacity 场景 | 在明确 mapping、剂量和不确定性下比较相对方向与资源代价 | 审核 baseline capacity、真实表达倍数或绝对产量 |
| 绝对 gene-level capacity | 使用经审核的 baseline capacity 执行绝对容量场景 | 当前没有合格锚点，因此保持 unavailable |
| curated phenotype | 特定 intervention/context 的已有表型方向 | 通用 mg/L 或跨条件保证 |
| 实验反馈 | 当前宿主、目标和条件下的观测结果 | 未测试条件的自动外推 |

`experiment_calibrated` 只表示存在严格匹配的高置信表型或内部实验支持，不表示能够预测绝对产量。

## 产品验收分层

### MVP 决策层

- 使用现有 KO、OE reaction proxy、相对 gene-capacity、证据分层和风险拦截生成候选优先级。
- 每条输出必须声明执行模式、校准状态、证据来源、不确定性和不可用原因。
- 该层可以支持研发候选选择，但不能声称绝对 capacity、mg/L、真实表达倍数或实验成功概率。
- 2026-07-20 在本层内深化四项免绝对数据的相对信号（见 [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）：影子价格瓶颈归因（深化主路径已有的 `analyze_target_protein_lp_attribution`，保留行级/复合体级与 `bound_type`，下界约束不报成 OE 瓶颈）、OE 剂量响应形状、排序对容量假设的稳健性标注（`ranking-insensitive-to-capacity`/`ranking-sensitive-to-capacity`，稳健性覆盖参数带宽+求解算法）、价值-of-information 实验优先级。四项均不产生绝对容量；稳健性的扫描带宽只用于稳健性分类，绝不断言为容量数值，绝对层保持 unavailable。

### 绝对容量研究层

- 继续沿用 ADR-001 的候选、审核和 promotion 门禁。
- 没有匹配的审核资产时返回 unavailable，不允许静默回退到 proxy、最优 flux、通用上界、固定 `1.0` 或 fixture。
- 新来源只有在开发前即具备明确单位、条件、版本、hash、license 和转换链时才允许启动接入。

### 实验校准层

- 实验反馈只校准候选排序、方向一致性和风险判断，不直接修改代谢矩阵、GPR 或正式容量资产。
- 没有真实实验数据时，软件阶段仍可通过脱敏 fixture、数值回归和边界测试验收；真实数据到来后再执行独立回填和校准。

分层决策见 [ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md)。

## 当前主要缺口

  1. 绝对 gene-level OE capacity 缺少可审核 baseline capacity，当前必须保持 unavailable。
  2. 分泌资源架构（`pcsec_pichia.secretory_resources`）已冻结全部七类的可执行契约，但完整机制求解、组合约束仍不完整；已核实转运/二硫键/糖基化/囊泡运输/folding-chaperone 五类有真实 kcat（前四类已在现有模型每次求解中生效），ER quality control/ERAD 有真实 kcat 但约束默认关闭（打开后对部分候选有 5%-14% 的真实排名影响），仅 target-specific 的目标蛋白降解速率（kdeg）真正缺数据。
  3. 目标蛋白自身的降解（区别于上一条的分泌通路 ERAD 机制）目前没有任何可执行路径：降解反应没有 GPR，任何基因敲除都无法影响它；已知有湿实验团队在做的液泡蛋白酶敲除（PEP4/PRB1/YPS1-3）里，PEP4/PRB1 在现有基因目录里的模型基因 ID 还被标注为低置信度待复核，YPS1-3 尚未进入基因目录。
  4. 大量模型外蛋白有序列和注释，但没有可执行 GPR，不能直接进入 KO/OE 求解。
  5. 筛查以单基因为主，缺少组合改造和上位性分析。
  6. 缺少跨碳源、生长率、氧供和参数扰动的完整稳定性排名；已用两次小范围检查（OE 候选跨条件排名稳定性、ERAD 约束开关敏感性）验证这类局部检查本身可行，尚未做成横向门禁。
  7. **目标函数是分泌通量/资源，不覆盖产物质量（糖链结构）与工艺（温度/pH/补料）**：糖基化人源化改造（改宿主糖基化通路）的目标是糖链**结构**，模型只数糖基化通量与供体消耗、判不了结构，相关敲除的细胞壁/生长风险也是结构性的、模型看不见；温度/pH/补料在放大培养带来约一个数量级增益但 FBA 看不到。已由 2026-07-13 hLF 周报佐证（**保密：具体基因/菌株/产量不入库**）。**同一份周报也验证了模型范围内的两条相对判断**——在已增强折叠/UPR 的本底上进一步增强折叠收益已饱和、某个运输因子过表达有毒——说明模型是“折叠/分泌容量这一片”的候选生成器，真实前沿（降解/调控/剂量/工艺/糖型质量）在其范围外。将来的 RNA-seq 拟用于表达约束的菌株特异建模（用表达量约束酶容量、把相对层贴向真实菌株），不解锁降解层建模。

研发发酵宽表适配和脱敏回放已经验收；尚未完成获批真实数据回填，但这不再是软件阶段阻塞。目前有一份周报级别的进度汇报，但数据量不构成这里说的"真实数据"。

## 架构边界

```text
app/ui
  -> app/services
    -> python_pichia
       -> loading / media / targets
       -> constraints / simulation / shadow_lp
       -> screens / analysis / reports
       -> services / homology / external_refs / experimental_feedback
       -> oe_capacity
```

- `python_pichia/`：核心科学逻辑、数据契约、约束、求解、筛查、证据和报告事实层。
- `app/services/`：facade、任务触发、缓存路径和错误汇总，不实现科学判断。
- `app/ui/`：Streamlit 展示和用户操作，不直接修改模型语义。
- `Code/`、`Model/`、`Enzymedata/`、`Results/`：legacy/reference 科学资产，默认只读；`Results/` 保留为 legacy MATLAB results 的只读参考，不是当前 Python/Streamlit 的输出目录。
- `local_runs/`：运行产物、缓存、报告和验证证据，默认 ignored。

### Phase 2 所有权边界

- `oe_capacity` 负责 gene-enzyme-reaction 映射、OE 剂量、参数不确定性、约束计划和 proxy 对照结果。
- `screens/gene_interventions.py` 继续拥有现有 GPR 解析与兼容 proxy 规划；Phase 2 调用它，不复制另一套 GPR 解释器。
- `core/pichia_enzymes.py` 和当前模型中的酶变量/稀释反应是本地执行语义的首要来源。
- `external_refs` 只提供外部 GEM、GPR、kcat、丰度和同源映射候选；未经当前模型 reaction/gene/enzyme 复核不得成为可执行约束。
- `shadow_lp`、SciPy HiGHS 和现有 reference solve 仍是求解路径；Phase 2 不引入新的默认 solver。
- gene-level capacity 与旧 reaction proxy 必须作为两个明确模式并行存在，不能通过改名掩盖降级。

### 外部容量候选与正式资产边界

研发组无法提供内部 baseline capacity 时，系统允许从外部来源建立候选，但外部候选不能直接成为可执行约束。来源优先级为：

1. 同宿主、同菌株且培养条件接近的 Pichia 定量蛋白组。
2. 经版本和许可审计的 iPichia/ecPichia 等酶约束模型。
3. Pichia 文献以及 BRENDA、SABIO-RK 等动力学来源与蛋白组的组合换算。
4. 通过 BLAST/RBH 确认的 S. cerevisiae 同源参数转移，仅作为低置信区间。

外部参数分为四种适用范围：

- `target_specific`：同宿主、同条件且针对当前目标蛋白获得的参数。
- `host_condition`：宿主基础容量，可在宿主、菌株、培养基、碳源和生长状态匹配时供 hLF/OPN 复用。
- `external_model_calibrated`：来自外部酶约束模型的校准参数，不等同于实测丰度。
- `homolog_transferred`：由其他物种同源转移，只能生成宽不确定性区间。

匹配优先级固定为 `target_specific > host_condition > external_model_calibrated > homolog_transferred`。候选必须完成单位换算、当前模型 gene/enzyme/formation 映射、条件匹配、来源版本/hash/license、冲突检查和人工审核后，才可提升到 `Enzymedata/oe_capacity_baseline_capacity.json`。`1000` 通用上界、baseline optimal flux、fixture 和未审核同源参数永远不能作为正式锚点。

正式容量 promotion 见 [ADR-001](adr/001-external-capacity-candidate-promotion.md)；相对决策层与绝对容量层的关系见 [ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md)。

## 数据与产物治理

本节吸收原独立数据治理策略的关键边界，作为需求与架构文档的一部分；完整历史细则归档于 `docs/archive/data_and_results_policy.md`。

### 目录职责

- `Data/`、`Model/`、`Enzymedata/`：人工确认、可追溯、可长期维护的稳定科学资产（curated 输入、GEM/pcSec 模型、酶容量与资源约束）。
- `Results/` 是 legacy MATLAB results，只读参考，不是当前 Python 或 Streamlit 的默认输出目录。
- `local_runs/` 是当前 Python、Streamlit、MATLAB harness 运行产物、临时缓存、外部下载、smoke、报告和待复核实验导入的统一落地目录，默认 ignored。
- `docs/archive/`：已完成计划、阶段验证和被 active 文档吸收的历史设计，默认不进入公开版本控制。

历史 `Results/` 迁移、Git LFS 改造或仓库历史瘦身不属于当前数据治理范围。

### 外部数据与容量参数

- BLAST/RBH、在线数据库响应、外部 GEM、MEMOTE 报告和 GPR mapping 默认进入 `local_runs/`，必须保留 source URL、version、retrieved_at、hash、license/provenance 和 warning。
- 外部名称、同源关系或 GPR 不得自动覆盖内部 `gene_id` 和当前模型规则。
- 正式 gene-level capacity 只能使用 `Enzymedata/oe_capacity_baseline_capacity.json` 中经审核的绝对 formation/dilution capacity；空资产或缺匹配记录时返回 `reviewed_baseline_capacity` 缺口并降级，不得从最优 flux、通用上界、fixture 或固定 `1.0` 推断。
- 缺失参数使用显式 low/nominal/high 区间；联网抓取与正式 screen 求解分离，求解只读取冻结的本地资产快照。相对 OE 与绝对容量分层按 ADR-002 执行，分层不降低正式容量门禁。

### 实验反馈数据

原始导入进入 `local_runs/experiment_feedback/inbox/`；标准化缓存进入 `.../validated/`；只有人工确认来源、权限、脱敏、单位和上下文后，才单独 checkpoint 提升到 `Data/` 正式目录。失败、无效、不可测和阴性结果必须保留，不能只收集成功实验。

### LLM 数据边界

LLM 只读取程序生成的 fact pack，不遍历运行目录；只有用户明确触发才调用外部 LLM API；LLM 不能修改模型、recommendation tier、curated evidence 或实验记录；最终报告必须同时通过程序 validator 和 Judge。

### 密钥与提交

`.env`、API key、认证文件和个人路径不得提交；日志与 manifest 不得打印密钥。新生成的 LP、solver output、BLAST cache、外部下载、Streamlit run、LLM report 和 smoke 默认进入 `local_runs/`。大文件进入 Git 前必须说明来源、可再生成性、license 和人工复核状态；仓库当前未启用 Git LFS，在没有独立评审前不得用它绕开这条规则直接提交大型二进制资产。`Data/Model/Enzymedata/Results` 的任何修改都必须明确声明为科学资产变更；每个涉及模型、筛查、证据、实验反馈或报告的任务结束前检查 `Code/Model/Enzymedata/Results` 与依赖声明无非预期 diff。

## 慢速测试网关

三类真实求解回归默认跳过，需要显式设置对应环境变量才会运行：

- `PCSEC_RUN_SLOW_PIPELINE_TESTS="1"`：单候选全流程求解（`run_pichia_secretion_simulation`），覆盖 `python_pichia/tests/test_pipeline_entrypoints.py`。
- `PCSEC_RUN_SLOW_SCREEN_TESTS="1"`：全模型 KO/OE 批量筛选，覆盖 `python_pichia/tests/test_screens_entrypoints.py`。
- `PCSEC_RUN_SLOW_PROBE_TESTS="1"`：probe 迁移回归，对照 MATLAB harness / baseline 生成的既有产物，覆盖 `python_pichia/tests/test_probe_migration.py`。

## 项目成功标准

项目是否更接近原始目标，不以新增多少字段或页面衡量，而以这些结果衡量：

- top-K 候选相对随机或人工经验基线是否有更高实验命中率。
- hLF 与 OPN 的候选排序是否能在重复实验中复现。
- 模型预测方向、实验方向和不确定性是否得到校准。
- 生长风险、不可执行候选和 proxy 结论是否被正确拦截。
- 每条实验建议是否能追溯到模型结果和证据来源。

当前范围以[项目级执行计划](EXECUTION_PLAN.md)为准，实际下一切片只从[当前 handoff](handoff.md)继续。
