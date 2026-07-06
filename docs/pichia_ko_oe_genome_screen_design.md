# Pichia 全基因组 KO/OE 分泌-生长权衡筛查设计

状态：active
最后更新：2026-07-02

## 目标

对 pcSecPichia 模型里全部约1025个基因，逐一评估敲除(KO)和过表达(OE)对目标蛋白分泌产量的影响，**同时**给出对细胞生长速率的影响——因为部分候选基因虽然能提升分泌产量，但会严重损害生长，这类候选需要额外的生物学补救手段才具备实验可行性，必须能被识别出来，而不是被产量这一个指标掩盖。

触发场景：研发同事在对某个基因做湿实验敲除，而模型侧对该基因目前只有策展目录里的 OE 数据，由此发现现有 KO/OE 覆盖存在系统性缺口。

## 现状核实

### 策展基因目录严重偏向 OE，是文献策展偏科，不是模型能力限制

`SECRETION_GENE_CATALOG`（`python_pichia/src/pcsec_pichia/services/gene_catalog.py`）共 33 条条目，其中仅 PEP4、PRB1 两条标注 `intervention="KO"`，其余全部是 `"OE"`。且多数 OE 条目只填了 `oe_reaction_id`（指向复合体形成反应，如 `sec_Arf1p_Sec3p_..._complex_formation`），没有独立 `gene_id`——一个反应经常对应多个基因（如 exocyst 复合体对应六个基因名）。

条目的 `evidence` 字段显示，这些基因是从"过表达提升分泌"的毕赤酵母分泌工程文献里摘录的，KO 类文献本来就少（PEP4/PRB1 这种降解蛋白酶敲除是例外）。不代表这些基因结构上不能做 KO，只是没人去查。

### 全模型路径已经能独立算任意基因，不依赖策展名单

`load_full_model_genes()` / `build_gene_capability_profile()` / `build_all_gene_capability_catalog()`（`screens/gene_interventions.py`）可以对模型里任意基因独立计算 KO/OE 支持状态，纯粹基于 GPR 规则结构和反应 ID 分类（`classify_secretory_process`），不需要跑 LP 求解，计算成本很低。

### 现成的批量求解入口已经是通用的，且自带预筛

`run_knockout_screen(genes: list[str], ...)` / `run_overexpression_screen(reactions: list[str], ...)`（`screens/__init__.py`）接受任意基因/反应列表，不写死策展名单，且已经有预筛优化：只对 `plan_gene_knockout` 判断会让某些反应失活的基因，才真正去求解 LP。

### 关键缺口：生长速率现在是锁死的输入，不是求解结果

`_prepare_screen_inputs` 里：

```python
fixed_model = build.model.with_bounds({"BIOMASS": (growth_rate, growth_rate)})
baseline, counts = solve_pcsec_maximize(fixed_model, exchange_reaction_id, ..., mu=growth_rate, ...)
```

生长速率被硬性锁定在同一个值（上下界相同），LP 只求解"在这个生长速率下，分泌产量最多能到多少"。这意味着"敲除后细胞实际能长多快"目前完全没有被计算出来，只能间接看这次求解是成功还是失败（在锁定的生长速率下是否可行）。

## 已排除的方向及原因

| 方向 | 结论 | 原因 |
| --- | --- | --- |
| 用 OE 结果反推 KO 结果（"分泌倒推敲除"） | 不采用 | OE 效果取决于反应是否原本就是限速步骤（不限速则 OE 基本无效，甚至因抢占蛋白质组预算而降低产量）；KO 效果取决于整条支路是否被切断。两者机制不同，对同一基因可能给出完全无关的答案，没有理论支撑做因果外推 |
| 引入外部数据库（UniProt/KEGG）扩大基因覆盖 | 暂不作为主要方向 | 已实测：`local_runs/gene_rule_evidence_cache/GENE_RULE_EVIDENCE_REPORT.md` 查了 PDI1/ERO1/ERV2/OCH1 四个基因，3 条记录 0 条达到高置信度、0 条变成可执行 overlay。常用名匹配不可靠，只有 locus tag 精确匹配才可信，产出率太低 |
| 跨物种（pcSecYeast/pcSecKmarx）kcat 借用扩大 OE 覆盖 | 暂缓 | 科学合理性待确认，需要研发组长评估三个物种分泌通路差异是否允许这样迁移 |

## 核心设计决策

### 决策1：范围为全部约1025个基因，不只筛分泌相关子集

**原因**：生长-产量权衡是核心关注点，而容易被漏掉的恰恰是"看起来和分泌无关，但敲除后连带影响生长"的基因（比如中心代谢基因）。只筛分泌相关子集会系统性漏掉这类候选。

### 决策2：两阶段筛查，只对有结构效应的基因做 LP 求解

第一阶段：对全部基因跑 `build_all_gene_capability_catalog`，零 LP 求解成本，拿到每个基因的 KO/OE 结构可行性 + 涉及的通路/反应分类。

第二阶段：只对第一阶段显示"KO 会让某些反应失活"或"OE 有可执行容量效应"的基因，才进入真正的 LP 求解——复用 `run_knockout_screen` 已有的预筛逻辑。

**取舍**：完全不影响任何反应的基因不会被求解，只会标记为"无结构效应"。这些基因即使求解了也会是零效应，跳过求解不会丢失有意义的信息。

**实测更新（2026-07-02）**：真跑了一遍全部1025个基因的第一阶段，耗时95秒。结果是**KO 这边预筛完全没起到收窄作用**——1025个基因全部显示"会让至少一个反应失活"（`ko_runnable_gpr_gene_deletion`），说明这个模型里同工酶冗余导致KO无效的情况很少；OE 这边收窄了14.2%（879个可执行、146个是复合体亚基只能解释不能算）。也就是说需要真正求解的候选是 1025(KO) + 879(OE) ≈ 1904 组，比最初设想的"筛掉大部分"要多得多，直接影响了第二阶段的实际耗时估算（见决策5）。

### 决策3：每个基因输出（生长影响，产量影响）两个数字，而不只是产量

对通过第二阶段预筛的每个基因：

1. 求解"敲除该基因后，最大可达生长速率"（目标函数 = BIOMASS，不锁定生长速率）
2. 在这个可达生长速率下（或原参考生长速率，取两者中可行的一个），求解目标蛋白的最大分泌产量

**原因**：只看产量会让"产量升高但生长严重受损"的候选和真正双赢的候选混在一起，无法区分哪些需要额外的生物学补救手段才具备实验可行性。

### 决策4：Fast / Precise 双模式 = 复用现成的 `run_pcsec_growth_tradeoff` 扫描一组固定 mu，而不是自己写迭代收敛

最初设想是"猜一个 mu → 求解 → 用结果更新 mu → 重新搭约束 → 再求解"这样的不动点迭代，因为蛋白质组预算约束依赖假设的生长速率、而"实际能达到的生长速率"恰恰是要求解的目标，两者存在循环依赖。

但排查代码时发现 `run_pcsec_growth_tradeoff`（`probe/_prototype.py`，已有 2 处调用）已经用另一种方式绕开了这个循环依赖：不去反解"最大可达生长速率"，而是**对一组给定的固定 mu 值逐一求解**（每个 mu 都和现有 `_prepare_screen_inputs` 一样，把 BIOMASS 上下界锁死在该值），返回每个 mu 点的可行性和该 mu 下的最大分泌产量（含 `secretion_per_biomass` 归一化产率）。这样完全不需要迭代收敛，而且天然给出的是一条完整的"生长速率-产量"关系曲线，而不是单个点，信息量更大。

- **Fast 模式**：`mu_points` 取一组稀疏的值（比如野生型参考生长速率的 100%/75%/50%/25%/10%），每个基因求解次数少。
- **Precise 模式**：`mu_points` 取更密的一组值（比如 20 个均匀分布的点），曲线分辨率更高，但求解次数更多。

两个模式复用同一个 `run_pcsec_growth_tradeoff`，区别只是传入的 `mu_points` 列表疏密不同；具体取点方式做成可配置参数，不写死。判断"敲除后是否严重影响生长"，看这条曲线在多高的 mu 处开始不可行（`success=False`）即可，不需要额外的收敛逻辑。

### 决策5：单次求解约9秒，必须做多进程并行才能在一晚内跑完全量

实测（HiGHS 求解 23024×29068 规模的LP）：单次 `solve_pcsec_maximize` 约9秒（搭约束矩阵1.1秒 + `linprog` 求解7.5秒，瓶颈在求解本身，不是搭约束，所以没有"复用约束矩阵"这种优化空间）。

按此推算，串行跑全部1025基因(KO)+879基因(OE)、Fast模式(3个mu点)、2个靶点(hLF+OPN_ALPHA_FULL_PROJECT)，大约要**48小时**，完全不可能在一晚内跑完。这是个天然可并行的问题（每个基因的求解互相独立），机器有8个逻辑核心，用 `concurrent.futures.ProcessPoolExecutor` 开6个worker并行，预计能把总耗时压到8小时左右，能在"下班到第二天上班"这个窗口内跑完。

**实现方式**：每个worker进程启动时独立加载一次模型（避免跨进程传递不可pickle的模型对象），并在worker内部按 `(target_id, mode)` 缓存"目标蛋白特异性的模型准备工作"（`build_supported_target_model`/`build_target_enzymedata`/野生型基线），避免同一个worker处理同一个靶点的多个基因时重复计算这部分。

## 实现计划（已完成）

- [python_pichia/src/pcsec_pichia/screens/genome_wide_tradeoff.py](../python_pichia/src/pcsec_pichia/screens/genome_wide_tradeoff.py)：核心逻辑，`gene_ko_tradeoff`/`gene_oe_tradeoff`/`run_genome_wide_tradeoff_screen`，复用 `run_pcsec_growth_tradeoff`/`run_pcsec_oe_screen` 的现成扰动机制
- [python_pichia/tools/run_genome_wide_ko_oe_screen.py](../python_pichia/tools/run_genome_wide_ko_oe_screen.py)：串行版命令行入口，适合小规模试跑
- [python_pichia/tools/run_genome_wide_ko_oe_screen_parallel.py](../python_pichia/tools/run_genome_wide_ko_oe_screen_parallel.py)：并行版命令行入口，全量跑用这个
- 两者都输出到 `local_runs/<run_name>/`：`gene_tradeoff_rows.csv`（每行一个基因+扰动类型，含生长影响、产量影响、涉及通路/反应、置信度）+ `SUMMARY.md`
- OE 侧的复合体级条目（有多个基因共享一个 `oe_reaction_id`）本轮按复合体整体处理，不拆解到单基因粒度（拆解是另一个独立问题，见下方待定事项）

## 验证发现：PEP4 的策展 gene_id 和模型实际 GPR 对不上

试跑验证阶段，特意测了策展名单里标注"敲除 PEP4 已知能提升分泌"的基因（`gene_id="PAS_chr2-2_0107"`），结果这个 ID 在模型里实际关联的是脂肪酸酯酶反应（`FACOAE140/160/180_no_1_fwd`），不是预期的液泡蛋白酶反应。这是一个具体的、可复现的"命名和对应关系有问题"的实例——印证了本文档"已排除的方向"里对外部数据库命名匹配的顾虑不是空穴来风，策展数据里的 `gene_id` 字段也不能直接当作可信来源，本轮采用"全模型路径独立算"而不是照抄策展 `gene_id` 的做法被这个发现进一步证实是对的。

## 产品化：Streamlit 页面 + LLM 报告（2026-07-03）

把上面这套筛查从脚本固化成了产品功能，任何人打开 Streamlit 都能直接用，不需要再找人临时跑一遍分析。关键设计决策：

- **后台任务用独立子进程，不用 Streamlit 内部线程**：现有的 `app/services/pichia_background_tasks.py` 用的是守护线程 + 10分钟无更新判定 stale，适合1-2分钟的仿真；这个筛查是小时级任务，Streamlit 重启/自动重载会直接杀掉线程、丢掉几个小时的进度，所以改用 `subprocess.Popen`（Windows 下用 `DETACHED_PROCESS` 标志）+ 文件心跳（`status.json`，30分钟无更新才判定 stale）。
- **并发处理**：`app/services/genome_wide_screen_registry.py` 维护运行记录；点击"启动"前先查是否有活跃任务，有的话弹出选择（取消/排队/强制并发），不是固定策略。排队的实现比较轻量——队列请求存在 session state 里，页面下次加载时如果没有活跃任务就自动启动，不是一个常驻的watcher进程，如果没人再打开这个页面，排队的任务不会自动触发（这是当前设计的已知局限）。
- **LLM 报告**：`app/services/llm_report_service.py` 是可插拔接口（`ReportGenerator` protocol），默认实现接 OpenAI（`OPENAI_API_KEY` 环境变量，需要用户自行配置，代码里不写死）。
- **5个分析维度**：`app/services/genome_wide_screen_analysis.py` 把之前临时脚本里的分析逻辑产品化——必需基因、产量升高但生长受损的KO候选、零代价KO候选、产量升高的OE候选、跨靶点差异。

新增/修改的文件：`app/services/genome_wide_screen_{registry,service,analysis}.py`、`app/services/llm_report_service.py`、`app/ui/views/genome_wide_screen.py`、`app/ui/common.py`（导航）、`app/ui/streamlit_app.py`（路由）、`python_pichia/tools/run_genome_wide_ko_oe_screen_parallel.py`（加了状态心跳写入）。用 Streamlit 自带的 `AppTest` 无浏览器验证过页面渲染无异常（浏览器自动化工具在这个环境里侧边栏渲染有问题，改用这个更可靠的方式验证）。

## 信号肽/引导肽对分泌效率的实证结论（2026-07-03）

- 问题：三段式构建是不是必要的？信号肽、引导肽本身会不会影响分泌效率，需不需要为此加设置项？
- 方法：直接跑了7个OPN变体的正面对比（用户提供的完整 leader 构建 + 6种候选 signal peptide/leader 组合），固定 mu=0.10，看目标蛋白分泌通量的最优解。
- 结果：`objective_value` 和 leader/全长的氨基酸长度直接负相关（越长代价越高），没有观察到"某种信号肽本身更高效"这种效率差异——当前模型不区分信号肽身份，只按它贡献的氨基酸数/PTM负担计入蛋白质组预算成本，是一个纯长度/成本模型，不是效率模型。
- 结论：引入信号肽特异性的效率参数（给每个信号肽一个独立的效率乘数）在实现上可行，但目前没有任何实验数据支撑具体取值，贸然加参数等于编造数据，不做。三段式构建保留为可选功能（用户原有需求：支持自定义信号肽/引导肽/成熟蛋白组合），但价值定位是"低摩擦切换组合、看清长度代价"，不是"帮用户挑更高效的信号肽"——UI 上也不需要为此加"推荐信号肽"之类的智能推荐。
- 因此原"待定事项"里"三段式目标蛋白构建"这一条已确认为成熟功能且效率问题已有实证结论，不再是待定项。

## UI重组：筛查结果到仿真验证的核实跳转（2026-07-03）

- 问题：全基因组筛查如果最后只筛出几个候选，是不是就不需要再跑一次完整仿真去确认了？另外当前UI把筛查阶段结果和仿真验证历史界面杂糅在一起，需要重新组织。
- 决策：不是"要不要跑仿真"的取舍——少数候选恰恰应该逐个跑仿真确认，真正要解决的是"怎么让这一步零摩擦"，而不是要不要做这一步。做法：在3个候选表格（KO产量升高有生长代价 / KO零代价 / OE产量升高）里加"选中一行 + 在仿真验证中核实"的跳转，点击后自动把仿真验证页面的靶点、构建模式、KO/OE基因输入填好并跳转过去，不用手动抄基因ID再重新选靶点。
- 实现：`app/ui/common.py` 的侧边栏导航 radio 加了显式 `key`（`NAV_RADIO_KEY` = `"app_page_nav"`），使其可以被其它页面编程式设置来跨页面跳转；`genome_wide_screen.py` 新增 `_render_verifiable_table`（用 `st.dataframe(on_select="rerun", selection_mode="single-row")` 做行选中）和 `_apply_verify_prefill`（选中后点按钮，写入仿真验证页面既有的 session_state key —— `pichia_tab_selector`/`pichia_draft_build_mode`/`pichia_template`/`pichia_draft_ko_genes`/`pichia_draft_oe_genes` 等 —— 再 `st.rerun()`）。
- 预填是替换而不是追加：点击核实会清空另一侧（KO/OE）的草稿字段，保证仿真验证的是且仅是这一个候选，不会混入之前探索时残留的其它基因输入。
- 侧边栏"推荐演示顺序"文案相应调整：全基因组筛查（探索）排在仿真验证（核实）之前，反映新的使用流程，并注明可以从筛查候选行直接跳转。
- 验证方式：浏览器自动化工具在这个环境里对本项目侧边栏渲染有已知问题（见产品化章节），继续用 Streamlit 的 `AppTest` 分两层验证——(1) 单独调用 `_apply_verify_prefill` 检查各 session_state key 写入的值是否正确；(2) 用完整 `streamlit_app.py` 预置这些 session_state 值（模拟"刚从筛查页跳转过来"的状态），确认仿真验证页面的 radio/selectbox/text_area 控件真的读到了预填值。另外对 `overnight_hLF_full`（1025基因真实产出）做了一遍页面渲染回归，无异常。

## 策展反应级 KO/OE 筛查：第二个分析入口 + 生长代价结论（2026-07-03）

- 问题：PDI1/ERO1/ERV2（hLF）、PMT1/2/4-6（OPN）这两个反应级 OE 候选有没有生长代价？另外，反应级 OE 测试要不要产品化成页面里独立于全基因组扫描之外的第二个分析入口？
- 决策：两个问题合并成一个功能解决——策展名单里的候选本来就是有限的一批（`SECRETION_GENE_CATALOG` 里去重后共 32 个唯一的 `(intervention_type, reaction_id)` 组合，多个基因常共享同一个复合体反应，比如 PDI1/ERO1/ERV2 三个条目都指向 `sec_PDI1_ERV2_Ero1p_complex_formation`），跑一遍这些候选的 mu 扫描既能拿到生长代价数据，又正好就是"策展反应级筛查"这个新入口本身，不需要分两步做。
- 实现：新增反应级版本的权衡函数 `reaction_ko_tradeoff`/`reaction_oe_tradeoff`（`genome_wide_tradeoff.py`）——和基因级的 `gene_ko_tradeoff`/`gene_oe_tradeoff` 区别在于跳过"基因→GPR→反应"解析这一步，直接对 `SecretionGeneEntry.oe_reaction_id`/`ko_reaction_id` 里的反应 ID 设边界/调用 `run_pcsec_oe_screen`，因为这些策展条目大多数本来就没有可解析的单一 `gene_id`（是复合体级 MATLAB 伪反应）。`catalog_reaction_candidates()` 负责去重。CLI 脚本 (`run_genome_wide_ko_oe_screen_parallel.py`) 加了 `--scope {gene,catalog}`，复用同一套 worker/心跳/CSV 基础设施，只是 catalog 模式下任务列表换成这 32 个候选而不是全部 1025 个基因。输出 CSV 新增 `common_name`/`candidate_kind` 两列（旧的纯基因级 CSV 没有这两列，`load_gene_tradeoff_csv` 会自动补默认值，不需要迁移旧数据）。UI 层加了"筛查范围"单选（全基因组 / 策展复合体反应对照表），运行记录表加了"范围"列，结果查看的 CSV 路径改成优先读 `status.json` 里的 `csv_path`（catalog 和 gene 两种范围输出文件名不同）。"在仿真验证中核实"按钮现在按 `candidate_kind` 分流：`gene` 填基因输入框（走 GPR 解析），`catalog_reaction` 填反应输入框（直接用反应 ID，不需要 GPR）。
- **实测结果（`local_runs/catalog_reaction_screen_hlf_opn`，32候选×2靶点=64任务，6 worker 并行实测 10.9 分钟）**：PDI1/ERO1/ERV2 对 hLF 分泌提升 8.15%（`secretion_ratio_vs_wildtype=1.0815`），PMT1/PMT2/PMT4-6 对 OPN 提升 16.96%（`1.1696`）——两者在整个 mu 扫描范围内 `growth_retention_ratio` 都恰好是 1.0，**零生长代价，是干净的win，不需要额外的生物学补救策略**。附带发现：hLF 上 OCH1/CWH41/BIP_NEFS 的 OE 也有 1-2% 的提升（同样零代价，量级较小）；`sec_Och1p_complex_formation` 的 KO 在两个靶点上都直接不可行（`max_feasible_mu=NaN`），提示这个复合体反应被建模成整体开关而不是可调节强度，一旦关闭对模型来说太剧烈；HRD1/CDC48/DOA10 等 ERAD 相关 KO 在模型里的分泌提升趋近于零（≈0.9999），和策展目录里"敲除 HRD1 可减少 ERAD、提升外源蛋白积累"的文献描述不一致——这是继 PEP4 之后又一处策展证据和模型实际行为对不上的案例，同样不采信策展描述，以模型直接求解结果为准。

## 仿真验证（单次运行）页面变慢排查 + 并行化/缓存修复（2026-07-03）

- 问题：用户反馈"仿真验证"页面单次运行（不是全基因组筛查）变得非常慢。要求先用 codebase-memory-mcp 通读架构和数据流、找出主要问题，再动代码。
- 数据流：按钮 → `pichia_background_tasks._run_background_task`（后台线程）→ `run_pichia_secretion_draft` → `pcsec_pichia.pipeline.run_pichia_secretion_simulation`，这个函数**串行**依次跑：基线求解（1次LP）→ `run_growth_tradeoff`（默认1个mu点，1次LP）→（可选）`run_protein_cost_slope_compatibility`（默认 2 mu × 5 ratio = 10次LP，嵌套循环）→ KO基因/KO反应/OE基因代理/OE反应 4 个screen（各自0~`screen_candidate_limit`默认20个候选，各1次LP）。
- 找到的4个主要问题（按影响力排序）：
  1. 全程零并行——每个候选的LP求解都是单线程逐个跑（对比全基因组筛查已经用了`ProcessPoolExecutor`）。
  2. `screen_candidate_limit`默认20、横跨4个输入框，容易不知不觉叠加到 4×20=80 次求解。
  3. "蛋白成本斜率对比"（`enable_cost_slope_compatibility`）是隐藏的10次求解乘数——2个生长速率×5个分泌比例的嵌套循环，是整个`simulation`模块里`transitive_loop_depth`最高（6）的函数；默认关闭、UI标注"较慢"。
  4. `load_pcsec_pichia_inputs()`（模型加载）完全没有跨请求缓存，每次点击"运行"都重新从磁盘解析`.mat`模型。
- 已实施的修复（用户选择先做 1 和 4）：
  - **并行化**：新增 `pcsec_pichia/core/concurrency.py` 的 `parallel_map()`（`ThreadPoolExecutor`），应用到 `run_pcsec_oe_screen`/`run_pcsec_reaction_ko_screen`（`probe/_prototype.py`）、`run_knockout_screen`的基因求解循环（`screens/__init__.py`）、`run_growth_tradeoff`/`run_protein_cost_slope_compatibility`（`simulation/__init__.py`）。选择线程而不是进程：HiGHS(`scipy.optimize.linprog(method="highs")`)的原生求解会释放GIL，且在已经运行的 Streamlit 服务进程内嵌 `ProcessPoolExecutor` 有 spawn 语义风险（Windows 上 spawn 会重新 import 当时的 `__main__`，在 `streamlit run` 场景下这个模块是什么并不确定）——这也是为什么全基因组筛查的并行化是走独立子进程脚本（`run_genome_wide_ko_oe_screen_parallel.py`），而不是在 Streamlit 进程内部直接嵌入进程池。
  - **实测速度提升比预期小**：6基因KO screen 串行60.1s → 并行(6线程)43.2s，约1.4倍，不是理想情况下的6倍。排查发现单次求解里约86%时间在`linprog`本身（9.0s vs 构建约束矩阵1.48s），理论上应该有更好的并行度；尝试用HiGHS的`threads`选项限制单次求解内部线程数以避免和外层线程池抢核，但这个选项在scipy里是"unrecognized option, passed through verbatim"，行为不可靠（`threads=2`直接导致求解静默失败，status=4）——判定这条路不安全，放弃，保留现有约1.4倍的线程池方案。后续如果要拿到更大提升，需要复用全基因组筛查那套独立子进程架构，是一个更大的、有独立风险的改造，本轮不做。
  - **模型缓存**：`loading/__init__.py`新增`_cached_base_pcsec_pichia_artifacts()`（`functools.lru_cache`），缓存`load_pcsec_pichia_model`/`load_aa_stoichiometry`/`load_metabolic_enzymedata`/`load_secretory_enzymedata`/`load_combined_enzymedata`这5个只依赖repo root、不依赖media_type/carbon_source_id的加载步骤（这些返回值全部是不可变对象，`with_*()`系列方法都是copy-then-return-new，多线程/多次请求共享同一份缓存实例是安全的）。**实测：冷启动6.48s → 缓存命中0.10s，约65倍**，且正确处理了media_type/carbon_source_id变化的情况（只有这两个依赖它们的轻量步骤会重新计算）。
  - 端到端验证：3基因KO候选的完整单次仿真，第一次60.2s（冷模型缓存），第二次50.7s（同进程内热缓存），两次都成功；全基因组筛查（复用了改动过的`run_pcsec_oe_screen`）跑了5基因烟雾测试确认没有破坏原有路径；`python_pichia/tests`全量150个测试通过（15 skip，都是环境限制导致的既有skip）；顶层`tests/`168个测试通过（另修复了1个无关的既有失败——`pichia_ko_oe_genome_screen_design.md`本身当初创建时漏加进了`test_docs_active_boundary.py`的`ACTIVE_DOCS`白名单，这轮顺手补上）。

## "目标蛋白成本分析"/"生长分析"定性还是定量？和MATLAB对比（2026-07-03）

- 问题：仿真结果里的"目标蛋白成本分析"和"生长分析"是不是没有真正计算，只是定性展示？如果是，MATLAB项目里是怎么做的？
- 结论：**两者情况不同**，且代码自己在警告文案里已经承认了这一点：
  - **"蛋白成本分析"（`analyze_target_protein_cost`，`analysis/__init__.py`）确认是纯定性的**——`_cost_items_from_target_and_plan()`里的"cost"是一组手工设定的经验权重（比如DSB数×18、OG数×12+min(ser_thr_count,og)、ER转运步骤数×20……），归一化成百分比"relative_score"用来排出"主要成本类别"。函数自带的警告原文："当前成本分析是 Python draft explanatory score，不代表真实发酵产量或湿实验成本" / "该分析不使用 LP shadow price，也不改变 corrected pipeline 的求解目标或约束"——不使用任何LP求解结果，纯粹是序列长度/PTM数量/路线步骤数的启发式打分。
  - **"生长分析"（`analyze_target_growth_impact`）是真正基于LP求解结果的**——用的是`run_growth_tradeoff`实际解出来的`secretion_flux`。但因为`growth_points`默认只有1个点（`(0.10,)`），`_growth_sensitivity()`在点数<2时会诚实地返回`"insufficient_points"`，不会编造趋势——默认配置下这个分析是"真实但退化"（数据点不够，说不出趋势），不是假的。
  - **MATLAB的真实做法**：读了`Code/pcSecPichia/Simulation/SimulateProteinCost.m`，做法是固定生长速率（`mulist=[0.05, 0.1]`）和目标蛋白分泌比例（`ratios=[5E-7,1E-6,5E-6,1E-5,2E-5]`），对每个组合最大化葡萄糖摄取（`Ex_glc_D`），从葡萄糖/核糖体通量随分泌比例变化的斜率里得到真实的"蛋白成本"。
  - **这个MATLAB方法已经被完整移植成Python**——就是`run_protein_cost_slope_compatibility()` + `_cost_slopes()`/`_linear_slope()`（`simulation/__init__.py`），默认参数`growth_rates=(0.05, 0.10)`、`secretion_ratios`空则退化为按当前分泌能力的5个比例——和MATLAB的`mulist`/`ratios`逐一对应，是真实的最小二乘斜率拟合，不是伪造的。但它被做成了一个默认关闭的"MATLAB兼容模式"选项（`enable_cost_slope_compatibility`），原因就是上面提到的"隐藏的10次求解乘数"性能问题——真实方法是有的，只是因为慢被藏起来了，默认展示给用户的反而是这个不用LP的定性替代品。
- 待定：现在KO/OE screen已经并行化，`run_protein_cost_slope_compatibility`的10次求解理论上会明显变快（虽然线程池只有约1.4倍提升，10次求解从~90s降到~65s左右，量级上还是不算快）——要不要把这个真实的蛋白成本斜率分析提升为默认展示（而不是默认关闭的兼容模式选项），这是一个会改变默认UI展示内容的产品决策，本轮没有做，需要单独确认。**（已在下一节实施）**

## 删除定性蛋白成本分析，生长分析补足数据点（2026-07-03）

- 决策（用户直接拍板）：删除`analyze_target_protein_cost`这个纯定性打分功能，"蛋白成本分析"这个概念以后完全等价于`run_protein_cost_slope_compatibility`那套基于LP的真实斜率分析——不勾选"启用蛋白成本分析"选项时，不展示任何替代内容（不再有不用LP结果的简化版本）。同时解决生长分析默认只有1个点、`_growth_sensitivity()`永远返回`insufficient_points`的问题。
- 实现：
  - 删除`ProteinCostItem`/`ProteinCostAnalysisResult`/`analyze_target_protein_cost`/`summarize_protein_cost_analysis`/`build_cost_item_table`/`_target_feature_payload`/`_cost_items_from_target_and_plan`/`_normalise_cost_scores`/`_dominant_categories`/`_cost_warnings`（`analysis/__init__.py`），确认这些函数除了互相调用之外没有被其它功能依赖，删除前用`search_graph`/`grep`核实过。保留`analyze_target_protein_lp_attribution`不变——它是基于基线求解的LP shadow price，本来就是真实的，且几乎零成本（不需要额外求解），只是之前被错误地嵌套展示在定性的"蛋白成本分析"区块里。
  - `pipeline.py`的`run_pichia_secretion_simulation`不再无条件计算`protein_cost`；`protein_cost_analysis`这个字段现在只在`request.enable_cost_slope_compatibility=True`时才是一个dict（包含`lp_attribution`+`cost_slope_compatibility`），否则是`None`——不再有回退到定性打分的分支。`reports/__init__.py`和`pipeline.py`各自的markdown报告生成函数（`_protein_cost_markdown_lines`/`_protein_cost_report_lines`，两处历史遗留的重复实现）和`simulation_results.py`的`_render_protein_cost_analysis`都同步做了裁剪，只保留LP归因和成本斜率两块。
  - `simulation_builder.py`的选项文案从"启用蛋白成本斜率对比（MATLAB 历史路线，可选，较慢）"改成"启用蛋白成本分析（固定生长率+分泌比例网格测算成本斜率，较慢）"，help文案明确写"这是目标蛋白成本分析功能本身"、"不勾选时不会展示任何蛋白成本分析——没有不使用LP结果的简化替代版本，因为那类替代版本不代表真实成本，容易造成误导"。
  - 生长分析数据点问题：根因是`growth_points`在3层schema（`app/services/pichia_secretion_schema.py`的`SecretionRunRequest`、`engines/base.py`的`PichiaSimulationRequest`、`app/api/pichia_secretion_api.py`的FastAPI模型）里都硬编码默认成单点`(0.10,)`/`[0.10]`。把三处默认值都改成空`()`/`[]`，然后在`pipeline.py`的`_growth_points()`里补上新逻辑：空输入时不再退化成`(fallback_mu,)`单点，而是用一个新常量`DEFAULT_GROWTH_TRADEOFF_MU_FRACTIONS=(1.0, 0.5, 0.1)`（复用全基因组筛查`FAST_MU_FRACTIONS`一样的100%/50%/10%取值习惯）生成围绕用户实际选定mu的3点网格。副作用：这也顺便修了一个之前没人注意到的潜在不一致——旧代码的单点默认值是字面量`0.10`，如果用户把mu改成别的值，生长权衡网格根本不会跟着变，新逻辑改成了`fallback_mu`的相对比例，跟随用户实际选择的mu走。
- 验证：4个`@slow_pipeline`门控的真实LP求解测试（`PCSEC_RUN_SLOW_PIPELINE_TESTS=1`跑）第一次暴露了一处必须修的既有测试断言（`assert len(tradeoff_rows) == 1`，因为默认网格从1点变3点），修完后连同新增的"`enable_cost_slope_compatibility=True`端到端应该产出真实`protein_cost_analysis`"测试一起，24个pipeline测试全部通过；顶层169个测试、`python_pichia`148个测试（另16个环境限制导致的既有skip）全部通过；用AppTest确认新选项文案正确渲染、5个页面都无异常。
- 记录一个方法论教训：这几个`@slow_pipeline`测试因为默认被环境变量门控跳过，在这次改动之前很可能从未被真正跑起来验证过——如果不主动加`PCSEC_RUN_SLOW_PIPELINE_TESTS=1`跑一次，会误以为"144 passed, 15 skipped"就代表安全，但skip不等于verified。以后改动任何被`@slow_pipeline`/`@slow_test`之类标记跳过的代码路径，必须主动打开环境变量跑一次真实求解，不能只看默认的快速测试结果。

## 补齐策展反应KO覆盖 + "KO降低分泌"维度 + OE能否补救的排查（2026-07-03）

- 问题链条（用户观察到"绕了一圈回到起点"引出的）：策展库32个反应里只有4个测过KO、28个从没人测过KO——这是策展遗漏不是技术限制；1025基因全量筛查里"必需基因"和"复合体亚基"这两类会不会才是真正影响分泌的地方；这些基因KO会降低分泌，那OE能不能补救；如果补救不了，缺什么信息、去哪能获得。
- **补齐策展反应KO覆盖**：`catalog_reaction_candidates()`（`genome_wide_tradeoff.py`）改成不论策展条目原本只填了`oe_reaction_id`还是`ko_reaction_id`，同一个反应ID都生成KO和OE两个候选（之前是哪个字段有值就只测哪个方向）。29个唯一反应→58个候选（之前32个）。重跑了`catalog_reaction_screen_hlf_opn_full_ko_oe`（116任务=58候选×2靶点，15.8分钟，0错误）。**新补的KO结果很有生物学意义**：SEC61、SRP受体、SPC、PDI1/ERO1/ERV2、OST复合体、PMT、COPII/COPI/exocyst等几乎所有核心ER/Golgi通路复合体，完全敲除后都直接不可行——和真实酵母生物学一致（这些确实是分泌通路里公认必需的核心机器），只有KAR2/HRD1/CDC48/DOA10/蛋白酶体/MNN2/EMP24/SEC12等少数几个复合体敲除后细胞仍可行且分泌基本不变。
- **必需基因(115) vs 复合体亚基(146)哪个才是"藏着的分泌效应"**：精确拆解后发现两者性质完全不同。115个必需基因里只有9个是复合体亚基GPR角色，而且这9个全部是"代谢或其它反应"分类（SERPT鞘脂合成、TRE6PP/TRE6PS海藻糖代谢、ACOATA/FA120ACPH脂肪酸合成复合体、PHETRS氨酰tRNA合成酶）——不是分泌通路本身，是核心代谢必需，只是GPR恰好是AND关系。146个复合体亚基基因（OE被跳过、只解释不执行）里，65个(44.5%)的KO有真实可测效应，**且全部是降低分泌，0个提升**——因为这些也是普通代谢复合体，敲一个亚基削弱整体功能，间接拖累和分泌共享的资源预算。
- **新增"KO降低分泌"维度**：这65个（以及所有类似的、非复合体亚基的降低案例）之前完全不出现在任何一张产品化的结果表里——现有三个维度表（必需基因/KO提升有代价/KO零代价提升）只捕捉"提升"和"完全不可行"两种情况，"可行但比野生型差"的情况两头不沾，等于筛查白跑了这部分。`genome_wide_screen_analysis.py`新增`SECRETION_DOWN_THRESHOLD=0.99`和`ko_yield_down`维度（按降低程度从重到轻排序，附带`gpr_role`列方便判断是不是复合体拖累），UI在结果页加了对应折叠区。这个模块之前是全会话唯一没写单元测试的分析模块，这次一并补了`tests/test_genome_wide_screen_analysis.py`（5个用例，用合成数据不依赖真实模型）。用hLF真实数据验证：`ko_yield_down`命中86个基因（之前完全不可见）。
- **"KO降低分泌的基因，OE能不能补救"——精确查了86个基因（hLF）**：56个(65%)有OE数据（`gpr_role=single_gene`），**全部效果≈1.0，没有一个OE后真的提升分泌**——说明这些基因是"必要但不是限速步骤"：拿掉会出问题，但当前2倍过表达倍数下多加一点也没用，通路瓶颈不在这里。剩下30个(35%)是复合体亚基，OE数据缺失是设计如此（避免在没有亚基化学计量/限速证据的情况下虚构容量提升）。
  - **缺什么**：这30个基因所在复合体的亚基化学计量比例、以及哪个亚基是真正的限速亚基。这个信息模型内部没有——这些是普通代谢反应（不是策展目录里`sec_*`/`Mach_*`那种有专门kcat参数的分泌复合体伪反应），模型里没有对应的"整体复合体产能"旋钮可调。
  - **去哪能获得**：(1) 模型里目前拿不到，只有策展的`sec_*`/`Mach_*`分泌复合体反应才有这种可调的整体kcat（就是这次一直在用的反应级OE代理机制）；(2) 真实亚基化学计量/丰度数据需要外部文献或蛋白质组学数据库（比如PaxDb，或者Pichia/酿酒酵母同源文献），这不是现有代码库里能查到的；(3) **更便宜、不需要外部新数据的折中方案**：对这30个基因各自的`affected_reactions`直接做"假设性整体过表达"测试——复用已经在用的反应级kcat乘数/bound放宽机制，跳过"单基因是否代表整个复合体"这层GPR判断，明确标注这是"如果能整体共表达这个复合体会怎样"的假设性结果，不代表单基因过表达就能做到。这条路技术上立即可做，本轮没有实施，需要用户确认是否要做。

## 复合体亚基"假设性整体过表达"测试：实施与结果（2026-07-03）

- 背景：上一节排查发现，86个"KO降低hLF分泌"的基因里，30个是复合体亚基、没有单基因OE数据（`plan_gene_overexpression`原则性地拒绝对复合体亚基做单基因OE，因为没有证据支持单基因过表达能按比例提升整个复合体的产能）。用户在看到"缺什么信息、去哪能获得"的排查结论后，明确要求："那就做测试吧，不过要明确标注基于哪些假设完成"——即做，但假设必须显式、不能悄悄冒充已验证结论。
- **假设声明**（`COMPLEX_OE_HYPOTHESIS_ASSUMPTION`，`genome_wide_tradeoff.py`）：数值上假设该反应涉及的复合体所有亚基能按比例协同过表达，用和策展分泌复合体OE同一套kcat/反应上限乘数机制（factor=2.0x）代表整体产能提升；**不代表**对任何单个基因做过表达就能达到此效果——单基因过表达默认不会让复合体产能跟着涨，这正是这些基因原本被跳过OE测试的原因。亚基化学计量比例、哪个亚基是真正限速步骤，这两项信息模型里没有、代码库里也没有，需要外部文献/蛋白质组学数据（比如亚基丰度数据库）才能确认。
- **实现**（`genome_wide_tradeoff.py`）：
  - `resolve_complex_subunit_oe_hypothesis_candidates(model, gene_ids, complex_subunits)`：对给定基因列表逐个调用`plan_gene_overexpression`（不信任CSV里的`affected_reactions`字段——早期跑的`overnight_hLF_full`该字段对跳过行是空的——而是直接从活模型重新解析GPR），取`explain_only_reactions`，按反应ID去重（多个亚基基因共享同一反应时合并到`common_name`）。
  - `run_complex_subunit_oe_hypothesis_screen(...)`：单靶点入口，仿照`run_catalog_reaction_tradeoff_screen`的结构，对每个去重后的反应调用`reaction_oe_tradeoff`（`candidate_kind="complex_oe_hypothesis"`，`hypothesis_note=COMPLEX_OE_HYPOTHESIS_ASSUMPTION`）。只测OE，不测KO（KO数据已经有了，这就是筛选这些基因的依据）。
  - `run_genome_wide_ko_oe_screen_parallel.py`新增`--scope complex_hypothesis`+`--source-run <run_name>`：从一个已完成的gene-scope运行的CSV里，用`genome_wide_screen_analysis.complex_subunit_oe_hypothesis_candidates(frame, target_id)`按（KO可行且显著降低分泌 + `gpr_role=complex_subunit`）筛出候选基因列表，再对每个靶点起一个worker任务（因为解析基因→反应需要建好的模型，不能像策展库那样在建模型之前就切分任务；规模够小，单任务内顺序跑该靶点的全部反应也就是分钟级）。CSV新增`hypothesis_note`列（`load_gene_tradeoff_csv`向后兼容补默认空字符串）。
- **真实运行结果**：hLF源自`overnight_hLF_full`，30个候选基因→去重后**6个**唯一反应，3.6分钟跑完；OPN源自`overnight_OPN_full`，35个候选基因→去重后**10个**唯一反应，5.5分钟跑完（`local_runs/complex_hypothesis_hlf`、`local_runs/complex_hypothesis_opn`）。**16个反应全部成功求解，比值全部落在0.999998~1.000005之间——没有一个显示出任何补救效果**。涉及的都是核心能量代谢复合体：ATP合酶（ATPS3m）、呼吸链复合体III/IV（CYOR/CYOOm）、磷酸果糖激酶（PFK）、丙酮酸脱氢酶（PDHa/PDHb，仅OPN候选里出现）。
- **结论更新**：加上这16个假设性测试结果，"KO降低分泌"的复合体亚基基因里，**不管是有真实单基因OE数据的56个，还是靠假设性整体过表达补测的这16个反应（覆盖剩下30个基因），没有一个显示出补救信号**——"必要但不限速"这个判断现在覆盖了整个复合体亚基降低分泌的候选集合，不只是有直接数据的那一部分。即使按最宽松的假设（整个复合体同步双倍表达），核心能量代谢复合体也不是当前操作点上的分泌瓶颈。
- **产品化**：`genome_wide_screen_analysis.py`新增`complex_oe_hypothesis`维度（`DimensionalResults`新字段）——这类结果比值全部接近1.0，如果沿用`oe_yield_up`的`>1.01`阈值过滤，会导致这16行在所有维度表里都不可见（重演`ko_yield_down`之前的"筛查白跑"问题），所以这个维度不做阈值过滤，全量展示。UI（`genome_wide_screen.py`的`_render_dimension_tables`）新增对应区块，用`st.warning`（不是`st.caption`）醒目展示`COMPLEX_OE_HYPOTHESIS_ASSUMPTION`全文，明确标注"这是假设性结果，不是已验证的过表达方案"，满足"必须明确标注基于哪些假设"的要求。
- **顺手修的一个路由 bug**：核实"仿真验证"跳转逻辑（`_apply_verify_prefill`）时发现它原来用`candidate_kind != "catalog_reaction"`判断"是不是基因"——这是黑名单式判断，新加的`complex_oe_hypothesis`（`gene_id`字段实际存的是反应ID，和`catalog_reaction`一样）会被错误地识别成"基因"，跳转后填进基因输入框而不是反应输入框，导致GPR解析静默失败。改成白名单`candidate_kind == "gene"`，只有真正是基因的情况才走基因输入框，任何其他/未来新增的候选类型默认走反应输入框（fail-safe）。顺手把判断逻辑抽成纯函数`_verify_prefill_field_values`，新增`tests/test_genome_wide_screen_view.py`（4个用例，包含这个具体的回归场景）。
- 验证：`tests/test_genome_wide_screen_analysis.py`新增复合体假设维度相关用例；用真实hLF/OPN数据端到端跑通（候选筛选→反应解析→LP求解→UI渲染），AppTest确认新区块正确显示、无异常。顶层184个测试、python_pichia148个测试全部通过（16个跳过，是既有的`@slow_pipeline`门控用例，不是本轮引入）。

## 策展目录扩充：模型已支持但从未收录的32个反应 + 2处修复（2026-07-06）

- 背景：用户问"下一步该如何寻找KO/OE基因位点"。排查发现全基因组盲扫(1025基因)和策展目录(37条)基本都挖干净了——目前仅有的两个零代价正向候选（PDI1/ERO1/ERV2、PMT）都来自策展目录而不是盲扫，暗示策展目录这条路更值得深挖。用户一开始不清楚"扩充目录"具体怎么做，讨论后选择"先列候选清单再逐条核实文献"的路线。
- **关键发现（比预想更有价值）**：在列候选清单之前，先直接查了模型里全部`sec_*`/`Mach_*`分泌机器反应——模型总共有**61个**，策展目录（37条）此前只引用了**29个**，**32个反应模型已经完整支持、从未被收录测试**。这批候选比"从文献找新基因"风险低得多：反应ID直接来自活模型，不存在PEP4那种命名对不上的问题，用户确认"32个全部测"。
- **顺手揪出2处已有条目"申报不完整"**（性质和PEP4/HRD1的坑类似——写的名字比实际测的范围大）：
  - "MPOLI/MPoLII"条目`oe_reaction_id`只指向`sec_MPOLI_complex_formation`，MPoLII对应的`sec_MPoLII_complex_formation`从未被测过——已拆成"MPOLI"和"MPoLII"两条独立条目
  - "MNN2"条目只指向`sec_Mnn2pA_complex_formation`，`Mnn2pB`/`Mnn2pC`两个反应完全没覆盖——已拆成"MNN2（A/B/C）"三条独立条目
- **32个新条目分组**（`gene_catalog.py`的`SECRETION_GENE_CATALOG`，37→69条，61个唯一反应全部通过模型`model.rxns`校验存在，零遗漏零typo）：
  - 明星靶点的平行反应（6）：KAR2/PDI1/SEC61/SRP已证明是高价值靶点，模型里还有额外变体反应从未测过——KAR2（辅助型）、PDI1（单独）、PDI1/ERO1（无ERV2）、MNL1/PDI1、SEC61C、SRC
  - GET通路（3，新分组）：尾锚定蛋白ER插入的独立路径，不同于已覆盖的SEC61/SRP共翻译路径——GET1/GET2、GET3、SGT2/GET4/GET5
  - ERAD穿梭因子+变体（3）：DSK2/RAD23/PNG1/UBA1、DSK2/RAD23/UBA1、HRD核心复合体（无YOS9/USA1）
  - COPII补充变体（4）：含SEC23/24旁系同源SHL23/LST1的3个平行路径 + 1个不含EMP24/ERV29货物受体的核心简化版
  - COPI补充反应（1）：RER1/RET2/COP1（变体）
  - **GPI锚定加工（新独立类别`CAT_GPI`，5个新+1个从`CAT_GENERAL`移入）**：BST1、GUP1、PER1、CWH43/LAS21/MCD4、TED1——这5个用PubMed查证了真实文献支持，链条完整（Bst1p脱去肌醇脂肪酸→Gup1p/Per1p做sn-2脂肪酸C18→C26重塑→Cwh43进一步换成神经酰胺→Ted1p监控重塑状态决定ER出口），引用：Fujita et al. 2005 *Mol Biol Cell* [DOI](https://doi.org/10.1091/mbc.e05-05-0443)、Ghugtyal et al. 2007 *Mol Microbiol* [DOI](https://doi.org/10.1111/j.1365-2958.2007.05883.x)、Rodriguez-Gallardo et al. 2022 *Cell Reports* [DOI](https://doi.org/10.1016/j.celrep.2022.110768)
  - **液泡/内体分选（新独立类别`CAT_VACUOLAR_SORTING`，5个）**：AP-3/AP-1/GGA衔接蛋白复合体、VPS1/CHC1/CLC1、VPS4/VPS27（ESCRT）——把货物从Golgi/内体分流到液泡降解，逻辑上和PEP4/PRB1类似但作用在更早的分选步骤，建议方向标为KO（减少目标蛋白被错误分流走）
  - 核糖体/翻译机器（2，收进`CAT_GENERAL`）：核糖体、核糖体装配因子——不是分泌通路特异性，更像整体产能杠杆
- **诚实度约定**：只对真正查证过的条目（GPI这5个）写具体文献引用；其余27个新条目的`evidence`字段只写"模型XX复合体"，不编造"已报道过表达提升分泌"这类没有验证过的强声明——和目录里大多数既有条目的严谨度保持一致，不虚高。
- 验证：新增/修改后重新对照`model.rxns`核实全部61个唯一反应ID真实存在（零缺失）；顶层184个测试、python_pichia148个测试全部通过（改动不影响任何现有测试，"MNN2"重命名前先确认代码库里没有依赖该精确字符串的引用）。
- 已启动`catalog_reaction_screen_hlf_opn_expanded61`后台运行（61个唯一反应×KO+OE×2靶点=244任务，预计30分钟量级，参照此前58候选×2靶点×15.8分钟的速率外推）——结果待补充。

## 事故排查：5个任务真实卡死2小时+，根因是HiGHS无超时保护（2026-07-06）

- **现象**：上面那轮244任务的筛查，进度显示239/244后完全停滞，状态文件超过2小时没有任何更新，但6个worker子进程里有4个的CPU时间还在持续增长——不是挂起等待（那样CPU会归零），是真的在算，只是算不完。
- **诊断过程**（遵循diagnose流程，先建立可复现、带硬超时的反馈循环，再假设排序，再单点验证）：
  1. 对比"预期任务集合"（`catalog_reaction_candidates()`）和"CSV里已完成的行"，精确定位到5个卡死的`(靶点, 反应, KO/OE)`组合，全部是本轮新加的32个反应之一：`AP-1衔接蛋白复合体`KO（hLF+OPN各一次）、`VPS1/CHC1/CLC1`KO（hLF+OPN各一次）、`核糖体装配因子`OE（仅OPN）。
  2. 用单进程+`timeout`硬超时重新单独跑其中一个（`sec_Vps1p_Chc1p_Clc1p_complex_formation`KO，hLF），完全隔离掉多进程调度和并发任务的干扰——30多分钟后依然卡在同一行代码，CPU持续在涨，排除了"多进程协调bug"和"资源争抢导致的假象"两种假设。
  3. 逐段插桩发现：卡点精确定位在`solve_pcsec_maximize`内部的`scipy.optimize.linprog(method="highs")`调用本身——在它之前的`build_pcsec_constraint_matrices`（构造约束矩阵，纯确定性的稀疏矩阵拼接，无循环依赖求解结果）只需1.3秒，问题完全在HiGHS求解器这一步。
  4. **根因确认**：`solve_pcsec_maximize`调用`linprog`时`options={"presolve": True, "disp": False}`，**没有设置`time_limit`**。这个函数在全代码库有10+处调用（全基因组筛查、策展筛查、假设性OE测试、单次仿真验证……），此前从未出问题，是因为"正常"的KO/OE扰动很少会撞上让HiGHS进退两难的退化LP结构；这次新加的32个反应里，有4个恰好和已经测过的反应共享网格蛋白（Chc1p/Clc1p）这类被多处引用的酶预算组分，把其中一个强制归零后产生了HiGHS既解不出可行解、也证明不了不可行的退化结构，会无限期占用一个worker且没有任何自我恢复机制。给同一个反应加`time_limit=30`重跑，30.06秒后干净返回`HiGHS Status 13: Time limit reached`，直接实锤。
  5. 额外发现：5个卡死任务里，`核糖体装配因子`OE（OPN）单独跑其实只要31秒就能正常求解成功（3个mu点全部收敛）——它本身不是退化LP，大概率是被另外4个真正卡死的任务长期占满worker资源池后的连带受害者，不是独立的求解器问题。
- **修复**（`python_pichia/src/pcsec_pichia/probe/_prototype.py`）：给`solve_pcsec_maximize`加`time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS`参数（新增常量`DEFAULT_SOLVER_TIME_LIMIT_SECONDS = 120.0`），透传进`linprog`的`options`。这是全代码库唯一的LP求解入口，加一次保护全部调用方都受益，不需要在每个screen脚本外面再包一层进程级超时。所有现有调用方不传这个新参数也能拿到120秒的默认保护，正常求解（几秒级）完全不受影响。
- **回归测试**：`test_probe_migration.py`新增`test_solve_pcsec_maximize_honors_time_limit_instead_of_hanging_forever`——不复现具体的退化LP（那依赖这次新加的策展数据，不适合固化进快速回归用例），而是验证保护机制本身：给一次正常求解传极短的`time_limit_seconds=0.01`，断言它快速返回`success=False`且message含"time limit"，而不是跑到正常完成或挂起。原有5个该文件的测试（含MATLAB对齐的精确数值断言）全部照常通过，确认新参数的默认值不改变任何现有数值结果。
- **修复验证**：用修复后的代码，通过真实的`reaction_ko_tradeoff`/`reaction_oe_tradeoff`（生产代码路径，不是简化复现脚本）把5个原本卡死的任务全部重新单独跑了一遍，全部在合理时间内（207~284秒）返回，不再无限挂起。观察到一个一致的模式：4个真正退化的任务里，卡住的都是`mu=0.1`（本轮测试的最高生长速率）这个点，`mu=0.01`/`0.05`大多能正常求解——推测是在最高生长速率下模型有更多资源约束同时逼近上限，强制其中一个网格蛋白相关反应归零后产生了更极端的退化顶点结构，这和产能受限代谢模型（ecGEM/pcFBA）文献里报道的"通量退化"现象一致，但没有进一步深挖具体的数值机制（性价比不高，修复本身已经解决了实际问题）。
- 这类"筛查了几百个此前没测过的扰动，其中几个撞上求解器退化"的情况，只要继续做全模型/全策展规模的扫描式筛查，以后大概率还会再遇到（不限于这5个反应）——现在有了通用超时保护，未来任何类似情况都会在120秒内以内清晰的失败状态收场，不会再无限占用worker。

## 扩容后244任务完整结果分析（2026-07-06）

修复后重跑（`catalog_reaction_screen_hlf_opn_expanded61_v2`，244任务，43.7分钟，0错误）的分析结果：

- **总体**：170/244求解成功，74个不可行——全部是KO方向（0个OE不可行）。新增32反应贡献96个求解成功+32个KO不可行；原有29反应贡献74个求解成功+42个KO不可行。
- **两个新的"干净赢家"OE候选**（分泌提升、生长零代价，两个靶点都成立）：
  - **AP-1衔接蛋白复合体 OE**：OPN +5.46%、hLF +5.15%
  - **VPS1/CHC1/CLC1 OE**：OPN +2.66%、hLF +2.54%
  - 按提升幅度把所有已知干净赢家排序：PMT(OPN,+16.96%) > PDI1/ERO1/ERV2(hLF,+8.15%) > **AP-1(两靶点,+5.15~5.46%)** > **VPS1/CHC1/CLC1(两靶点,+2.54~2.66%)** > OCH1(hLF,+2.10%) > CWH41(hLF,+2.04%) > BIP/NEFS(hLF,+1.36%)。这两个新候选排第3、4位，而且是仅有的**两个靶点都成立**的候选（PMT只对OPN有效，PDI1/ERO1/ERV2只对hLF有效）——提示它们影响的是更底层、不针对特定目标蛋白修饰需求的分泌通路环节。
  - **重要保留意见**：这两个反应KO方向的结果（生长速率降到野生型50%、分泌趋近于0）**不能直接采信为"KO会严重损害生长和分泌"**——查了原始的每个mu点求解状态，最高的mu=0.1这个点，4条记录全部是`status=1`（HiGHS"时间限制到达"），不是`status=2`（真正的不可行证明）。也就是说，求解器在120秒内没能对mu=0.1给出确定答案，不代表真的不可行——当前CSV/分析层的`max_feasible_mu`统计口径把"求解超时"和"证明不可行"混为一谈（这是这次加超时保护后才第一次暴露出来的口径缺口，之前从来没有"超时"这种结果类型，所以没人注意到）。如果要拿到确定结论，需要对这4条（KO, mu=0.1）单独用更长的超时或者别的求解器设置重新跑，本轮没有做，只影响这2个新反应的KO侧解读，不影响已经确认的OE侧+5%/+2.6%发现。
  - 生物学上这个OE结果本身合理：AP-1和VPS1/CHC1/CLC1都涉及Golgi-内体间的囊泡分选/网格蛋白装配，过表达可能提升整体囊泡通量或改善分泌相关组分的回收效率；反而是我们最初把它们归到"液泡/内体分选，建议KO以减少竞争性分流"的生物学假设（见上文目录扩充部分）被数据推翻——这是"为什么要做经验性筛查而不是单纯按通路逻辑推理"的一个正面案例。
  - `KO`方向确认的效果：0个新反应KO显示分泌提升（新旧候选在这一点上一致：过去也是0个KO候选提升分泌）。
- **必需性确认**：GPI锚定加工整条链路（BST1/GUP1/PER1/CWH43/TED1，PubMed查证过的那5个）和GET通路（GET1/GET2/GET3/SGT2-GET4-GET5）全部KO不可行——和之前"核心复合体反应KO几乎都是全有全无致死开关"的结论一致；OE方向全部趋近1.0（0.003%~0.008%变化），符合"必要但不限速"模式，OE不能补救也不需要补救。
- 结论：这轮策展目录扩容除了补全数据完整性和修掉2处历史遗漏之外，**新增了2个真实、可复现、跨靶点通用的干净分泌提升候选**（AP-1、VPS1/CHC1/CLC1过表达），是这次深挖过程里除PDI1/ERO1/ERV2和PMT之外新增的实质性发现。

## 待定事项

- 复合体级 OE 反应拆解到具体基因（wet-lab 可操作性问题）——本轮不做，后续单独立项
- 跨物种 kcat 借用——待用户请教研发组长后决定
- 策展名单（`SECRETION_GENE_CATALOG`）里的 `gene_id`/文献描述可信度存疑（PEP4 的 `gene_id` 错配、HRD1 的 KO 效应和文献描述不符，见上文两处案例），后续如果要清理策展数据，需要逐条用全模型路径重新核实
- `sec_Och1p_complex_formation` KO 为什么在两个靶点上都直接不可行——目前只是记录了现象，还没深入排查是模型对这个复合体的建模方式导致的（全有全无开关），还是有其他结构性原因（**部分解答**：这轮批量补测发现几乎所有核心ER/Golgi复合体敲除都直接不可行，不是OCH1特例，是这种"整体复合体反应"建模方式本身的通性——全有全无开关，没有部分抑制的中间状态）
- 复合体亚基化学计量比例/限速亚基的真实数据（PaxDb等蛋白质组学数据库或同源文献）——如果要把"假设性整体过表达"升级成有真实依据的结论，需要这类外部数据，目前代码库和模型里都没有
- **"求解超时"和"证明不可行"目前在`max_feasible_mu`统计口径里无法区分**——`AP-1衔接蛋白复合体`/`VPS1/CHC1/CLC1`KO在mu=0.1这个点的`status=1`（超时）被当成失败处理，和真正的`status=2`（不可行）产生同样的`max_feasible_mu`降级效果，可能低估这两个新OE候选对应KO方向的真实生长/分泌上限。要拿到确定结论需要单独用更长超时重跑这4条记录；如果这类超时案例以后随着筛查规模增长变多，值得考虑把每个mu点的真实HiGHS status存进CSV，而不是只存聚合后的`max_feasible_mu`。

## 候选库简化：把"发现"职责交给筛查页（2026-07-03）

- 问题：仿真验证页面的"候选库"混乱——已验证分泌工程候选库（37条）、全模型GPR基因库（1025个基因，分页+多重筛选器）、反应级代理（sec_*/Mach_*）、外部证据GPR overlay，四块分散的浏览器面板，同一件事（找到要测试的基因/反应ID）要在4个地方找。
- 核实结论：先核实了全基因组筛查的实际覆盖范围——gene-scope已经把全部1025个基因（"全模型GPR基因库"那一批）的KO/OE都测完了（`overnight_hLF_full`/`overnight_OPN_full`各1025/1025）；catalog-scope已经把策展库去重后的32个反应级候选（"已验证候选库"+"反应级代理"背后同一份`SECRETION_GENE_CATALOG`数据）也测完了（`catalog_reaction_screen_hlf_opn`64/64）。这意味着候选库里"已验证候选库"和"全模型基因库"这两个大浏览器面板本质上是在重复"筛查已经做完的发现工作"，而且筛查结果+"在仿真验证中核实"跳转比手动翻页/搜索更直接。**唯一没被筛查覆盖的是外部证据GPR overlay**（两种scope都没有接入`enable_gene_rule_overlay`）。
- 决策（用户认可后拍板）：把候选库的"发现"职责交给筛查页，候选库本身退化成一个轻量的"直接输入/搜索已知ID"工具——合并"已验证候选库"和"全模型基因库"成一个搜索框（因为`list_verified_secretion_gene_library`已经带有`ko_reaction_id`/`oe_reaction_id`字段，和"反应级代理"面板读的是同一份底层数据，天然可以合并），保留组合测试能力（多选+一键添加到KO/OE输入），保留外部证据GPR overlay作为独立的"高级"区（因为它是唯一不冗余的部分）。
- 实现（`app/ui/views/simulation_gene_catalog.py`）：
  - 新的`render_gene_lookup_panel()`只有一个搜索框，非空查询时调用`_collect_search_rows(query)`合并`list_verified_secretion_gene_library(query)`（策展库，本身已含反应ID字段）和`load_pichia_full_model_gene_catalog()`（全模型1025基因，用已有的`_filter_full_model_gene_rows`按query过滤，截断显示前30条避免长列表）；按`(kind, id)`去重后统一渲染成一张表（来源/名称/类型/ID/可用于/说明），多选后"添加到敲除输入"/"添加到过表达输入"两个按钮通过`_partition_selection_by_kind()`把选中项按基因/反应类型分流到`pichia_draft_{ko,oe}_{genes,reactions}`四个输入框之一。
  - 删除了`_render_verified_secretion_gene_library`/`_render_full_model_gene_lookup`/`_render_reaction_proxy_lookup`/`_paginate_full_model_gene_rows`/`_page_input_widget_key`/`_add_curated_knockout_selection`/`_add_curated_oe_reaction_selection`，以及只服务于被删大浏览器的展示格式化函数（`_gene_action_label`/`_ko_status_label`/`_oe_status_label`/`_gpr_role_label`/`_process_label`/`_external_id_summary`/`_wet_lab_readiness_label`/`_curated_mapping_status_label`/`_curated_recommended_use_label`/`_reaction_proxy_evidence_label`）。顺手删了一个已经死代码的`_render_matlab_gene_target_lookup`（定义了但从未被调用）。
  - 外部证据GPR overlay面板（`_render_gene_rule_overlay_lookup`）原样保留；两个原本分散在不同面板的缓存刷新按钮（策展证据缓存、全模型湿实验注释缓存）合并到这个"高级：外部证据GPR overlay/候选库维护"区，作为唯一的维护入口。
- 验证：新增`_partition_selection_by_kind`纯函数单元测试和`_collect_search_rows`的monkeypatch合并去重测试；重写了原本断言"分页/筛选器/多个独立面板"存在的既有测试（改成断言"这些UI元素不存在"+"合并搜索的关键行为存在"）；用AppTest对真实数据（不是mock）跑了一遍完整链路——搜索"PDI1"→找到`sec_PDI1_ERV2_Ero1p_complex_formation`→选中→点"添加到过表达输入"→确认`pichia_draft_oe_reactions`被正确写入，全程无异常。顶层170个测试、python_pichia148个测试全部通过。
