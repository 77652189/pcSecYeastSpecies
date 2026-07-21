# pcSecPichia 项目级执行计划

状态：active
最后更新：2026-07-20

## 决策摘要

项目目标保持为：围绕 hLF、OPN 等目标蛋白，持续生成、解释和校准 KO、OE、分泌通路及组合改造候选，降低实验候选选择成本。项目不是绝对产量预测器。

湿实验没有固定预算和固定批次，不作为软件阶段启动条件。系统必须在没有新实验数据时继续完成科学契约、模型、数值回归和可解释输出；实验数据不定期到来后，通过独立校准层回填，不直接修改代谢矩阵、GPR 或正式科学资产。

当前保留五个研发方向：

1. 实验反馈闭环。
2. gene-level OE capacity。
3. 分泌通路机制约束。
4. 组合改造设计。
5. 条件鲁棒性筛查。

EVO2、GPU 推理和云端模型集成不属于当前范围，不创建占位接口。

### 2026-07-20：相对信号深化（ADR-004）

绝对 OE 表达量与酶容量数据经三次独立调研确认在公开来源中永久缺失（见 [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）。据此把方向 2-3 的后续投入从“找绝对数据”转向“深化不需要绝对数值也能算的相对信号”，在 [ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md) 相对决策层内新增四个方向，全部不产生绝对容量、绝对层保持 unavailable：

1. **影子价格瓶颈归因**（方向 3）：深化并透出主路径上已有的 `analyze_target_protein_lp_attribution`（非新建、不接 shadow_lp 休眠路径），保留行级/复合体级归因与 `bound_type`，下界约束不报成 OE 瓶颈。
2. **OE 剂量响应形状**（方向 2）：倍数扫描替代固定 `2.0×`，输出形状类别（平坦/单调上升/快速饱和）。
3. **排序对容量假设的稳健性标注**（方向 2-3）：`ranking-insensitive-to-capacity` / `ranking-sensitive-to-capacity`（不用 `capacity-robust` 一类会被误读的名字），稳健性同时覆盖参数带宽与求解算法，扫描带宽绝不断言为容量值。
4. **价值-of-information 实验优先级**：对 top 候选按“哪次最小测量最能消解排序歧义”排序。

## 1. 距离原始 MVP 的主要缺口

| 方向 | 当前完成度估算 | 主要缺口 |
| --- | ---: | --- |
| 实验反馈闭环 | 软件闭环已完成 | 研发发酵宽表适配和脱敏回放已验收（[ADR-003](adr/003-fermentation-feedback-minimal-fields.md) 已实现）；等待获批真实数据首次回填——已有一份周报级别的进度汇报，但数据量不构成这里说的"真实数据"，不触发回填。Streamlit 界面本地化（标签页/指标/表格列名中英文一致化、结论先行摘要）也已完成，详见 `pichia_current_architecture_and_requirements.md`"当前主要缺口"第 7 条 |
| gene-level OE capacity | 产品分层已完成 | 相对未校准决策层已可用；绝对层等待可审核 baseline capacity |
| 分泌机制约束 | Round 0 与 ERAD 激活决定均已完成 | 全部七类资源的架构与可执行契约已冻结，并核实过哪些类别已有真实 kcat 数据：转运/二硫键/糖基化/囊泡运输/folding-chaperone 五类确认有真实数据（前四类已无条件参与现有模型每次求解）；ER quality control/ERAD 确认有真实数据，hLF/OPN 均验证求解可行，已知 MATLAB 兼容性差异已确认是复核过的修正，**最终决定：约束保持可选，不改默认值**（理由见第 6 节成功条件）；仅 target-specific 的目标蛋白自身降解速率（kdeg）真正缺数据，且该降解反应目前没有 GPR，**已明确决定当前阶段不建模**（PEP4/PRB1 基因身份本身低置信度，YPS1-3 未入基因目录，且没有合理动力学建模路径；等真实湿实验结果） |
| 组合改造设计 | 约 10% | 没有组合契约、上位性计算或组合搜索后端；仍依赖方向 3 形成可靠单基因评分，ERAD 激活即使完成也只覆盖一小部分基因，未达到这个门槛 |
| 条件鲁棒性筛查 | 局部验收维度已验证 | 未做成独立产品，按第 7 节范围只作为方向 2/3 的局部验收维度使用；已用两次真实检查验证这个用法本身可行（OE 候选跨碳源/生长率排名稳定性检查、ERAD 约束开关敏感性检查） |

原始 MVP 的近期硬缺口——独立分泌资源层的架构冻结——已经完成。方向 1-3 完成交付后，用户决定继续推进方向 3 的 ERAD 约束验证与激活一轮，已完成；目标蛋白降解通路建模明确不做。真实数据回填和绝对容量证据仍是等待中的条件，不构成当前阻塞；一份周报级别的进度汇报不构成方向 1 的"真实数据回填"条件。

## 2. 当前路线的预期决策价值

### 方向 1：实验反馈闭环

已有 schema、研发发酵宽表 CSV/XLSX/JSONL 导入、质量校验、预测关联、对照匹配、top-K/evidence tier 指标、service、Streamlit 和脱敏回放证据。软件闭环已验收；获批真实数据到来后进入独立回填 checkpoint，不重新开启实现阶段。

### 方向 2：gene-level OE capacity

mapping、剂量、参数区间、复合体语义、约束、求解、报告和 UI 已形成。产品分层已验收：reaction proxy、relative uncalibrated、absolute unavailable 和 not executable 由核心层统一判定；绝对容量等待真正合格的新证据。

### 方向 3-5

方向 3 能补上传统 GEM 对真实分泌瓶颈覆盖不足的最大科学缺口，下一阶段价值最高。方向 4 依赖可信的单基因与资源层评分；方向 5 应逐步成为方向 2-4 的横向门禁，而不是最后一次性追加。

## 3. 已完成的交付

仓库没有工时或财务台账，历史投入不做量化统计，也不追踪人日消耗。当前已完成的交付：

- **实验反馈闭环**：研发发酵宽表回填收口（CSV/XLSX/JSONL 导入、质量校验、预测关联、对照匹配、离线回放）；随后按 [ADR-003](adr/003-fermentation-feedback-minimal-fields.md) 完成模板对齐现场真实 14 字段结构的重构；Streamlit 界面本地化（标签页/指标/表格列名中英文一致化、新增结论先行摘要）。
- **gene-level OE capacity**：产品层级与门禁收口，reaction proxy、相对未校准、绝对 unavailable 和 not-executable 状态清晰分层。
- **secretory resource Round 0**：冻结资源池、单位、来源、开关、边界和验收，全部七类；核实过每类是否有真实 kcat 数据。
- **ERAD 约束验证与激活决定**：hLF 可行性验证已补齐（新增 `test_pipeline_runs_builtin_hlf_with_optional_constraints`，与 OPN 版本一起通过）、已知 MATLAB 兼容性差异确认是复核过的修正、决定保持可选并把理由写进 `engines/base.py` 字段注释。
- **相对信号深化 R1（影子价格瓶颈归因 + 求解器稳健性）**（[ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）：深化既有 `analyze_target_protein_lp_attribution`，新增 bound_type 分离的 OE 可缓解瓶颈派生（下界 floor 永不进 OE 清单）与 opt-in 的求解器稳健性重解比对；核心加 `solver_method` 参数 + 常量 `DEFAULT_SOLVER_METHOD`（**默认求解算法从 highs 改为确定性的 highs-ds**，见下条）、开关全链路镜像 cost-slope、Streamlit 展示。R3-R4 仍前瞻。
- **相对信号深化 R2（OE 剂量响应形状）**（[ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）：固定 `2.0×` 单点升级为对扫描 factor 网格的分泌响应做形状分类（`saturating`/`linear`/`threshold`/`flat_no_response`/`non_monotonic_numerical_artifact`）；纯分类器 `classify_oe_dose_response_shape`（AUC 判凹凸）+ 编排 `run_oe_dose_response_sweep`；opt-in `enable_oe_dose_response` 全链路镜像 cost-slope + Streamlit 面板 + 工具 flag。相对形状信号，绝不输出真实表达倍数或 mg/L；OE 单调不减是硬约束，越噪声下降判数值假象。
- **求解器确定性（Fix A，已修复）**：`test_simulation_entrypoints.py` 目标值断言此前脆弱，根因是 `method="highs"` 在退化最优面上落到不唯一顶点（run-to-run 抖 ~2e-6，偶超 `pytest.approx` 1e-6）。**修复**：把默认求解算法钉为 `highs-ds`（对偶单纯形，天然确定，稳定落到被标定顶点）。**验证**：`python_pichia/tests/` 全量隔离 551 passed / 0 failed，计数与基线逐字一致（不挪动任何期望值）。此前"跨文件污染 + threads=1 放大"的推断已证伪，详见 `docs/handoff.md`。

当前验证基线（2026-07-20）：

- 仓库根 `tests/`：319 passed，0 failed。
- `python_pichia/tests/`：549 passed，20 skipped，0 failed。

"OE 产品层级与门禁收口"交付后新发现一处需要纠正：该切片当时产出的 Streamlit 页面（基因级 OE 容量对照）把人工输入倍率和全局统一假设的不确定性区间包装成看起来比实际更严谨的样子，已整体移除，只保留 `python_pichia.oe_capacity` 库和对应 service 供内部脚本调用。这不改变本切片的验收结论——分层判定逻辑（相对/绝对/proxy/not-executable）本身没变，仍满足第 6 节方向 2 成功条件；纠正的是产品呈现方式，不是重新开工。

## 4. 路线比较

| 路线 | 决策价值 | 风险 | 结论 |
| --- | --- | --- | --- |
分泌资源 Round 0 与 ERAD 约束验证 | 高，已有 5%-14% 真实排名影响作证据 | 范围小、已独立验证（hLF 可行性 + 默认值决定） | **已完成** |
| 建模目标蛋白降解通路（PEP4/PRB1/YPS） | 低 | 基因身份低置信度、无合理动力学路径 | 明确不做，等真实湿实验结果 |
| 继续扩大绝对容量来源搜索 | 低到中 | 重复获得不闭合单位链 | 降级为机会驱动研究 |
| 直接实现组合搜索 | 低 | 单基因和资源评分尚未收口 | 不做 |
| 全面并行五个方向 | 低 | 接口和科学语义同时漂移 | 不做 |
| 因暂无实验数据而冻结方向 2-5 | 低 | 被不确定实验时间长期阻塞 | 不采用 |

## 5. 推荐路线

推荐顺序（全部已完成，见下一节范围边界）：

1. 方向 1 的真实模板回填入口已完成；真实数据到来后走独立回填 checkpoint（周报级别的进度汇报不算）。
2. 方向 2 已完成：reaction proxy、相对未校准场景、绝对 unavailable 和不可执行状态已经分层，绝对容量继续等待证据。
3. 方向 3 的独立 secretory resource layer Round 0 已完成。
4. 方向 3 的 ERAD 约束验证与激活决定已完成：保持可选，不改默认值。
5. 目标蛋白降解通路（PEP4/PRB1/YPS）建模明确不做，未进入本阶段范围。
6. 方向 3 形成可靠单基因评分后，再决定是否推进方向 4；ERAD 激活覆盖的基因范围不构成这个门槛，方向 4 仍不做。
7. 方向 5 从方向 3 开始逐步作为横向验收门禁接入，已用两次真实检查验证可行。

每一项已产生独立可验证结果，完成一项不自动把范围扩大到下一个方向。下一阶段范围由用户决定何时推进。

## 6. 成功、失败与停止条件

### 方向 1 成功条件

- 保留现有发酵模板字段，只新增四项：改造方案（含对应基因）、数据状态、亲本对照组编号、重复编号。
- 目标蛋白和批次可作为表单级元数据，不要求逐行填写。
- 数据状态至少区分正常、污染、培养失败、检测失败和其他排除。
- 污染、失败和排除数据保留原值但不进入校准。
- 同一克隆的独立培养可用重复编号区分。
- 没有真实数据时可用脱敏 fixture 完成导入、关联、校准资格和离线回放验证。

以上是当时已验收交付的历史条件。[ADR-003](adr/003-fermentation-feedback-minimal-fields.md) 已正式接受一次后续修订并已完成实现：模板以研发组实际在用的现场 14 字段结构为基准（五条最小必要信息只是验证这 14 字段够不够用的透镜，不是替换它们的独立精简设计），补齐了唯一缺失的改造确认字段。实现细节见 `docs/handoff.md`。

### 方向 2 成功条件

- 相对 gene-capacity、reaction proxy 和 not-executable 状态不混淆。
- 没有审核 baseline capacity 时保持 unavailable。
- 不使用通用上界、最优 flux、固定 1.0 或 fixture 伪造正式容量。
- 现有 mapping、剂量、约束、求解和回归能力继续通过验证。

### 方向 3 Round 0 成功条件

- 明确代谢层、protein resource、secretory resource 和实验校准层边界。
- 每类资源池有单位、来源、适用条件、不确定性、开关和基线回归要求。
- 文献基因名单不能直接变成可执行约束。
- hLF/OPN 的目标特异成本和不可用状态有明确契约。

### 方向 3 ERAD 约束验证与激活成功条件

- hLF 和 OPN 都必须有真实求解可行性验证，不能只覆盖一个 target。
- 现有 `aligned_except_known_matlab_compatibility_differences` 这条已知差异必须先理解清楚是良性的记录差异还是真实建模问题，不能带着未解释的差异做默认开启决定。
- 默认开启还是保持可选，必须有明确记录的决定和理由，不能通过默认值静默改变方向 2 现有候选排名的语义。
- 现有 hLF/OPN 的既有回归（product tiering、relative smoke、genome-wide screen 等）必须继续通过。
- 目标蛋白降解通路（PEP4/PRB1/YPS）不在这轮范围内，不得借这轮工作顺带实现。

### 相对信号深化 R1-R4 成功条件（ADR-004）

- R1：扩展主路径已有的 `analyze_target_protein_lp_attribution`（不接 `analysis/shadow_lp` 休眠的第二份对偶提取），对 hLF/OPN 产出行级/复合体级相对 binding 贡献（可选再附层级汇总）；每条归因保留 `bound_type`，下界约束不报成 OE 瓶颈；proxy 边界下附带宽+求解器双重稳健性说明，不作绝对瓶颈断言。
- R2（已实现 2026-07-21）：OE 倍数扫描输出形状类别（saturating/linear/threshold/flat/非单调假象），固定 `2.0×` 保留为曲线上一个点，opt-in feature-off 时不改默认求解；只输出相对增益与形状，不输出真实表达倍数或 mg/L；OE 单调不减为硬约束、越噪声下降判数值假象。
- R3：测试覆盖“跨带宽稳健”“跨带宽翻转”与“跨求解器翻转”三类；各类下绝对状态均 unavailable，扫描带宽不写正式资产、不作 promotion 依据；标签名不含 `capacity-robust` 一类措辞。
- R4：价值-of-information 清单可回溯到候选、排序歧义和建议测量；不含绝对产量预测，不自动提升候选为 `experiment_calibrated` 或绝对可执行。
- service/UI 只透传与展示上述判断，不重新实现科学降级或容量推断逻辑。

### 停止条件

- 投入明显增加仍无法形成可验证交付：停止扩展并报告真实阻塞。
- 新来源在开发前不能证明含明确单位、条件、版本/hash/license 和转换链：不启动 adapter。
- 方向 3 尚未形成可靠单基因评分时，不进入组合搜索；ERAD 激活覆盖的基因范围不构成这个门槛。
- 任一路线要求用伪参数或未审核实验数据修改科学资产：立即停止。
- 目标蛋白降解通路建模没有可靠基因身份或合理动力学路径前，不得启动：当前判断是这个条件在可预见时间内不会满足，等真实湿实验结果。

## 7. 当前范围边界

### 范围内

第 3 节列出的交付（发酵模板回填、OE 产品层级、secretory resource Round 0、ERAD 约束验证与激活决定、实验反馈闭环 UI 本地化）均已完成，不需要重新实现。

2026-07-20 新增的相对信号深化四个方向（见决策摘要与 [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)），全部限定在 ADR-002 相对决策层内。**R1（影子价格瓶颈归因）已实现**（见第 3 节），扩展主路径已有的 `analyze_target_protein_lp_attribution`（不接 `secretory_resources` 休眠层，也不接 shadow_lp 休眠对偶路径）。R2（OE 剂量响应形状，相对独立）进行中；R3（容量带宽稳健性）、R4（价值-of-information，消费 R1-R3 输出）前瞻。2026-07-20 用户授权在其离开期间自主推进 R2-R4（R2 保底），每方向本地提交待其回来 review，不 push。

### 会在以下情况发生时处理

- 新实验数据到来后执行导入、质量检查、预测关联和核对更新；一份周报级别的进度汇报不满足这个条件。
- 新绝对容量来源只有在开发前满足最小 provenance 条件时，才启动接入。

### 明确不做

- 固定湿实验 pilot 或预留湿实验预算。
- EVO2、GPU、云端推理或相关占位接口。
- 方向 4 的组合搜索实现。
- 方向 5 的完整跨条件排名产品；当前只允许作为方向 2-3 的局部验收维度。
- 新增默认 solver、重写核心模型或修改受保护科学资产。
- 目标蛋白降解通路（PEP4/PRB1/YPS 家族）建模：基因身份本身低置信度待复核、没有合理动力学建模路径，明确不做，不是"以后再看"；等真实湿实验结果。

下一次重新评估范围的时机：出现合格绝对容量新证据，真实实验数据（非周报级别）首次完成回填，或有新方向/新范围需要决定是否投入。

## 文档治理

- 需求文档定义五个方向、范围和验收。
- 架构文档定义模型层、证据层、校准层和用户界面的边界。
- 长期高影响权衡进入独立 ADR。
- handoff 只记录当前状态、下一步、必读材料和验证方式。
- 提示词只负责指向当前文档和阶段，不复制完整项目设计。
