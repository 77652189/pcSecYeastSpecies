# ADR-004：数据永久缺失下的相对信号深化

状态：accepted
日期：2026-07-20
关系：补充 ADR-001、ADR-002 与 ADR-003，不取代任何一方的门禁与不变量

## 背景

对 Komagataella（Pichia）的绝对 OE 表达量与酶容量定量数据，已经过三次独立调研确认在当前公开来源中不存在：ecPichia/GECKO provenance closure（A0c，出口为 `architecture_decision_required`）、PRIDE `PXD055501`（仅相对 iBAQ、条件不匹配）、以及对 PRIDE 全部 Komagataella 项目的系统普查（全部为相对定量方法）。因此“补上绝对容量/表达量数据”不再是一个可达目标，而是一个永久性缺口。

在这个前提下，继续把绝对容量当作 OE 能力的唯一入口，或用 proxy、默认值、同源参数填充容量，都无法服务项目的原始目标（降低 hLF/OPN 候选选择成本）。需要一个明确决策：把工程投入转向那些**不需要绝对数值也能计算**的相对信号。

两点使这个转向可行：其一，当前 KO 判定基于化学计量矩阵与模型 GPR，是确切数据，结论可信；不可信的只有 OE 与酶容量，这是模型固有属性。其二，LP 的对偶（影子价格）本来就能回答“当前解被哪个约束卡住”，而**这一能力已经在主求解路径上实现并无条件运行**：`analysis/__init__.py` 的 `analyze_target_protein_lp_attribution` 读取主 pipeline 求解产出的 `lp_sensitivity`（来自 `probe/_prototype.py::_linprog_sensitivity`），已按命名约束块归因（stoichiometric / metabolic_coupling / secretory_coupling / protein_mass / proteasome / ribosome_assembly / ribosome_translation / misfolding / mitochondrial），并已在其 `warnings` 里记录“下界 marginal 不是 OE 线索”这条教训。因此本 ADR 的方向 1 是**深化并透出这个既有能力**，不是新建。

（注意：`analysis/shadow_lp` 后端另有一份 `_extract_duals` 对偶提取，但 `include_duals` 全项目无一处打开、`lp_sensitivity` 在该路径写死为 `None`，是休眠的第二份重复实现。方向 1 不接这条休眠路径，避免维护两份对偶提取逻辑。）

## 决策

在既有相对 OE 决策层（ADR-002 第 1 层）内，授权四项**免绝对数据**的相对信号能力：

1. **影子价格瓶颈归因（深化既有能力，非新建）**：扩展并透出主路径上已有的 `analyze_target_protein_lp_attribution`，回答“hLF/OPN 分泌当前被哪个约束卡住”。约束三点：
   - **不接休眠的 shadow_lp `_extract_duals` 路径**，只用主 pipeline 已算好的 `lp_sensitivity`；不新增第二份对偶提取。
   - **保留按约束行/按复合体的粒度，不过度聚合成粗资源层**。现有 `top_constraint_marginals` 是按行保留的，一行对应一个具体复合体（能看出是 PDI/ERO1/ERV2 还是 OST 在卡），这直接对应可操作的 OE 靶点；聚合成一个笼统的“secretory_coupling 层”反而丢掉这个可操作细节。要提供的是行级/复合体级归因，可选地再给一个层级汇总，但层级汇总不替代行级明细。
   - **`bound_type`（下界/上限）必须一路带到任何汇总层级**：现有块聚合器用绝对值求和会丢符号，若把变量下界 marginal 折进“限制层”会把一个 OE 根本动不了的下界约束（floor）报成瓶颈——这是已确认的 PDI1/核糖体假阳性（见 `analysis/__init__.py:166` 告警）。汇总时下界约束不得计入“OE 可缓解的瓶颈”。
   - 澄清两套分类不是一回事：LP 约束块分类（secretory_coupling/protein_mass/…）与 `secretory_resources` 的生物过程分类（转运/折叠/糖基化/囊泡…）是不同粒度；`secretory_coupling` 一个块里其实糊合了折叠/糖基化/转运/囊泡多个生物步骤。若未来要生物过程级分辨率（区分“折叠瓶颈”vs“糖基化瓶颈”），需按反应 stage 标签（`core/target_protein_plan.py` 已有 stage 标注）拆分，那是额外工作，不在本方向的“读既有归因”范围内。
   这是相对排序信号，不依赖任何绝对容量数值。
2. **OE 剂量响应形状**：用一组倍数（factor sweep）替代当前固定 `2.0×`，报告响应曲线的形状类别（平坦=非瓶颈 / 单调上升=候选 / 快速饱和=弱表达即足够）。输出是形状分类，对“真实表达量未知”稳健。
3. **排序对容量假设的稳健性标注**：把不确定的容量在一个合理带宽内作为**敏感性扫描输入**，只输出排序是否翻转的**稳健性分类**。标签命名为 `ranking-insensitive-to-capacity` / `ranking-sensitive-to-capacity`（排序对容量假设稳健/敏感）——**刻意不使用 `capacity-robust` 这类带 `capacity` 的名字**，避免被读成“容量已知且稳健”；标签只描述排序稳定性，绝不描述容量本身。稳健性检查必须同时覆盖两个维度：**参数带宽**（容量在带宽内扫）和**求解算法**（换 highs-ds/highs-ipm 重解）——因为影子价格/对偶解在退化最优处本身不唯一，只扫带宽不换求解器会把数值假象当成稳健。带宽永远不被断言为容量数值，绝对状态保持 unavailable。
4. **价值-of-information 实验优先级**：对相对层筛出的 top 候选，按“哪一次最小测量最能消解排序歧义”排序，产出最小湿实验清单，对接 ADR-001 的 `target_specific` 优先级。

## 不变量

ADR-001、ADR-002 与 ADR-003 的全部不变量继续有效，本 ADR 额外固定：

- 上述四项能力全部位于相对、未校准决策层内，不产生、不提升绝对 gene-capacity；绝对层在缺审核 baseline capacity 时仍为 unavailable。
- **排序稳健性（第 3 项）的扫描带宽是不确定性分析输入，不是被断言的容量**：不得写入正式容量资产，不得作为 promotion 依据，不得在任何输出中呈现为容量数值或 mg/L。`ranking-insensitive-to-capacity` 标签只表示“相对排序对容量假设不敏感”，绝不表示“绝对容量已知或可用”；标签名中不得出现会被误读为绝对断言的 `capacity-robust` 一类措辞。
- 影子价格归因基于当前模型内的约束边界；若这些边界含 proxy 值，归因结论必须一并报告稳健性，不得当作绝对瓶颈断言。稳健性检查必须同时覆盖参数带宽与求解算法（换 highs-ds/highs-ipm 重解）——对偶解在退化最优处不唯一，只扫带宽不够。
- **`bound_type`（下界/上限）必须在任何归因汇总中保留**：OE 放宽的是上限产能，对下界（最低要求类约束）无效；下界 marginal 再大也不得计入“OE 可缓解的瓶颈”，也不得在层级聚合时被绝对值求和吞掉方向信息（已确认的 PDI1/核糖体假阳性，见 `analysis/__init__.py:166`）。
- OE 剂量响应形状是相对 flux 对倍数的响应，不输出真实表达倍数、拷贝数或 mg/L；固定 `2.0×` 旧路径作为兼容对照保留。
- 价值-of-information 只对测量做优先级排序，不预测测量结果，也不自动把任何候选提升为 `experiment_calibrated` 或绝对可执行。

## 验收边界

- 影子价格归因对 hLF/OPN 分别产出**行级/复合体级**相对 binding 贡献（可选再附层级汇总），扩展既有 `analyze_target_protein_lp_attribution`、不新增第二份对偶提取；每条归因保留 `bound_type`，下界约束不被报成 OE 瓶颈；proxy 边界存在时附带宽+求解器双重稳健性说明。
- OE 剂量响应对每个候选产出形状类别与曲线；固定 `2.0×` 旧路径 feature-off 回归通过。
- 排序稳健性标注有测试覆盖“跨带宽稳健”“跨带宽翻转”**以及“跨求解器翻转”**三类，且断言绝对状态在各类下都保持 unavailable、不写正式资产；标签名不含 `capacity-robust` 一类措辞。
- 价值-of-information 输出可回溯到具体候选、具体排序歧义与建议测量；不包含绝对产量预测。
- service/UI 只透传与展示上述判断，不重新实现科学降级或容量推断逻辑。

## 影响

优点：项目可以在承认绝对数据永久缺失的前提下继续深化相对决策层，把“无法回答的绝对问题”转化为“可回答的相对与稳健性问题”，并让湿实验投入集中在真正能改变决策的测量上。代价：产品与报告必须长期维护“相对信号 vs 绝对断言”的清晰边界，尤其是容量稳健性标注绝不能被读成绝对容量；四项能力都必须保留来源、扫描带宽、稳健性依据和不可解释范围。
