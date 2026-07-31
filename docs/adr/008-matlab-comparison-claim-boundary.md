# ADR-008：与旧 MATLAB 实现「对照」的声称边界

状态：accepted（追认 2026-06 已在代码中生效的边界）  
日期：2026-07-31（决策发生于 2026-06-23 ~ 06-29，本文为补记）  
关联：与 [ADR-002](002-relative-oe-and-absolute-capacity-layers.md)（绝对容量恒 `unavailable`）同类——都是**声称边界**，限制系统可以说什么，而非改变它算什么。**不替代**任何既有 ADR。

## 背景：一条代码正在强制、却没写下来的约束

2026-06 的 MATLAB 对齐工作收口时定下了一组「什么能叫已对齐」的规矩。这些规矩**只落进了代码和一条报告字符串**，没有进任何文档。后果在 2026-07-29 显形：项目负责人看到界面显示「待对齐」，判断为「明明已经完成了，显示错了」——即**这条边界连决策者本人都不知道它还在生效**。

这正是 ADR 存在的理由：它是唯一不能从仓库推导的东西——代码能告诉你 `corrected_condition` 映射到 `pending`，但说不出**为什么不许改这一行**。

## 决策

### 1. 修正条件下不作任何对照声称

默认运行条件是 `compatibility_mode="corrected"`（服务层各入口的默认值）。该模式下 `python_result_status` 恒为 `corrected_condition`，而分类函数对它**无条件**返回 `pending`——不看目标值差多少、不看约束是否逐行匹配。

理由：修正后的培养基条件与旧 MATLAB 基线**不是同一个条件**。跨条件比对即使数值接近也不构成「对齐」，只会制造一个看起来有依据、实际无意义的结论。

**因此界面显示「待处理 / 与历史实现对照：待对齐」是设计意图，不是缺陷。** 它表达的是「这次运行没有可比对的参照物」，与结果好坏无关。任何把这条映射「修好」的改动都是在移除一条诚实性边界。

### 2. 七个状态，各自只说一件事

`ALIGNMENT_STATUSES` 是封闭集合：

| 状态 | 含义 |
| --- | --- |
| `baseline_missing` | 没有基线产物可比 |
| `matlab_failed` | 基线侧自己就没跑成功 |
| `python_draft` | Python 侧还是草稿状态 |
| `pending` | **不具备比对前提**（含 corrected 条件） |
| `aligned` | 目标值在容差内且约束逐行匹配 |
| `not_aligned` | 有基线、比过了、不匹配 |
| `aligned_except_known_matlab_compatibility_differences` | 除**已登记**的兼容性差异外一致 |

最后一个状态只在差异**已被逐条登记**（`KNOWN_OPN_...` / `KNOWN_HLF_PROJECT_710_MATLAB_COMPATIBILITY_EXCEPTIONS`）时给出。未登记的差异不得归入此类——否则这个状态会退化成「凡是对不上的都算已知差异」。

### 3. 原始 hLF 恒为 `matlab_failed`；harness 归一化产物不等于原始目标已对齐

原始 hLF 目标在基线侧从未跑成功。为使比对能进行而做的 harness 归一化产物（`hLF_CLEAN`）使用的目标蛋白定义与 Python 内置定义**不同**，它的对齐结论**不可以**转述成「原始 hLF 已对齐」。

### 4. 否掉的方案

- **「目标值相对差 ≤ 1% 即判 aligned」**：只看目标值会掩盖约束结构差异——两套不同的 LP 完全可能给出接近的最优值。
- **「行级 LP diff 收敛到个位数即判 aligned」**：2026-06 的探针确实把差异压到了个位数行，但残留差异的**来源**未定（基线侧漏写约束 / Python 侧多写耦合项 / 纯 row-mapping 假象）。来源不明就判定对齐，等于把未解释的差异改名成已解决。
- **「corrected 条件下也给出对照结论」**：见决策 1。

## 不变量（守卫已焊）

- `classify_alignment_status(python_result_status="corrected_condition", ...)` 恒返回 `pending`，与其余参数无关。
- 原始 hLF 恒为 `matlab_failed`。
- 默认管线的对照状态**不得**为 `aligned`。

守卫见 `python_pichia/tests/test_alignment_entrypoints.py`。改动这些行为必须先写一条取代本 ADR 的新 ADR。

## 后果

- 正面：一条已在生效、却只有代码知道的边界被写下来了；界面上的「待对齐」不再被误读为「还有一步没做」；将来有人想「修好」那行映射时，会先撞上守卫和本文。
- 代价：默认条件下永远看不到 `aligned`。这是有意的——想要对照结论必须显式切到可比对的条件。
- 风险边界：本 ADR 不改变任何数值行为，只约束**可以说什么**。不影响 glucose 的 `corrected_reference` 基准，不触碰绝对容量的 `unavailable` 状态。
