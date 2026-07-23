# ADR-006：碳源条件标定与分档（数据契约）

状态：accepted（2026-07-23）
关联：[ADR-002](002-relative-oe-and-absolute-capacity-layers.md)（相对/绝对分层）、[ADR-005](005-rnaseq-expression-constrained-enzyme-capacity.md)（同为"数据契约"型）、[ADR-003](003-fermentation-feedback-minimal-fields.md)（发酵反馈数据入口）

## 背景

模型的碳源条件轴目前有五个：`glucose / glycerol / methanol / glucose_glycerol / glycerol_methanol`（`loading` 已接好摄取/生长反应的 bound 切换，`simulation` 已按 formulation 选生长反应）。但只有 **glucose 标为 `corrected_reference`**（对齐过 MATLAB 基准）；其余在 `media` 里的 formulation_status 是 `draft_carbon_source_boundary` / `draft_induction_boundary_requires_calibration`——**边界能打开、能求解，但没有条件专属标定**，且代码自带 warning 说明。

代码核查发现的具体缺口（非移植 bug，是有意的"glucose 先行、其余标 draft"分档，标注诚实）：
- **蛋白含量硬编码 0.37（glucose 值）**：`total_protein_content` 在 constraints/simulation/lp_writer/shadow_lp 都是默认 0.37 参数，无碳源条件化；原始 MATLAB 甲醇条件用 0.40。
- **蛋白成本约束锚定 `"BIOMASS"`**（`lp_writer` 中 `modeled_protein_fraction("BIOMASS", …)`）：而甲醇/甘油走的是 `BIOMASS_meoh` / `BIOMASS_glyc`。
- **甲醇诱导/摄取/启动子调控（AOX 主导的蛋白负担）未标定**；代码注明原始 MATLAB artifact 不是可靠的同条件甲醇基准。

真实工艺需要这些条件可信：hLF 走**甘油生长 → 甘油-葡萄糖过渡 → 葡萄糖生产**（组成型启动子，不涉甲醇）；OPN 走**甲醇**（小规模、验证少）。

## 决策

1. **机械标定（无需外部数据）**：碳源条件化蛋白含量（甲醇取 0.40）、蛋白成本/生长约束改认 formulation 选定的生长反应（**glucose 的 `corrected_reference` 行为必须逐字不变**，回归锁定）、核实并补齐条件专属生物量组成。
2. **引入三档 formulation 状态（诚实分档，不假装 corrected）**：
   - `corrected_reference`：对齐过可信基准（当前仅 glucose）。
   - `internally_calibrated`：机械标定后内部一致、量级合理，但**未对齐外部实测**（A 完成后 glycerol/methanol/mixed 的目标档）。
   - `draft_boundary`：仅接了 bound、未标定。
3. **数据契约——把某条件升到 `corrected_reference` 需要**（相对信号；绝对容量无论如何恒 `unavailable`）：
   - 各碳源的比生长速率 μ、比摄取速率；
   - per-condition 生物量蛋白含量；
   - 条件专属酶负担（甲醇：AOX 主导那部分占总蛋白比例）；
   - 一个验证锚点（某确定条件下实测分泌量/产量 + 该条件参数）。
4. **在手数据先用**：已有的 hLF 发酵验证数据（本地私有区、不入库）可提供 hLF 葡萄糖/甘油工艺的 μ 与 titer 验证锚点，现在即可用于内部标定验证；甲醇/OPN 的定量数据仍缺。

## 后果

- 解锁方向5 的跨条件稳健性：稳健性只在**可信档**条件间做，并如实标注各条件的标定档（避免拿未标定条件比排序＝garbage-in）。
- 保持诚实：不给缺数据的条件贴 `corrected`；甲醇在拿到 AOX 负担 + 验证锚点前停在 `internally_calibrated`。
- 与 ADR-005 并列为"数据契约"族；湿实验数据经 [ADR-003](003-fermentation-feedback-minimal-fields.md) 入口/本地私有区读取。
- **保密边界**：一切湿实验专有数据（基因型/位点/产量）仅存仓库外本地私有区，提交进仓库的产物只含机制层抽象。

## 被否决的备选

- **不给数据就直接升 corrected**：否决——不诚实，把未标定当已验证。
- **其余条件全留 draft**：否决——堵死方向5 与真实工艺对齐。
- **照搬 MATLAB 甲醇特例脚本当基准**：否决——代码已核实该 artifact 不是可靠同条件基准。
