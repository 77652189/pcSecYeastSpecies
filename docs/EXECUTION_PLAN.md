# pcSecPichia 项目级执行计划

状态：active
最后更新：2026-07-23

> 本文件只列**当前与待启动的阶段及任务**。已完成阶段不在此保留——其形成的能力见 [需求与架构文档](pichia_current_architecture_and_requirements.md)，逐次历史见 git。五方向的完整定义与状态、分层架构、数据与产物治理同样以需求与架构文档为准；技术计划不得绕过本文件扩大范围。

## 项目目标

围绕 hLF、OPN 等目标蛋白，持续生成、解释、校准 KO/OE/分泌通路候选，降低实验候选选择成本。**不是绝对产量预测器**；绝对 OE/容量数据经三次独立调研确认公开来源永久缺失，主力放在相对决策层（[ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md) / [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）。湿实验数据不定期到来，经独立校准层回填，不直接改代谢矩阵 / GPR / 正式科学资产。

## 阶段① 碳源条件标定 + 跨条件稳健性（当前 · 进行中）

方向5 从"方向 2-3 的局部验收维度"**有界升级**为一个 slice。目标：把真实工艺碳源条件标定到可信档，并给短名单排序加"跨条件稳健性 + 湿实验一致性"标注。依据见 [ADR-006](adr/006-carbon-source-condition-calibration.md)。

- [x] **A1** 碳源条件化蛋白含量（甲醇 0.40 vs 葡萄糖 0.37），从 formulation 喂到求解 — 已完成（`media.protein_content` + 4 处求解/装配点接线）
- [x] **A2** 蛋白成本 / 生长约束改认 formulation 选定的生长反应（去掉 `"BIOMASS"` 硬编码）；**glucose 的 `corrected_reference` 结果逐字不变，回归锁定** — 已完成（`biomass_reaction_id` 贯穿；护栏 `test_simulation_entrypoints` 10 passed）
- [x] **A3** 核实并按需补齐甲醇条件专属生物量组成（`*_meoh`）— 已完成（核实 `BIOMASS`/`BIOMASS_glyc`/`BIOMASS_meoh` 各自 PROTEIN[c] 系数 0.8798/0.5903/0.4977，无需补齐）
- [x] **A4** 内部标定验证：各条件可行求解、生长 / 分泌量级合理、蛋白预算正确 — 已完成（五条件全部 success/status=0；量级 glucose>glycerol≫methanol；混合条件数值≈其选定单碳源，即已记录的单一生物量近似）
- [x] **A5** 三档 formulation 状态（`corrected_reference` / `internally_calibrated` / `draft_boundary`）+ UI 诚实标注 — 已完成（非葡萄糖=`internally_calibrated`；结果页"碳源标定档"行 + i18n 三档标签）
- [x] **B1** 求解结果缓存层（内容寻址 key、默认读缓存、显式重算）+ 预热常见条件 — 已完成（`pcsec_pichia.solve_cache`：`SecretionSolveCacheKey` 含碳源/μ/flags/solver + KO/OE 模型变体指纹 + schema 版本，只缓存 success、不改求解语义，缓存存 `local_runs/solve_cache/` gitignored；`tools/prewarm_secretion_solve_cache.py` 预热常见 (目标×碳源×μ)；实测 ~11s 求解 → 二次 0.0s 命中）
- [x] **B2** 短名单跑真实工艺条件矩阵 — 已完成（`tools/run_shortlist_condition_matrix.py`：复用 `_oe_shortlist` + `enable_oe_dose_response`，干净单碳源 hLF{甘油,葡萄糖}/OPN{甲醇}、统一 μ=0.10，缓存 `{target}_condition_matrix.json`：per-条件 baseline + reaction_shapes + 跨条件 per-reaction 视图；**只产矩阵数据**，分类留 B3、面板留 B4；混合/过渡条件因单一生物量近似不进排序比较。验证 hLF×{甘油,葡萄糖} 2/2 成功）
- [ ] **B3** 噪声门控（**当前条件集 0 翻转、暂不建**）：真实工艺矩阵（hLF 甘油/葡萄糖 top-15、OPN 甲醇）跑出 **0 个条件敏感候选**——hLF 短名单 15/15 两碳源均 `saturating`、排序逐字不变（分泌机器瓶颈杠杆与碳源无关，合 R1 folding-limited），无表观敏感可甄别。换 highs-ds/highs-ipm 分真敏感 vs 数值假象（复用 R3 `compare_ranking_robustness`）门控在未来出现表观敏感（更多靶点/条件）时再建。
- [ ] **B4** ① 短名单读出加"跨条件稳健性"+"湿实验一致性"标注（后者读本地发酵数据，输出相对 / 抽象结论）
  - [x] 跨条件稳健性：服务读 B2 条件矩阵缓存、按 reaction 附稳健性判定（`cross_condition_stable`/`_sensitive`/`_single`）；面板"跨条件稳健性"列 + note（未扫描优雅降级；真敏感 vs 数值假象留 B3）
  - [ ] 湿实验一致性：待 B5 私有数据读取护栏 + validated titer
- [ ] **B5** 在手发酵数据本地接线：仓库外私有区（`CursorProject/pcSec_wetlab_private/`）+ gitignored 本地路径配置 + 提交护栏；方向1 本地摄入，提 μ + titer 验证锚点 + UPR×折叠一致性交叉验证
  - [x] hLF 发酵 μ 数据已核验（14 个 μ 从 OD 重算吻合、量级合理）：**模型默认 μ=0.10 与 hLF 甘油生长相实测一致**（跨温度/pH 稳健），确认默认合理；葡萄糖生产相实测慢（限量补料工艺操作点）但按决定模型仍用生长相 μ、不单独锚生产相 μ。原始数据归档私有区
  - [x] gitignored 私有路径配置 + 运行时读取护栏：`pcsec_pichia/wetlab_private.py`（env `PCSEC_WETLAB_PRIVATE_DIR` / 约定同级回退；**拒绝仓库树内路径**防误入 git + 拒绝越界穿越；未配置优雅降级；`.env.example` 文档化，7 单测）
  - [x] UPR×折叠交叉验证：模型 hLF 折叠/二硫键受限（PDI/ERO1 头号杠杆）论断有**独立湿实验佐证**（湿实验聚焦折叠/UPR 轴 + 发酵 UPR 信号，两条独立路径收敛）；机制层结论进架构方向1、细节留私有区
  - [ ] 余：titer 验证锚点（待验证数据）、方向1 本地摄入（经护栏读私有区）
- [ ] 测试 + 重跑文档锚点 / 契约回归

**范围边界**：无界完整跨条件排名产品**仍不做**；只产相对信号、不产绝对容量；不动 glucose 的 `corrected_reference` 基准；不换默认 solver；保密湿实验数据只存仓库外本地私有区、提交产物只含机制层抽象。

## 阶段② RNA-seq 表达约束建模（待启动 · 数据门控）

触发：拿到生产菌株 RNA-seq。契约与方法已定，见 [ADR-005](adr/005-rnaseq-expression-constrained-enzyme-capacity.md)。数据到位后：transcript→酶丰度上界 → 经 curated 基因→复合体映射触达分泌 binding 层 → 相对 / opt-in → 验证。绝对容量恒 `unavailable`，**不**解锁降解层建模。

## 收尾条件：碳源条件升 corrected（非独立阶段）

阶段① 把非葡萄糖条件做到 `internally_calibrated`；升到 `corrected_reference` 需研发 / 文献补齐（[ADR-006](adr/006-carbon-source-condition-calibration.md) 数据契约）：各碳源比生长 / 比摄取速率、per-condition 蛋白含量、甲醇 AOX 酶负担、一个验证锚点。**hLF 的 μ 与 titer 锚点已在手**（本地私有区）先用；甲醇 / OPN 定量仍缺。

## 明确不做

- **方向4 组合 / 多基因搜索**（含 GA / SA / MILP）：真实组合改造在模型范围之外，模型内搜遗传组合低价值；留到有稳定性标注的可信排序后，仅在 top 短名单做有界两两上位性（O(k²)、小 k、仍相对）。
- **目标蛋白降解通路（PEP4/PRB1/YPS 家族）建模**：基因身份低置信度待复核、无合理动力学路径；等真实湿实验结果，不是"以后再看"。
- 伪造绝对容量（通用上界 / 最优 flux / 固定 1.0 / fixture）、改写受保护科学资产、引入新默认 solver。
- 固定湿实验 pilot / 预留湿实验预算；EVO2 / GPU / 云端推理或占位接口。

## 文档治理

- **本执行计划**：只列当前 / 待启动阶段与任务，剔除已完成阶段。
- **[需求与架构文档](pichia_current_architecture_and_requirements.md)**：五方向完整定义与状态、分层架构、能力边界、数据与产物治理。
- **ADR**：长期高影响决策（分层 ADR-002、相对信号深化 ADR-004、RNA-seq 数据契约 ADR-005、碳源标定数据契约 ADR-006 等）。
- **handoff**：当前目标、下一步、必读材料、验证方式。
