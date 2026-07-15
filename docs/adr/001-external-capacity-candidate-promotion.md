# ADR-001：外部容量候选与正式资产分层

状态：accepted
日期：2026-07-14

## 背景

Phase 2 gene-level OE capacity 已具备 mapping、参数、约束、求解、报告和正式验收门禁，但当前仓库和研发组都无法提供 hLF/OPN 的审核 baseline capacity。使用通用 reaction upper bound 或 baseline optimal flux 会造成循环定义或无效约束，不能作为长期方案。

## 决策

引入独立的外部容量候选层。联网来源、外部模型和人工导入首先进入 ignored `local_runs/oe_capacity/`，经过当前模型映射、单位转换、条件匹配、冲突检查和人工审核后，才可提升到 `Enzymedata/oe_capacity_baseline_capacity.json`。

候选适用范围固定为：

- `target_specific`
- `host_condition`
- `external_model_calibrated`
- `homolog_transferred`

匹配优先级固定为：

```text
target_specific
> host_condition
> external_model_calibrated
> homolog_transferred
```

同源转移不能单独使 Phase 2 正式验收通过。正式求解不联网，只读取冻结的本地容量资产快照。

## 不变量

- 外部名称、GPR 或同源关系不能覆盖当前模型 gene/enzyme/reaction identity。
- BLAST/RBH 只负责映射证据，不生成 capacity 数值。
- 每个正式锚点必须能回溯到原始值、转换公式、单位链、条件、版本、hash、license 和 reviewer。
- 不使用 `1000` 通用上界、baseline optimal flux、固定 `1.0` 或 smoke fixture 补齐正式容量。
- 没有合格候选时保持不可执行和 `reviewed_baseline_capacity` 缺口。

## 真实来源完成标准

以下能力只能证明候选工作流存在，不能证明外部容量来源已经接入：

- UniProt、NCBI 或其他数据库只确认 gene/protein identity。
- 人工 CSV/TSV/JSON 导入器及其 smoke fixture。
- 尚未包含定量值的 source inventory 或 `manual_import_required` 记录。
- 只验证公式、单位和 promotion 流程的合成候选。

Round 6A 的真实来源 checkpoint 至少需要一个公开定量来源完成可重复获取或正式文件解析，并产出带原始值、版本、hash、license、条件和单位链的记录。优先来源为 Pichia 定量蛋白组或 iPichia/ecPichia；BRENDA/SABIO-RK 只能提供动力学部分，仍需可追溯的丰度或直接容量来源。

在尚未完成上述来源接入时，不得把 Round 6A 标记为 `awaiting_candidate_review`，也不得重新把提供容量数值的责任默认交回研发组。此时状态保持 `in_progress`。

## 证据闭合与停止决策

接入一个公开定量来源只证明 source workflow 可用，不等于已经获得正式 baseline capacity。对同一候选连续获得相对丰度、外部模型字段或间接动力学证据后，必须进入一次有边界的 provenance closure，不得无限增加同类 adapter。

provenance closure 必须逐项回答：

- abundance/concentration 的物理单位、biomass basis、归一化方法和原始测量来源是什么；
- kcat 是否能追溯到明确的 enzyme、reaction/substrate、宿主或转移关系、条件和原始记录；
- 外部模型 coefficient 的定义和生成公式是什么，是否能由公开字段独立复算；
- 外部 gene/enzyme/reaction 如何映射到当前模型 formation/dilution handle；
- 从原始量到 `model_flux` capacity 的每一步单位是否闭合；
- 菌株、培养基、碳源和生长率差异如何进入 applicability、uncertainty 和 warning；
- source artifact、version、hash、license 和 reviewer 是否完整。

closure 只有两个合法出口：

1. **evidence chain closed**：生成可重复的候选、换算 trace 和 promotion preview；仍需明确人工批准后才能更新正式资产。
2. **evidence chain unresolved**：生成结构化 gap report，并将当前 checkpoint 标记为 `architecture_decision_required`。此后不得继续接入只提供相同类型间接证据的来源，必须先决定是否调整 Phase 2 的产品验收层级。

`architecture_decision_required` 不代表允许降低科学门禁。若后续拆分验收等级，绝对容量校准必须继续保持独立状态；相对 OE 场景也必须明确标为相对、未校准且不能用于绝对产量解释。

该产品分层决策已由 [ADR-002](002-relative-oe-and-absolute-capacity-layers.md) 接受。ADR-001 继续拥有绝对容量候选与 promotion 门禁，未被取代。

## 模块和调用边界

- `external_refs/capacity_sources.py` 拥有联网获取、许可元数据和原始 cache。
- `oe_capacity/external_candidate_schema.py` 只拥有候选契约和基础校验。
- `oe_capacity/external_candidate_io.py` 只拥有本地序列化、人工导入和离线回放。
- `oe_capacity/external_candidate_evaluation.py` 拥有 current-model binding、单位换算、冲突和状态推导。
- `oe_capacity/external_candidate_promotion.py` 拥有 preview、asset hash、merge、validation 和原子替换。
- `oe_capacity/external_candidate_audit.py` 提供 CLI 和 service 共用的公开 orchestration API。
- `app/services` 只调用上述公开 API；`python_pichia/tools` 不得导入 `app.services` 或任何私有 `_...` 函数。

可以暂时保留 `external_candidates.py` 作为兼容 facade，但不得继续把新职责加入该文件。

## 影响

优点是研发组无需自行产生容量数值，项目仍可利用公开 Pichia 蛋白组、酶约束模型和动力学数据推进；同时保留正式求解的可追溯性。代价是需要维护独立的来源、候选评估和 promotion 边界，并接受部分候选只能保持低置信、不可正式执行。
