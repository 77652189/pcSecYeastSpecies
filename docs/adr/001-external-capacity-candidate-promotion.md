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

## 影响

优点是研发组无需自行产生容量数值，项目仍可利用公开 Pichia 蛋白组、酶约束模型和动力学数据推进；同时保留正式求解的可追溯性。代价是需要新增在线候选获取、单位转换、审核 UI 和 promotion 工作流，并接受部分候选只能保持低置信、不可正式执行。
