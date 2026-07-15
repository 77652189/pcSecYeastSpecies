# pcSecPichia Handoff

状态：active
最后更新：2026-07-15

## 当前执行位置

```yaml
current_program: mvp_directions_1_to_3
current_slice: direction_3_secretory_resource_round_0
slice_status: ready
direction_1_status: accepted_waiting_for_real_data_replay
direction_2_status: accepted_product_tiering_closed
relative_oe_status: available_uncalibrated_with_independent_solver_path
absolute_capacity_status: unavailable_waiting_for_qualified_evidence
```

## 当前状态

- A0c 已完成，现有 PRIDE/ecPichia 证据不能形成审核后的绝对 baseline capacity；正式 registry 未修改。
- ADR-002 已接受相对 OE 决策层与绝对容量研究层分离。绝对容量继续保持 unavailable，不再扩大同类低信息来源接入。
- 方向 1 已通过验收：研发发酵宽表 CSV/XLSX/JSONL 已接入 canonical validation、cache、prediction linkage、calibration eligibility 和报告链路。
- 脱敏回放覆盖正常、污染、培养失败、检测失败、其他排除、亲本对照、独立培养重复和阴性结果；失败/排除原值保留且不进入校准。
- 尚未读取获批真实研发数据；真实数据到来后只执行独立回填 checkpoint，不重新开启方向 1 开发。
- 方向 2 已验收：核心层统一判定 reaction proxy、relative uncalibrated、absolute unavailable 和 not executable；report、service 与 Streamlit 只透传和展示。
- relative uncalibrated 使用独立的 current-model enzyme-coupling 成对求解，不读取 formation 通用上界或 baseline optimal flux；absolute 公共入口反查 runtime 审核 catalog、asset hash/version 和 baseline 数值。
- hLF/OPN G6PDH2 新鲜 smoke 均为 `relative_uncalibrated`，绝对模式稳定为 `absolute_unavailable`；正式绝对容量 acceptance 仍为 `passed=false`，正式 registry 未修改。

## 已授权切片

只执行 `direction_3_secretory_resource_round_0`：冻结 secretory resource layer 的资源池、单位、来源、适用条件、不确定性、开关、基线回归和 hLF/OPN 不可用状态契约。

本切片不得实现完整 secretory mechanism 求解、进入组合搜索或完整跨条件排名，也不得用文献基因名单直接生成约束。

## 必读材料

1. [项目级执行与预算计划：方向 3 Round 0 成功条件与授权边界](EXECUTION_PLAN.md#方向-3-round-0-成功条件)
2. [当前架构：实验校准层与产品验收分层](pichia_current_architecture_and_requirements.md#产品验收分层)
3. [Phase 3：分泌资源与蛋白稳态约束](pichia_next_plan.md#phase-3分泌资源与蛋白稳态约束)
4. [ADR-002：相对 OE 与绝对容量分层](adr/002-relative-oe-and-absolute-capacity-layers.md)
5. [数据与结果治理策略](data_and_results_policy.md)

## 验收与停止线

- Round 0 只冻结可执行契约和边界，必须明确代谢层、protein resource、secretory resource 与实验校准层的所有权。
- 每类资源池必须声明单位、来源、适用条件、不确定性、开关和 baseline/feature-off 回归。
- 文献基因名单、外部注释或同源关系不能直接提升为可执行约束；无合格参数时必须 unavailable/not executable。
- hLF/OPN 的目标特异成本与不可用状态必须分开表达。
- 完成 review/fix/verify 后更新状态并停止；不得自动进入方向 3 的机制实现或方向 4。
