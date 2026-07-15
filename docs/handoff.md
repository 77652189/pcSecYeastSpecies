# pcSecPichia Handoff

状态：active
最后更新：2026-07-15

## 当前执行位置

```yaml
current_program: mvp_directions_1_to_3
current_slice: direction_2_oe_product_tiering_closure
slice_status: ready
direction_1_status: accepted_waiting_for_real_data_replay
relative_oe_status: authorized_ready_for_product_closure
absolute_capacity_status: unavailable_waiting_for_qualified_evidence
```

## 当前状态

- A0c 已完成，现有 PRIDE/ecPichia 证据不能形成审核后的绝对 baseline capacity；正式 registry 未修改。
- ADR-002 已接受相对 OE 决策层与绝对容量研究层分离。绝对容量继续保持 unavailable，不再扩大同类低信息来源接入。
- 方向 1 已通过验收：研发发酵宽表 CSV/XLSX/JSONL 已接入 canonical validation、cache、prediction linkage、calibration eligibility 和报告链路。
- 脱敏回放覆盖正常、污染、培养失败、检测失败、其他排除、亲本对照、独立培养重复和阴性结果；失败/排除原值保留且不进入校准。
- 尚未读取获批真实研发数据；真实数据到来后只执行独立回填 checkpoint，不重新开启方向 1 开发。
- 当前进入方向 2 产品收口。现有 PRIDE/ecPichia 证据不能形成审核后的绝对 baseline capacity，正式 registry 保持不变。
- ADR-002 已接受 reaction proxy、相对未校准 gene-capacity 与绝对 unavailable 三种产品状态必须分离。

## 已授权切片

只执行 `direction_2_oe_product_tiering_closure`：让已有 OE 能力按照 ADR-002 形成一致、可审计的核心状态、报告、service 和 Streamlit 产品门禁。

本切片不得新增外部容量来源、伪造绝对容量、进入方向 3 secretory resource layer、组合搜索或完整跨条件排名。

## 必读材料

1. [项目级执行与预算计划：方向 2 成功条件与授权边界](EXECUTION_PLAN.md#方向-2-成功条件)
2. [当前架构：实验校准层与产品验收分层](pichia_current_architecture_and_requirements.md#产品验收分层)
3. [ADR-002：相对 OE 与绝对容量分层](adr/002-relative-oe-and-absolute-capacity-layers.md)
4. [Phase 2 既有实现与边界](pichia_next_plan.md#phase-2gene-level-oe-与酶容量)
5. [数据与结果治理策略](data_and_results_policy.md)

## 验收与停止线

- 必须产生实际产品状态或门禁收口及 focused tests，不能只修改文档或重复外部来源审计。
- reaction proxy、relative uncalibrated、absolute unavailable 和 not executable 必须由核心层判定并在报告、service、UI 一致透传。
- 缺少审核 baseline capacity 时不得调用绝对容量求解或生成 nominal capacity；不得用 proxy、最优 flux、通用上界、固定值或 fixture 替代。
- 保留旧 proxy 数值和 feature-off 回归；hLF/OPN 必须保持 target/context/model version 与风险说明隔离。
- `app/services` 和 Streamlit 只做 facade/展示，不重新实现模式判断或科学降级逻辑。
- 完成 review/fix/verify 后更新状态并停止；不得自动进入方向 3。
- 验证 focused tests、hLF/OPN smoke、相关 service/UI contract、`compileall`、保护目录、依赖、密钥和 ignore 边界。
