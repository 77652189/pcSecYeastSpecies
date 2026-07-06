# COBRApy 引入前 Phase 0 基线评估

日期：2026-07-06

范围：`pcSecYeastSpecies / python_pichia`

## 结论

当前项目没有使用 COBRApy。模型数据来自 MATLAB/COBRA 风格的 `.mat` 结构，但 Python 运行路径使用项目内部轻量模型容器和 `scipy.optimize.linprog(method="highs")` 求解。

建议进入 Phase 1，但仅限于 opt-in shadow mode：

- 不替换当前 `PichiaModel` / `CobraModel` 容器。
- 不替换当前 SciPy HiGHS / MATLAB 对齐路径。
- 不默认安装或启用 COBRApy。
- 不让 COBRApy 结果驱动 UI 推荐、报告结论或 KO/OE 排名。
- 只在小范围 baseline FBA 对照中比较 COBRApy 转换结果与当前 SciPy HiGHS 结果。

Phase 1 的目标应是验证“COBRApy 是否能忠实复现基础 GEM FBA”，不是迁移 pcSec 约束、目标蛋白分泌能力、KO/OE 推荐或产量结论。

## 当前基线

### 模型容器

当前存在两套内部轻量模型结构：

- `pcsec_pichia.core.pichia_model.PichiaModel`
  - 字段包括 `rxns`、`mets`、`genes`、`lb`、`ub`、`c`、`b`、`s_matrix`、`rules`、`gr_rules`、`rxn_gene_mat`。
  - 提供 `change_rxn_bounds()`、`with_reaction_bounds()`、`set_bound()`、`set_objective()`、`add_reaction()` 等不可变式操作。
- `pcsec_pichia.probe._prototype.CobraModel`
  - 名字里有 `Cobra`，但不是 COBRApy `cobra.Model`。
  - 是迁移期 prototype 使用的内部 dataclass，字段包括 `rxns`、`mets`、`genes`、`lb`、`ub`、`b`、`s_matrix`、`rules`、`gr_rules`。

这个命名容易误导。Phase 1 如引入 COBRApy，建议明确使用类似 `CobraPyModelAdapter` / `CobraPyShadowSolver` 的命名，避免把内部 `CobraModel` 误认为 COBRApy 对象。

### 模型加载

模型加载当前由 `MatStructLoader.load_pcsec_pichia_model()` 读取 MATLAB `.mat`：

- 读取 `model` 结构。
- 提取 `modelID`、`rxns`、`mets`、`genes`、`lb`、`ub`、`c`、`b`、`S`、`rules`、`grRules`、`rxnGeneMat`。
- `S` 被转换为 SciPy sparse matrix。

prototype 路径也直接从 `Model/pcSecPichia.mat` 读取同类字段，转换为内部 `CobraModel`。

这说明当前基线是“MATLAB 模型结构 -> 内部 Python 容器 -> SciPy LP”，不是“SBML/COBRApy 模型 -> COBRApy optimize”。

### 基础 FBA 求解

基础 LP/FBA 求解由两条路径承担：

- `pcsec_pichia.adapters.lp_solver.ScipyHiGHSSolver.solve()`
- `pcsec_pichia.probe._prototype.solve_maximize()`

核心形式一致：

- 目标：最大化目标反应时将目标向量取负，因为 `linprog` 默认最小化。
- 等式约束：`A_eq = S`，`b_eq = b`。
- 变量边界：`bounds = zip(lb, ub)`。
- 求解器：`scipy.optimize.linprog(..., method="highs")`。

当前测试 `tests/test_lp_solver.py` 将 `Ex_glc_D` 作为标准 FBA 目标之一，验证 glucose uptake 为负。

### pcSec 分泌能力求解

目标蛋白分泌能力入口是 `solve_secretion_capacity()`：

- 先准备目标蛋白相关反应和 secretory/combined enzyme data。
- 固定 growth reaction bounds。
- 调用 `solve_pcsec_maximize()` 最大化目标 exchange/secretion reaction。

`solve_pcsec_maximize()` 不是单纯 GEM FBA。它调用 `build_pcsec_constraint_matrices()` 追加 pcSec 约束后再用 SciPy HiGHS 求解。

当前 pcSec 约束包括：

- stoichiometric constraints
- metabolic coupling
- secretory coupling
- protein mass
- proteasome
- ribosome assembly
- optional ribosome translation
- optional misfolding
- mitochondrial inequality constraints

这些约束是项目科学核心，不应在 Phase 1 直接交给 COBRApy 替代。

### KO/OE screen

当前 screen 边界包括：

- `run_knockout_screen()`
  - 通过 `plan_gene_knockout()` 从 GPR 解释 gene KO 会失活哪些 reactions。
  - 对可执行 KO，将相关 reaction bounds 设为 `(0, 0)`。
  - 调用 `solve_pcsec_maximize()` 与 baseline 对比。
- `run_reaction_knockout_screen()`
  - 对 reaction-level KO 直接运行 pcSec KO screen。
- `run_overexpression_screen()`
  - 当前 OE 是 reaction capacity proxy。
  - 通过放大 reaction bound/capacity 近似评估，不是完整基因表达调控模型。
  - 调用 prototype pcSec OE screen，再标准化输出。
- `plan_gene_overexpression()`
  - 判断 gene 关联 reactions 哪些可作为 reaction-level proxy，哪些只能 explain-only。
  - complex-subunit 或 mixed GPR 情况保留人工解释/警告。

因此，KO/OE 逻辑包含两类边界：

- 可执行 LP 边界：reaction bounds / capacity proxy / pcSec solve。
- 解释和证据边界：GPR 可解析性、gene mapping、OE proxy warning、phenotype evidence tier。

COBRApy 若引入，只能先对第一类中的基础 GEM FBA 做 shadow 对照；不能把 OE proxy、GPR planning、phenotype evidence 或 recommendation tier 直接升级为 COBRApy 结论。

## 适合 COBRApy 替代或对照的部分

Phase 1 可考虑用 COBRApy shadow mode 对照以下基础能力：

- 从当前内部模型转换为 COBRApy `cobra.Model` 的最小 adapter。
- 基础 reaction bounds 和 objective 设置。
- 基础 GEM FBA：`S·v=b`、`lb/ub`、单目标 reaction maximize/minimize。
- 小样本关键 flux 对照：
  - `BIOMASS`
  - `Ex_glc_D`
  - `Ex_glyc`
  - `Ex_meoh`
  - `Ex_o2`
  - 当前目标 exchange reaction
- missing objective / infeasible / solver failure 的状态映射。

这些对照应只用于验证转换正确性，不应改变当前生产路径。

## 不应由 COBRApy 替代的部分

Phase 1 不应替代以下部分：

- MATLAB `.mat` 作为当前权威模型输入的读取路径。
- `PichiaModel` / prototype `CobraModel` 作为当前 engine 内部数据结构。
- `solve_pcsec_maximize()` 中的 pcSec constraint matrix 构建。
- protein mass、secretory coupling、metabolic enzyme coupling、ribosome/proteasome/mitochondrial constraints。
- target protein reaction construction、signal peptide / mature protein / PTM burden 相关准备逻辑。
- gene KO/OE planning 的 GPR 解释逻辑。
- OE reaction proxy 语义和 warning。
- phenotype evidence / recommendation tier。
- Streamlit UI 行为、报告结论、用户推荐排序。
- 任何 mg/L 或绝对产量预测。

这些是项目当前核心科学语义或用户解释层。COBRApy 可以成为对照工具，不应在 Phase 1 成为默认事实来源。

## Shadow Mode 必须对齐的 Baseline 指标

Phase 1 若建立 `CobraPyShadowSolver`，至少需要记录并比较：

1. 模型结构
   - reaction count
   - metabolite count
   - gene count
   - stoichiometric matrix shape
   - nonzero stoichiometric entries
   - lower/upper bound checksum 或差异列表
2. 目标设置
   - objective reaction id
   - objective sense
   - objective coefficient
3. Solver 状态
   - success
   - status code/category
   - message/category
   - infeasible/unbounded/missing objective 映射
4. 目标值
   - SciPy HiGHS objective value
   - COBRApy objective value
   - absolute difference
   - relative difference
5. 关键 flux
   - `BIOMASS`
   - selected carbon source exchanges
   - oxygen exchange
   - target exchange reaction, if present
6. KO/OE 边界
   - KO disabled reactions must match current GPR planning output.
   - reaction KO bound changes must match current bound mutations.
   - OE proxy must retain reaction-level proxy warning.
7. pcSec 边界
   - Phase 1 不要求 COBRApy 复现 pcSec constraints。
   - 若仅比较 base GEM FBA，文档和输出必须明确说明不含 pcSec protein/secretion constraints。

建议数值容忍：

- objective 和 flux 对照先使用 `abs(diff) <= 1e-7` 或 `rel(diff) <= 1e-6`。
- 若存在 solver convention 差异，应记录具体 reaction、bound、objective，而不是放宽为模糊通过。

## 推荐最小 Phase 1 切片

建议的 Phase 1 是一个 opt-in、不可影响主路径的 shadow adapter：

1. 新增可选依赖检测
   - 不在默认 requirements 中强制加入 COBRApy。
   - 如果 `cobra` 不可 import，shadow mode 返回 `unavailable`。
2. 新增内部转换器
   - 输入：当前 `PichiaModel` 或 prototype `CobraModel`。
   - 输出：COBRApy `cobra.Model`。
   - 只转换 metabolites、reactions、stoichiometry、bounds、objective。
   - 不在第一步承诺完整 gene/reaction rule 语义。
3. 新增 shadow solve
   - 只跑基础 GEM FBA。
   - 不跑 pcSec constraints。
   - 不参与 UI 推荐和报告排序。
4. 新增 focused tests
   - COBRApy 未安装时行为稳定：skip 或 `unavailable`。
   - 安装 COBRApy 的环境中，比对 2-3 个小型基础 FBA cases。
   - 验证 no-default-enable。
5. 新增文档标记
   - COBRApy shadow result 是 solver parity check。
   - 不是实验产量预测。
   - 不是 pcSec 约束替代。

推荐 Phase 1 checkpoint 名称：

`feat(pichia): add optional cobrapy shadow fba baseline`

只有当 shadow mode 完全 opt-in 且默认路径不变时才适合进入该 checkpoint。

## 风险清单

- 命名混淆：内部 `CobraModel` 不是 COBRApy model。
- solver 差异：SciPy HiGHS 与 COBRApy 后端可能在状态码、objective sign、容忍度上不同。
- 约束差异：COBRApy 标准 FBA 不能直接表达当前 pcSec 追加矩阵的全部语义。
- GPR 差异：当前 KO planning 使用项目自有 rule token 解析；COBRApy gene deletion 语义可能与当前 MATLAB/rxnGeneMat 对齐方式不同。
- OE 误读：当前 OE 是 reaction-level proxy，不是 gene expression simulation。
- 依赖风险：`cobra` 常带 `optlang`、SBML、solver backend 依赖；默认安装可能扩大环境复杂度。
- 用户解释风险：shadow mode 若暴露不清，容易被误解为绝对产量、mg/L 或实验成功率预测。
- 工作区风险：当前工作区已有大量与 Phase 0 无关的未提交改动，后续 checkpoint 必须 patch-level 分离。

## Phase 0 验收标准

本 Phase 0 完成标准：

- 明确当前是否使用 COBRApy。
- 明确当前 FBA/pcSec LP 是如何构建和求解的。
- 明确 COBRApy 可影响和不可影响的边界。
- 明确 Phase 1 只能作为 opt-in shadow mode。
- 明确 baseline 指标和推荐 focused tests。
- 不修改 `Code/`、`Model/`、`Enzymedata/`、`Results/`。
- 不安装依赖。
- 不提交、不 push。
- 不改变 Streamlit UI 或核心求解行为。

## 已运行命令与结果

仓库状态：

```powershell
git status --short --branch
```

结果摘要：

- 当前分支：`main...origin/main`
- 工作区已有大量既有 modified/deleted/untracked 文件。
- 本文档新增前，工作区已包含 docs 清理、UI、pipeline、screen、genome-wide 相关未提交改动。

最近提交：

```powershell
git log -5 --oneline
```

结果：

```text
c679889 docs: refresh bilingual readmes
30e7708 chore(streamlit): preserve app icon asset
35e5b29 test(pichia): skip missing local scope artifacts
2437840 chore(dev): add pcsec orphan process cleanup script
47f6edb feat(pichia): add gene evidence catalog cache tools
```

保护目录检查：

```powershell
git diff --name-only -- Code Model Enzymedata Results
```

结果：空输出。

依赖安装检查：

```powershell
python -c "import importlib.util; names=['cobra','optlang','swiglpk','libsbml']; [print(f'{name}_installed={importlib.util.find_spec(name) is not None}') for name in names]"
```

结果：

```text
cobra_installed=False
optlang_installed=False
swiglpk_installed=False
libsbml_installed=False
```

依赖声明检查：

```powershell
Select-String -Path requirements.txt, python_pichia\pyproject.toml -Pattern 'cobra|cobrapy|optlang|swiglpk|python-libsbml|libsbml' -CaseSensitive:$false
```

结果：空输出。

Focused tests：

```powershell
python -m pytest -q tests\test_mat_loader.py tests\test_lp_solver.py
```

结果：

```text
10 passed in 57.69s
```

```powershell
python -m pytest -q python_pichia\tests\test_simulation_entrypoints.py python_pichia\tests\test_screens_entrypoints.py
```

结果：

```text
29 passed, 5 skipped in 146.21s (0:02:26)
```

语法检查：

```powershell
python -m compileall -q python_pichia\src\pcsec_pichia app
```

结果：通过，空输出。

## 建议后续 focused tests

后续 Phase 1 前，建议继续保留以下 baseline tests：

```powershell
python -m pytest -q tests\test_mat_loader.py tests\test_lp_solver.py
python -m pytest -q python_pichia\tests\test_simulation_entrypoints.py python_pichia\tests\test_screens_entrypoints.py
python -m compileall -q python_pichia\src\pcsec_pichia app
```

如果 Phase 1 添加 COBRApy shadow adapter，应新增独立测试文件，且在 COBRApy 未安装时不失败默认 CI。

## Phase 1 状态：Opt-in COBRApy Shadow FBA

日期：2026-07-06

Phase 1 已按最小切片设计为完全可选的 shadow FBA adapter：

- 新增 `pcsec_pichia.adapters.cobrapy_shadow`。
- 不在默认 `requirements.txt` 或 `python_pichia/pyproject.toml` 中加入 COBRApy、optlang、swiglpk 或 libSBML。
- 模块导入时不强制 import `cobra`；只有调用转换/求解函数时才检测可选依赖。
- 当前环境未安装 COBRApy 时，shadow API 返回 `available=False`、`status="unavailable"`，message 明确包含 `COBRApy is not installed`。
- 默认 pipeline、simulation、screen、report、UI、service contract 均不 import 或调用该 shadow adapter。

Phase 1 adapter 的转换范围严格限定为基础 GEM FBA parity check：

- metabolites
- reactions
- stoichiometry
- lower / upper bounds
- objective reaction
- objective sense

Phase 1 不转换或不承诺：

- pcSec protein / secretion constraints
- metabolic / secretory coupling rows
- protein mass、proteasome、ribosome、misfolding、mitochondrial constraints
- gene reaction rule 语义
- KO/OE GPR planning
- OE reaction proxy 语义
- phenotype evidence / recommendation tier
- Streamlit UI、报告结论、候选排序或实验建议
- mg/L、绝对产量或实验成功率预测

可用 API：

```python
from pcsec_pichia.adapters.cobrapy_shadow import (
    cobrapy_available,
    convert_to_cobrapy_model,
    build_cobrapy_shadow_model,
    solve_cobrapy_shadow_fba,
    compare_shadow_fba,
)
```

当前推荐使用方式：

1. 主路径仍先运行当前 SciPy HiGHS 基础 FBA 或 pcSec 求解。
2. 只有在开发者显式调用 `solve_cobrapy_shadow_fba(...)` 时，才运行 COBRApy shadow FBA。
3. 若 COBRApy 未安装，记录 `unavailable`，不让默认测试或用户工作流失败。
4. 若 COBRApy 已安装，仅将结果用于基础 GEM solver parity 对照，不进入 UI 推荐或报告排序。

Phase 1 focused tests：

```powershell
python -m pytest -q python_pichia\tests\test_cobrapy_shadow_adapter.py
```

该测试覆盖：

- COBRApy 未安装时不抛 `ImportError`。
- 缺失 objective reaction 的稳定状态返回。
- 默认 pipeline/simulation/screen/report/UI/service 不 import 或调用 shadow adapter。
- shadow comparison 在 unavailable 时稳定返回不可比较状态。
- 若 COBRApy 已安装，才运行 tiny model parity；否则 skip。

## Phase 2 状态：真实 pcSecPichia 模型 Shadow FBA Parity Harness

日期：2026-07-06

Phase 2 新增一个仅供开发者显式调用的真实模型 parity harness：

- 新增 `pcsec_pichia.analysis.cobrapy_shadow_baseline`。
- 可通过 `run_cobrapy_shadow_baseline(...)` 或 `python -m pcsec_pichia.analysis.cobrapy_shadow_baseline` 显式运行。
- 加载当前 `Model/pcSecPichia.mat` 到既有 `PichiaModel`。
- 只对基础 GEM FBA 做 SciPy/HiGHS 与 COBRApy shadow FBA 对照。
- objective cases 只在模型中存在时纳入，默认候选为 `BIOMASS`、`Ex_glc_D`、`Ex_glyc`、`Ex_meoh`、`Ex_o2`。
- 输出 JSON 与 Markdown 到 ignored `local_runs/cobrapy_shadow_baseline/`。

Phase 2 harness 记录每个 case 的：

- `case_id`
- `objective_reaction`
- `sense`
- `shadow_available`
- `current_success`
- `shadow_success`
- `current_objective_value`
- `shadow_objective_value`
- `objective_abs_diff`
- `objective_rel_diff`
- `within_tolerance`
- `key_flux_diffs`
- `model_summary`
- `warnings`

Phase 2 仍保持以下边界：

- COBRApy 仍是 optional dependency。
- 不在默认 `requirements.txt` 或 `python_pichia/pyproject.toml` 中加入 COBRApy、optlang、swiglpk 或 libSBML。
- COBRApy 未安装或安装残缺时，结果为 `shadow_available=false`、`status="unavailable"`，harness 不失败。
- 默认 pipeline、simulation、screen、report、UI、FastAPI/service 均不 import 或调用 Phase 2 harness。
- harness 不转换 pcSec protein/secretion constraints。
- harness 不作为 KO/OE recommendation、phenotype tier、报告排序、mg/L、绝对产量或实验成功率的依据。
- harness 输出是本地验证 artifact，默认不提交。

Phase 2 focused tests：

```powershell
python -m pytest -q python_pichia\tests\test_cobrapy_shadow_baseline.py
```

该测试覆盖：

- COBRApy spec 存在但 import 失败时，真实模型 harness 返回 unavailable，不抛异常。
- artifact 输出路径必须位于 `local_runs/` 下。
- harness 拒绝写入 `Results/` 等受保护科学资产目录。
- 默认 pipeline/simulation/screen/report/UI 不 import Phase 2 harness。
