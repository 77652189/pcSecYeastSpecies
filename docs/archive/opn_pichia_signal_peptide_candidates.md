# 人骨桥蛋白 OPN 的毕赤酵母候选信号肽

## 目标

本清单用于把用户提供的成熟人骨桥蛋白 OPN/SPP1 序列接到不同分泌 leader 后，在 `pcSecPichia` 模型中做低成本、可追踪的候选比较。这里的结果只能回答“在模型约束下，不同 leader 对序列长度、分泌路径和求解可行性的影响”，不能直接证明工业表达产量。

## 当前候选

| 候选 ID | leader / signal peptide | 分类 | 主要用途 |
|---|---:|---|---|
| `OPN_ALPHA_FULL_PROJECT` | 89 aa / 19 aa | 项目基线 | 与项目内现有 alpha-factor 建模方式保持一致，是默认对照。 |
| `OPN_ALPHA_PRE_ONLY` | 19 aa / 19 aa | 酵母短信号肽 | 去掉 alpha pro 区，观察只保留 signal peptide 时的模型成本。 |
| `OPN_NATIVE_SPP1` | 16 aa / 16 aa | OPN 人源天然信号肽 | 作为生物学参考组，比较天然信号肽和酵母信号肽。 |
| `OPN_OST1N23_ALPHA_PRO` | 93 aa / 23 aa | 酵母混合 leader | 用 OST1 N 端 pre 序列替换 alpha-factor pre 区，保留 alpha pro 区。 |
| `OPN_PPA_DDDK18` | 18 aa / 18 aa | 毕赤酵母来源短信号肽 | 避免 alpha pro/Kex2 路线，适合作为 OPN 的重点候选之一。 |
| `OPN_PPA_PASCHR3_0030` | 20 aa / 20 aa | 毕赤酵母来源短信号肽 | Pichia-native 短 leader 候选，用于筛选比较。 |
| `OPN_PPA_EPX1_SA` | 20 aa / 20 aa | 毕赤酵母来源短信号肽 | 另一个 Pichia-native 短 leader 候选，用于筛选比较。 |

模型输入文件：

- `Data/pcSecPichia/TargetProtein_OPN.csv`：单个项目基线 OPN 目标，保持向后兼容。
- `Data/pcSecPichia/TargetProtein_OPN_candidates.csv`：多个候选目标，可直接被 `local_opn_pichia_glc.m` 读取。
- `Data/pcSecPichia/TargetProtein_OPN_candidates_meta.csv`：候选来源、理由和注意事项。

## 重要生物学注意点

用户提供的是去除信号肽后的成熟 OPN，长度 298 aa。成熟 OPN 内部包含二碱性位点，例如 `RR` 和 `KR`。如果使用 alpha-factor prepro leader，理论上会经过 Kex2/Ste13 类加工路线，因此需要在真实构建设计中关注内部二碱性位点带来的异常切割风险。

当前 pcSec 模型会把 leader 长度、氨基酸组成、ER 分泌路径、翻译/折叠/分泌资源约束写入 LP，但不会模拟 Kex2 是否真的切割 OPN 内部位点，也不会预测糖基化异质性、蛋白降解、宿主蛋白酶、培养工艺或实际滴度。

## 外部依据

- 人 OPN/SPP1 天然信号肽：UniProt P10451。
- 酵母 OST1 N 端序列：UniProt P41543。
- DDDK 18 aa 信号肽：`A new signal sequence for recombinant protein secretion in Pichia pastoris`，PubMed 24317479。
- PAS_chr3_0030 信号肽：`Potential of the Signal Peptide Derived from the PAS_chr3_0030 Gene Product for Secretory Expression of Valuable Enzymes in Pichia pastoris`，PMCID PMC9088389。

## 推荐下一步

1. 先跑 `OPN_ALPHA_FULL_PROJECT`，作为项目基线。
2. 再跑 `OPN_PPA_DDDK18` 和 `OPN_PPA_PASCHR3_0030`，优先比较无 alpha pro/Kex2 路线的 Pichia 短 leader。
3. 如果模型结果都可行，再把这些候选交给实验侧做构建设计、密码子优化、切割位点检查和小规模表达验证。

示例 MATLAB 调用：

```matlab
opts = struct('mediaType', 4, 'writeMisfoldingConstraints', false, 'writeRibosomeConstraint', false);
local_opn_pichia_glc(0.10, 1e-8, [], opts, 'OPN_PPA_DDDK18');
```
