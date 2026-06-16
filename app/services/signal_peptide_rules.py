from __future__ import annotations

import re
from dataclasses import dataclass


KYTE_DOOLITTLE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}
HYDROPHOBIC = set("AILMFWVYC")
SMALL_NEUTRAL = set("ACSTVG")
POSITIVE = set("KRH")
NEGATIVE = set("DE")
AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


@dataclass(frozen=True)
class RuleScreeningResult:
    score: int
    passed: bool
    tier: str
    reasons: list[str]
    risks: list[str]
    length: int
    n_region_positive_count: int
    h_region_max_hydrophobicity: float
    c_region_small_neutral_rule: bool


def score_signal_peptide(sequence: str) -> RuleScreeningResult:
    seq = sequence.strip().upper()
    reasons: list[str] = []
    risks: list[str] = []
    score = 0

    length = len(seq)
    if 15 <= length <= 35:
        score += 25
        reasons.append("长度在常见可切割信号肽范围内")
    elif 12 <= length <= 45:
        score += 12
        risks.append("长度略偏离常见范围，需要人工复核")
    else:
        risks.append("长度明显偏离常见信号肽范围")

    if not AA_PATTERN.fullmatch(seq):
        risks.append("序列含非标准氨基酸字符")

    n_region = seq[:8]
    positive_count = sum(1 for aa in n_region if aa in POSITIVE)
    negative_count = sum(1 for aa in n_region if aa in NEGATIVE)
    if positive_count >= 1:
        score += 15
        reasons.append("N 端含正电残基，符合许多分泌信号肽特征")
    elif negative_count == 0:
        score += 8
        risks.append("N 端缺少正电残基，但没有明显负电性")
    else:
        risks.append("N 端负电性偏强，可能不利于典型分泌信号肽")

    max_hydrophobicity = max_window_hydrophobicity(seq, window=8)
    middle = seq[5 : min(length, 22)]
    hydrophobic_count = sum(1 for aa in middle if aa in HYDROPHOBIC)
    if max_hydrophobicity >= 1.8 and hydrophobic_count >= 6:
        score += 30
        reasons.append("中段疏水核心明显")
    elif max_hydrophobicity >= 1.3 and hydrophobic_count >= 5:
        score += 20
        risks.append("中段疏水核心可接受但不强")
    else:
        risks.append("中段疏水核心偏弱")

    c_rule = c_region_small_neutral_rule(seq)
    if c_rule:
        score += 20
        reasons.append("切割位点附近符合 small-neutral residue 偏好")
    else:
        risks.append("切割位点附近 small-neutral residue 特征不明显")

    if _has_low_complexity_run(seq):
        risks.append("存在较长重复残基，可能是低复杂度风险")
    else:
        score += 10

    hard_fail = length < 10 or not AA_PATTERN.fullmatch(seq) or max_hydrophobicity < 0.8
    passed = score >= 65 and not hard_fail
    tier = "规则推荐" if passed else ("需人工复核" if score >= 50 and not hard_fail else "规则不推荐")
    return RuleScreeningResult(
        score=min(score, 100),
        passed=passed,
        tier=tier,
        reasons=reasons,
        risks=risks,
        length=length,
        n_region_positive_count=positive_count,
        h_region_max_hydrophobicity=round(max_hydrophobicity, 3),
        c_region_small_neutral_rule=c_rule,
    )


def max_window_hydrophobicity(sequence: str, window: int = 8) -> float:
    seq = sequence.strip().upper()
    if not seq:
        return 0.0
    if len(seq) < window:
        values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq]
        return sum(values) / len(values)
    best = -99.0
    for index in range(0, len(seq) - window + 1):
        values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in seq[index : index + window]]
        best = max(best, sum(values) / window)
    return best


def c_region_small_neutral_rule(sequence: str) -> bool:
    seq = sequence.strip().upper()
    if len(seq) < 4:
        return False
    minus_three = seq[-3]
    minus_one = seq[-1]
    if minus_three == "P" or minus_one == "P":
        return False
    return minus_three in SMALL_NEUTRAL and minus_one in SMALL_NEUTRAL


def _has_low_complexity_run(sequence: str) -> bool:
    return bool(re.search(r"([A-Z])\1{6,}", sequence.upper()))
