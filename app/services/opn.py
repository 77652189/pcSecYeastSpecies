from __future__ import annotations

import csv
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from app.adapters.matlab import MatlabAdapter
from pcsec_pichia.adapters.process_runner import ProcessRunner
from app.adapters.soplex_parser import parse_soplex_file
from app.core.models import OpnCandidate, OpnCandidateRank, OpnConstructDesign, OpnSimulationResult
from app.core.opn_inputs import BuiltinOpnInputProvider
from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.core.target_inputs import target_input_set_from_opn_input_set
from pcsec_pichia.core.target_protein_plan import build_target_protein_plan
from app.engines.matlab_pichia_engine import MatlabPichiaEngine
from pcsec_pichia.engines.base import PichiaSimulationRequest
from app.services.simulation import MU_MAX, MU_MIN
from pcsec_pichia.services.target_simulation import PichiaTargetSimulationService as PythonPichiaSimulationService
from pcsec_pichia.services.target_simulation import PichiaTargetSimulationService


DEFAULT_OPN_CANDIDATE = "OPN_PPA_DDDK18"
DEFAULT_OPN_MEDIA_TYPE = 4
DEFAULT_OPN_PRODUCTION_RATIO = 1e-8
OPN_SHORTLIST = ("OPN_PPA_PASCHR3_0030", "OPN_PPA_DDDK18", "OPN_ALPHA_FULL_PROJECT")
ALPHA_PRO_MARKER = "APVNTTTEDETAQIPAEAVIGY"
OpnEngineMode = Literal["matlab", "python"]


OPN_CATEGORY_LABELS = {
    "project_baseline": "项目基线",
    "yeast_signal_only": "酵母短信号肽",
    "target_native_signal": "OPN 天然信号肽",
    "hybrid_yeast_leader": "酵母混合 leader",
    "pichia_native_signal": "毕赤酵母来源信号肽",
}

OPN_RECOMMENDATION_PROFILES = {
    "OPN_PPA_PASCHR3_0030": {
        "rank": 1,
        "recommendation": "首选实验候选",
        "experimental_role": "优先构建 1",
        "evidence_level": "较强：Pichia 来源，公开研究报道",
        "risk_level": "较低：signal peptidase 路线，避开 alpha pro/Kex2",
        "reason": "Pichia-native 短 leader，最近一轮模型成本最低，并且不经过 alpha-factor pro 区。",
    },
    "OPN_PPA_DDDK18": {
        "rank": 2,
        "recommendation": "并行实验候选",
        "experimental_role": "优先构建 2",
        "evidence_level": "较强：Pichia DDDK 短信号肽报道",
        "risk_level": "较低：signal peptidase 路线，避开 alpha pro/Kex2",
        "reason": "虽然模型成本不是最低，但它是 Pichia 相关短信号肽，适合与 PAS_chr3_0030 并行小试。",
    },
    "OPN_ALPHA_FULL_PROJECT": {
        "rank": 3,
        "recommendation": "阳性对照/项目基线",
        "experimental_role": "必须保留的对照",
        "evidence_level": "强：Pichia 常用 alpha-factor 分泌 leader",
        "risk_level": "中等：alpha pro/Kex2 路线，OPN 内部有 RR/KR 位点",
        "reason": "它是最常用、最容易对照的分泌 leader，但 OPN 内部二碱性位点使异常切割风险更值得关注。",
    },
    "OPN_PPA_EPX1_SA": {
        "rank": 4,
        "recommendation": "备选筛选候选",
        "experimental_role": "预算允许时加入",
        "evidence_level": "中等：Pichia-native 筛选候选",
        "risk_level": "中等：signal-anchor 行为可能依赖下游蛋白",
        "reason": "模型成本居中，可作为扩大筛选的备选，但证据不如 PAS_chr3_0030/DDDK18 清晰。",
    },
    "OPN_NATIVE_SPP1": {
        "rank": 5,
        "recommendation": "生物学参考组",
        "experimental_role": "参考，不建议作为首选",
        "evidence_level": "中等：OPN 天然信号肽，但非 Pichia 优化",
        "risk_level": "中等：哺乳动物信号肽在 Pichia 中表现不确定",
        "reason": "模型成本较好，但宿主适配性证据弱于 Pichia-native 候选。",
    },
    "OPN_OST1N23_ALPHA_PRO": {
        "rank": 6,
        "recommendation": "工程化参考组",
        "experimental_role": "可选对照",
        "evidence_level": "中等：酵母 Ost1-alpha 混合思路",
        "risk_level": "中等：仍保留 alpha pro/Kex2 路线",
        "reason": "模型表现接近 alpha-factor 基线，但没有解决 OPN 内部二碱性位点风险。",
    },
    "OPN_ALPHA_PRE_ONLY": {
        "rank": 7,
        "recommendation": "不建议首轮构建",
        "experimental_role": "模型对照",
        "evidence_level": "弱：去掉 pro 区后分泌成熟证据不足",
        "risk_level": "较低：不走 alpha pro/Kex2，但分泌效率不确定",
        "reason": "主要用于拆分 alpha pre/pro 区的模型成本，不适合作为首轮工业表达候选。",
    },
}


def opn_category_label(category: str) -> str:
    return OPN_CATEGORY_LABELS.get(category, category)


@dataclass
class OpnCandidateCatalog:
    paths: ProjectPaths

    def list_candidates(self) -> list[OpnCandidate]:
        if not self.paths.opn_candidate_meta_csv.exists():
            return []
        with self.paths.opn_candidate_meta_csv.open(newline="", encoding="utf-8") as handle:
            return [OpnCandidate.model_validate(row) for row in csv.DictReader(handle)]

    def get_candidate(self, candidate_id: str) -> OpnCandidate:
        for candidate in self.list_candidates():
            if candidate.candidate_id == candidate_id:
                return candidate
        raise KeyError(f"未找到 OPN 候选：{candidate_id}")

    def construct_designs(self) -> list[OpnConstructDesign]:
        rows_by_id = self._load_model_rows()
        candidates = {candidate.candidate_id: candidate for candidate in self.list_candidates()}
        rankings = {ranking.candidate_id: ranking for ranking in self.rank_candidates()}
        designs: list[OpnConstructDesign] = []
        for candidate_id, row in rows_by_id.items():
            candidate = candidates.get(candidate_id)
            if candidate is None:
                continue
            ranking = rankings.get(candidate_id)
            full_sequence = row["sequence"]
            leader_sequence = candidate.leader_sequence
            mature_sequence = full_sequence[len(leader_sequence) :] if full_sequence.startswith(leader_sequence) else ""
            contains_alpha_pro = ALPHA_PRO_MARKER in leader_sequence
            designs.append(
                OpnConstructDesign(
                    candidate_id=candidate_id,
                    experimental_role=ranking.experimental_role if ranking else "待定",
                    recommendation=ranking.recommendation if ranking else "待评估",
                    leader_sequence=leader_sequence,
                    signal_peptide_sequence=candidate.signal_peptide_sequence,
                    mature_opn_sequence=mature_sequence,
                    full_protein_sequence=full_sequence,
                    leader_length=len(leader_sequence),
                    signal_peptide_length=len(candidate.signal_peptide_sequence),
                    mature_opn_length=len(mature_sequence),
                    full_protein_length=len(full_sequence),
                    contains_alpha_pro_region=contains_alpha_pro,
                    processing_route=candidate.processing_route,
                    kex2_risk=(
                        "需要重点检查：alpha pro/Kex2 加工路线，成熟 OPN 内部含 RR/KR 位点"
                        if contains_alpha_pro
                        else "相对较低：不含 alpha pro 区，主要按 signal peptidase 路线考虑"
                    ),
                    codon_optimization_next="需要：用 PichiaCLM 或同等工具为毕赤酵母生成 CDS",
                    note=ranking.reason if ranking else candidate.rationale,
                )
            )
        return sorted(designs, key=lambda design: rankings.get(design.candidate_id).rank if rankings.get(design.candidate_id) else 99)

    def rank_candidates(self) -> list[OpnCandidateRank]:
        candidates = self.list_candidates()
        objective_by_id = latest_opn_objectives(self.paths)
        solved = [(cid, value) for cid, (value, _text, optimal, _path) in objective_by_id.items() if optimal and value is not None]
        model_order = {
            cid: index + 1
            for index, (cid, _value) in enumerate(sorted(solved, key=lambda item: item[1], reverse=True))
        }
        best_objective = max((value for _cid, value in solved), default=None)
        rows: list[OpnCandidateRank] = []
        for candidate in candidates:
            value, text, optimal, output_file = objective_by_id.get(candidate.candidate_id, (None, None, False, None))
            profile = OPN_RECOMMENDATION_PROFILES.get(
                candidate.candidate_id,
                {
                    "rank": 99,
                    "recommendation": "待评估",
                    "experimental_role": "待定",
                    "evidence_level": "待补充",
                    "risk_level": "待补充",
                    "reason": candidate.rationale,
                },
            )
            delta = None
            if best_objective is not None and value is not None and best_objective != 0:
                delta = abs((best_objective - value) / best_objective) * 100
            rows.append(
                OpnCandidateRank(
                    candidate_id=candidate.candidate_id,
                    category=candidate.category,
                    category_label=opn_category_label(candidate.category),
                    recommendation=profile["recommendation"],
                    experimental_role=profile["experimental_role"],
                    rank=int(profile["rank"]),
                    model_rank=model_order.get(candidate.candidate_id),
                    objective_value=value,
                    objective_text=text,
                    objective_delta_percent=round(delta, 4) if delta is not None else None,
                    optimal=optimal,
                    leader_length=candidate.leader_length,
                    construct_length=candidate.construct_length,
                    processing_route=candidate.processing_route,
                    evidence_level=profile["evidence_level"],
                    risk_level=profile["risk_level"],
                    reason=profile["reason"],
                    output_file=output_file,
                )
            )
        return sorted(rows, key=lambda row: row.rank)

    def _load_model_rows(self) -> dict[str, dict[str, str]]:
        if not self.paths.opn_candidate_csv.exists():
            return {}
        with self.paths.opn_candidate_csv.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        return {row["Protein name"]: row for row in rows}


@dataclass
class OpnSimulationService:
    paths: ProjectPaths
    matlab: MatlabAdapter
    runner: ProcessRunner = field(default_factory=ProcessRunner)
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def run_candidate_smoke(
        self,
        candidate_id: str = DEFAULT_OPN_CANDIDATE,
        mu: float = 0.10,
        production_ratio: float = DEFAULT_OPN_PRODUCTION_RATIO,
        timeout_seconds: int = 300,
        engine_mode: OpnEngineMode = "matlab",
        output_dir: Path | None = None,
    ) -> OpnSimulationResult:
        if not MU_MIN <= mu <= MU_MAX:
            raise ValueError("生长速率 mu 必须在 0.01 到 0.44 h^-1 之间。")
        if engine_mode not in ("matlab", "python"):
            raise ValueError(f"不支持的 OPN 仿真引擎：{engine_mode}")
        catalog = OpnCandidateCatalog(self.paths)
        catalog.get_candidate(candidate_id)
        if not self._lock.acquire(blocking=False):
            return OpnSimulationResult(
                success=False,
                mu=mu,
                candidate_id=candidate_id,
                production_ratio=production_ratio,
                media_type=DEFAULT_OPN_MEDIA_TYPE,
                message="已有 OPN 仿真任务正在运行，请稍后再试。",
            )
        try:
            if engine_mode == "python":
                plan = _builtin_opn_plan(candidate_id)
                engine_result = PythonPichiaSimulationService(self.paths).run_glucose_smoke(
                    plan,
                    output_dir=output_dir or self.paths.local_runs_dir / "pichia_python" / "opn_service",
                    mu=mu,
                    production_ratio=production_ratio,
                    media_type=DEFAULT_OPN_MEDIA_TYPE,
                )
            else:
                engine = MatlabPichiaEngine(self.paths, self.matlab, self.runner)
                engine_result = engine.run_target_smoke(
                    PichiaSimulationRequest(
                        target_id="OPN",
                        candidate_id=candidate_id,
                        mu=mu,
                        production_ratio=production_ratio,
                        media_type=DEFAULT_OPN_MEDIA_TYPE,
                        timeout_seconds=timeout_seconds,
                    )
                )
            return OpnSimulationResult(
                success=engine_result.success,
                mu=engine_result.mu,
                candidate_id=candidate_id,
                production_ratio=engine_result.production_ratio,
                media_type=engine_result.media_type,
                message=engine_result.message,
                lp_file=engine_result.lp_file,
                output_file=engine_result.output_file,
                objective_value=engine_result.objective_value,
                command_output=engine_result.command_output,
            )
        finally:
            self._lock.release()

    def latest_candidate_result(self, candidate_id: str | None = None) -> tuple[Path | None, object | None]:
        pattern = "*.lp.float.out" if candidate_id is None else f"*{candidate_id}*.lp.float.out"
        outputs = list(self.paths.opn_run_dir.glob(pattern)) if self.paths.opn_run_dir.exists() else []
        if not outputs:
            return None, None
        latest = max(outputs, key=lambda path: path.stat().st_mtime)
        return latest, parse_soplex_file(latest)


def _builtin_opn_plan(candidate_id: str):
    generic = target_input_set_from_opn_input_set(BuiltinOpnInputProvider().load_input_set())
    leader = next((item for item in generic.leaders if item.candidate_id == candidate_id), None)
    if leader is None:
        raise KeyError(f"当前 Python engine 只支持内置 OPN 候选，未找到：{candidate_id}")
    return build_target_protein_plan(generic.target, leader)


def latest_opn_objectives(paths: ProjectPaths) -> dict[str, tuple[float | None, str | None, bool, Path | None]]:
    candidates = OpnCandidateCatalog(paths).list_candidates()
    result: dict[str, tuple[float | None, str | None, bool, Path | None]] = {}
    if not paths.opn_run_dir.exists():
        return result
    for candidate in candidates:
        outputs = list(paths.opn_run_dir.glob(f"*{candidate.candidate_id}*.lp.float.out"))
        if not outputs:
            result[candidate.candidate_id] = (None, None, False, None)
            continue
        latest = max(outputs, key=lambda path: path.stat().st_mtime)
        summary = parse_soplex_file(latest)
        value = None
        if summary.objective_value is not None:
            try:
                value = float(summary.objective_value)
            except ValueError:
                value = None
        result[candidate.candidate_id] = (value, summary.objective_value, summary.optimal, latest)
    return result
