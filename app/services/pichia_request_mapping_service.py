from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from pcsec_pichia.core.paths import ProjectPaths

from app.services.pichia_secretion_schema import (
    NormalizationMode,
    SecretionRunRequest,
    TerminalStopPolicy,
)

_ALLOWED_SEQUENCE_SYMBOLS = set("ACDEFGHIKLMNPQRSTVWY*")


def resolve_output_dir(request: SecretionRunRequest, paths: ProjectPaths) -> Path:
    if request.output_dir is not None:
        return request.output_dir
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_target = _safe_id(request.target_id or "target")
    return paths.local_runs_dir / "streamlit_pichia_runs" / f"{stamp}_{safe_target}"


def engine_target_id(request: SecretionRunRequest) -> str:
    if request.target_source == "builtin":
        aliases = {
            "OPN": "OPN_ALPHA_FULL_PROJECT",
            "OPN_ALPHA_FULL_PROJECT": "OPN_ALPHA_FULL_PROJECT",
            "hLF": "hLF",
        }
        return aliases.get(request.target_id, request.target_id)
    return request.target_id or _safe_id(request.target_name or "custom_target")


def target_input_payload(request: SecretionRunRequest) -> Any | None:
    if request.target_source == "builtin":
        return None
    if request.target_source == "custom_json":
        return request.custom_json_path
    sequence = _normalize_sequence(request.sequence, request.normalization_mode, request.terminal_stop_policy)
    if not sequence:
        raise ValueError("选择自定义序列时必须填写蛋白序列。")
    target_id = engine_target_id(request)
    return {
        "target_id": target_id,
        "protein_id": target_id,
        "mature_sequence": sequence,
        "leader_sequence": _normalize_sequence(
            request.leader_sequence,
            request.normalization_mode,
            request.terminal_stop_policy,
        ),
        "signal_peptide_sequence": _normalize_sequence(
            request.signal_peptide_sequence,
            request.normalization_mode,
            request.terminal_stop_policy,
        ),
        "through_er": True,
        "localization": "e",
        "disulfide_sites": int(request.disulfide_sites),
        "n_glycosylation_sites": int(request.n_glycosylation_sites),
        "o_glycosylation_sites": int(request.o_glycosylation_sites),
        "transmembrane": 0,
        "gpi_sites": 0,
        "cotranslation": 0,
    }


def sequence_contract_for_engine(request: SecretionRunRequest) -> dict[str, Any]:
    if request.target_source != "custom_sequence":
        return {}
    mature_raw = _compact_sequence(request.sequence)
    leader_raw = _compact_sequence(request.leader_sequence)
    signal_raw = _compact_sequence(request.signal_peptide_sequence)
    mature_norm = _normalize_sequence(request.sequence, request.normalization_mode, request.terminal_stop_policy)
    leader_norm = _normalize_sequence(request.leader_sequence, request.normalization_mode, request.terminal_stop_policy)
    signal_norm = _normalize_sequence(
        request.signal_peptide_sequence,
        request.normalization_mode,
        request.terminal_stop_policy,
    )
    contains_signal = bool(signal_norm) if request.contains_signal_peptide is None else bool(request.contains_signal_peptide)
    contains_leader = bool(leader_norm) if request.contains_leader is None else bool(request.contains_leader)
    return {
        "sequence_role": request.sequence_role,
        "normalization_mode": request.normalization_mode,
        "contains_signal_peptide": contains_signal,
        "contains_leader": contains_leader,
        "terminal_stop_policy": request.terminal_stop_policy,
        "original_sequence_length": len(mature_raw),
        "normalized_sequence_length": len(mature_norm),
        "original_full_sequence_length": len(leader_raw) + len(mature_raw),
        "normalized_full_sequence_length": len(leader_norm) + len(mature_norm),
        "original_leader_sequence_length": len(leader_raw),
        "normalized_leader_sequence_length": len(leader_norm),
        "original_signal_peptide_length": len(signal_raw),
        "normalized_signal_peptide_length": len(signal_norm),
        "terminal_stop_present": mature_raw.endswith("*"),
        "terminal_stop_removed": mature_raw.endswith("*") and mature_norm == mature_raw[:-1],
    }


def request_warnings(request: SecretionRunRequest) -> list[str]:
    warnings: list[str] = []
    resolved_engine_target_id = engine_target_id(request)
    if request.target_source == "builtin" and resolved_engine_target_id == "hLF":
        warnings.append(
            "hLF 使用用户提供的 710aa 项目序列；当前 Python hLF 对应 MATLAB artifact target "
            "`hLF_PROJECT_710`，可报告为 `aligned_except_known_matlab_compatibility_differences`，但不是 fully aligned。"
        )
        warnings.append("旧 MATLAB `hLF` baseline 是 historical matlab_failed，不代表当前项目 hLF 710aa 失败。")
    if request.target_source == "custom_sequence":
        warnings.append("自定义序列需要手动提供信号肽、引导肽序列和 DSB/NG/OG 计数，模型不做智能推断。")
        if request.sequence_role == "unknown":
            warnings.append("序列角色为「未知」，建议明确声明前体蛋白(preprotein)或成熟分泌蛋白(mature_secreted)。")
        warnings.extend(_sequence_input_warnings(request))
        if not any((request.disulfide_sites, request.n_glycosylation_sites, request.o_glycosylation_sites)):
            warnings.append("当前 DSB/NG/OG 均为 0；如果目标蛋白存在二硫键或糖基化位点，需要显式填写。")
        if request.sequence.strip().endswith("*") and request.terminal_stop_policy == "strip":
            warnings.append("序列末尾包含终止符 *，已按 terminal_stop_policy=strip 移除后进入模型。")
        elif request.sequence.strip().endswith("*") and request.normalization_mode == "remove_terminal_stop":
            warnings.append("序列末尾包含终止符 *，已按「移除末尾终止符」规范化模式处理。")
        elif request.sequence.strip().endswith("*"):
            warnings.append("序列末尾包含终止符 *，如需移除请选择「移除末尾终止符」规范化模式。")
    if request.oe_gene_ids:
        warnings.append("过表达基因会解析到模型反应，按 reaction-level OE proxy 运行；这不是完整的基因表达调控模拟。")
    if request.growth_points != (0.10,):
        warnings.append("生长权衡使用指定的小网格生长速率，为快速 smoke 而非完整扫描。")
    return warnings


def _normalize_sequence(
    sequence: str,
    mode: NormalizationMode,
    terminal_stop_policy: TerminalStopPolicy = "allow_for_record_only",
) -> str:
    text = "".join(str(sequence or "").split()).upper()
    if terminal_stop_policy == "reject" and text.endswith("*"):
        raise ValueError("序列末尾包含终止符 *；terminal_stop_policy=reject 时不能进入模型。")
    if (mode == "remove_terminal_stop" or terminal_stop_policy == "strip") and text.endswith("*"):
        return text[:-1]
    return text


def _compact_sequence(sequence: str) -> str:
    return "".join(str(sequence or "").split()).upper()


def _sequence_input_warnings(request: SecretionRunRequest) -> list[str]:
    raw_parts = {
        "成熟蛋白序列": request.sequence,
        "引导肽序列": request.leader_sequence,
        "信号肽序列": request.signal_peptide_sequence,
    }
    warnings: list[str] = []
    for label, raw in raw_parts.items():
        text = str(raw or "")
        if not text:
            continue
        compact = _compact_sequence(text)
        if any(ch.isspace() for ch in text):
            warnings.append(f"{label}包含空白字符；进入模型前会移除空白并转为大写。")
        invalid = sorted({ch for ch in compact if ch not in _ALLOWED_SEQUENCE_SYMBOLS})
        if invalid:
            warnings.append(
                f"{label}包含非标准氨基酸字符：{''.join(invalid)}；当前不会智能修正，请人工确认序列。"
            )
    return warnings


def _safe_id(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value).strip("._")
    return cleaned or "target"


__all__ = [
    "engine_target_id",
    "request_warnings",
    "resolve_output_dir",
    "sequence_contract_for_engine",
    "target_input_payload",
]
