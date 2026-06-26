from __future__ import annotations

from pathlib import Path

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.core.paths import ProjectPaths

from app.services.pichia_secretion_schema import BuiltinTargetTemplate


def list_builtin_target_templates(paths: ProjectPaths | None = None) -> list[BuiltinTargetTemplate]:
    """Return formal engine builtin targets for Streamlit/FastAPI facades.

    This is presentation metadata only. The simulation still resolves targets
    through ``pcsec_pichia.pipeline`` and does not use these labels for model
    construction.
    """
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    ensure_python_pichia_on_path()

    from pcsec_pichia.targets import list_supported_builtin_targets, load_builtin_targets

    specs_by_id = {target.target_id: target for target in load_builtin_targets(resolved_paths.repo_root)}
    templates: list[BuiltinTargetTemplate] = []
    for summary in list_supported_builtin_targets():
        spec = specs_by_id.get(summary.target_id)
        if spec is None:
            continue
        semantics = _builtin_target_semantics(summary.target_id)
        templates.append(
            BuiltinTargetTemplate(
                target_id=summary.target_id,
                label=_builtin_target_label(summary.target_id),
                parameter_status=summary.parameter_status,
                source=summary.source,
                leader_length=len(spec.leader_sequence),
                signal_peptide_length=len(spec.signal_peptide_sequence),
                full_sequence_length=len(spec.full_sequence),
                alignment_target_kind=semantics["alignment_target_kind"],
                sequence_role=semantics["sequence_role"],
                normalization_mode=semantics["normalization_mode"],
                mature_sequence_length=len(spec.mature_sequence),
                target_warning=semantics["target_warning"],
                note=_builtin_target_note(summary.target_id),
            )
        )
    return templates


def _builtin_target_label(target_id: str) -> str:
    labels = {
        "OPN_ALPHA_FULL_PROJECT": "OPN alpha-factor 项目基线",
        "OPN_PPA_DDDK18": "OPN Pichia DDDK18 短信号肽",
        "OPN_PPA_PASCHR3_0030": "OPN Pichia PAS_chr3_0030 短 leader",
        "OPN_ALPHA_PRE_ONLY": "OPN alpha-factor pre-only 对照",
        "OPN_NATIVE_SPP1": "OPN 天然 SPP1 信号肽参考",
        "OPN_OST1N23_ALPHA_PRO": "OPN Ost1N23 + alpha pro 工程参考",
        "OPN_PPA_EPX1_SA": "OPN Pichia EPX1 signal-anchor 备选",
        "hLF": "hLF 项目目标（人源天然信号肽 + 成熟 hLF）",
    }
    return labels.get(target_id, target_id)


def _builtin_target_note(target_id: str) -> str:
    if target_id.startswith("OPN_"):
        return "OPN leader 比较是模型约束下的相对 smoke 结果，不等于真实发酵产量。"
    if target_id == "hLF":
        return (
            "hLF 使用用户提供的 710aa native-signal 项目序列；对齐 artifact 为 hLF_PROJECT_710，"
            "状态可为 aligned_except_known_matlab_compatibility_differences，但不是旧 MATLAB hLF fully aligned。"
        )
    return ""


def _builtin_target_semantics(target_id: str) -> dict[str, str]:
    if target_id == "hLF":
        return {
            "alignment_target_kind": "project_defined_hLF",
            "sequence_role": "native_signal_plus_mature_hLF",
            "normalization_mode": "user_provided_as_provided",
            "target_warning": (
                "hLF 使用用户提供的 710aa 项目序列；当前对齐 artifact target 是 hLF_PROJECT_710，"
                "状态可为 aligned_except_known_matlab_compatibility_differences，但不是 fully aligned；"
                "旧 MATLAB hLF baseline 仍是 historical matlab_failed。"
            ),
        }
    if target_id.startswith("OPN_"):
        return {
            "alignment_target_kind": "opn_leader_candidate",
            "sequence_role": "mature_secreted_with_leader_candidate",
            "normalization_mode": "as_provided",
            "target_warning": "OPN corrected 结果不能标记为旧 MATLAB fully aligned。",
        }
    return {
        "alignment_target_kind": "custom",
        "sequence_role": "custom_user_sequence",
        "normalization_mode": "as_provided",
        "target_warning": "",
    }


__all__ = [
    "list_builtin_target_templates",
    "_builtin_target_semantics",
]
