from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

from pcsec_pichia.experimental_feedback.schema import (
    ASSAY_TYPE_OD600,
    ASSAY_TYPE_TITER,
    ASSAY_TYPE_UPR,
    CANONICAL_UNITS,
    COMPARTMENT_EXTRACELLULAR,
    COMPARTMENT_INTRACELLULAR,
    COMPARTMENT_WHOLE_CULTURE,
    ConditionContext,
    ExperimentRecord,
    FermentationDataStatus,
    HostContext,
    InterventionRecord,
    InterventionType,
    MeasurementRecord,
    MeasurementStatus,
    ModificationConfirmationStatus,
    QualityStatus,
    SchemaValidationError,
)


FERMENTATION_TEMPLATE_ADAPTER_ID = "pcsec_pichia.fermentation_template.v2"


@dataclass(frozen=True)
class FermentationTemplateImport:
    records: tuple[tuple[str, object], ...]
    warnings: tuple[str, ...]
    metadata: Mapping[str, object]


_ALIASES: Mapping[str, tuple[str, ...]] = {
    "experiment_id": ("experiment_id", "实验编号"),
    "clone_id": ("clone_id", "克隆编号", "单克隆编号", "克隆id", "克隆ID"),
    "target_id": ("target_id", "目标蛋白", "目标蛋白ID"),
    "target_name": ("target_name", "目标蛋白名称"),
    "batch_id": ("batch_id", "批次", "发酵批次", "批次编号"),
    "context_id": ("context_id", "预测条件ID", "条件ID"),
    "host_species": ("host_species", "宿主物种", "物种"),
    "host_strain": ("host_strain", "发酵菌株", "宿主菌株", "菌株"),
    "parent_strain": ("parent_strain", "本底菌株", "亲本菌株"),
    "condition_description": ("condition_description", "发酵条件"),
    "sampling_time_h": ("sampling_time_h", "取样时间_h", "取样时间(h)", "取样时间"),
    "od600_value": ("od600_value", "72h-OD600", "OD600"),
    "upr_value": ("upr_value", "72h-UPR", "UPR"),
    "extracellular_titer_value": (
        "extracellular_titer_value",
        "72h-胞外产量mg/L",
        "胞外产量mg/L",
        "胞外产量",
    ),
    "extracellular_dilution_factor": (
        "extracellular_dilution_factor",
        "胞外ELISA稀释倍数",
        "胞外稀释倍数",
    ),
    "extracellular_above_range": ("extracellular_above_range", "胞外是否超标曲"),
    "intracellular_titer_value": (
        "intracellular_titer_value",
        "72h-胞内产量mg/L",
        "胞内产量mg/L",
        "胞内产量",
    ),
    "intracellular_dilution_factor": (
        "intracellular_dilution_factor",
        "胞内ELISA稀释倍数",
        "胞内稀释倍数",
    ),
    "intracellular_above_range": ("intracellular_above_range", "胞内是否超标曲"),
    "notes": ("notes", "备注"),
    "assay_method": ("assay_method", "检测方法"),
    "modification_plan": ("modification_plan", "改造方案", "改造方案（含对应基因）"),
    "intervention_type": ("intervention_type", "改造类型"),
    "gene_id": ("gene_id", "对应基因", "基因ID", "gene"),
    "data_status": ("data_status", "数据状态"),
    "parent_control_group_id": (
        "parent_control_group_id",
        "亲本对照组编号",
        "对照组编号",
    ),
    "replicate_id": ("replicate_id", "重复编号", "生物学重复编号"),
    "prediction_run_id": ("prediction_run_id", "预测Run编号", "预测run编号"),
    "evidence_id": ("evidence_id", "证据编号"),
    "construction_method": ("construction_method", "构建方法"),
    "construct_id": ("construct_id", "构建编号"),
    "promoter": ("promoter", "启动子"),
    "induction_mode": ("induction_mode", "诱导方式"),
    "copy_number": ("copy_number", "拷贝数"),
    "confirmation_status": ("confirmation_status", "改造确认"),
    "confirmation_method": ("confirmation_method", "确认方式", "改造确认方式"),
    "status_reason": ("status_reason", "状态原因"),
    "exclusion_reason": ("exclusion_reason", "排除原因"),
}

_STATUS_ALIASES = {
    "normal": FermentationDataStatus.NORMAL,
    "正常": FermentationDataStatus.NORMAL,
    "contamination": FermentationDataStatus.CONTAMINATION,
    "污染": FermentationDataStatus.CONTAMINATION,
    "culture_failed": FermentationDataStatus.CULTURE_FAILED,
    "培养失败": FermentationDataStatus.CULTURE_FAILED,
    "assay_failed": FermentationDataStatus.ASSAY_FAILED,
    "检测失败": FermentationDataStatus.ASSAY_FAILED,
    "other_excluded": FermentationDataStatus.OTHER_EXCLUDED,
    "其他排除": FermentationDataStatus.OTHER_EXCLUDED,
}

_CONFIRMATION_STATUS_ALIASES = {
    "": ModificationConfirmationStatus.NOT_CONFIRMED,
    "未确认": ModificationConfirmationStatus.NOT_CONFIRMED,
    "not_confirmed": ModificationConfirmationStatus.NOT_CONFIRMED,
    "已确认成功": ModificationConfirmationStatus.CONFIRMED_SUCCESS,
    "确认成功": ModificationConfirmationStatus.CONFIRMED_SUCCESS,
    "confirmed_success": ModificationConfirmationStatus.CONFIRMED_SUCCESS,
    "已确认失败": ModificationConfirmationStatus.CONFIRMED_FAILED,
    "确认失败": ModificationConfirmationStatus.CONFIRMED_FAILED,
    "confirmed_failed": ModificationConfirmationStatus.CONFIRMED_FAILED,
}

_TRUE_FLAG_VALUES = {"是", "true", "yes", "y", "1"}
_FALSE_FLAG_VALUES = {"否", "false", "no", "n", "0", ""}

_DEFAULT_ASSAY_METHODS = {
    ASSAY_TYPE_OD600: "OD600",
    ASSAY_TYPE_UPR: "UPR_marker",
    ASSAY_TYPE_TITER: "ELISA",
}


def map_fermentation_template_rows(
    rows: Iterable[tuple[int, Mapping[str, object]]],
    *,
    metadata: Mapping[str, object] | None = None,
    source_sheet: str = "",
) -> FermentationTemplateImport:
    supplied_metadata = _normalize_mapping(metadata or {})
    canonical_records: list[tuple[str, object]] = []
    warnings: list[str] = []
    seen_headers: set[str] = set()
    row_count = 0
    for row_number, raw_row in rows:
        if not any(value not in (None, "") for value in raw_row.values()):
            continue
        row_count += 1
        normalized_row, unknown_columns = _canonicalize_row(raw_row)
        seen_headers.update(normalized_row)
        warnings.extend(f"unmapped_template_column:{column}" for column in unknown_columns)
        merged = _merge_metadata(normalized_row, supplied_metadata, row_number=row_number)
        canonical_records.extend(
            _records_from_row(
                merged,
                raw_row=raw_row,
                row_number=row_number,
                source_sheet=source_sheet,
            )
        )
    if row_count == 0:
        raise SchemaValidationError("fermentation template contains no data rows.")
    required_template_fields = {
        "clone_id",
        "modification_plan",
    }
    missing_headers = sorted(required_template_fields - seen_headers - set(supplied_metadata))
    if missing_headers:
        raise SchemaValidationError(
            "fermentation template is missing required direction-1 fields: "
            + ", ".join(missing_headers)
        )
    return FermentationTemplateImport(
        records=tuple(canonical_records),
        warnings=tuple(dict.fromkeys(warnings)),
        metadata=dict(supplied_metadata),
    )


def merge_fermentation_template_metadata(
    *sources: Mapping[str, object] | None,
) -> dict[str, object]:
    merged: dict[str, object] = {}
    for source in sources:
        normalized = _normalize_mapping(source or {})
        for key, value in normalized.items():
            if key in merged and value not in (None, "") and merged[key] not in (None, ""):
                if _text(merged[key]) != _text(value):
                    raise SchemaValidationError(
                        f"fermentation import metadata conflict: {key}"
                    )
            elif value not in (None, ""):
                merged[key] = value
    return merged


def _records_from_row(
    row: Mapping[str, object],
    *,
    raw_row: Mapping[str, object],
    row_number: int,
    source_sheet: str,
) -> tuple[tuple[str, object], ...]:
    clone_id = _required_text(row, "clone_id", row_number)
    target_id = _required_text(row, "target_id", row_number)
    batch_id = _required_text(row, "batch_id", row_number)
    replicate_id = _text(row.get("replicate_id")) or f"row{row_number}"
    control_group_id = _text(row.get("parent_control_group_id"))
    raw_data_status = _text(row.get("data_status"))
    data_status = (
        _data_status(raw_data_status, row_number=row_number)
        if raw_data_status
        else FermentationDataStatus.NORMAL
    )
    intervention_type, gene_id, design_label = _intervention(row, row_number=row_number)
    condition_description = _text(row.get("condition_description")) or "unknown"
    experiment_id = _text(row.get("experiment_id")) or _generated_experiment_id(
        target_id=target_id,
        batch_id=batch_id,
        clone_id=clone_id,
        replicate_id=replicate_id,
        identity_context=(
            control_group_id,
            _text(row.get("context_id")),
            _text(row.get("host_species")),
            _text(row.get("host_strain")),
            _text(row.get("parent_strain")),
            condition_description,
        ),
    )
    sampling_time_h = _optional_float(row.get("sampling_time_h"), "sampling_time_h", row_number)
    condition = ConditionContext(
        condition_description=condition_description,
        sampling_time_h=72.0 if sampling_time_h is None else sampling_time_h,
    )
    quality_status, quality_reason = _quality_status(data_status, row)
    experiment = ExperimentRecord(
        experiment_id=experiment_id,
        target_id=target_id,
        target_name=_text(row.get("target_name")),
        host=HostContext(
            species=_text(row.get("host_species")) or "unknown",
            strain=_text(row.get("host_strain")) or "unknown",
            parent_strain=_text(row.get("parent_strain")) or "unknown",
        ),
        batch_id=batch_id,
        condition=condition,
        context_id=_text(row.get("context_id")),
        biological_replicate_id=replicate_id,
        clone_id=clone_id,
        parent_control_group_id=control_group_id,
        notes=_text(row.get("notes")),
        fermentation_data_status=data_status,
        quality_status=quality_status,
        quality_reason=quality_reason,
    )
    warnings: tuple[str, ...] = ()
    if intervention_type is InterventionType.OE and row.get("copy_number") in (None, ""):
        warnings = ("copy_number_unknown",)
    confirmation_status, confirmation_method = _confirmation(row, row_number=row_number)
    intervention = InterventionRecord(
        experiment_id=experiment_id,
        intervention_id="CONTROL-1" if intervention_type is InterventionType.CONTROL else "INTERVENTION-1",
        component_index=1,
        intervention_type=intervention_type,
        gene_id=gene_id,
        construction_method=(
            _text(row.get("construction_method")) or design_label or "unknown"
        ),
        construct_id=_text(row.get("construct_id")) or (
            design_label if intervention_type is InterventionType.OE else ""
        ),
        promoter=_text(row.get("promoter")) or (
            "unknown" if intervention_type is InterventionType.OE else ""
        ),
        induction_mode=_text(row.get("induction_mode")) or (
            "unknown" if intervention_type is InterventionType.OE else ""
        ),
        copy_number=_optional_float(row.get("copy_number"), "copy_number", row_number),
        prediction_run_id=_text(row.get("prediction_run_id")),
        evidence_id=_text(row.get("evidence_id")),
        design_label=design_label,
        confirmation_status=confirmation_status,
        confirmation_method=confirmation_method,
        warnings=warnings,
    )
    raw_fields_json = json.dumps(
        {str(key): _json_scalar(value) for key, value in raw_row.items()},
        ensure_ascii=False,
        sort_keys=True,
    )
    measurements = tuple(
        _build_measurement(
            experiment_id=experiment_id,
            assay_type=assay_type,
            compartment=compartment,
            value_key=value_key,
            dilution_key=dilution_key,
            above_range_key=above_range_key,
            row=row,
            row_number=row_number,
            data_status=data_status,
            raw_fields_json=raw_fields_json,
            source_sheet=source_sheet,
        )
        for assay_type, compartment, value_key, dilution_key, above_range_key in (
            (ASSAY_TYPE_OD600, COMPARTMENT_WHOLE_CULTURE, "od600_value", None, None),
            (ASSAY_TYPE_UPR, COMPARTMENT_INTRACELLULAR, "upr_value", None, None),
            (
                ASSAY_TYPE_TITER,
                COMPARTMENT_EXTRACELLULAR,
                "extracellular_titer_value",
                "extracellular_dilution_factor",
                "extracellular_above_range",
            ),
            (
                ASSAY_TYPE_TITER,
                COMPARTMENT_INTRACELLULAR,
                "intracellular_titer_value",
                "intracellular_dilution_factor",
                "intracellular_above_range",
            ),
        )
    )
    return (
        ("experiment", experiment),
        ("intervention", intervention),
    ) + tuple(("measurement", measurement) for measurement in measurements)


def _build_measurement(
    *,
    experiment_id: str,
    assay_type: str,
    compartment: str,
    value_key: str,
    dilution_key: str | None,
    above_range_key: str | None,
    row: Mapping[str, object],
    row_number: int,
    data_status: FermentationDataStatus,
    raw_fields_json: str,
    source_sheet: str,
) -> MeasurementRecord:
    canonical_unit = CANONICAL_UNITS[assay_type]
    raw_value = _optional_float(row.get(value_key), value_key, row_number)
    dilution_factor = (
        _optional_float(row.get(dilution_key), dilution_key, row_number)
        if dilution_key is not None
        else None
    )
    above_range = (
        _parse_bool_flag(row.get(above_range_key), field_name=above_range_key, row_number=row_number)
        if above_range_key is not None
        else False
    )
    measurement_status, excluded, status_reason = _measurement_status(
        data_status,
        raw_value=raw_value,
        above_range=above_range,
        row=row,
    )
    return MeasurementRecord(
        experiment_id=experiment_id,
        measurement_id=f"{assay_type}-{compartment}",
        assay_type=assay_type,
        assay_method=_text(row.get("assay_method")) or _DEFAULT_ASSAY_METHODS[assay_type],
        compartment=compartment,
        raw_value=raw_value,
        raw_unit=canonical_unit,
        canonical_value=(raw_value if measurement_status is MeasurementStatus.VALID else None),
        canonical_unit=canonical_unit,
        status=measurement_status,
        dilution_factor=dilution_factor,
        technical_replicate_id="",
        status_reason=status_reason,
        excluded=excluded,
        exclusion_reason=(status_reason if excluded else ""),
        source_row_number=row_number,
        source_sheet=source_sheet,
        raw_fields_json=raw_fields_json,
    )


def _canonicalize_row(
    raw_row: Mapping[str, object],
) -> tuple[dict[str, object], tuple[str, ...]]:
    lookup = {
        _normalize_header(alias): canonical
        for canonical, aliases in _ALIASES.items()
        for alias in aliases
    }
    normalized: dict[str, object] = {}
    unknown: list[str] = []
    for raw_key, value in raw_row.items():
        key = str(raw_key or "").strip()
        if not key:
            continue
        canonical = lookup.get(_normalize_header(key))
        if canonical is None:
            unknown.append(key)
            continue
        if canonical in normalized and normalized[canonical] not in (None, "") and value not in (None, ""):
            raise SchemaValidationError(f"duplicate fermentation template column alias: {key}")
        if value not in (None, "") or canonical not in normalized:
            normalized[canonical] = value
    return normalized, tuple(unknown)


def _normalize_mapping(values: Mapping[str, object]) -> dict[str, object]:
    normalized, unknown = _canonicalize_row(values)
    if unknown:
        raise SchemaValidationError(
            "unsupported fermentation import metadata: " + ", ".join(unknown)
        )
    return normalized


def _merge_metadata(
    row: Mapping[str, object],
    metadata: Mapping[str, object],
    *,
    row_number: int,
) -> dict[str, object]:
    merged = dict(row)
    for key, metadata_value in metadata.items():
        row_value = merged.get(key)
        if row_value not in (None, "") and metadata_value not in (None, ""):
            if _text(row_value) != _text(metadata_value):
                raise SchemaValidationError(
                    f"fermentation template metadata conflict at row {row_number}: {key}"
                )
            continue
        if metadata_value not in (None, ""):
            merged[key] = metadata_value
    return merged


def _intervention(
    row: Mapping[str, object],
    *,
    row_number: int,
) -> tuple[InterventionType, str, str]:
    plan = _required_text(row, "modification_plan", row_number)
    explicit_type = _text(row.get("intervention_type")).upper()
    explicit_gene = _text(row.get("gene_id"))
    if plan.strip().lower() in {"control", "parent_control", "亲本对照", "对照"}:
        return InterventionType.CONTROL, "", plan
    match = re.match(r"^\s*(KO|OE)\s*[:：/\-\s]\s*([A-Za-z0-9_.-]+)\s*$", plan, re.I)
    if match:
        parsed_type = match.group(1).upper()
        parsed_gene = match.group(2)
        if explicit_type and explicit_type != parsed_type:
            raise SchemaValidationError(
                f"modification type conflict at fermentation template row {row_number}."
            )
        if explicit_gene and explicit_gene != parsed_gene:
            raise SchemaValidationError(
                f"gene_id conflict at fermentation template row {row_number}."
            )
        return InterventionType(parsed_type), explicit_gene or parsed_gene, plan
    incomplete_match = re.match(r"^\s*(KO|OE)\s*[:：/\-]?\s*$", plan, re.I)
    if incomplete_match:
        parsed_type = incomplete_match.group(1).upper()
        if explicit_type and explicit_type != parsed_type:
            raise SchemaValidationError(
                f"modification type conflict at fermentation template row {row_number}."
            )
        return InterventionType(parsed_type), explicit_gene, plan
    if explicit_type in {"KO", "OE"}:
        return InterventionType(explicit_type), explicit_gene, plan
    raise SchemaValidationError(
        f"invalid modification_plan at fermentation template row {row_number}: {plan}"
    )


def _confirmation(
    row: Mapping[str, object], *, row_number: int
) -> tuple[ModificationConfirmationStatus, str]:
    raw_status = _text(row.get("confirmation_status"))
    status = _CONFIRMATION_STATUS_ALIASES.get(raw_status) or _CONFIRMATION_STATUS_ALIASES.get(
        raw_status.lower()
    )
    if status is None:
        raise SchemaValidationError(
            f"invalid confirmation_status at fermentation template row {row_number}: {raw_status}"
        )
    method = _text(row.get("confirmation_method"))
    if status is not ModificationConfirmationStatus.NOT_CONFIRMED and not method:
        raise SchemaValidationError(
            f"confirmation_method is required at fermentation template row {row_number} "
            f"once confirmation_status is {status.value}."
        )
    return status, method


def _parse_bool_flag(value: object, *, field_name: str, row_number: int) -> bool:
    normalized = _text(value).lower()
    if normalized in _TRUE_FLAG_VALUES:
        return True
    if normalized in _FALSE_FLAG_VALUES:
        return False
    raise SchemaValidationError(
        f"{field_name} at fermentation template row {row_number} must be a yes/no flag, got: {value!r}"
    )


def _data_status(value: object, *, row_number: int) -> FermentationDataStatus:
    normalized = _text(value)
    status = _STATUS_ALIASES.get(normalized) or _STATUS_ALIASES.get(normalized.lower())
    if status is None:
        raise SchemaValidationError(
            f"invalid data_status at fermentation template row {row_number}: {normalized or '<empty>'}"
        )
    return status


def _quality_status(
    status: FermentationDataStatus,
    row: Mapping[str, object],
) -> tuple[QualityStatus, str]:
    reason = _text(row.get("status_reason")) or _text(row.get("exclusion_reason"))
    if status is FermentationDataStatus.NORMAL:
        return QualityStatus.VALID, reason
    if status is FermentationDataStatus.ASSAY_FAILED:
        return QualityStatus.WARNING, reason or status.value
    return QualityStatus.EXCLUDED, reason or status.value


def _measurement_status(
    status: FermentationDataStatus,
    *,
    raw_value: float | None,
    above_range: bool,
    row: Mapping[str, object],
) -> tuple[MeasurementStatus, bool, str]:
    explicit_reason = _text(row.get("status_reason")) or _text(row.get("exclusion_reason"))
    if status is FermentationDataStatus.ASSAY_FAILED:
        return MeasurementStatus.ASSAY_FAILED, False, explicit_reason or status.value
    if status in {
        FermentationDataStatus.CONTAMINATION,
        FermentationDataStatus.CULTURE_FAILED,
        FermentationDataStatus.OTHER_EXCLUDED,
    }:
        return MeasurementStatus.EXCLUDED, True, explicit_reason or status.value
    if above_range:
        return (
            MeasurementStatus.ABOVE_RANGE,
            False,
            explicit_reason or "超出ELISA标准曲线可定量范围（现场标记超标曲）",
        )
    if raw_value is None:
        return MeasurementStatus.MISSING, False, explicit_reason or "measurement_value_missing"
    return MeasurementStatus.VALID, False, ""


def _generated_experiment_id(
    *,
    target_id: str,
    batch_id: str,
    clone_id: str,
    replicate_id: str,
    identity_context: Sequence[str] = (),
) -> str:
    readable = "-".join(_safe_token(item) for item in (target_id, batch_id, clone_id, replicate_id))
    digest = hashlib.sha256(
        "\0".join(
            (target_id, batch_id, clone_id, replicate_id, *identity_context)
        ).encode("utf-8")
    ).hexdigest()[:10]
    return f"FERM-{readable}-{digest}"


def _safe_token(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", value.strip()).strip(".-")
    return cleaned or "item"


def _required_text(row: Mapping[str, object], key: str, row_number: int) -> str:
    value = _text(row.get(key))
    if not value:
        raise SchemaValidationError(
            f"fermentation template row {row_number} requires {key}."
        )
    return value


def _optional_float(value: object, field_name: str, row_number: int) -> float | None:
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        raise SchemaValidationError(
            f"{field_name} at fermentation template row {row_number} must be numeric."
        )
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise SchemaValidationError(
            f"{field_name} at fermentation template row {row_number} must be numeric."
        ) from exc


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _normalize_header(value: str) -> str:
    return re.sub(r"[\s_\-]+", "", str(value).strip().casefold())


def _json_scalar(value: object) -> object:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


__all__ = [
    "FERMENTATION_TEMPLATE_ADAPTER_ID",
    "FermentationTemplateImport",
    "map_fermentation_template_rows",
    "merge_fermentation_template_metadata",
]
