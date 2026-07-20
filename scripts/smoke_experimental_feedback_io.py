from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from enum import Enum
from pathlib import Path

from pcsec_pichia.experimental_feedback import (
    ConditionContext,
    ExperimentRecord,
    HostContext,
    InterventionRecord,
    InterventionType,
    MeasurementRecord,
    MeasurementStatus,
    load_experiment_bundle,
    validate_experiment_bundle,
    write_experiment_feedback_cache,
)


def _records() -> tuple[tuple[str, object], ...]:
    condition = ConditionContext("sanitized defined medium, methanol, shake flask, sanitized agitation setting", 72.0)
    host = HostContext("Komagataella phaffii", "sanitized-strain", "sanitized-parent")
    return (
        ("experiment", ExperimentRecord("IO-HLF-1", "hLF", host, "IO-B1", condition)),
        ("experiment", ExperimentRecord("IO-OPN-1", "OPN", host, "IO-B2", condition)),
        ("intervention", InterventionRecord("IO-HLF-1", "CONTROL-1", 1, InterventionType.CONTROL)),
        (
            "intervention",
            InterventionRecord(
                "IO-OPN-1",
                "OE-1",
                1,
                InterventionType.OE,
                gene_id="PAS_chr2-1_0140",
                construct_id="sanitized-construct",
                promoter="sanitized-promoter",
                induction_mode="constitutive",
                warnings=("copy_number_unknown",),
            ),
        ),
        (
            "measurement",
            MeasurementRecord(
                "IO-HLF-1",
                "TITER-T1",
                "titer",
                "sanitized-assay",
                "extracellular",
                10.0,
                "mg/L",
                10.0,
                "mg/L",
                MeasurementStatus.VALID,
            ),
        ),
        (
            "measurement",
            MeasurementRecord(
                "IO-OPN-1",
                "TITER-T1",
                "titer",
                "sanitized-assay",
                "extracellular",
                None,
                "mg/L",
                None,
                "mg/L",
                MeasurementStatus.ASSAY_FAILED,
                status_reason="sanitized assay failure",
            ),
        ),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("local_runs/experiment_feedback/round1_io"),
    )
    args = parser.parse_args()
    inbox = args.output_root / "inbox"
    inbox.mkdir(parents=True, exist_ok=True)
    csv_path = inbox / "sanitized_bundle.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["record_type", "payload_json"])
        writer.writeheader()
        for record_type, record in _records():
            writer.writerow(
                {
                    "record_type": record_type,
                    "payload_json": json.dumps(asdict(record), default=_enum_value),
                }
            )
    bundle = load_experiment_bundle(csv_path)
    validation = validate_experiment_bundle(bundle)
    outputs = write_experiment_feedback_cache(bundle, args.output_root / "validated")
    roundtrip = load_experiment_bundle(outputs.validated_records_path)
    targets = sorted(record.target_id for record in roundtrip.experiments)
    if not validation.is_valid or targets != ["OPN", "hLF"]:
        return 1
    print(outputs.manifest_path)
    return 0


def _enum_value(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    raise TypeError(type(value).__name__)


if __name__ == "__main__":
    raise SystemExit(main())
