from __future__ import annotations

import argparse
import json
from pathlib import Path

from pcsec_pichia.experimental_feedback import (
    ConditionContext,
    ExperimentBundle,
    ExperimentRecord,
    HostContext,
    InterventionRecord,
    InterventionType,
    validate_experiment_bundle,
)


def build_round0_smoke_bundle() -> ExperimentBundle:
    condition = ConditionContext(
        medium="sanitized_defined_medium",
        carbon_source="methanol",
        culture_mode="shake_flask",
        temperature_c=30.0,
        ph=6.0,
        oxygen_or_agitation="sanitized agitation setting",
        sampling_time_h=72.0,
    )
    host = HostContext(
        species="Komagataella phaffii",
        strain="sanitized-strain",
        parent_strain="sanitized-parent",
    )
    experiments = (
        ExperimentRecord("SMOKE-HLF-1", "hLF", host, "SMOKE-B1", condition),
        ExperimentRecord("SMOKE-OPN-1", "OPN", host, "SMOKE-B2", condition),
    )
    interventions = (
        InterventionRecord("SMOKE-HLF-1", "CONTROL-1", 1, InterventionType.CONTROL),
        InterventionRecord(
            "SMOKE-OPN-1",
            "OE-1",
            1,
            InterventionType.OE,
            gene_id="PAS_chr2-1_0140",
            construct_id="sanitized-construct",
            promoter="sanitized-promoter",
            induction_mode="constitutive",
            warnings=("copy_number_unknown",),
        ),
    )
    return ExperimentBundle(experiments=experiments, interventions=interventions)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("local_runs/experiment_feedback/round0_contract"),
    )
    args = parser.parse_args()
    result = validate_experiment_bundle(build_round0_smoke_bundle())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "contract_smoke_summary.json"
    output_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "is_valid": result.is_valid,
                "targets": [record.target_id for record in result.bundle.experiments],
                "error_count": len(result.errors),
                "warning_count": len(result.warnings),
                "uses_sanitized_fixture": True,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(output_path)
    return 0 if result.is_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
