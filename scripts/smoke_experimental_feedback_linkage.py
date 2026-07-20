from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from enum import Enum
from pathlib import Path

from pcsec_pichia.experimental_feedback import (
    ConditionContext,
    ExperimentBundle,
    ExperimentRecord,
    HostContext,
    InterventionRecord,
    InterventionType,
    build_prediction_index,
    link_experiments_to_predictions,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("local_runs/experiment_feedback/round2_linkage"),
    )
    args = parser.parse_args()
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "sanitized-screen-run",
                "evidence_items": [
                    _prediction("hLF-KO-0001", "hLF", "G1", "KO", "ctx-hlf"),
                    _prediction("OPN-OE-0001", "OPN", "G2", "OE", "ctx-opn-prediction"),
                    _prediction("hLF-KO-0002", "hLF", "G3", "KO", "ctx-hlf"),
                    _prediction("hLF-KO-0003", "hLF", "G3", "KO", "ctx-hlf"),
                ],
            },
        )
    )
    host = HostContext("Komagataella phaffii", "sanitized-strain", "sanitized-parent")
    condition = ConditionContext("sanitized-medium, methanol, shake_flask, sanitized agitation", 72.0)
    bundle = ExperimentBundle(
        experiments=(
            ExperimentRecord("LINK-HLF", "hLF", host, "B1", condition, context_id="ctx-hlf"),
            ExperimentRecord("LINK-OPN", "OPN", host, "B2", condition, context_id="ctx-opn-experiment"),
        ),
        interventions=(
            _ko("LINK-HLF", "MATCHED", "G1", evidence_id="hLF-KO-0001"),
            _ko("LINK-HLF", "AMBIGUOUS", "G3"),
            _ko("LINK-HLF", "MISSING", "G4"),
            InterventionRecord(
                "LINK-OPN",
                "CONTEXT",
                1,
                InterventionType.OE,
                gene_id="G2",
                construct_id="sanitized-construct",
                promoter="sanitized-promoter",
                induction_mode="constitutive",
                prediction_run_id="sanitized-screen-run",
                warnings=("copy_number_unknown",),
            ),
        ),
    )
    result = link_experiments_to_predictions(bundle, index)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "linkage_summary.json"
    output_path.write_text(
        json.dumps(
            {
                "matched_count": result.matched_count,
                "ambiguous_count": result.ambiguous_count,
                "missing_prediction_count": result.missing_prediction_count,
                "context_mismatch_count": result.context_mismatch_count,
                "links": [_json_ready(asdict(link)) for link in result.links],
                "uses_sanitized_fixture": True,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(output_path)
    return 0 if (result.matched_count, result.ambiguous_count, result.missing_prediction_count, result.context_mismatch_count) == (1, 1, 1, 1) else 1


def _prediction(evidence_id: str, target_id: str, gene_id: str, intervention_type: str, context_id: str) -> dict[str, object]:
    return {
        "evidence_id": evidence_id,
        "target_id": target_id,
        "gene_id": gene_id,
        "intervention_type": intervention_type,
        "context_id": context_id,
        "secretion_ratio_vs_wildtype": 1.1,
        "recommendation_tier": "model_supported",
    }


def _ko(experiment_id: str, intervention_id: str, gene_id: str, *, evidence_id: str = "") -> InterventionRecord:
    return InterventionRecord(
        experiment_id,
        intervention_id,
        1,
        InterventionType.KO,
        gene_id=gene_id,
        construction_method="CRISPR-Cas9",
        prediction_run_id="sanitized-screen-run",
        evidence_id=evidence_id,
    )


def _json_ready(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


if __name__ == "__main__":
    raise SystemExit(main())
