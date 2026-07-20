from __future__ import annotations

import pytest

from pcsec_pichia.experimental_feedback import (
    ConditionContext,
    ExperimentBundle,
    ExperimentRecord,
    HostContext,
    InterventionRecord,
    InterventionType,
    PredictionLinkStatus,
    SchemaValidationError,
    build_prediction_index,
    link_experiments_to_predictions,
)


def test_exact_prediction_link_requires_target_gene_type_run_and_evidence() -> None:
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "screen-run-001",
                "evidence_items": [
                    {
                        "evidence_id": "hLF-KO-0001",
                        "target_id": "hLF",
                        "gene_id": "PAS_chr1-1_0001",
                        "intervention_type": "KO",
                        "context_id": "methanol-bmmy-72h",
                        "rank": 1,
                        "recommendation_tier": "model_supported",
                    }
                ],
            },
        )
    )
    experiment = ExperimentRecord(
        experiment_id="HLF-LINK-1",
        target_id="hLF",
        host=HostContext("Komagataella phaffii", "X33", "X33"),
        batch_id="B01",
        condition=ConditionContext("BMMY, methanol, shake_flask, 250 rpm", 72.0),
        context_id="methanol-bmmy-72h",
    )
    intervention = InterventionRecord(
        experiment_id=experiment.experiment_id,
        intervention_id="KO-1",
        component_index=1,
        intervention_type=InterventionType.KO,
        gene_id="PAS_chr1-1_0001",
        construction_method="CRISPR-Cas9",
        prediction_run_id="screen-run-001",
        evidence_id="hLF-KO-0001",
    )

    result = link_experiments_to_predictions(
        ExperimentBundle(experiments=(experiment,), interventions=(intervention,)),
        index,
    )

    assert result.matched_count == 1
    assert result.links[0].status is PredictionLinkStatus.MATCHED
    assert result.links[0].evidence_id == "hLF-KO-0001"


def test_linkage_preserves_ambiguous_missing_and_context_mismatch_states() -> None:
    index = build_prediction_index(
        (
            {
                "prediction_run_id": "screen-run-002",
                "evidence_items": [
                    _prediction("OPN-OE-1", "OPN", "G2", "OE", "ctx-other"),
                    _prediction("HLF-KO-A", "hLF", "G3", "KO", "ctx-hlf"),
                    _prediction("HLF-KO-B", "hLF", "G3", "KO", "ctx-hlf"),
                    _prediction("HLF-KO-COMMON", "hLF", "G5", "KO", "ctx-hlf", common_name="KAR2"),
                ],
            },
        )
    )
    host = HostContext("Komagataella phaffii", "X33", "X33")
    condition = ConditionContext("BMMY, methanol, shake_flask, 250 rpm", 72.0)
    experiments = (
        ExperimentRecord("HLF-A", "hLF", host, "B1", condition, context_id="ctx-hlf"),
        ExperimentRecord("OPN-A", "OPN", host, "B2", condition, context_id="ctx-opn"),
        ExperimentRecord("HLF-TARGET", "hLF", host, "B3", condition, context_id="ctx-other"),
    )
    interventions = (
        _ko("HLF-A", "AMBIG", "G3"),
        _ko("HLF-A", "MISSING", "G4"),
        _ko("HLF-A", "COMMON", "", common_name="KAR2"),
        InterventionRecord(
            "OPN-A",
            "CONTEXT",
            1,
            InterventionType.OE,
            gene_id="G2",
            construct_id="construct",
            promoter="promoter",
            induction_mode="constitutive",
            prediction_run_id="screen-run-002",
            warnings=("copy_number_unknown",),
        ),
        InterventionRecord(
            "HLF-TARGET",
            "TARGET",
            1,
            InterventionType.OE,
            gene_id="G2",
            construct_id="construct",
            promoter="promoter",
            induction_mode="constitutive",
            prediction_run_id="screen-run-002",
            warnings=("copy_number_unknown",),
        ),
    )

    result = link_experiments_to_predictions(
        ExperimentBundle(experiments=experiments, interventions=interventions),
        index,
    )

    assert result.ambiguous_count == 2
    assert result.missing_prediction_count == 1
    assert result.context_mismatch_count == 2
    assert result.matched_count == 0
    assert next(link for link in result.links if link.intervention_id == "COMMON").reason == "common_name_only"
    assert [record.rank for record in index.records] == [None, None, None, None]


def test_prediction_index_rejects_duplicate_run_and_evidence_identity() -> None:
    duplicate = _prediction("DUPLICATE", "hLF", "G1", "KO", "ctx")

    with pytest.raises(SchemaValidationError, match="duplicate prediction identity"):
        build_prediction_index(
            (
                {
                    "prediction_run_id": "screen-run-duplicate",
                    "evidence_items": [duplicate, duplicate],
                },
            )
        )


def _prediction(
    evidence_id: str,
    target_id: str,
    gene_id: str,
    intervention_type: str,
    context_id: str,
    *,
    common_name: str = "",
) -> dict[str, object]:
    return {
        "evidence_id": evidence_id,
        "target_id": target_id,
        "gene_id": gene_id,
        "intervention_type": intervention_type,
        "context_id": context_id,
        "common_name": common_name,
    }


def _ko(experiment_id: str, intervention_id: str, gene_id: str, *, common_name: str = "") -> InterventionRecord:
    return InterventionRecord(
        experiment_id,
        intervention_id,
        1,
        InterventionType.KO,
        gene_id=gene_id,
        common_name=common_name,
        construction_method="CRISPR-Cas9",
        prediction_run_id="screen-run-002",
    )
