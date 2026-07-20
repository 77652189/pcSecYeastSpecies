from __future__ import annotations

from dataclasses import dataclass

from pcsec_pichia.secretory_resources.schema import (
    ExecutionStatus,
    ResourceCategory,
    SecretoryResourceCatalog,
)


@dataclass(frozen=True)
class SecretoryResourcePlanEntry:
    resource_id: str
    category: ResourceCategory
    status: ExecutionStatus
    action: str
    reason: str = ""


@dataclass(frozen=True)
class SecretoryResourcePlan:
    target_id: str
    backend: str
    entries: tuple[SecretoryResourcePlanEntry, ...]


def plan_secretory_resource_constraints(catalog: SecretoryResourceCatalog) -> SecretoryResourcePlan:
    """Backend-neutral plan. Never calls a solver, never mutates model assets.

    Round 0 only records which resources *would* be relative-comparison
    candidates for a future solving round, and which stay unplanned and why.
    """

    entries: list[SecretoryResourcePlanEntry] = []
    for resource in catalog.resources:
        if resource.status is ExecutionStatus.EXECUTABLE:
            entries.append(
                SecretoryResourcePlanEntry(
                    resource_id=resource.resource_id,
                    category=resource.category,
                    status=resource.status,
                    action="relative_comparison_ready",
                    reason="real model handles identified; no solve has run this round.",
                )
            )
        else:
            reason = resource.limitations[0] if resource.limitations else resource.status.value
            entries.append(
                SecretoryResourcePlanEntry(
                    resource_id=resource.resource_id,
                    category=resource.category,
                    status=resource.status,
                    action="not_planned",
                    reason=reason,
                )
            )
    return SecretoryResourcePlan(target_id=catalog.target_id, backend="none", entries=tuple(entries))


@dataclass(frozen=True)
class SecretoryResourceCoverageSummary:
    target_id: str
    feature_enabled: bool
    total: int
    executable: int
    evidence_only: int
    unavailable: int
    not_applicable: int
    conflict: int
    unavailable_or_conflict_explanations: tuple[str, ...]


def summarize_secretory_resource_coverage(catalog: SecretoryResourceCatalog) -> SecretoryResourceCoverageSummary:
    counts = {status: 0 for status in ExecutionStatus}
    explanations: list[str] = []
    for resource in catalog.resources:
        counts[resource.status] += 1
        if resource.status in (ExecutionStatus.UNAVAILABLE, ExecutionStatus.CONFLICT):
            if resource.limitations:
                reason = resource.limitations[0]
            elif resource.conflicts:
                reason = resource.conflicts[0]
            else:
                reason = resource.status.value
            explanations.append(f"{resource.resource_id}: {reason}")
    return SecretoryResourceCoverageSummary(
        target_id=catalog.target_id,
        feature_enabled=catalog.feature_enabled,
        total=len(catalog.resources),
        executable=counts[ExecutionStatus.EXECUTABLE],
        evidence_only=counts[ExecutionStatus.EVIDENCE_ONLY],
        unavailable=counts[ExecutionStatus.UNAVAILABLE],
        not_applicable=counts[ExecutionStatus.NOT_APPLICABLE],
        conflict=counts[ExecutionStatus.CONFLICT],
        unavailable_or_conflict_explanations=tuple(explanations),
    )


__all__ = [
    "SecretoryResourceCoverageSummary",
    "SecretoryResourcePlan",
    "SecretoryResourcePlanEntry",
    "plan_secretory_resource_constraints",
    "summarize_secretory_resource_coverage",
]
