from __future__ import annotations

from dataclasses import dataclass

from pcsec_pichia.secretory_resources.schema import (
    IN_SCOPE_CATEGORIES,
    ExecutionStatus,
    SecretoryResourceCatalog,
)


# Round 0 forbids the same fallback values ADR-001 forbids for OE capacity:
# no generic upper bound, no optimal-flux stand-in, no fixed 1.0, no fixture
# literal masquerading as a real model handle.
_FORBIDDEN_HANDLE_LITERALS = {"1000", "1.0", "1", "fixture", "optimal_flux", "default"}


@dataclass(frozen=True)
class SecretoryResourceValidationIssue:
    resource_id: str
    code: str
    message: str


@dataclass(frozen=True)
class SecretoryResourceValidationResult:
    catalog: SecretoryResourceCatalog
    issues: tuple[SecretoryResourceValidationIssue, ...]

    @property
    def is_valid(self) -> bool:
        return not self.issues


def validate_secretory_resource_catalog(catalog: SecretoryResourceCatalog) -> SecretoryResourceValidationResult:
    """Schema-level contract checks plus the Round 0 governance gates.

    Raises SecretoryResourceValidationError (via catalog.validate()) for
    structural contract violations. Softer governance concerns that don't
    corrupt the data model are returned as issues instead of raised, so a
    caller can inspect every problem in one pass.
    """

    catalog.validate()

    issues: list[SecretoryResourceValidationIssue] = []
    for resource in catalog.resources:
        if resource.category not in IN_SCOPE_CATEGORIES:
            issues.append(
                SecretoryResourceValidationIssue(
                    resource_id=resource.resource_id,
                    code="category_out_of_round0_scope",
                    message=(
                        f"{resource.category.value} is not authorized for this Round 0 sub-round; "
                        "it must stay deferred until the next sub-round."
                    ),
                )
            )
        if resource.status is ExecutionStatus.EVIDENCE_ONLY and resource.calibration_mode is not None:
            issues.append(
                SecretoryResourceValidationIssue(
                    resource_id=resource.resource_id,
                    code="evidence_only_must_not_carry_calibration_mode",
                    message="evidence_only must never be silently promoted toward an executable calibration mode.",
                )
            )
        for handle in resource.model_handles:
            if handle.strip().lower() in _FORBIDDEN_HANDLE_LITERALS:
                issues.append(
                    SecretoryResourceValidationIssue(
                        resource_id=resource.resource_id,
                        code="forbidden_handle_literal",
                        message=f"model_handle {handle!r} looks like a placeholder/fallback, not a real model handle.",
                    )
                )
    return SecretoryResourceValidationResult(catalog=catalog, issues=tuple(issues))


__all__ = [
    "SecretoryResourceValidationIssue",
    "SecretoryResourceValidationResult",
    "validate_secretory_resource_catalog",
]
