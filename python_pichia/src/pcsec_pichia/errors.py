"""Shared exception contracts that must not create package-layer cycles."""


class OECapacityError(RuntimeError):
    """Base error for gene-level OE capacity workflows."""


class OECapacityValidationError(OECapacityError, ValueError):
    """Raised when an OE capacity contract violates a frozen invariant."""


__all__ = ["OECapacityError", "OECapacityValidationError"]
