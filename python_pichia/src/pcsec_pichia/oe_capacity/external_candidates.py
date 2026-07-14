"""Compatibility facade for external OE-capacity candidate workflows.

New code should import the responsibility-specific modules directly.  This
module preserves the original public API for existing callers.
"""

from pcsec_pichia.external_refs.capacity_sources import (
    ExternalCapacitySource,
    ExternalCapacitySourceType,
    RetrievalMode,
    cache_uniprot_identity_source,
)
from pcsec_pichia.oe_capacity.external_candidate_evaluation import (
    build_capacity_candidate,
    build_capacity_model_binding,
)
from pcsec_pichia.oe_capacity.external_candidate_io import (
    ExternalCapacityCandidateOutputs,
    import_capacity_measurements,
    load_external_capacity_candidate_bundle,
    write_external_capacity_candidate_cache,
)
from pcsec_pichia.oe_capacity.external_candidate_promotion import (
    build_capacity_promotion_manifest,
    promote_capacity_candidates,
)
from pcsec_pichia.oe_capacity.external_candidate_schema import (
    CANDIDATE_MANIFEST_FILENAME,
    CANDIDATE_RECORDS_FILENAME,
    RAW_MEASUREMENTS_FILENAME,
    CapacityApplicabilityScope,
    CapacityCandidateStatus,
    CapacityConfidence,
    CapacityConversionStep,
    CapacityModelBinding,
    CapacityParameterKind,
    CapacityPromotionManifest,
    ExternalCapacityCandidate,
    ExternalCapacityCandidateBundle,
    HostCondition,
    PromotionDecision,
    RawCapacityMeasurement,
    validate_external_capacity_candidate_bundle,
)


__all__ = [
    "CANDIDATE_MANIFEST_FILENAME",
    "CANDIDATE_RECORDS_FILENAME",
    "CapacityApplicabilityScope",
    "CapacityCandidateStatus",
    "CapacityConfidence",
    "CapacityConversionStep",
    "CapacityModelBinding",
    "CapacityParameterKind",
    "CapacityPromotionManifest",
    "ExternalCapacityCandidate",
    "ExternalCapacityCandidateBundle",
    "ExternalCapacityCandidateOutputs",
    "ExternalCapacitySource",
    "ExternalCapacitySourceType",
    "HostCondition",
    "PromotionDecision",
    "RAW_MEASUREMENTS_FILENAME",
    "RawCapacityMeasurement",
    "RetrievalMode",
    "build_capacity_candidate",
    "build_capacity_model_binding",
    "build_capacity_promotion_manifest",
    "cache_uniprot_identity_source",
    "import_capacity_measurements",
    "load_external_capacity_candidate_bundle",
    "promote_capacity_candidates",
    "validate_external_capacity_candidate_bundle",
    "write_external_capacity_candidate_cache",
]
