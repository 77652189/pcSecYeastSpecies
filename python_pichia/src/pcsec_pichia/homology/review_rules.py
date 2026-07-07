from __future__ import annotations


MODEL_READY_RBH_HIGH_CONFIDENCE = "model_ready_rbh_high_confidence"
RBH_NOT_IN_MODEL = "rbh_not_in_model"
LOW_IDENTITY_REVIEW_REQUIRED = "low_identity_review_required"
COVERAGE_REVIEW_REQUIRED = "coverage_review_required"
PARALOG_RISK_REVIEW_REQUIRED = "paralog_risk_review_required"
NO_RECIPROCAL_HIT = "no_reciprocal_hit"
UNRESOLVED_QUERY_SYMBOL = "unresolved_query_symbol"
MANUAL_REVIEW_REQUIRED = "manual_review_required"

NAME_CONFIRMED_BY_RBH = "name_confirmed_by_rbh"
ALIAS_CONFIRMED_BY_RBH = "alias_confirmed_by_rbh"
SEQUENCE_NAME_CONFLICT = "sequence_name_conflict"
EXTERNAL_NAME_MISSING = "external_name_missing"
INTERNAL_NAME_MISSING = "internal_name_missing"


def classify_homology_review_status(
    *,
    is_rbh: bool,
    identity_pct: float | None,
    query_coverage: float | None,
    subject_coverage: float | None,
    in_model_gene_index: bool,
    min_identity: float = 30.0,
    min_coverage: float = 50.0,
    paralog_count: int = 0,
    unresolved_query: bool = False,
) -> tuple[str, tuple[str, ...]]:
    """Classify homolog evidence without turning it into phenotype evidence."""

    warnings: list[str] = []
    if unresolved_query:
        return UNRESOLVED_QUERY_SYMBOL, ("query symbol could not be resolved to an SCE sequence",)
    if not is_rbh:
        warnings.append("forward best hit is not reciprocal best hit")
        return NO_RECIPROCAL_HIT, tuple(warnings)
    if paralog_count > 1:
        warnings.append(f"{paralog_count} close forward hits indicate paralog risk")
        return PARALOG_RISK_REVIEW_REQUIRED, tuple(warnings)
    if identity_pct is None or identity_pct < min_identity:
        warnings.append(f"identity below threshold {min_identity}")
        return LOW_IDENTITY_REVIEW_REQUIRED, tuple(warnings)
    if (
        query_coverage is None
        or subject_coverage is None
        or query_coverage < min_coverage
        or subject_coverage < min_coverage
    ):
        warnings.append(f"coverage below threshold {min_coverage}")
        return COVERAGE_REVIEW_REQUIRED, tuple(warnings)
    if not in_model_gene_index:
        warnings.append("RBH candidate is not present in current Pichia GEM gene_index")
        return RBH_NOT_IN_MODEL, tuple(warnings)
    return MODEL_READY_RBH_HIGH_CONFIDENCE, ()


def classify_name_consistency(
    *,
    internal_common_name: str,
    external_gene_name: str | None,
    external_aliases: tuple[str, ...] = (),
    is_rbh: bool,
) -> str:
    """Classify name agreement separately from sequence-level homology."""

    internal = _normalize_name(internal_common_name)
    external = _normalize_name(external_gene_name or "")
    aliases = {_normalize_name(alias) for alias in external_aliases if alias}
    if not internal:
        return INTERNAL_NAME_MISSING
    if not external and not aliases:
        return EXTERNAL_NAME_MISSING
    if not is_rbh:
        return PARALOG_RISK_REVIEW_REQUIRED
    if external and (internal == external or internal in _split_name_tokens(external)):
        return NAME_CONFIRMED_BY_RBH
    if internal in aliases:
        return ALIAS_CONFIRMED_BY_RBH
    return SEQUENCE_NAME_CONFLICT


def _normalize_name(value: str) -> str:
    return value.strip().upper().replace(" ", "").replace("_", "").replace("-", "")


def _split_name_tokens(value: str) -> set[str]:
    normalized = value.replace("/", " ").replace(",", " ").replace(";", " ")
    return {_normalize_name(token) for token in normalized.split() if token}
