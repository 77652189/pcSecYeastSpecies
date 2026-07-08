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
PICHIA_LOCUS_CONFIRMED_BY_RBH = "pichia_locus_confirmed_by_rbh"
SEQUENCE_NAME_CONFLICT = "sequence_name_conflict"
EXTERNAL_NAME_MISSING = "external_name_missing"
INTERNAL_NAME_MISSING = "internal_name_missing"

RULE_TRANSFER_READY = "rule_transfer_ready"
RULE_TRANSFER_SUPPORTED_NOT_MODEL_OPERABLE = "rule_transfer_supported_not_model_operable"
RULE_TRANSFER_LOW_CONFIDENCE = "rule_transfer_low_confidence"
RULE_TRANSFER_PARALOG_RISK = "rule_transfer_paralog_risk"
RULE_TRANSFER_UNRESOLVED = "rule_transfer_unresolved"
RULE_TRANSFER_NOT_SUPPORTED = "rule_transfer_not_supported"

EXTERNAL_CROSSCHECK_NOT_AVAILABLE = "not_available"
EXTERNAL_MATCH_CONFIRMED = "external_match_confirmed"
EXTERNAL_ALIAS_CONFIRMED = "external_alias_confirmed"
EXTERNAL_LOCUS_CONFIRMED = "external_locus_confirmed"
EXTERNAL_CONFLICT = "external_conflict"
EXTERNAL_REFERENCE_INCOMPLETE = "external_reference_incomplete"


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
    if any(_looks_like_stable_pichia_identifier(value) for value in (external_gene_name, *external_aliases)):
        return PICHIA_LOCUS_CONFIRMED_BY_RBH
    return SEQUENCE_NAME_CONFLICT


def classify_rule_transfer_status(
    *,
    homology_review_status: str,
    is_rbh: bool,
    in_model_gene_index: bool,
) -> tuple[str, tuple[str, ...]]:
    """Classify whether an SCE rule can be transferred as review evidence."""

    if homology_review_status == UNRESOLVED_QUERY_SYMBOL:
        return RULE_TRANSFER_UNRESOLVED, ("query symbol is unresolved",)
    if homology_review_status == PARALOG_RISK_REVIEW_REQUIRED:
        return RULE_TRANSFER_PARALOG_RISK, ("paralog risk requires manual review",)
    if homology_review_status in {LOW_IDENTITY_REVIEW_REQUIRED, COVERAGE_REVIEW_REQUIRED, MANUAL_REVIEW_REQUIRED}:
        return RULE_TRANSFER_LOW_CONFIDENCE, (f"{homology_review_status} prevents high-confidence transfer",)
    if not is_rbh or homology_review_status == NO_RECIPROCAL_HIT:
        return RULE_TRANSFER_NOT_SUPPORTED, ("no reciprocal best hit support",)
    if is_rbh and not in_model_gene_index:
        return RULE_TRANSFER_SUPPORTED_NOT_MODEL_OPERABLE, (
            "RBH supports homology, but candidate is not present in current Pichia GEM gene_index",
        )
    if homology_review_status == MODEL_READY_RBH_HIGH_CONFIDENCE and in_model_gene_index:
        return RULE_TRANSFER_READY, ()
    return RULE_TRANSFER_NOT_SUPPORTED, (f"unsupported homology status: {homology_review_status}",)


def classify_external_name_crosscheck(
    *,
    internal_common_name: str,
    external_gene_name: str | None,
    external_locus_tag: str | None,
    external_aliases: tuple[str, ...] = (),
    reference_gene_name: str | None,
    reference_locus_tag: str | None,
    reference_aliases: tuple[str, ...] = (),
) -> tuple[str, tuple[str, ...]]:
    """Classify offline external database agreement without changing RBH facts."""

    candidate_names = _candidate_name_tokens(
        internal_common_name=internal_common_name,
        external_gene_name=external_gene_name,
        external_aliases=external_aliases,
    )
    reference_name = _normalize_name(reference_gene_name or "")
    reference_alias_names = {_normalize_name(alias) for alias in reference_aliases if alias}
    candidate_locus = _normalize_name(external_locus_tag or "")
    reference_locus = _normalize_name(reference_locus_tag or "")

    if reference_name and reference_name in candidate_names:
        return EXTERNAL_MATCH_CONFIRMED, ()
    if candidate_names & reference_alias_names:
        return EXTERNAL_ALIAS_CONFIRMED, ()
    if candidate_locus and reference_locus and candidate_locus == reference_locus:
        return EXTERNAL_LOCUS_CONFIRMED, ()
    if reference_name or reference_alias_names:
        return EXTERNAL_CONFLICT, (
            "external reference name or alias does not match current RBH-derived name fields",
        )
    return EXTERNAL_REFERENCE_INCOMPLETE, ("external reference lacks comparable name, alias, or locus data",)


def _normalize_name(value: str) -> str:
    return value.strip().upper().replace(" ", "").replace("_", "").replace("-", "")


def _split_name_tokens(value: str) -> set[str]:
    normalized = value.replace("/", " ").replace(",", " ").replace(";", " ")
    return {_normalize_name(token) for token in normalized.split() if token}


def _looks_like_stable_pichia_identifier(value: object) -> bool:
    text = str(value or "").strip().upper().replace(" ", "")
    if not text:
        return False
    if text.startswith(("PAS_", "PAS-", "PASCHR", "AT250_GQ_")):
        return True
    return len(text) >= 6 and text.startswith("C4") and text[:6].isalnum()


def _candidate_name_tokens(
    *,
    internal_common_name: str,
    external_gene_name: str | None,
    external_aliases: tuple[str, ...],
) -> set[str]:
    names = _split_name_tokens(internal_common_name)
    external = _normalize_name(external_gene_name or "")
    if external:
        names.add(external)
    names.update(_normalize_name(alias) for alias in external_aliases if alias)
    return {name for name in names if name}
