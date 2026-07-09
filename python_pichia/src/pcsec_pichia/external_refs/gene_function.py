from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping

from pcsec_pichia.external_refs.name_resolution import (
    EXTERNAL_REFERENCE_MISSING,
    classify_external_name_consistency,
    select_external_records_for_name_audit_row,
)
from pcsec_pichia.external_refs.schema import (
    ExternalGeneFunctionEvidence,
    ExternalGprCandidateEvidence,
    ExternalReferenceRecord,
)


@dataclass(frozen=True)
class KoOeExternalGeneEvidence:
    pichia_gene_id: str
    standard_name: str | None
    external_name_status: str
    function_evidence: tuple[ExternalGeneFunctionEvidence, ...]
    gpr_candidates: tuple[ExternalGprCandidateEvidence, ...] = ()
    model_executable_gene_id: str | None = None
    model_gpr_executable: bool = False
    manual_review_reasons: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["function_evidence"] = [asdict(item) for item in self.function_evidence]
        payload["gpr_candidates"] = [asdict(item) for item in self.gpr_candidates]
        return payload


def build_gene_function_evidence(
    *,
    internal_gene_id: str,
    external_records: Iterable[ExternalReferenceRecord],
) -> tuple[ExternalGeneFunctionEvidence, ...]:
    """Extract annotation-only gene function evidence from external references."""

    matched = _records_for_gene(internal_gene_id, external_records)
    evidence: list[ExternalGeneFunctionEvidence] = []
    for record in matched:
        if not _has_function_annotation(record):
            continue
        item = ExternalGeneFunctionEvidence(
            provenance=record.provenance,
            gene_id=internal_gene_id,
            protein_name=record.protein_name,
            function_description=record.function_description,
            ec_numbers=record.ec_numbers,
            go_terms=record.go_terms,
            pathways=record.pathways,
            orthology=record.orthology,
            reviewed=record.reviewed,
            evidence_scope="annotation_only",
        )
        confidence, _warnings = classify_gene_function_confidence(item)
        item = ExternalGeneFunctionEvidence(
            provenance=item.provenance,
            gene_id=item.gene_id,
            protein_name=item.protein_name,
            function_description=item.function_description,
            ec_numbers=item.ec_numbers,
            go_terms=item.go_terms,
            pathways=item.pathways,
            orthology=item.orthology,
            reviewed=item.reviewed,
            evidence_scope=confidence,
        )
        evidence.append(item)
    return tuple(sorted(evidence, key=lambda item: (item.provenance.source_database, item.provenance.source_query)))


def classify_gene_function_confidence(
    evidence: ExternalGeneFunctionEvidence,
) -> tuple[str, tuple[str, ...]]:
    """Classify annotation evidence while keeping it separate from phenotype evidence."""

    warnings: list[str] = []
    has_function = bool(evidence.function_description or evidence.protein_name)
    has_structured_terms = bool(evidence.ec_numbers or evidence.go_terms or evidence.pathways or evidence.orthology)
    if evidence.reviewed and has_structured_terms:
        return "reviewed_structured_annotation", (
            "external annotation is not phenotype evidence and must not calibrate recommendation_tier",
        )
    if evidence.reviewed and has_function:
        return "reviewed_name_annotation", (
            "external annotation is not phenotype evidence and must not calibrate recommendation_tier",
        )
    if has_structured_terms:
        return "structured_annotation", (
            "external annotation is not phenotype evidence and must not calibrate recommendation_tier",
        )
    if has_function:
        return "name_annotation", (
            "external annotation is not phenotype evidence and must not calibrate recommendation_tier",
        )
    warnings.append("external record lacks protein, function, EC, GO, pathway, or orthology annotation")
    return "manual_review_required", tuple(warnings)


def build_ko_oe_external_gene_evidence(
    *,
    pichia_gene_id: str,
    external_records: Iterable[ExternalReferenceRecord],
    standard_name: str | None = None,
    aliases: Iterable[str] = (),
    model_executable_gene_id: str | None = None,
    model_gpr_executable: bool = False,
) -> KoOeExternalGeneEvidence:
    matched = _records_for_gene(pichia_gene_id, external_records)
    name_candidate = classify_external_name_consistency(
        internal_gene_id=pichia_gene_id,
        internal_common_name=standard_name,
        internal_aliases=aliases,
        external_records=matched,
    )
    function_evidence = build_gene_function_evidence(
        internal_gene_id=pichia_gene_id,
        external_records=matched,
    )
    manual_review_reasons = tuple(
        dict.fromkeys(
            (
                *name_candidate.manual_review_reasons,
                *(
                    ("external reference has no function annotation",)
                    if matched and not function_evidence
                    else ()
                ),
            )
        )
    )
    return KoOeExternalGeneEvidence(
        pichia_gene_id=pichia_gene_id,
        standard_name=standard_name,
        external_name_status=name_candidate.external_name_status,
        function_evidence=function_evidence,
        model_executable_gene_id=model_executable_gene_id,
        model_gpr_executable=model_gpr_executable,
        manual_review_reasons=manual_review_reasons,
    )


def attach_ko_oe_external_gene_evidence(
    candidate_rows: Iterable[Mapping[str, Any]],
    external_records: Iterable[ExternalReferenceRecord],
) -> tuple[Mapping[str, Any], ...]:
    """Attach external annotation fields to KO/OE rows without changing tiers."""

    records = tuple(external_records)
    rows: list[Mapping[str, Any]] = []
    for row in candidate_rows:
        payload = dict(row)
        gene_id = _candidate_gene_id(payload)
        matched = _matched_records_for_candidate(payload, records)
        evidence = build_ko_oe_external_gene_evidence(
            pichia_gene_id=gene_id,
            standard_name=_optional_text(
                payload.get("standard_gene_symbol")
                or payload.get("common_name")
                or payload.get("internal_common_name")
                or payload.get("display_name")
                or payload.get("input_gene_id")
            ),
            aliases=_candidate_aliases(payload),
            external_records=matched,
            model_executable_gene_id=_optional_text(payload.get("canonical_gene_id") or payload.get("gene_id")),
            model_gpr_executable=bool(payload.get("model_gpr_executable") or payload.get("ko_support_status") == "ko_runnable_gpr_gene_deletion"),
        )
        rows.append(
            {
                **payload,
                "external_gene_function_evidence": tuple(asdict(item) for item in evidence.function_evidence),
                "ko_oe_external_gene_evidence": evidence.to_dict(),
                "external_gene_function_confidence": tuple(item.evidence_scope for item in evidence.function_evidence),
                "external_gene_function_sources": tuple(
                    dict.fromkeys(item.provenance.source_database for item in evidence.function_evidence)
                ),
                "external_gene_function_warnings": tuple(
                    dict.fromkeys(
                        warning
                        for item in evidence.function_evidence
                        for warning in classify_gene_function_confidence(item)[1]
                    )
                ),
            }
        )
    return tuple(rows)


def _matched_records_for_candidate(
    row: Mapping[str, Any],
    records: tuple[ExternalReferenceRecord, ...],
) -> tuple[ExternalReferenceRecord, ...]:
    selected = select_external_records_for_name_audit_row(row, records)
    if selected:
        return selected
    gene_id = _candidate_gene_id(row)
    return _records_for_gene(gene_id, records)


def _records_for_gene(
    internal_gene_id: str,
    records: Iterable[ExternalReferenceRecord],
) -> tuple[ExternalReferenceRecord, ...]:
    gene_tokens = _tokens(internal_gene_id)
    matched = [
        record
        for record in records
        if gene_tokens & _record_tokens(record)
    ]
    return tuple(sorted(matched, key=lambda record: record.cache_key))


def _has_function_annotation(record: ExternalReferenceRecord) -> bool:
    return bool(
        record.protein_name
        or record.function_description
        or record.ec_numbers
        or record.go_terms
        or record.pathways
        or record.orthology
    )


def _candidate_gene_id(row: Mapping[str, Any]) -> str:
    for key in ("canonical_gene_id", "gene_id", "internal_gene_id", "mapped_model_gene_id", "input_gene_id", "pichia_gene_id"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    return ""


def _candidate_aliases(row: Mapping[str, Any]) -> tuple[str, ...]:
    aliases: list[str] = []
    for key in ("standard_gene_symbol", "common_name", "internal_common_name", "display_name", "input_gene_id", "external_aliases"):
        aliases.extend(_tuple_values(row.get(key)))
    return tuple(dict.fromkeys(aliases))


def _record_tokens(record: ExternalReferenceRecord) -> set[str]:
    return _tokens(
        record.primary_accession,
        record.gene_id,
        record.gene_name,
        record.locus_tag,
        record.provenance.source_query,
        *record.aliases,
    )


def _tokens(*values: object) -> set[str]:
    result: set[str] = set()
    for value in values:
        text = str(value or "")
        for sep in ("/", ",", ";", "|"):
            text = text.replace(sep, " ")
        for part in (str(value or ""), *text.split()):
            normalized = part.strip().upper().replace(" ", "").replace("_", "").replace("-", "")
            if normalized:
                result.add(normalized)
    return result


def _tuple_values(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    text = str(value).strip()
    if not text:
        return ()
    return tuple(part.strip() for part in text.split(";") if part.strip())


def _optional_text(value: object) -> str | None:
    text = str(value or "").strip()
    return text or None


__all__ = [
    "KoOeExternalGeneEvidence",
    "attach_ko_oe_external_gene_evidence",
    "build_gene_function_evidence",
    "build_ko_oe_external_gene_evidence",
    "classify_gene_function_confidence",
]
