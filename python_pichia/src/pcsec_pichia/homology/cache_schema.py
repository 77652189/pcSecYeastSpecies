from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class ProteinRecord:
    organism: str
    gene_id: str
    sequence: str
    symbol: str | None = None
    aliases: tuple[str, ...] = ()
    accession: str | None = None
    source: str = ""


@dataclass(frozen=True)
class BlastConfig:
    max_evalue: float = 1e-10
    min_identity: float = 30.0
    min_coverage: float = 50.0
    max_target_seqs: int = 5
    threads: int = 1
    blast_bin: Path | None = None


@dataclass(frozen=True)
class BlastHit:
    query_id: str
    subject_id: str
    identity_pct: float
    alignment_length: int
    query_length: int
    subject_length: int
    evalue: float
    bitscore: float
    query_coverage: float
    subject_coverage: float

    @classmethod
    def from_lengths(
        cls,
        *,
        query_id: str,
        subject_id: str,
        identity_pct: float,
        alignment_length: int,
        query_length: int,
        subject_length: int,
        evalue: float,
        bitscore: float,
    ) -> "BlastHit":
        query_coverage = _coverage(alignment_length, query_length)
        subject_coverage = _coverage(alignment_length, subject_length)
        return cls(
            query_id=query_id,
            subject_id=subject_id,
            identity_pct=identity_pct,
            alignment_length=alignment_length,
            query_length=query_length,
            subject_length=subject_length,
            evalue=evalue,
            bitscore=bitscore,
            query_coverage=query_coverage,
            subject_coverage=subject_coverage,
        )


@dataclass(frozen=True)
class ReciprocalBestHit:
    query_id: str
    subject_id: str | None
    is_rbh: bool
    forward_hit: BlastHit | None
    reverse_hit: BlastHit | None
    failure_reason: str | None = None


@dataclass(frozen=True)
class CatalogHomologyQuery:
    internal_common_name: str
    query_symbol: str
    aliases: tuple[str, ...] = ()
    internal_gene_id: str = ""
    source: str = "secretion_gene_catalog"


@dataclass(frozen=True)
class HomologyCrosswalkRow:
    internal_common_name: str
    query_symbol: str
    sce_orf: str | None
    pichia_gene_id: str | None
    pichia_model_gene_id: str | None
    is_rbh: bool
    identity_pct: float | None
    evalue: float | None
    query_coverage: float | None
    subject_coverage: float | None
    in_model_gene_index: bool
    review_status: str
    warnings: tuple[str, ...] = ()
    external_accession: str = ""
    external_gene_name: str = ""
    external_locus_tag: str = ""
    external_aliases: tuple[str, ...] = ()


@dataclass(frozen=True)
class ExternalNameReference:
    source_database: str
    source_version: str
    taxon: str
    accession: str = ""
    gene_name: str = ""
    locus_tag: str = ""
    aliases: tuple[str, ...] = field(default_factory=tuple)
    retrieved_at: str = ""
    warnings: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class NameAuditRow:
    internal_gene_id: str
    internal_common_name: str
    internal_sequence_id: str
    external_accession: str
    external_gene_name: str
    external_locus_tag: str
    external_aliases: tuple[str, ...]
    identity_pct: float | None
    query_coverage: float | None
    subject_coverage: float | None
    evalue: float | None
    is_rbh: bool
    in_model_gene_index: bool
    name_consistency_status: str
    review_status: str
    external_crosscheck_status: str = "not_available"
    external_crosscheck_sources: tuple[str, ...] = field(default_factory=tuple)
    external_crosscheck_warnings: tuple[str, ...] = field(default_factory=tuple)
    warnings: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class ExternalDatabaseCrosscheck:
    internal_gene_id: str = ""
    internal_common_name: str = ""
    internal_sequence_id: str = ""
    external_accession: str = ""
    source_database: str = ""
    source_version: str = ""
    taxon: str = ""
    accession: str = ""
    gene_name: str = ""
    locus_tag: str = ""
    aliases: tuple[str, ...] = field(default_factory=tuple)
    retrieved_at: str = ""
    match_status: str = "not_available"
    warnings: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class RuleTransferAuditRow:
    internal_common_name: str
    query_symbol: str
    sce_orf: str
    pichia_gene_id: str
    pichia_model_gene_id: str
    is_rbh: bool
    in_model_gene_index: bool
    identity_pct: float | None
    query_coverage: float | None
    subject_coverage: float | None
    evalue: float | None
    homology_review_status: str
    rule_transfer_status: str
    warnings: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class HomologyAuditSummary:
    blast_status: str
    homology_row_count: int
    name_audit_row_count: int
    rule_transfer_row_count: int
    homology_review_status_counts: dict[str, int]
    name_consistency_status_counts: dict[str, int]
    rule_transfer_status_counts: dict[str, int]


@dataclass(frozen=True)
class CacheWriteResult:
    jsonl_path: Path
    tsv_path: Path
    row_count: int


def _coverage(alignment_length: int, sequence_length: int) -> float:
    if sequence_length <= 0:
        return 0.0
    return min(100.0, round(100.0 * alignment_length / sequence_length, 6))
