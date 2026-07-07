from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from pcsec_pichia.homology.crosswalk import load_name_audit_cache, load_rule_transfer_audit_cache


DEFAULT_HOMOLOGY_AUDIT_CACHE_DIR = Path("local_runs") / "pichia_homology_cache"
NAME_AUDIT_JSONL = "sce_to_pichia_name_audit.jsonl"
RULE_TRANSFER_JSONL = "sce_to_pichia_rule_transfer_audit.jsonl"


@dataclass(frozen=True)
class GeneHomologyEvidence:
    gene_id: str
    internal_common_name: str = ""
    query_symbol: str = ""
    sce_orf: str = ""
    pichia_gene_id: str = ""
    pichia_model_gene_id: str = ""
    is_rbh: bool = False
    in_model_gene_index: bool = False
    identity_pct: float | None = None
    query_coverage: float | None = None
    subject_coverage: float | None = None
    evalue: float | None = None
    homology_review_status: str = ""
    rule_transfer_status: str = ""
    name_consistency_status: str = ""
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["warnings"] = list(self.warnings)
        return payload


def load_homology_evidence_cache(
    cache_dir: Path | str | None = None,
) -> dict[str, GeneHomologyEvidence]:
    root = Path(cache_dir) if cache_dir is not None else DEFAULT_HOMOLOGY_AUDIT_CACHE_DIR
    run_dir = _latest_valid_run(root)
    if run_dir is None:
        return {}
    try:
        name_rows = load_name_audit_cache(run_dir / NAME_AUDIT_JSONL)
        rule_rows = load_rule_transfer_audit_cache(run_dir / RULE_TRANSFER_JSONL)
    except (FileNotFoundError, OSError, ValueError, TypeError):
        return {}
    return build_homology_evidence_map(name_rows=name_rows, rule_rows=rule_rows)


def build_homology_evidence_map(
    *,
    name_rows: tuple[Any, ...] = (),
    rule_rows: tuple[Any, ...] = (),
) -> dict[str, GeneHomologyEvidence]:
    name_status_by_key = {
        _name_key(getattr(row, "internal_common_name", ""), getattr(row, "internal_sequence_id", "")): getattr(
            row, "name_consistency_status", ""
        )
        for row in name_rows
    }
    evidence_by_key: dict[str, GeneHomologyEvidence] = {}
    for row in rule_rows:
        pichia_model_gene_id = str(getattr(row, "pichia_model_gene_id", "") or "")
        pichia_gene_id = str(getattr(row, "pichia_gene_id", "") or "")
        gene_id = pichia_model_gene_id or pichia_gene_id
        if not gene_id:
            continue
        name_status = name_status_by_key.get(
            _name_key(getattr(row, "internal_common_name", ""), getattr(row, "sce_orf", "")),
            "",
        )
        evidence = GeneHomologyEvidence(
            gene_id=gene_id,
            internal_common_name=str(getattr(row, "internal_common_name", "") or ""),
            query_symbol=str(getattr(row, "query_symbol", "") or ""),
            sce_orf=str(getattr(row, "sce_orf", "") or ""),
            pichia_gene_id=pichia_gene_id,
            pichia_model_gene_id=pichia_model_gene_id,
            is_rbh=bool(getattr(row, "is_rbh", False)),
            in_model_gene_index=bool(getattr(row, "in_model_gene_index", False)),
            identity_pct=getattr(row, "identity_pct", None),
            query_coverage=getattr(row, "query_coverage", None),
            subject_coverage=getattr(row, "subject_coverage", None),
            evalue=getattr(row, "evalue", None),
            homology_review_status=str(getattr(row, "homology_review_status", "") or ""),
            rule_transfer_status=str(getattr(row, "rule_transfer_status", "") or ""),
            name_consistency_status=name_status,
            warnings=tuple(str(item) for item in getattr(row, "warnings", ()) or ()),
        )
        for key in _evidence_index_keys(evidence):
            evidence_by_key.setdefault(key, evidence)
    return evidence_by_key


def homology_evidence_for_gene(
    gene_id: str,
    evidence_by_gene: dict[str, GeneHomologyEvidence] | None,
    aliases: tuple[str, ...] = (),
) -> GeneHomologyEvidence | None:
    if not evidence_by_gene:
        return None
    for key in (gene_id, *aliases):
        text = str(key or "").strip()
        normalized = _normalize_key(key)
        if text and text in evidence_by_gene:
            return evidence_by_gene[text]
        if normalized and normalized in evidence_by_gene:
            return evidence_by_gene[normalized]
    return None


def _latest_valid_run(root: Path) -> Path | None:
    if (root / NAME_AUDIT_JSONL).exists() and (root / RULE_TRANSFER_JSONL).exists():
        return root
    if not root.exists():
        return None
    candidates = [
        path
        for path in root.iterdir()
        if path.is_dir() and (path / NAME_AUDIT_JSONL).exists() and (path / RULE_TRANSFER_JSONL).exists()
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _name_key(common_name: object, sce_orf: object) -> tuple[str, str]:
    return (_normalize_key(common_name), _normalize_key(sce_orf))


def _evidence_index_keys(evidence: GeneHomologyEvidence) -> tuple[str, ...]:
    keys = (
        evidence.gene_id,
        evidence.pichia_model_gene_id,
        evidence.pichia_gene_id,
        evidence.internal_common_name,
        evidence.query_symbol,
        evidence.sce_orf,
    )
    return tuple(dict.fromkeys(key for key in (_normalize_key(item) for item in keys) if key))


def _normalize_key(value: object) -> str:
    return str(value or "").strip().lower()


__all__ = [
    "DEFAULT_HOMOLOGY_AUDIT_CACHE_DIR",
    "GeneHomologyEvidence",
    "build_homology_evidence_map",
    "homology_evidence_for_gene",
    "load_homology_evidence_cache",
]
