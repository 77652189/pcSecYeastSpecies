from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


INVENTORY_JSONL_FILENAME = "external_model_inventory.jsonl"
INVENTORY_TSV_FILENAME = "external_model_inventory.tsv"
INVENTORY_REPORT_FILENAME = "external_model_inventory_report.md"


@dataclass(frozen=True)
class ExternalModelInventoryRecord:
    model_id: str
    model_name: str
    organism: str
    source_database_or_repository: str
    source_url: str
    publication_url: str
    license: str
    available_artifact_types: tuple[str, ...]
    download_status: str
    local_path: str
    checksum_sha256: str
    has_gpr: bool
    has_gene_ids: bool
    has_reaction_ids: bool
    has_sbml: bool
    notes: str
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        for field_name in (
            "model_id",
            "model_name",
            "organism",
            "source_database_or_repository",
            "download_status",
        ):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"{field_name} must be non-empty.")
        if not self.source_url and not self.publication_url:
            raise ValueError("source_url or publication_url must be present.")
        if self.local_path and not self.checksum_sha256:
            raise ValueError("checksum_sha256 is required when local_path is present.")
        if self.checksum_sha256 and (
            len(self.checksum_sha256) != 64
            or any(char not in "0123456789abcdef" for char in self.checksum_sha256.lower())
        ):
            raise ValueError("checksum_sha256 must be a 64-character hex digest.")

    def to_dict(self) -> dict[str, object]:
        self.validate()
        return _json_ready(asdict(self))


@dataclass(frozen=True)
class ExternalModelInventoryOutputs:
    jsonl_path: Path
    tsv_path: Path
    report_path: Path
    record_count: int


def default_external_model_inventory_records() -> tuple[ExternalModelInventoryRecord, ...]:
    """Curated Round A audit seed for external GEM/GPR source prioritization.

    These records intentionally describe availability and provenance only. They
    do not download artifacts or import external GPR rules into the Pichia GEM.
    """

    return (
        ExternalModelInventoryRecord(
            model_id="iPichia",
            model_name="iPichia",
            organism="Komagataella phaffii / Pichia pastoris",
            source_database_or_repository="publication",
            source_url="",
            publication_url="https://doi.org/10.1016/j.bej.2025.109940",
            license="needs_manual_review",
            available_artifact_types=("publication",),
            download_status="needs_manual_access",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=False,
            notes=(
                "Prioritized as a Pichia-specific GEM/GPR source when an author "
                "or supplementary model artifact is obtained."
            ),
            warnings=("public_model_artifact_not_found_in_round_a",),
        ),
        ExternalModelInventoryRecord(
            model_id="ecPichia",
            model_name="ecPichia enzyme-constrained model",
            organism="Komagataella phaffii / Pichia pastoris",
            source_database_or_repository="publication",
            source_url="",
            publication_url="https://doi.org/10.1016/j.bej.2025.109940",
            license="needs_manual_review",
            available_artifact_types=("publication",),
            download_status="needs_manual_access",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=False,
            notes=(
                "Pichia-specific enzyme-constrained evidence source; artifact "
                "access needs manual confirmation before parsing."
            ),
            warnings=("public_model_artifact_not_found_in_round_a",),
        ),
        ExternalModelInventoryRecord(
            model_id="Kp.1.0",
            model_name="Kp.1.0 genome-scale model",
            organism="Komagataella phaffii",
            source_database_or_repository="Cambridge Apollo / BioModels",
            source_url="https://www.repository.cam.ac.uk/items/02da7483-3966-4d96-b90d-eda1e890e104",
            publication_url="https://doi.org/10.1002/bit.26380",
            license="CC BY 4.0",
            available_artifact_types=("SBML", "MATLAB", "model archive"),
            download_status="downloadable",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=True,
            notes="Pichia-specific GEM source; high priority for Round B artifact cache.",
            warnings=("not_downloaded_in_round_a",),
        ),
        ExternalModelInventoryRecord(
            model_id="iAUKM",
            model_name="iAUKM",
            organism="Komagataella phaffii / Pichia pastoris",
            source_database_or_repository="publication_or_author_repository",
            source_url="",
            publication_url="https://pubmed.ncbi.nlm.nih.gov/37597025/",
            license="needs_manual_review",
            available_artifact_types=("publication",),
            download_status="needs_manual_access",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=False,
            notes="Pichia-specific model candidate; public artifact location needs manual confirmation.",
            warnings=("public_model_artifact_not_found_in_round_a",),
        ),
        ExternalModelInventoryRecord(
            model_id="Yeast8_Yeast9",
            model_name="Yeast8 / Yeast9 yeast-GEM",
            organism="Saccharomyces cerevisiae",
            source_database_or_repository="GitHub SysBioChalmers/yeast-GEM",
            source_url="https://github.com/SysBioChalmers/yeast-GEM",
            publication_url="https://doi.org/10.1038/s41467-019-11581-3",
            license="MIT",
            available_artifact_types=("SBML", "YAML", "MATLAB", "repository"),
            download_status="repository_available",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=True,
            notes="Use after Pichia-specific sources, through RBH/homology rule-transfer review.",
            warnings=("not_downloaded_in_round_a", "cross_species_mapping_required"),
        ),
        ExternalModelInventoryRecord(
            model_id="BioModels_Kp.1.0_MODEL1703150000",
            model_name="BioModels Kp.1.0 entry",
            organism="Komagataella phaffii",
            source_database_or_repository="BioModels",
            source_url="https://www.ebi.ac.uk/biomodels/MODEL1703150000",
            publication_url="https://doi.org/10.1002/bit.26380",
            license="model_license_needs_manual_review",
            available_artifact_types=("SBML", "BioModels entry"),
            download_status="downloadable",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=True,
            notes="BioModels mirror/entry for Kp.1.0-style Pichia GEM audit.",
            warnings=("not_downloaded_in_round_a",),
        ),
        ExternalModelInventoryRecord(
            model_id="BioModels_iMT1026",
            model_name="iMT1026 Pichia pastoris GEM",
            organism="Pichia pastoris / Komagataella phaffii",
            source_database_or_repository="BioModels / publication",
            source_url="https://www.ebi.ac.uk/biomodels/",
            publication_url="https://doi.org/10.1186/1752-0509-6-24",
            license="model_license_needs_manual_review",
            available_artifact_types=("SBML_or_supplementary_model",),
            download_status="needs_manual_access",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=False,
            notes="Relevant historical Pichia GEM candidate; exact BioModels artifact should be confirmed before parsing.",
            warnings=("exact_biomodels_accession_needs_manual_confirmation",),
        ),
        ExternalModelInventoryRecord(
            model_id="GPRuler",
            model_name="GPRuler",
            organism="multi-organism tool",
            source_database_or_repository="publication / GitHub",
            source_url="https://github.com/qLSLab/GPRuler",
            publication_url="https://doi.org/10.1371/journal.pcbi.1009550",
            license="repository_license_needs_manual_review",
            available_artifact_types=("software", "publication", "repository"),
            download_status="tool_only_not_primary_gem",
            local_path="",
            checksum_sha256="",
            has_gpr=False,
            has_gene_ids=False,
            has_reaction_ids=False,
            has_sbml=False,
            notes=(
                "Supplemental automatic rule-generation tool only; not a primary "
                "GEM/GPR evidence source for current prioritization."
            ),
            warnings=("supplemental_gpr_tool_only", "do_not_treat_as_model_gpr_executable"),
        ),
    )


def write_external_model_inventory(
    records: Iterable[ExternalModelInventoryRecord],
    output_dir: Path,
) -> ExternalModelInventoryOutputs:
    resolved = tuple(records)
    for record in resolved:
        record.validate()
    output_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = output_dir / INVENTORY_JSONL_FILENAME
    tsv_path = output_dir / INVENTORY_TSV_FILENAME
    report_path = output_dir / INVENTORY_REPORT_FILENAME
    _write_jsonl(resolved, jsonl_path)
    _write_tsv(resolved, tsv_path)
    report_path.write_text(render_external_model_inventory_report(resolved), encoding="utf-8")
    return ExternalModelInventoryOutputs(
        jsonl_path=jsonl_path,
        tsv_path=tsv_path,
        report_path=report_path,
        record_count=len(resolved),
    )


def load_external_model_inventory(path: Path) -> tuple[ExternalModelInventoryRecord, ...]:
    """Load inventory records from a JSONL file or an inventory output directory."""

    resolved_path = path / INVENTORY_JSONL_FILENAME if path.is_dir() else path
    records: list[ExternalModelInventoryRecord] = []
    with resolved_path.open("r", encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                payload = json.loads(line)
                record = _record_from_json_payload(payload)
                record.validate()
            except Exception as exc:
                raise ValueError(
                    f"Invalid external model inventory record at {resolved_path}:{row_number}: {exc}"
                ) from exc
            records.append(record)
    return tuple(records)


def render_external_model_inventory_report(records: Iterable[ExternalModelInventoryRecord]) -> str:
    resolved = tuple(records)
    counts: dict[str, int] = {}
    for record in resolved:
        counts[record.download_status] = counts.get(record.download_status, 0) + 1
    lines = [
        "# External GEM / GPR Resource Inventory",
        "",
        "This audit records external model resources as evidence sources only; it does not import external GPR rules into the current Pichia GEM.",
        "",
        "## Summary",
        "",
        f"- record_count: {len(resolved)}",
    ]
    for status, count in sorted(counts.items()):
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Records",
            "",
            "| model_id | organism | source | status | has_gpr | has_sbml | warnings |",
            "|---|---|---|---|---:|---:|---|",
        ]
    )
    for record in resolved:
        warnings = "; ".join(record.warnings)
        lines.append(
            f"| {record.model_id} | {record.organism} | {record.source_database_or_repository} | "
            f"{record.download_status} | {record.has_gpr} | {record.has_sbml} | {warnings} |"
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- External GEM evidence remains `external_gpr_candidate` until reaction and gene mapping to the current Pichia GEM is proven.",
            "- UniProt/NCBI/SGD annotation and automatic tools such as GPRuler are supplemental evidence, not phenotype evidence and not `experiment_calibrated`.",
            "- Missing public artifacts are recorded as `needs_manual_access`; no URL, local file, or checksum is fabricated.",
        ]
    )
    return "\n".join(lines) + "\n"


def _write_jsonl(records: tuple[ExternalModelInventoryRecord, ...], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record.to_dict(), ensure_ascii=False, sort_keys=True))
            handle.write("\n")


def _write_tsv(records: tuple[ExternalModelInventoryRecord, ...], path: Path) -> None:
    fieldnames = tuple(ExternalModelInventoryRecord.__dataclass_fields__)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({key: _tsv_value(value) for key, value in record.to_dict().items()})


def _tsv_value(value: object) -> str:
    if isinstance(value, list):
        return "|".join(str(item) for item in value)
    return str(value)


def _record_from_json_payload(payload: object) -> ExternalModelInventoryRecord:
    if not isinstance(payload, dict):
        raise ValueError("record must be a JSON object.")
    field_names = set(ExternalModelInventoryRecord.__dataclass_fields__)
    unknown_fields = sorted(str(key) for key in payload if key not in field_names)
    if unknown_fields:
        raise ValueError(f"unexpected fields: {', '.join(unknown_fields)}")
    values = dict(payload)
    values["available_artifact_types"] = _tuple_from_json_list(values.get("available_artifact_types", ()))
    values["warnings"] = _tuple_from_json_list(values.get("warnings", ()))
    try:
        return ExternalModelInventoryRecord(**values)
    except TypeError as exc:
        raise ValueError(str(exc)) from exc


def _tuple_from_json_list(value: object) -> tuple[str, ...]:
    if isinstance(value, tuple):
        return tuple(str(item) for item in value)
    if isinstance(value, list):
        return tuple(str(item) for item in value)
    raise ValueError("tuple fields must be encoded as JSON lists.")


def _json_ready(value: object) -> object:
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "INVENTORY_JSONL_FILENAME",
    "INVENTORY_REPORT_FILENAME",
    "INVENTORY_TSV_FILENAME",
    "ExternalModelInventoryOutputs",
    "ExternalModelInventoryRecord",
    "default_external_model_inventory_records",
    "load_external_model_inventory",
    "render_external_model_inventory_report",
    "write_external_model_inventory",
]
