from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from typing import Iterable

import pandas as pd

from app.adapters.uniprot import UniProtSignalPeptideSource
from app.core.signal_peptides import (
    AA_PATTERN,
    CandidateDiscoveryResult,
    SignalPeptideCandidate,
    UniProtCandidateLibraryResult,
)


REQUIRED_IMPORT_COLUMNS = [
    "candidate_id",
    "leader_sequence",
    "signal_peptide_sequence",
    "category",
    "processing_route",
    "source_note",
    "rationale",
    "caution",
]


@dataclass(frozen=True)
class CandidateImportValidation:
    valid: bool
    rows: list[dict[str, object]]
    errors: list[str]


@dataclass
class SignalPeptideLibraryService:
    candidates: Iterable[SignalPeptideCandidate] = field(default_factory=list)

    def __post_init__(self) -> None:
        self._candidates = list(self.candidates)

    def list_candidates(self) -> list[SignalPeptideCandidate]:
        return list(self._candidates)

    def library_rows(self) -> list[dict[str, object]]:
        return [candidate.as_row() for candidate in self._candidates]

    def template_csv(self) -> bytes:
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=REQUIRED_IMPORT_COLUMNS)
        writer.writeheader()
        writer.writerow(
            {
                "candidate_id": "OPN_NEW_SIGNAL_001",
                "leader_sequence": "MKFAISTLLIILQAAAVFAA",
                "signal_peptide_sequence": "MKFAISTLLIILQAAAVFAA",
                "category": "pichia_native_signal",
                "processing_route": "signal peptidase only",
                "source_note": "请填写 UniProt/论文/外部工具/内部实验来源",
                "rationale": "请填写为什么值得进入候选库",
                "caution": "请填写待验证风险，例如切割位点、宿主适配性或文献证据不足",
            }
        )
        return output.getvalue().encode("utf-8-sig")

    def discover_uniprot_candidates(
        self,
        *,
        taxon_id: int = 4922,
        size: int = 25,
        reviewed_only: bool = False,
    ) -> CandidateDiscoveryResult:
        library = self.discover_uniprot_candidate_library(
            taxon_id=taxon_id,
            max_records=size,
            reviewed_only=reviewed_only,
            exclude_existing=True,
        )
        return CandidateDiscoveryResult(library.rows, library.source_url, library.errors)

    def discover_uniprot_candidate_library(
        self,
        *,
        taxon_id: int = 4922,
        max_records: int = 300,
        reviewed_only: bool = False,
        exclude_existing: bool = False,
        page_size: int = 100,
    ) -> UniProtCandidateLibraryResult:
        return UniProtSignalPeptideSource(self._candidates).discover(
            taxon_id=taxon_id,
            max_records=max_records,
            reviewed_only=reviewed_only,
            exclude_existing=exclude_existing,
            page_size=page_size,
        )

    def rows_from_uniprot_payload(self, payload: dict) -> tuple[list[dict[str, object]], list[str]]:
        return UniProtSignalPeptideSource(self._candidates).rows_from_payload(payload)

    def validate_import_csv(self, content: bytes) -> CandidateImportValidation:
        try:
            frame = pd.read_csv(io.BytesIO(content))
        except Exception as exc:
            return CandidateImportValidation(False, [], [f"CSV 读取失败：{exc}"])

        errors = self._missing_column_errors(frame.columns)
        rows: list[dict[str, object]] = []
        existing_ids = {row["candidate_id"] for row in self.library_rows()}
        seen_ids: set[str] = set()
        if not errors:
            for index, raw in frame.iterrows():
                row_number = index + 2
                row = {column: _clean(raw.get(column)) for column in REQUIRED_IMPORT_COLUMNS}
                row_errors = self._validate_row(row, row_number, existing_ids, seen_ids)
                errors.extend(row_errors)
                if row_errors:
                    continue
                seen_ids.add(str(row["candidate_id"]))
                rows.append(
                    {
                        **row,
                        "leader_length": len(str(row["leader_sequence"])),
                        "signal_peptide_length": len(str(row["signal_peptide_sequence"])),
                        "library_stage": "待审核",
                        "category_label": str(row["category"]),
                    }
                )
        return CandidateImportValidation(valid=not errors, rows=rows, errors=errors)

    def merged_draft_csv(self, new_rows: Iterable[dict[str, object]]) -> bytes:
        current = pd.DataFrame(self.library_rows())
        additions = pd.DataFrame(list(new_rows))
        merged = pd.concat([current, additions], ignore_index=True, sort=False)
        return merged.to_csv(index=False).encode("utf-8-sig")

    def _missing_column_errors(self, columns) -> list[str]:
        missing = [column for column in REQUIRED_IMPORT_COLUMNS if column not in columns]
        if not missing:
            return []
        return [f"缺少必填列：{', '.join(missing)}"]

    def _validate_row(
        self,
        row: dict[str, object],
        row_number: int,
        existing_ids: set[object],
        seen_ids: set[str],
    ) -> list[str]:
        errors = []
        candidate_id = str(row["candidate_id"])
        leader = str(row["leader_sequence"])
        signal = str(row["signal_peptide_sequence"])
        if not candidate_id:
            errors.append(f"第 {row_number} 行：candidate_id 不能为空")
        if candidate_id in existing_ids:
            errors.append(f"第 {row_number} 行：candidate_id 已存在：{candidate_id}")
        if candidate_id in seen_ids:
            errors.append(f"第 {row_number} 行：candidate_id 在导入文件中重复：{candidate_id}")
        if not AA_PATTERN.fullmatch(leader):
            errors.append(f"第 {row_number} 行：leader_sequence 只能包含标准氨基酸单字母代码")
        if not AA_PATTERN.fullmatch(signal):
            errors.append(f"第 {row_number} 行：signal_peptide_sequence 只能包含标准氨基酸单字母代码")
        if signal and leader and signal not in leader:
            errors.append(f"第 {row_number} 行：signal_peptide_sequence 应包含在 leader_sequence 中")
        return errors


def _clean(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()
