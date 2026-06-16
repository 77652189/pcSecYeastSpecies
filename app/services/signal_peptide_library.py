from __future__ import annotations

import csv
import io
import json
import re
from dataclasses import dataclass
from typing import Iterable
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

from app.services.opn import OPN_SHORTLIST, OpnCandidateCatalog, opn_category_label


AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
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


@dataclass(frozen=True)
class CandidateDiscoveryResult:
    rows: list[dict[str, object]]
    source_url: str
    errors: list[str]


@dataclass(frozen=True)
class UniProtCandidateLibraryResult:
    rows: list[dict[str, object]]
    source_url: str
    errors: list[str]
    initial_hit_count: int
    fetched_record_count: int
    extracted_signal_count: int
    deduplicated_count: int


@dataclass
class SignalPeptideLibraryService:
    catalog: OpnCandidateCatalog

    def library_rows(self) -> list[dict[str, object]]:
        rows = []
        for candidate in self.catalog.list_candidates():
            stage = self._library_stage(candidate.candidate_id, candidate.category)
            rows.append(
                {
                    "candidate_id": candidate.candidate_id,
                    "category": candidate.category,
                    "category_label": opn_category_label(candidate.category),
                    "library_stage": stage,
                    "leader_length": candidate.leader_length,
                    "signal_peptide_length": len(candidate.signal_peptide_sequence),
                    "construct_length": candidate.construct_length,
                    "processing_route": candidate.processing_route,
                    "source_type": self._source_type(candidate.source_note),
                    "source_note": candidate.source_note,
                    "rationale": candidate.rationale,
                    "caution": candidate.caution,
                    "leader_sequence": candidate.leader_sequence,
                    "signal_peptide_sequence": candidate.signal_peptide_sequence,
                }
            )
        return rows

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
        query = f"(organism_id:{taxon_id}) AND (ft_signal:*)"
        if reviewed_only:
            query += " AND (reviewed:true)"
        safe_page_size = max(1, min(int(page_size), 500))
        params = {
            "query": query,
            "format": "json",
            "size": str(safe_page_size),
            "fields": "accession,id,protein_name,organism_name,ft_signal,sequence",
        }
        url = "https://rest.uniprot.org/uniprotkb/search?" + urlencode(params)
        next_url: str | None = url
        items: list[dict] = []
        initial_hit_count = 0
        errors: list[str] = []
        try:
            while next_url and len(items) < max_records:
                request = Request(next_url, headers={"User-Agent": "pcSecYeastSpecies-local/1.0"})
                with urlopen(request, timeout=30) as response:
                    if not initial_hit_count:
                        initial_hit_count = _safe_int(response.headers.get("x-total-results"))
                    payload = json.loads(response.read().decode("utf-8"))
                    items.extend(payload.get("results", []))
                    next_url = _next_link(response.headers.get("Link"))
        except Exception as exc:
            return UniProtCandidateLibraryResult(
                rows=[],
                source_url=url,
                errors=[f"UniProt API 请求失败：{exc}"],
                initial_hit_count=initial_hit_count,
                fetched_record_count=len(items),
                extracted_signal_count=0,
                deduplicated_count=0,
            )
        limited_items = items[:max_records]
        rows, row_errors, extracted_signal_count = self._rows_from_uniprot_items(
            limited_items,
            exclude_existing=exclude_existing,
        )
        errors.extend(row_errors)
        if not initial_hit_count:
            initial_hit_count = len(items)
        return UniProtCandidateLibraryResult(
            rows=rows,
            source_url=url,
            errors=errors,
            initial_hit_count=initial_hit_count,
            fetched_record_count=len(limited_items),
            extracted_signal_count=extracted_signal_count,
            deduplicated_count=len(rows),
        )

    def rows_from_uniprot_payload(self, payload: dict) -> tuple[list[dict[str, object]], list[str]]:
        rows, errors, _extracted_count = self._rows_from_uniprot_items(
            payload.get("results", []),
            exclude_existing=True,
        )
        return rows, errors

    def _rows_from_uniprot_items(
        self,
        items: Iterable[dict],
        *,
        exclude_existing: bool,
    ) -> tuple[list[dict[str, object]], list[str], int]:
        existing_ids = {row["candidate_id"] for row in self.library_rows()}
        existing_leaders = {row["leader_sequence"] for row in self.library_rows()}
        rows: list[dict[str, object]] = []
        errors: list[str] = []
        seen_sequences: set[str] = set()
        extracted_signal_count = 0
        for item in items:
            accession = item.get("primaryAccession", "")
            sequence = item.get("sequence", {}).get("value", "")
            signal, signal_start, signal_end = _signal_feature(item, sequence)
            if not accession or not signal:
                continue
            extracted_signal_count += 1
            candidate_id = f"OPN_UNIPROT_{_safe_id(accession)}"
            already_in_formal_library = candidate_id in existing_ids or signal in existing_leaders
            if exclude_existing and candidate_id in existing_ids:
                continue
            if exclude_existing and signal in existing_leaders:
                continue
            if signal in seen_sequences:
                continue
            seen_sequences.add(signal)
            protein_name = _protein_name(item)
            organism_name = item.get("organism", {}).get("scientificName", "")
            uniprot_id = item.get("uniProtkbId", "")
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "accession": accession,
                    "uniprot_id": uniprot_id,
                    "protein_name": protein_name,
                    "organism_name": organism_name,
                    "protein_sequence": sequence,
                    "protein_length": len(sequence),
                    "uniprot_signal_start": signal_start,
                    "uniprot_signal_end": signal_end,
                    "leader_sequence": signal,
                    "signal_peptide_sequence": signal,
                    "category": "pichia_native_signal",
                    "category_label": opn_category_label("pichia_native_signal"),
                    "processing_route": "signal peptidase only",
                    "source_note": f"UniProt {accession} {uniprot_id}; {organism_name}; {protein_name}",
                    "rationale": "UniProt 中带 signal peptide 注释的 Komagataella/Pichia 蛋白，适合作为外部候选草案。",
                    "caution": "来自数据库自动发现；进入模型前需要人工确认切割位点、文献证据和是否适合 OPN。",
                    "leader_length": len(signal),
                    "signal_peptide_length": len(signal),
                    "library_stage": "外部发现草案",
                    "source_type": "UniProt",
                    "already_in_formal_library": already_in_formal_library,
                }
            )
        if not rows:
            errors.append("没有发现可加入草案的新 signal peptide；可能都已存在或查询结果没有明确序列。")
        return rows, errors, extracted_signal_count

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
                        "category_label": opn_category_label(str(row["category"])),
                    }
                )
        return CandidateImportValidation(valid=not errors, rows=rows, errors=errors)

    def merged_draft_csv(self, new_rows: Iterable[dict[str, object]]) -> bytes:
        current = pd.DataFrame(self.library_rows())
        additions = pd.DataFrame(list(new_rows))
        merged = pd.concat([current, additions], ignore_index=True, sort=False)
        return merged.to_csv(index=False).encode("utf-8-sig")

    def _library_stage(self, candidate_id: str, category: str) -> str:
        if candidate_id in OPN_SHORTLIST:
            return "首轮推荐"
        if category == "project_baseline":
            return "对照基线"
        return "候选库"

    def _source_type(self, source_note: str) -> str:
        lowered = source_note.lower()
        if "uniprot" in lowered:
            return "UniProt"
        if "reported" in lowered or "paper" in lowered or "pmcid" in lowered:
            return "文献"
        if "project" in lowered:
            return "项目基线"
        return "待补充"

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


def _signal_sequence(item: dict, sequence: str) -> str:
    signal, _start, _end = _signal_feature(item, sequence)
    return signal


def _signal_feature(item: dict, sequence: str) -> tuple[str, int | None, int | None]:
    for feature in item.get("features", []):
        if feature.get("type") != "Signal":
            continue
        location = feature.get("location", {})
        start = location.get("start", {}).get("value")
        end = location.get("end", {}).get("value")
        if not start or not end:
            continue
        start_int = int(start)
        end_int = int(end)
        signal = sequence[start_int - 1 : end_int]
        if AA_PATTERN.fullmatch(signal):
            return signal, start_int, end_int
    return "", None, None


def _protein_name(item: dict) -> str:
    recommended = item.get("proteinDescription", {}).get("recommendedName", {})
    full_name = recommended.get("fullName", {})
    if full_name.get("value"):
        return str(full_name["value"])
    submission_names = item.get("proteinDescription", {}).get("submissionNames", [])
    if submission_names:
        return str(submission_names[0].get("fullName", {}).get("value", ""))
    return ""


def _safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").upper()


def _safe_int(value: object) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return 0


def _next_link(link_header: str | None) -> str | None:
    if not link_header:
        return None
    match = re.search(r"<([^>]+)>\s*;\s*rel=\"next\"", link_header)
    return match.group(1) if match else None
