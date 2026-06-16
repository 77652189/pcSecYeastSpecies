from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Iterable
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from app.core.signal_peptides import AA_PATTERN, SignalPeptideCandidate, UniProtCandidateLibraryResult


PICHIA_NATIVE_SIGNAL_LABEL = "毕赤酵母来源信号肽"


@dataclass
class UniProtSignalPeptideSource:
    existing_candidates: Iterable[SignalPeptideCandidate] = field(default_factory=list)

    def discover(
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
                duplicate_count=0,
                duplicate_rows=[],
            )

        limited_items = items[:max_records]
        rows, row_errors, extracted_signal_count, duplicate_rows = self.rows_from_items(
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
            duplicate_count=len(duplicate_rows),
            duplicate_rows=duplicate_rows,
        )

    def rows_from_payload(
        self,
        payload: dict,
        *,
        exclude_existing: bool = True,
    ) -> tuple[list[dict[str, object]], list[str]]:
        rows, errors, _extracted_count, _duplicate_rows = self.rows_from_items(
            payload.get("results", []),
            exclude_existing=exclude_existing,
        )
        return rows, errors

    def rows_from_items(
        self,
        items: Iterable[dict],
        *,
        exclude_existing: bool,
    ) -> tuple[list[dict[str, object]], list[str], int, list[dict[str, object]]]:
        existing_ids, existing_leaders = self._existing_sets()
        rows: list[dict[str, object]] = []
        duplicate_rows: list[dict[str, object]] = []
        errors: list[str] = []
        seen_sequences: dict[str, str] = {}
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
            candidate = SignalPeptideCandidate(
                candidate_id=candidate_id,
                accession=accession,
                uniprot_id=item.get("uniProtkbId", ""),
                protein_name=_protein_name(item),
                organism_name=item.get("organism", {}).get("scientificName", ""),
                protein_sequence=sequence,
                protein_length=len(sequence),
                uniprot_signal_start=signal_start,
                uniprot_signal_end=signal_end,
                leader_sequence=signal,
                signal_peptide_sequence=signal,
                category="pichia_native_signal",
                category_label=PICHIA_NATIVE_SIGNAL_LABEL,
                processing_route="signal peptidase only",
                source_note="",
                rationale="UniProt 中带 signal peptide 注释的 Komagataella/Pichia 蛋白，适合作为外部候选草案。",
                caution="来自数据库自动发现；进入模型前需要人工确认切割位点、文献证据和是否适合 OPN。",
                library_stage="外部发现草案",
                source_type="UniProt",
                already_in_formal_library=already_in_formal_library,
            )
            row = candidate.as_row()
            row["source_note"] = (
                f"UniProt {row['accession']} {row['uniprot_id']}; "
                f"{row['organism_name']}; {row['protein_name']}"
            )
            if already_in_formal_library:
                duplicate_rows.append(
                    {
                        **row,
                        "duplicate_reason": "与正式候选库已有序列重复",
                        "duplicate_of": "formal_library",
                    }
                )
                if exclude_existing:
                    continue
            if signal in seen_sequences:
                duplicate_rows.append(
                    {
                        **row,
                        "duplicate_reason": "UniProt 结果中信号肽序列重复",
                        "duplicate_of": seen_sequences[signal],
                    }
                )
                continue
            seen_sequences[signal] = candidate_id
            rows.append(row)
        if not rows:
            errors.append("没有发现可加入草案的新 signal peptide；可能都已存在或查询结果没有明确序列。")
        return rows, errors, extracted_signal_count, duplicate_rows

    def _existing_sets(self) -> tuple[set[str], set[str]]:
        ids = set()
        leaders = set()
        for candidate in self.existing_candidates:
            ids.add(candidate.candidate_id)
            leaders.add(candidate.leader_sequence)
        return ids, leaders


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
