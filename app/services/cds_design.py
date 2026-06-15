from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from app.adapters.pichia_clm import PichiaClmAdapter
from app.core.models import CdsDesignRecord, CdsDesignResult
from app.core.paths import ProjectPaths
from app.services.opn import OPN_SHORTLIST, OpnCandidateCatalog


DEFAULT_FORBIDDEN_MOTIFS = ("AATAAA",)
DEFAULT_RESTRICTION_SITES = ("EcoRI=GAATTC", "NotI=GCGGCCGC")
CONSTRUCT_PLAN_CSV = "OPN_Pichia_construct_plan.csv"
CONSTRUCT_PLAN_XLSX = "OPN_Pichia_construct_plan.xlsx"
CONSTRUCT_PLAN_FASTA = "OPN_Pichia_construct_sequences.fasta"

CODON_TO_AA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


@dataclass
class CdsDesignService:
    paths: ProjectPaths
    adapter: PichiaClmAdapter | None = None

    def __post_init__(self) -> None:
        if self.adapter is None:
            self.adapter = PichiaClmAdapter(self.paths.pichia_clm_repo)

    def design_opn_shortlist(
        self,
        candidate_ids: list[str] | None = None,
        *,
        cds_candidates_per_construct: int = 3,
        seed: int = 42,
        write_files: bool = True,
    ) -> CdsDesignResult:
        assert self.adapter is not None
        ok, message = self.adapter.is_available()
        if not ok:
            return CdsDesignResult(available=False, message=message)

        designs = {
            design.candidate_id: design
            for design in OpnCandidateCatalog(self.paths).construct_designs()
        }
        requested = candidate_ids or list(OPN_SHORTLIST)
        records: list[CdsDesignRecord] = []
        errors: list[str] = []
        for construct_id in requested:
            design = designs.get(construct_id)
            if design is None:
                errors.append(f"{construct_id}: 未找到 OPN 构建设计")
                continue
            try:
                candidate_set = self.adapter.predict_candidates(
                    design.full_protein_sequence,
                    num_candidates=cds_candidates_per_construct,
                    subset_size=min(3, cds_candidates_per_construct),
                    seed=seed,
                    motifs=DEFAULT_FORBIDDEN_MOTIFS,
                    custom_restriction_sites=DEFAULT_RESTRICTION_SITES,
                )
            except Exception as exc:
                errors.append(f"{construct_id}: {exc}")
                continue
            selected_ranks = set()
            if candidate_set.recommended_subset is not None:
                selected_ranks = set(candidate_set.recommended_subset.selected_ranks)
            for candidate in candidate_set.candidates:
                records.append(
                    self._record_from_candidate(
                        design=design,
                        candidate=candidate,
                        expected_amino_acids=candidate_set.amino_acids,
                        recommended_subset=candidate.rank in selected_ranks,
                    )
                )
        if records:
            suffix = f"；部分构建失败：{'；'.join(errors)}" if errors else ""
            csv_file = xlsx_file = fasta_file = None
            if write_files:
                csv_file, xlsx_file, fasta_file = self.write_exports(records)
            return CdsDesignResult(
                available=True,
                message=f"已生成 {len(records)} 条毕赤酵母 CDS 候选{suffix}",
                records=records,
                csv_file=csv_file,
                xlsx_file=xlsx_file,
                fasta_file=fasta_file,
            )
        return CdsDesignResult(
            available=True,
            message="未生成 CDS 候选。" + ("；".join(errors) if errors else ""),
        )

    def write_exports(self, records: list[CdsDesignRecord]) -> tuple[Path, Path, Path]:
        self.paths.opn_design_dir.mkdir(parents=True, exist_ok=True)
        csv_file = self.paths.opn_design_dir / CONSTRUCT_PLAN_CSV
        xlsx_file = self.paths.opn_design_dir / CONSTRUCT_PLAN_XLSX
        fasta_file = self.paths.opn_design_dir / CONSTRUCT_PLAN_FASTA

        frame = pd.DataFrame([self._export_row(record) for record in records])
        frame.to_csv(csv_file, index=False, encoding="utf-8-sig")
        frame.to_excel(xlsx_file, index=False)
        fasta_file.write_text(self._format_fasta(records), encoding="utf-8")
        return csv_file, xlsx_file, fasta_file

    def _record_from_candidate(
        self,
        *,
        design,
        candidate: Any,
        expected_amino_acids: str,
        recommended_subset: bool,
    ) -> CdsDesignRecord:
        translation = _check_translation(candidate.cds, expected_amino_acids)
        analysis = candidate.analysis
        cai = getattr(analysis, "cai", None)
        internal_stop_codons = _as_int_list(
            getattr(analysis, "internal_stop_codons", translation["internal_stop_codons"])
        )
        if not internal_stop_codons:
            internal_stop_codons = translation["internal_stop_codons"]
        return CdsDesignRecord(
            construct_id=design.candidate_id,
            experimental_role=design.experimental_role,
            recommendation=design.recommendation,
            cds_candidate_rank=int(candidate.rank),
            recommended_subset=recommended_subset,
            leader_sequence=design.leader_sequence,
            mature_opn_sequence=design.mature_opn_sequence,
            full_protein_sequence=design.full_protein_sequence,
            aa_length=len(expected_amino_acids),
            cds_length=int(getattr(analysis, "cds_length", len(candidate.cds))),
            gc_percent=getattr(analysis, "gc_percent", None),
            gc_status=str(getattr(analysis, "gc_status", "")),
            cai_training=getattr(cai, "training", None),
            cai_public=getattr(cai, "public", None),
            quality_status=str(getattr(candidate.quality, "status", "")),
            warnings=int(getattr(candidate.quality, "warnings", 0)),
            restriction_sites=len(getattr(analysis, "restriction_sites", [])),
            motif_hits=len(getattr(analysis, "motif_hits", [])),
            length_multiple_of_three=translation["length_multiple_of_three"],
            translation_matches_input=translation["translation_matches_input"],
            internal_stop_codons=internal_stop_codons,
            kex2_risk=design.kex2_risk,
            risk_note=design.note,
            source=str(candidate.source),
            cds=candidate.cds,
        )

    def _export_row(self, record: CdsDesignRecord) -> dict[str, object]:
        return {
            "候选ID": record.construct_id,
            "实验角色": record.experimental_role,
            "推荐级别": record.recommendation,
            "CDS候选序号": record.cds_candidate_rank,
            "推荐保留": _yes_no(record.recommended_subset),
            "leader序列": record.leader_sequence,
            "成熟OPN序列": record.mature_opn_sequence,
            "完整蛋白序列": record.full_protein_sequence,
            "CDS序列": record.cds,
            "氨基酸长度": record.aa_length,
            "CDS长度bp": record.cds_length,
            "GC百分比": record.gc_percent,
            "GC状态": record.gc_status,
            "CAI训练参考": record.cai_training,
            "CAI公开毕赤参考": record.cai_public,
            "质控状态": record.quality_status,
            "提醒数": record.warnings,
            "限制性位点数": record.restriction_sites,
            "不期望motif数": record.motif_hits,
            "长度是否为3倍数": _yes_no(record.length_multiple_of_three),
            "回译是否匹配蛋白序列": _yes_no(record.translation_matches_input),
            "内部终止密码子位置": ",".join(str(item) for item in record.internal_stop_codons),
            "Kex2风险": record.kex2_risk,
            "风险说明": record.risk_note,
            "CDS来源": record.source,
        }

    def _format_fasta(self, records: list[CdsDesignRecord]) -> str:
        chunks: list[str] = []
        for record in records:
            header = (
                f">{record.construct_id}|cds_candidate_rank={record.cds_candidate_rank}|"
                f"role={record.experimental_role}|recommended={_yes_no(record.recommended_subset)}"
            )
            chunks.append(header)
            chunks.extend(_wrap_sequence(record.cds))
        return "\n".join(chunks) + "\n"


def _check_translation(cds: str, expected_amino_acids: str) -> dict[str, Any]:
    normalized = "".join(cds.split()).upper().replace("U", "T")
    codons = [normalized[index : index + 3] for index in range(0, len(normalized) - 2, 3)]
    translated = "".join(CODON_TO_AA.get(codon, "X") for codon in codons)
    internal_stops = [
        index
        for index, codon in enumerate(codons, start=1)
        if CODON_TO_AA.get(codon) == "*" and index < len(codons)
    ]
    return {
        "length_multiple_of_three": len(normalized) % 3 == 0,
        "translated_amino_acids": translated,
        "translation_matches_input": translated == expected_amino_acids,
        "internal_stop_codons": internal_stops,
    }


def _as_int_list(values: Any) -> list[int]:
    if values is None:
        return []
    return [int(value) for value in values]


def _yes_no(value: bool) -> str:
    return "是" if value else "否"


def _wrap_sequence(sequence: str, width: int = 70) -> list[str]:
    return [sequence[index : index + width] for index in range(0, len(sequence), width)]
