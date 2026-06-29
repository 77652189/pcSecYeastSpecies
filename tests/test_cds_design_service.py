from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import pandas as pd

from pcsec_pichia.core.paths import ProjectPaths
from app.services.cds_design import CODON_TO_AA, CdsDesignService


AA_TO_CODON = {}
for codon, amino_acid in CODON_TO_AA.items():
    if amino_acid != "*" and amino_acid not in AA_TO_CODON:
        AA_TO_CODON[amino_acid] = codon


@dataclass
class FakePichiaClmAdapter:
    def is_available(self) -> tuple[bool, str]:
        return True, "ok"

    def predict_candidates(self, amino_acids: str, **kwargs):
        num_candidates = int(kwargs.get("num_candidates", 1))
        candidates = []
        for rank in range(1, num_candidates + 1):
            cds = _cds_for_amino_acids(amino_acids)
            analysis = SimpleNamespace(
                cds_length=len(cds),
                gc_percent=50.0,
                gc_status="ok",
                cai=SimpleNamespace(training=0.71, public=0.69),
                restriction_sites=[] if rank == 1 else ["EcoRI"],
                motif_hits=[],
                internal_stop_codons=[],
            )
            candidates.append(
                SimpleNamespace(
                    rank=rank,
                    source="reference" if rank == 1 else "kazusa_diverse",
                    cds=cds,
                    analysis=analysis,
                    quality=SimpleNamespace(status="pass", warnings=0),
                )
            )
        return SimpleNamespace(
            amino_acids=amino_acids,
            candidates=candidates,
            recommended_subset=SimpleNamespace(selected_ranks=[1, 2, 3]),
        )


@dataclass
class MissingPichiaClmAdapter:
    def is_available(self) -> tuple[bool, str]:
        return False, "未找到 PichiaCLM 仓库"


def test_cds_design_service_generates_records_and_exports_for_opn_shortlist() -> None:
    paths = ProjectPaths.discover()
    result = CdsDesignService(paths, FakePichiaClmAdapter()).design_opn_shortlist(
        cds_candidates_per_construct=3,
    )

    assert result.available is True
    assert len(result.records) == 9
    assert {record.construct_id for record in result.records} == {
        "OPN_PPA_PASCHR3_0030",
        "OPN_PPA_DDDK18",
        "OPN_ALPHA_FULL_PROJECT",
    }
    record = result.records[0]
    assert record.recommended_subset is True
    assert record.cds_length == record.aa_length * 3
    assert len(record.cds) % 3 == 0
    assert record.length_multiple_of_three is True
    assert record.translation_matches_input is True
    assert record.internal_stop_codons == []
    assert record.experimental_role
    assert record.leader_sequence
    assert record.full_protein_sequence
    assert result.csv_file is not None and result.csv_file.exists()
    assert result.xlsx_file is not None and result.xlsx_file.exists()
    assert result.fasta_file is not None and result.fasta_file.exists()
    csv = pd.read_csv(result.csv_file)
    assert "候选ID" in csv.columns
    assert "CDS序列" in csv.columns
    assert result.fasta_file.read_text(encoding="utf-8").startswith(">OPN_PPA_PASCHR3_0030|cds_candidate_rank=1")


def test_cds_design_service_returns_chinese_error_when_pichia_clm_missing() -> None:
    paths = ProjectPaths.discover()
    result = CdsDesignService(paths, MissingPichiaClmAdapter()).design_opn_shortlist()

    assert result.available is False
    assert "未找到 PichiaCLM" in result.message


def _cds_for_amino_acids(amino_acids: str) -> str:
    return "".join(AA_TO_CODON[amino_acid] for amino_acid in amino_acids)
