from __future__ import annotations

from pathlib import Path

from app.adapters.uspnet import USPNetPrediction, USPNetRunResult
from app.core.paths import ProjectPaths
from app.services.signal_peptide_library import UniProtCandidateLibraryResult
from app.services.signal_peptide_screening import SignalPeptideScreeningService


def test_signal_peptide_screening_compares_rules_and_uspnet(tmp_path: Path) -> None:
    service = SignalPeptideScreeningService(
        ProjectPaths(tmp_path),
        library_service=FakeLibraryService(),
        uspnet_adapter=FakeUSPNetAdapter(),
    )

    result = service.screen_uniprot_candidates(max_records=2)

    assert result.success is True
    assert result.summary["uniprot_initial_hits"] == 2
    assert result.summary["deduplicated_candidates"] == 2
    assert result.summary["uniprot_duplicate_count"] == 0
    assert result.summary["rules_passed"] == 1
    assert result.summary["rules_high_priority"] == 1
    assert result.summary["uspnet_completed"] == 2
    assert result.summary["uspnet_passed"] == 1
    assert result.summary["consensus_passed"] == 1
    assert result.uniprot_csv and result.uniprot_csv.exists()
    assert result.duplicate_csv and result.duplicate_csv.exists()
    assert result.input_fasta and result.input_fasta.exists()
    assert result.comparison_csv and result.comparison_csv.exists()
    assert result.recommended_fasta and result.recommended_fasta.exists()
    assert "OPN_UNIPROT_X12345" in result.recommended_fasta.read_text(encoding="utf-8")


def test_signal_peptide_screening_reuses_persisted_uniprot_candidates(tmp_path: Path) -> None:
    library_service = CountingFakeLibraryService()
    service = SignalPeptideScreeningService(
        ProjectPaths(tmp_path),
        library_service=library_service,
        uspnet_adapter=FakeMissingUSPNetAdapter(),
    )

    service.discover_and_persist_uniprot_candidates(max_records=2)
    result = service.screen_uniprot_candidates(max_records=2)

    assert library_service.calls == 1
    assert result.summary["uniprot_reused_from_disk"] is True
    assert result.summary["uniprot_candidate_source"] == "已复用本地保存的 UniProt 候选"
    assert result.uniprot_csv and result.uniprot_csv.exists()


def test_signal_peptide_screening_loads_saved_result(tmp_path: Path) -> None:
    service = SignalPeptideScreeningService(
        ProjectPaths(tmp_path),
        library_service=FakeLibraryService(),
        uspnet_adapter=FakeUSPNetAdapter(),
    )

    result = service.screen_uniprot_candidates(max_records=2)
    loaded = service.load_persisted_screening_result()

    assert loaded is not None
    assert loaded.summary["consensus_passed"] == result.summary["consensus_passed"]
    assert loaded.rows[0]["rules_score_note"]
    assert loaded.rows[0]["uspnet_prediction_label"] == "SP：USPNet 判断为信号肽"


def test_signal_peptide_screening_keeps_rule_prefilter_when_uspnet_missing(tmp_path: Path) -> None:
    service = SignalPeptideScreeningService(
        ProjectPaths(tmp_path),
        library_service=FakeLibraryService(),
        uspnet_adapter=FakeMissingUSPNetAdapter(),
    )

    result = service.screen_uniprot_candidates(max_records=2)

    assert result.available is False
    assert result.success is True
    assert "USPNet 尚未安装" in result.message
    assert result.summary["rules_passed"] == 1
    assert result.summary["rules_high_priority"] == 1
    assert result.summary["uspnet_completed"] == 0
    assert result.summary["consensus_passed"] == 0
    assert result.summary["needs_external_review"] == 1
    assert result.comparison_csv and result.comparison_csv.exists()
    assert result.recommended_fasta and result.recommended_fasta.exists()


def test_signal_peptide_screening_keeps_rule_prefilter_when_uspnet_run_fails(tmp_path: Path) -> None:
    service = SignalPeptideScreeningService(
        ProjectPaths(tmp_path),
        library_service=FakeLibraryService(),
        uspnet_adapter=FakeFailingUSPNetAdapter(),
    )

    result = service.screen_uniprot_candidates(max_records=2)

    assert result.available is True
    assert result.success is True
    assert "USPNet 已检测到，但本次运行未完成" in result.message
    assert result.summary["uspnet_available"] is True
    assert result.summary["uspnet_success"] is False
    assert result.summary["uspnet_completed"] == 0
    assert result.summary["consensus_passed"] == 0
    assert result.summary["needs_external_review"] == 1

    by_id = {row["candidate_id"]: row for row in result.rows}
    high_priority = by_id["OPN_UNIPROT_X12345"]
    assert high_priority["screening_status"] == "规则高优先级，待 USPNet 复核"
    assert high_priority["recommended_for_draft_library"] is True
    assert result.recommended_fasta and "OPN_UNIPROT_X12345" in result.recommended_fasta.read_text(encoding="utf-8")


class FakeLibraryService:
    def discover_uniprot_candidate_library(self, **_kwargs):
        rows = [
            {
                "candidate_id": "OPN_UNIPROT_X12345",
                "accession": "X12345",
                "uniprot_id": "TEST1_PICPA",
                "protein_name": "Secreted test protein",
                "organism_name": "Komagataella phaffii",
                "protein_sequence": "MKALLLALLALAAASAGAQREST",
                "protein_length": 24,
                "uniprot_signal_start": 1,
                "uniprot_signal_end": 18,
                "leader_sequence": "MKALLLALLALAAASAGA",
                "signal_peptide_sequence": "MKALLLALLALAAASAGA",
                "category": "pichia_native_signal",
                "processing_route": "signal peptidase only",
                "source_note": "UniProt X12345",
                "rationale": "fixture",
                "caution": "fixture",
                "leader_length": 18,
                "signal_peptide_length": 18,
                "library_stage": "外部发现草案",
                "source_type": "UniProt",
                "already_in_formal_library": False,
            },
            {
                "candidate_id": "OPN_UNIPROT_Y12345",
                "accession": "Y12345",
                "uniprot_id": "TEST2_PICPA",
                "protein_name": "Low complexity test protein",
                "organism_name": "Komagataella phaffii",
                "protein_sequence": "MNNNNNNNNNNNNNNNNNNNN",
                "protein_length": 20,
                "uniprot_signal_start": 1,
                "uniprot_signal_end": 18,
                "leader_sequence": "MNNNNNNNNNNNNNNNNN",
                "signal_peptide_sequence": "MNNNNNNNNNNNNNNNNN",
                "category": "pichia_native_signal",
                "processing_route": "signal peptidase only",
                "source_note": "UniProt Y12345",
                "rationale": "fixture",
                "caution": "fixture",
                "leader_length": 18,
                "signal_peptide_length": 18,
                "library_stage": "外部发现草案",
                "source_type": "UniProt",
                "already_in_formal_library": False,
            },
        ]
        return UniProtCandidateLibraryResult(
            rows=rows,
            source_url="https://rest.uniprot.org/fixture",
            errors=[],
            initial_hit_count=2,
            fetched_record_count=2,
            extracted_signal_count=2,
            deduplicated_count=2,
        )


class CountingFakeLibraryService(FakeLibraryService):
    def __init__(self) -> None:
        self.calls = 0

    def discover_uniprot_candidate_library(self, **kwargs):
        self.calls += 1
        return super().discover_uniprot_candidate_library(**kwargs)


class FakeUSPNetAdapter:
    def run(self, _fasta_file: Path, output_dir: Path, **_kwargs):
        output_dir.mkdir(parents=True, exist_ok=True)
        return USPNetRunResult(
            available=True,
            success=True,
            message="USPNet fixture",
            output_dir=output_dir,
            predictions=[
                USPNetPrediction(
                    candidate_id="OPN_UNIPROT_X12345",
                    predicted_type="SP",
                    predicted_cleavage="MKALLLALLALAAASAGA",
                    passed=True,
                    raw_sequence="MKALLLALLALAAASAGAQREST",
                ),
                USPNetPrediction(
                    candidate_id="OPN_UNIPROT_Y12345",
                    predicted_type="NO_SP",
                    predicted_cleavage="",
                    passed=False,
                    raw_sequence="MNNNNNNNNNNNNNNNNNNNN",
                ),
            ],
        )


class FakeMissingUSPNetAdapter:
    def run(self, _fasta_file: Path, output_dir: Path, **_kwargs):
        output_dir.mkdir(parents=True, exist_ok=True)
        return USPNetRunResult(
            available=False,
            success=False,
            message="未检测到 USPNet 本地仓库。",
            output_dir=output_dir,
            predictions=[],
        )


class FakeFailingUSPNetAdapter:
    def run(self, _fasta_file: Path, output_dir: Path, **_kwargs):
        output_dir.mkdir(parents=True, exist_ok=True)
        return USPNetRunResult(
            available=True,
            success=False,
            message="USPNet runtime failed",
            output_dir=output_dir,
            predictions=[],
        )
