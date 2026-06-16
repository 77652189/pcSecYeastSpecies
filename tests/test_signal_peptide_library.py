from __future__ import annotations

from app.core.paths import ProjectPaths
from app.services.opn import OpnCandidateCatalog
from app.services.signal_peptide_library import SignalPeptideLibraryService, _next_link


def test_signal_peptide_library_labels_current_candidates() -> None:
    service = SignalPeptideLibraryService(OpnCandidateCatalog(ProjectPaths.discover()))
    rows = service.library_rows()
    by_id = {row["candidate_id"]: row for row in rows}

    assert len(rows) >= 7
    assert by_id["OPN_PPA_PASCHR3_0030"]["library_stage"] == "首轮推荐"
    assert by_id["OPN_ALPHA_FULL_PROJECT"]["library_stage"] == "首轮推荐"
    assert by_id["OPN_PPA_DDDK18"]["source_type"] == "文献"


def test_signal_peptide_library_validates_import_csv() -> None:
    service = SignalPeptideLibraryService(OpnCandidateCatalog(ProjectPaths.discover()))
    content = (
        "candidate_id,leader_sequence,signal_peptide_sequence,category,processing_route,"
        "source_note,rationale,caution\n"
        "OPN_NEW_SIGNAL_001,MKFAISTLLIILQAAAVFAA,MKFAISTLLIILQAAAVFAA,"
        "pichia_native_signal,signal peptidase only,UniProt candidate,"
        "adds a new Pichia-native candidate,needs external confirmation\n"
    ).encode("utf-8")

    result = service.validate_import_csv(content)

    assert result.valid is True
    assert result.rows[0]["candidate_id"] == "OPN_NEW_SIGNAL_001"
    assert result.rows[0]["leader_length"] == 20
    assert service.merged_draft_csv(result.rows).startswith(b"\xef\xbb\xbf")


def test_signal_peptide_library_rejects_duplicate_and_bad_sequence() -> None:
    service = SignalPeptideLibraryService(OpnCandidateCatalog(ProjectPaths.discover()))
    content = (
        "candidate_id,leader_sequence,signal_peptide_sequence,category,processing_route,"
        "source_note,rationale,caution\n"
        "OPN_PPA_DDDK18,MKFAI123,MKFAI,pichia_native_signal,signal peptidase only,"
        "duplicate,bad sequence,bad\n"
    ).encode("utf-8")

    result = service.validate_import_csv(content)

    assert result.valid is False
    assert any("已存在" in error for error in result.errors)
    assert any("标准氨基酸" in error for error in result.errors)


def test_signal_peptide_library_extracts_uniprot_signal_features() -> None:
    service = SignalPeptideLibraryService(OpnCandidateCatalog(ProjectPaths.discover()))
    payload = {
        "results": [
            {
                "primaryAccession": "X12345",
                "uniProtkbId": "TEST_PICPA",
                "organism": {"scientificName": "Komagataella pastoris"},
                "proteinDescription": {
                    "recommendedName": {"fullName": {"value": "Secreted test protein"}}
                },
                "sequence": {"value": "MKTLLALALALAAPAAQREST"},
                "features": [
                    {
                        "type": "Signal",
                        "location": {"start": {"value": 1}, "end": {"value": 16}},
                    }
                ],
            }
        ]
    }

    rows, errors = service.rows_from_uniprot_payload(payload)

    assert errors == []
    assert rows[0]["candidate_id"] == "OPN_UNIPROT_X12345"
    assert rows[0]["accession"] == "X12345"
    assert rows[0]["protein_sequence"] == "MKTLLALALALAAPAAQREST"
    assert rows[0]["uniprot_signal_end"] == 16
    assert rows[0]["signal_peptide_sequence"] == "MKTLLALALALAAPAA"
    assert rows[0]["source_type"] == "UniProt"
    assert rows[0]["library_stage"] == "外部发现草案"


def test_uniprot_next_link_parser_handles_commas_inside_url() -> None:
    header = (
        '<https://rest.uniprot.org/uniprotkb/search?format=json&'
        'fields=accession,id,protein_name,organism_name,ft_signal,sequence&'
        'cursor=abc&size=100>; rel="next"'
    )

    assert _next_link(header) == (
        "https://rest.uniprot.org/uniprotkb/search?format=json&"
        "fields=accession,id,protein_name,organism_name,ft_signal,sequence&"
        "cursor=abc&size=100"
    )
