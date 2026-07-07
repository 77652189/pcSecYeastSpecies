from __future__ import annotations

from pathlib import Path

from pcsec_pichia.homology.blast_runner import find_blastp_executable, parse_blast_tsv, run_blastp
from pcsec_pichia.homology.cache_schema import BlastConfig


def test_parse_blast_tsv_computes_coverage_from_lengths(tmp_path: Path) -> None:
    path = tmp_path / "blast.tsv"
    path.write_text(
        "YJL034W\tPAS_chr2-1_0140\t75.53\t658\t682\t678\t1e-200\t900.5\n",
        encoding="ascii",
    )

    hits = parse_blast_tsv(path)

    assert len(hits) == 1
    hit = hits[0]
    assert hit.query_id == "YJL034W"
    assert hit.subject_id == "PAS_chr2-1_0140"
    assert hit.identity_pct == 75.53
    assert hit.evalue == 1e-200
    assert hit.bitscore == 900.5
    assert hit.query_coverage == round(100 * 658 / 682, 6)
    assert hit.subject_coverage == round(100 * 658 / 678, 6)


def test_parse_blast_tsv_rejects_unexpected_column_count(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("query\tsubject\t99\n", encoding="ascii")

    try:
        parse_blast_tsv(path)
    except ValueError as exc:
        assert "expected 8 BLAST columns" in str(exc)
    else:
        raise AssertionError("parse_blast_tsv should reject malformed rows")


def test_run_blastp_reports_unavailable_without_failing(tmp_path: Path) -> None:
    config = BlastConfig(blast_bin=tmp_path / "missing_blastp")

    result = run_blastp(tmp_path / "query.fasta", tmp_path / "db", tmp_path / "out.tsv", config)

    assert result.available is False
    assert result.message == "blastp unavailable"


def test_find_blastp_executable_accepts_configured_file(tmp_path: Path) -> None:
    exe = tmp_path / "blastp.exe"
    exe.write_text("", encoding="ascii")

    assert find_blastp_executable(BlastConfig(blast_bin=exe)) == exe
