from __future__ import annotations

from pathlib import Path

from scipy.io import savemat

from pcsec_pichia.homology.cache_schema import ProteinRecord
from pcsec_pichia.homology.sequence_sources import load_protein_sequences_from_mat, write_fasta


def test_load_protein_sequences_from_mat_reads_id_and_sequence(tmp_path: Path) -> None:
    path = tmp_path / "Protein_Sequence.mat"
    savemat(
        path,
        {
            "ProteinSequence": {
                "id": ["YJL034W", "YCL043C"],
                "seq": ["MKAR2*", "MPDI1"],
                "fullseq": ["MKAR2*", "MPDI1"],
            }
        },
    )

    records = load_protein_sequences_from_mat(path, "sce")

    assert records == (
        ProteinRecord(organism="sce", gene_id="YJL034W", sequence="MKAR2", source=str(path)),
        ProteinRecord(organism="sce", gene_id="YCL043C", sequence="MPDI1", source=str(path)),
    )


def test_write_fasta_is_deterministic_and_wraps_sequences(tmp_path: Path) -> None:
    path = tmp_path / "records.fasta"
    records = [
        ProteinRecord(organism="pichia", gene_id="b", sequence="M" * 81),
        ProteinRecord(organism="pichia", gene_id="a", sequence="ABC"),
    ]

    write_fasta(records, path)

    assert path.read_text(encoding="ascii") == ">a\nABC\n>b\n" + "M" * 80 + "\nM\n"
