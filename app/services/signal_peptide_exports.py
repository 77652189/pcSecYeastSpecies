from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(rows)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def write_candidate_fasta(path: Path, rows: Iterable[dict[str, object]]) -> None:
    records = []
    for row in rows:
        sequence = str(row.get("protein_sequence", ""))
        if not sequence:
            continue
        header = f"{row.get('candidate_id')}|accession={row.get('accession')}|source=UniProt"
        records.append((header, sequence))
    write_fasta(path, records)


def write_signal_peptide_fasta(path: Path, rows: Iterable[dict[str, object]]) -> None:
    records = []
    for row in rows:
        sequence = str(row.get("signal_peptide_sequence", ""))
        if not sequence:
            continue
        header = (
            f"{row.get('candidate_id')}|accession={row.get('accession')}|"
            f"status={row.get('screening_status')}|rules={row.get('rules_score')}"
        )
        records.append((header, sequence))
    write_fasta(path, records)


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
