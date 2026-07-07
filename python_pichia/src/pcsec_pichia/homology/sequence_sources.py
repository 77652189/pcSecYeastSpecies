from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat

from pcsec_pichia.homology.cache_schema import ProteinRecord
from pcsec_pichia.loading import load_pcsec_pichia_inputs, repo_root


def load_protein_sequences_from_mat(path: Path, organism: str) -> tuple[ProteinRecord, ...]:
    """Load protein records from local pcSec Protein_Sequence.mat assets."""

    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "ProteinSequence" not in data:
        available = ", ".join(sorted(key for key in data if not key.startswith("__")))
        raise KeyError(f"{path} does not contain MATLAB variable 'ProteinSequence'. Available: {available}")
    protein_sequence = data["ProteinSequence"]
    ids = _string_items(getattr(protein_sequence, "id"))
    sequences = _string_items(getattr(protein_sequence, "seq", getattr(protein_sequence, "fullseq", [])))
    records: list[ProteinRecord] = []
    for gene_id, sequence in zip(ids, sequences):
        cleaned = _clean_sequence(sequence)
        if not gene_id or not cleaned:
            continue
        records.append(
            ProteinRecord(
                organism=organism,
                gene_id=gene_id,
                sequence=cleaned,
                source=str(path),
            )
        )
    return tuple(records)


def load_pichia_model_gene_index(root: Path | None = None) -> set[str]:
    """Return gene ids currently present in the Pichia GEM gene index."""

    inputs = load_pcsec_pichia_inputs(root or repo_root())
    return set(inputs.model.gene_index)


def write_fasta(records: Iterable[ProteinRecord], path: Path) -> Path:
    """Write deterministic FASTA used by local BLAST+ runs."""

    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records, key=lambda record: (record.organism, record.gene_id))
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for record in ordered:
            handle.write(f">{record.gene_id}\n")
            for start in range(0, len(record.sequence), 80):
                handle.write(record.sequence[start : start + 80] + "\n")
    return path


def _string_items(value: Any) -> list[str]:
    array = np.asarray(value, dtype=object)
    if array.shape == ():
        return [_string_scalar(array.item())]
    return [_string_scalar(item) for item in array.reshape(-1)]


def _string_scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _string_scalar(value.item())
        if value.dtype.kind in {"U", "S"}:
            return "".join(str(item) for item in value.reshape(-1)).strip()
        return " ".join(_string_scalar(item) for item in value.reshape(-1)).strip()
    return str(value).strip()


def _clean_sequence(sequence: str) -> str:
    return "".join(sequence.replace("*", "").split()).upper()
