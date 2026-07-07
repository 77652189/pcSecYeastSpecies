from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from pcsec_pichia.homology.cache_schema import BlastConfig, BlastHit


@dataclass(frozen=True)
class BlastCommandResult:
    available: bool
    command: tuple[str, ...]
    output_path: Path | None = None
    returncode: int | None = None
    message: str = ""


def find_blastp_executable(config: BlastConfig | None = None) -> Path | None:
    """Locate local blastp; return None when unavailable."""

    candidates: list[Path] = []
    if config and config.blast_bin:
        blast_bin = Path(config.blast_bin)
        candidates.append(blast_bin / _exe_name("blastp"))
        candidates.append(blast_bin)
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        return None
    which = shutil.which("blastp")
    if which:
        candidates.append(Path(which))
    candidates.extend(_local_blast_candidates())
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def make_blast_db(fasta_path: Path, db_prefix: Path, blast_bin: Path | None = None) -> BlastCommandResult:
    makeblastdb = _sibling_tool(blast_bin, "makeblastdb") if blast_bin else _find_tool("makeblastdb")
    command = (str(makeblastdb), "-in", str(fasta_path), "-dbtype", "prot", "-out", str(db_prefix)) if makeblastdb else ()
    if makeblastdb is None:
        return BlastCommandResult(False, command, db_prefix, message="makeblastdb unavailable")
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    return BlastCommandResult(
        completed.returncode == 0,
        command,
        db_prefix,
        completed.returncode,
        (completed.stderr or completed.stdout).strip(),
    )


def run_blastp(
    query_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    config: BlastConfig,
) -> BlastCommandResult:
    blastp = find_blastp_executable(config)
    command = (
        str(blastp),
        "-query",
        str(query_fasta),
        "-db",
        str(db_prefix),
        "-out",
        str(out_tsv),
        "-outfmt",
        "6 qseqid sseqid pident length qlen slen evalue bitscore",
        "-evalue",
        str(config.max_evalue),
        "-max_target_seqs",
        str(config.max_target_seqs),
        "-num_threads",
        str(config.threads),
    ) if blastp else ()
    if blastp is None:
        return BlastCommandResult(False, command, out_tsv, message="blastp unavailable")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    return BlastCommandResult(
        completed.returncode == 0,
        command,
        out_tsv,
        completed.returncode,
        (completed.stderr or completed.stdout).strip(),
    )


def parse_blast_tsv(path: Path) -> tuple[BlastHit, ...]:
    """Parse BLAST outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore."""

    hits: list[BlastHit] = []
    with path.open("r", encoding="ascii", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            if len(parts) != 8:
                raise ValueError(f"{path}:{line_number} expected 8 BLAST columns, found {len(parts)}")
            hits.append(
                BlastHit.from_lengths(
                    query_id=parts[0],
                    subject_id=parts[1],
                    identity_pct=float(parts[2]),
                    alignment_length=int(parts[3]),
                    query_length=int(parts[4]),
                    subject_length=int(parts[5]),
                    evalue=float(parts[6]),
                    bitscore=float(parts[7]),
                )
            )
    return tuple(hits)


def _find_tool(name: str) -> Path | None:
    which = shutil.which(name)
    if which:
        return Path(which)
    for blastp in _local_blast_candidates():
        sibling = _sibling_tool(blastp, name)
        if sibling:
            return sibling
    return None


def _sibling_tool(blast_bin: Path | None, name: str) -> Path | None:
    if blast_bin is None:
        return None
    path = Path(blast_bin)
    if path.is_dir():
        candidate = path / _exe_name(name)
    else:
        candidate = path.with_name(_exe_name(name))
    return candidate if candidate.is_file() else None


def _local_blast_candidates() -> list[Path]:
    root = Path(__file__).resolve().parents[4]
    return sorted(root.glob(f"local_runs/blast_homolog_feasibility/bin/**/{_exe_name('blastp')}"))


def _exe_name(name: str) -> str:
    return f"{name}.exe" if "\\" in str(Path.cwd()) else name
