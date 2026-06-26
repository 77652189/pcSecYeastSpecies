"""pcSecYeastSpecies local web application package."""

from __future__ import annotations

import sys
from pathlib import Path


def ensure_python_pichia_on_path() -> None:
    """Make the local python_pichia package importable from any app entrypoint."""
    repo_root = Path(__file__).resolve().parents[1]
    src = repo_root / "python_pichia" / "src"
    if src.exists():
        src_text = str(src)
        if src_text not in sys.path:
            sys.path.insert(0, src_text)


ensure_python_pichia_on_path()


__all__ = ["ensure_python_pichia_on_path"]
