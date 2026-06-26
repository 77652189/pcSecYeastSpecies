from __future__ import annotations

import sys
from pathlib import Path


PYTHON_PICHIA_SRC = Path(__file__).resolve().parents[1] / "python_pichia" / "src"

if str(PYTHON_PICHIA_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_PICHIA_SRC))
