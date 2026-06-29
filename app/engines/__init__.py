"""Legacy MATLAB engine facades used by old app services.

The current pcSecPichia Python draft workbench must enter through
``app.services.pichia_secretion_service`` and the formal ``python_pichia``
pipeline. Keep this package as a reference/compatibility boundary for legacy
MATLAB OPN smoke paths only.
"""

from pcsec_pichia.engines.base import PichiaEngine, PichiaSimulationRequest, PichiaSimulationRunResult
from app.engines.matlab_pichia_engine import MatlabPichiaEngine

__all__ = [
    "MatlabPichiaEngine",
    "PichiaEngine",
    "PichiaSimulationRequest",
    "PichiaSimulationRunResult",
]
