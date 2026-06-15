from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PichiaClmAdapter:
    repo_path: Path
    device: str = "cpu"

    def is_available(self) -> tuple[bool, str]:
        if not self.repo_path.exists():
            return False, f"未找到 PichiaCLM 仓库：{self.repo_path}"
        weights = (
            self.repo_path
            / "Model_PichiaCLM"
            / "Training"
            / "PichiaData"
            / "2Target_AllData"
            / "Arch1-0404.weights.pt"
        )
        if not weights.exists():
            return False, f"未找到 PichiaCLM 权重文件：{weights}"
        try:
            self._ensure_import_path()
            import torch  # noqa: F401
            from Model_PichiaCLM.core.predictor import PichiaCLMPredictor  # noqa: F401
        except Exception as exc:
            return False, f"PichiaCLM 依赖不可用：{exc}"
        return True, "PichiaCLM 可用"

    def predict_candidates(
        self,
        amino_acids: str,
        *,
        num_candidates: int = 3,
        subset_size: int | None = 3,
        temperature: float = 0.8,
        seed: int | None = 42,
        motifs: Iterable[str] | None = None,
        custom_restriction_sites: Iterable[str] | None = None,
    ):
        self._ensure_import_path()
        from Model_PichiaCLM.core.predictor import PichiaCLMPredictor

        predictor = PichiaCLMPredictor(device=self.device)
        return predictor.predict_candidates(
            amino_acids,
            num_candidates=num_candidates,
            subset_size=subset_size,
            temperature=temperature,
            seed=seed,
            motifs=motifs,
            custom_restriction_sites=custom_restriction_sites,
        )

    def _ensure_import_path(self) -> None:
        repo = str(self.repo_path)
        if repo not in sys.path:
            sys.path.insert(0, repo)
