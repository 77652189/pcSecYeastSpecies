from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

from app.adapters.process_runner import CommandResult, ProcessRunner


USPNET_REPO_URL = "https://github.com/ml4bio/USPNet"
USPNET_MODEL_FILE = "USPNet_fast_no_group_info.pth"


@dataclass(frozen=True)
class USPNetStatus:
    available: bool
    repo_dir: Path | None = None
    model_dir: Path | None = None
    message: str = ""


@dataclass(frozen=True)
class USPNetPrediction:
    candidate_id: str
    predicted_type: str
    predicted_cleavage: str
    passed: bool
    raw_sequence: str


@dataclass(frozen=True)
class USPNetRunResult:
    available: bool
    success: bool
    message: str
    output_dir: Path
    predictions: list[USPNetPrediction]
    results_csv: Path | None = None
    command_outputs: list[CommandResult] | None = None


class USPNetAdapter:
    def __init__(
        self,
        repo_dir: str | Path | None = None,
        model_dir: str | Path | None = None,
        python_bin: str | Path | None = None,
        runner: ProcessRunner | None = None,
    ) -> None:
        self.repo_dir = Path(repo_dir) if repo_dir is not None else None
        self.model_dir = Path(model_dir) if model_dir is not None else None
        self.python_bin = str(python_bin) if python_bin is not None else os.environ.get("USPNET_PYTHON")
        self.runner = runner or ProcessRunner()

    def status(self) -> USPNetStatus:
        repo_dir = self._repo_dir()
        model_dir = self._model_dir(repo_dir)
        if repo_dir is None or not (repo_dir / "predict_fast.py").exists() or not (repo_dir / "data_processing.py").exists():
            return USPNetStatus(
                available=False,
                repo_dir=repo_dir,
                model_dir=model_dir,
                message=(
                    "未检测到 USPNet 本地仓库。USPNet 是 MIT License，可作为商业友好的外部复核工具；"
                    f"请克隆 {USPNET_REPO_URL}，并按 README 下载 USPNet-fast no_group_info 模型权重。"
                ),
            )
        if model_dir is None or not (model_dir / USPNET_MODEL_FILE).exists():
            return USPNetStatus(
                available=False,
                repo_dir=repo_dir,
                model_dir=model_dir,
                message=(
                    f"已检测到 USPNet 仓库，但缺少模型权重 {USPNET_MODEL_FILE}。"
                    "请把权重放在 USPNet 根目录，或设置 USPNET_MODEL_DIR。"
                ),
            )
        return USPNetStatus(
            available=True,
            repo_dir=repo_dir,
            model_dir=model_dir,
            message=f"已检测到 USPNet-fast：{repo_dir}",
        )

    def run(
        self,
        fasta_file: Path,
        output_dir: Path,
        *,
        timeout_seconds: int = 3600,
    ) -> USPNetRunResult:
        status = self.status()
        output_dir.mkdir(parents=True, exist_ok=True)
        if not status.available or status.repo_dir is None or status.model_dir is None:
            return USPNetRunResult(
                available=False,
                success=False,
                message=status.message,
                output_dir=output_dir,
                predictions=[],
            )

        data_dir = output_dir / "data_processed"
        process_result = self.runner.run(
            [
                self._python_bin(status.repo_dir),
                "data_processing.py",
                "--fasta_file",
                str(fasta_file),
                "--data_processed_dir",
                str(data_dir),
            ],
            cwd=status.repo_dir,
            timeout_seconds=timeout_seconds,
        )
        if process_result.returncode != 0:
            return USPNetRunResult(
                available=True,
                success=False,
                message=f"USPNet 数据预处理失败：{process_result.combined_output}",
                output_dir=output_dir,
                predictions=[],
                command_outputs=[process_result],
            )

        predict_result = self.runner.run(
            [
                self._python_bin(status.repo_dir),
                "predict_fast.py",
                "--data_dir",
                str(data_dir),
                "--group_info",
                "no_group_info",
                "--model_dir",
                str(status.model_dir),
            ],
            cwd=status.repo_dir,
            timeout_seconds=timeout_seconds,
        )
        results_csv = data_dir / "results.csv"
        predictions = parse_uspnet_results(results_csv, _fasta_ids(fasta_file)) if results_csv.exists() else []
        success = predict_result.returncode == 0 and bool(predictions)
        message = (
            f"USPNet-fast 分析完成，解析到 {len(predictions)} 条预测结果。"
            if success
            else f"USPNet-fast 未生成可解析结果：{predict_result.combined_output}"
        )
        return USPNetRunResult(
            available=True,
            success=success,
            message=message,
            output_dir=output_dir,
            predictions=predictions,
            results_csv=results_csv if results_csv.exists() else None,
            command_outputs=[process_result, predict_result],
        )

    def _repo_dir(self) -> Path | None:
        configured = self.repo_dir or _env_path("USPNET_REPO")
        return configured if configured is not None else None

    def _model_dir(self, repo_dir: Path | None) -> Path | None:
        configured = self.model_dir or _env_path("USPNET_MODEL_DIR")
        if configured is not None:
            return configured
        return repo_dir

    def _python_bin(self, repo_dir: Path) -> str:
        if self.python_bin:
            return self.python_bin
        venv_python = repo_dir / ".venv" / "Scripts" / "python.exe"
        if venv_python.exists():
            return str(venv_python)
        return "python"


def parse_uspnet_results(results_csv: Path, candidate_ids: list[str]) -> list[USPNetPrediction]:
    with results_csv.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    predictions: list[USPNetPrediction] = []
    for index, row in enumerate(rows):
        candidate_id = candidate_ids[index] if index < len(candidate_ids) else f"row_{index + 1}"
        predicted_type = str(row.get("predicted_type", "")).strip()
        predictions.append(
            USPNetPrediction(
                candidate_id=candidate_id,
                predicted_type=predicted_type,
                predicted_cleavage=str(row.get("predicted_cleavage", "")).strip(),
                passed=predicted_type == "SP",
                raw_sequence=str(row.get("sequence", "")).strip(),
            )
        )
    return predictions


def _fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            ids.append(header.split("|")[0])
    return ids


def _env_path(name: str) -> Path | None:
    value = os.environ.get(name)
    if not value:
        return None
    return Path(value)
