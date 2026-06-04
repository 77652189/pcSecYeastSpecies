from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px
from scipy.io import loadmat

from app.core.i18n import category_label, species_label
from app.core.models import DatasetInfo, LoadedDataset, SpeciesCode
from app.core.paths import ProjectPaths


SPECIES_HINTS: dict[str, SpeciesCode] = {
    "SCE": "SCE",
    "Sce": "SCE",
    "Yeast": "SCE",
    "PP": "PPA",
    "Ppa": "PPA",
    "Pichia": "PPA",
    "KMX": "KMX",
    "Kma": "KMX",
    "KM": "KMX",
    "Kmarx": "KMX",
}


def infer_species(path: Path) -> SpeciesCode:
    name = path.stem
    for hint, species in SPECIES_HINTS.items():
        if hint in name:
            return species
    return "Unknown"


@dataclass
class ResultCatalog:
    paths: ProjectPaths

    def list_datasets(self) -> list[DatasetInfo]:
        datasets: list[DatasetInfo] = []
        for path in sorted(self.paths.results_dir.rglob("*")):
            if path.suffix.lower() not in {".xlsx", ".mat"} or not path.is_file():
                continue
            rel = path.relative_to(self.paths.repo_root)
            datasets.append(
                DatasetInfo(
                    id=rel.as_posix(),
                    name=path.stem,
                    path=path,
                    category=path.parent.name,
                    suffix=path.suffix.lower(),
                    species=infer_species(path),
                    size_bytes=path.stat().st_size,
                    modified_at=pd.Timestamp(path.stat().st_mtime, unit="s").strftime("%Y-%m-%d %H:%M"),
                )
            )
        return datasets

    def get_dataset(self, dataset_id: str) -> DatasetInfo:
        for dataset in self.list_datasets():
            if dataset.id == dataset_id:
                return dataset
        raise KeyError(f"未知数据集：{dataset_id}")


class ResultLoader:
    def load_dataset(self, dataset: DatasetInfo) -> LoadedDataset:
        if dataset.suffix == ".xlsx":
            return self._load_xlsx(dataset)
        if dataset.suffix == ".mat":
            return self._load_mat(dataset)
        raise ValueError(f"暂不支持这种结果文件类型：{dataset.suffix}")

    def _load_xlsx(self, dataset: DatasetInfo) -> LoadedDataset:
        workbook = pd.ExcelFile(dataset.path)
        tables = {}
        for sheet in workbook.sheet_names:
            frame = pd.read_excel(workbook, sheet_name=sheet)
            tables[sheet] = self._records(frame)
        return LoadedDataset(info=dataset, kind="table", tables=tables)

    def _load_mat(self, dataset: DatasetInfo) -> LoadedDataset:
        data = loadmat(dataset.path, squeeze_me=True, struct_as_record=False)
        summary: list[dict[str, Any]] = []
        tables: dict[str, list[dict[str, Any]]] = {}
        for key, value in data.items():
            if key.startswith("__"):
                continue
            shape = getattr(value, "shape", None)
            dtype = getattr(value, "dtype", None)
            summary.append(
                {
                    "变量": key,
                    "类型": type(value).__name__,
                    "形状": str(shape) if shape is not None else "结构体/标量",
                    "数据类型": str(dtype) if dtype is not None else "",
                }
            )
            frame = self._array_preview(value)
            if frame is not None:
                tables[key] = self._records(frame)
        return LoadedDataset(info=dataset, kind="mat", tables=tables, variable_summary=summary)

    def _array_preview(self, value: Any) -> pd.DataFrame | None:
        if not hasattr(value, "ndim"):
            return None
        if value.ndim == 1:
            return pd.DataFrame({"序号": range(len(value)), "数值": value}).head(500)
        if value.ndim == 2:
            frame = pd.DataFrame(value)
            frame.index.name = "行号"
            return frame.reset_index().head(500)
        return None

    def _records(self, frame: pd.DataFrame) -> list[dict[str, Any]]:
        clean = frame.copy()
        clean.columns = [str(column) for column in clean.columns]
        return clean.where(pd.notna(clean), None).to_dict(orient="records")


class PlotBuilder:
    def build_chart(self, loaded: LoadedDataset, table_name: str | None = None):
        table_name = table_name or next(iter(loaded.tables), None)
        if not table_name:
            return None
        frame = pd.DataFrame(loaded.tables[table_name])
        if frame.empty:
            return None
        numeric = frame.select_dtypes(include="number")
        title_prefix = f"{category_label(loaded.info.category)} - {species_label(loaded.info.species)}"
        if loaded.info.category == "CSource" and {"mu-in-vivo (h-1)", "mu-in-silico(h-1)"}.issubset(frame.columns):
            figure = px.scatter(
                frame,
                x="mu-in-vivo (h-1)",
                y="mu-in-silico(h-1)",
                color="Carbon-source" if "Carbon-source" in frame.columns else None,
                title=f"{title_prefix}：实验生长速率与模型预测对比",
                labels={
                    "mu-in-vivo (h-1)": "实验测得生长速率 mu（h^-1）",
                    "mu-in-silico(h-1)": "模型预测生长速率 mu（h^-1）",
                    "Carbon-source": "碳源",
                },
            )
            figure.update_layout(legend_title_text="碳源")
            return figure
        if numeric.shape[1] >= 2:
            return px.line(
                numeric.reset_index(),
                x="index",
                y=list(numeric.columns[: min(5, numeric.shape[1])]),
                title=f"{title_prefix}：数值变量趋势预览",
                labels={"index": "行号", "value": "数值", "variable": "变量"},
            )
        if numeric.shape[1] == 1:
            column = numeric.columns[0]
            return px.histogram(frame, x=column, title=f"{title_prefix}：{column} 分布", labels={column: str(column)})
        return None
