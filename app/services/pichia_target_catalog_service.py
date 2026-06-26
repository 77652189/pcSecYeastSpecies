from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from pcsec_pichia.core.paths import ProjectPaths


@dataclass(frozen=True)
class BuiltinTarget:
    """A predefined target configuration used in the Streamlit catalog."""

    signal_peptide_id: str
    leader_id: str
    mature_id: str
    label: str
    disulfide_sites: int = 0
    n_glycosylation_sites: int = 0
    o_glycosylation_sites: int = 0


def _load_catalog_from_targets() -> dict[str, dict[str, object]]:
    """Load known signal peptides, leaders, and mature sequences from formal targets."""
    from pcsec_pichia.targets import load_builtin_targets

    targets = {target.target_id: target for target in load_builtin_targets()}
    opn = targets["OPN_ALPHA_FULL_PROJECT"]
    hlf = targets["hLF"]
    opn_alpha_pro = opn.leader_sequence[len(opn.signal_peptide_sequence) :]
    return {
        "signal_peptides": {
            "alpha-factor_MRFPS": {
                "id": "alpha-factor_MRFPS",
                "label": "酵母 alpha-factor 信号肽 (MRFPSIFTAVLFAASSALA)",
                "sequence": opn.signal_peptide_sequence,
                "length": len(opn.signal_peptide_sequence),
                "source": "S. cerevisiae alpha-factor",
            },
            "native_hLF": {
                "id": "native_hLF",
                "label": "人源 hLF 天然信号肽 (MKLVFLVLLFLGALGLCLA)",
                "sequence": "MKLVFLVLLFLGALGLCLA",
                "length": 19,
                "source": "Homo sapiens lactoferrin",
            },
            "none": {
                "id": "none",
                "label": "无信号肽",
                "sequence": "",
                "length": 0,
                "source": "",
            },
        },
        "leaders": {
            "OPN_alpha-pro": {
                "id": "OPN_alpha-pro",
                "label": f"OPN alpha-factor pro 引导肽 ({len(opn_alpha_pro)}aa)",
                "sequence": opn_alpha_pro,
                "length": len(opn_alpha_pro),
                "source": "OPN prepro leader",
            },
            "none": {
                "id": "none",
                "label": "无引导肽",
                "sequence": "",
                "length": 0,
                "source": "",
            },
        },
        "mature_proteins": {
            "OPN_ALPHA_FULL_PROJECT": {
                "id": "OPN_ALPHA_FULL_PROJECT",
                "label": f"OPN (骨桥蛋白, {len(opn.mature_sequence)}aa)",
                "sequence": opn.mature_sequence,
                "length": len(opn.mature_sequence),
                "disulfide_sites": opn.disulfide_sites,
                "n_glycosylation_sites": opn.n_glycosylation_sites,
                "o_glycosylation_sites": opn.o_glycosylation_sites,
                "source": "Data/pcSecPichia/TargetProtein_OPN_candidates.csv",
            },
            "hLF": {
                "id": "hLF",
                "label": f"hLF (人乳铁蛋白, {len(hlf.mature_sequence)}aa)",
                "sequence": hlf.mature_sequence,
                "length": len(hlf.mature_sequence),
                "disulfide_sites": hlf.disulfide_sites,
                "n_glycosylation_sites": hlf.n_glycosylation_sites,
                "o_glycosylation_sites": hlf.o_glycosylation_sites,
                "source": hlf.source,
            },
        },
        "builtin_targets": {
            "OPN_ALPHA_FULL_PROJECT": BuiltinTarget(
                signal_peptide_id="alpha-factor_MRFPS",
                leader_id="OPN_alpha-pro",
                mature_id="OPN_ALPHA_FULL_PROJECT",
                label="OPN（酵母 alpha-factor 信号肽 + pro 引导肽）",
                disulfide_sites=opn.disulfide_sites,
                n_glycosylation_sites=opn.n_glycosylation_sites,
                o_glycosylation_sites=opn.o_glycosylation_sites,
            ),
            "hLF": BuiltinTarget(
                signal_peptide_id="native_hLF",
                leader_id="none",
                mature_id="hLF",
                label="hLF（人源天然信号肽 + hLF 成熟序列）",
                disulfide_sites=hlf.disulfide_sites,
                n_glycosylation_sites=hlf.n_glycosylation_sites,
                o_glycosylation_sites=hlf.o_glycosylation_sites,
            ),
        },
    }


_CATALOG: dict[str, dict[str, object]] | None = None


def _get_catalog() -> dict[str, dict[str, object]]:
    global _CATALOG
    if _CATALOG is None:
        _CATALOG = _load_catalog_from_targets()
    return _CATALOG


def known_signal_peptides() -> dict[str, dict[str, object]]:
    """Return {id: info} for all known signal peptides."""
    return dict(_get_catalog()["signal_peptides"])


def known_leaders() -> dict[str, dict[str, object]]:
    """Return {id: info} for all known leaders."""
    return dict(_get_catalog()["leaders"])


def known_mature_proteins() -> dict[str, dict[str, object]]:
    """Return {id: info} for all known mature proteins."""
    return dict(_get_catalog()["mature_proteins"])


@dataclass
class TargetBuildConfig:
    """A user's saved target configuration."""

    name: str
    signal_peptide_id: str
    signal_peptide_sequence: str
    leader_id: str
    leader_sequence: str
    mature_id: str
    mature_sequence: str
    disulfide_sites: int = 0
    n_glycosylation_sites: int = 0
    o_glycosylation_sites: int = 0
    enable_ribosome: bool = True
    enable_misfolding: bool = True


def _safe_id(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value).strip("._")
    return cleaned or "target"


def _config_dir(paths: ProjectPaths) -> Path:
    directory = paths.local_runs_dir / "streamlit_pichia_runs" / ".user_configs"
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def save_target_config(config: TargetBuildConfig, paths: ProjectPaths) -> None:
    """Save a user target configuration to a local JSON file."""
    safe = _safe_id(config.name) or f"config_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    path = _config_dir(paths) / f"{safe}.json"
    data = {
        "name": config.name,
        "signal_peptide_id": config.signal_peptide_id,
        "signal_peptide_sequence": config.signal_peptide_sequence,
        "leader_id": config.leader_id,
        "leader_sequence": config.leader_sequence,
        "mature_id": config.mature_id,
        "mature_sequence": config.mature_sequence,
        "disulfide_sites": config.disulfide_sites,
        "n_glycosylation_sites": config.n_glycosylation_sites,
        "o_glycosylation_sites": config.o_glycosylation_sites,
        "enable_ribosome": config.enable_ribosome,
        "enable_misfolding": config.enable_misfolding,
        "saved_at": datetime.now().isoformat(),
    }
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def load_target_config(name: str, paths: ProjectPaths) -> TargetBuildConfig | None:
    """Load a saved user target configuration by name."""
    safe = _safe_id(name)
    path = _config_dir(paths) / f"{safe}.json"
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return TargetBuildConfig(
            name=data.get("name", name),
            signal_peptide_id=data.get("signal_peptide_id", "custom"),
            signal_peptide_sequence=data.get("signal_peptide_sequence", ""),
            leader_id=data.get("leader_id", "custom"),
            leader_sequence=data.get("leader_sequence", ""),
            mature_id=data.get("mature_id", "custom"),
            mature_sequence=data.get("mature_sequence", ""),
            disulfide_sites=int(data.get("disulfide_sites", 0)),
            n_glycosylation_sites=int(data.get("n_glycosylation_sites", 0)),
            o_glycosylation_sites=int(data.get("o_glycosylation_sites", 0)),
            enable_ribosome=bool(data.get("enable_ribosome", True)),
            enable_misfolding=bool(data.get("enable_misfolding", True)),
        )
    except (json.JSONDecodeError, KeyError, TypeError, ValueError):
        return None


def list_saved_configs(paths: ProjectPaths) -> list[dict[str, object]]:
    """List all saved user configurations."""
    configs: list[dict[str, object]] = []
    for path in sorted(_config_dir(paths).glob("*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            configs.append(
                {
                    "name": data.get("name", path.stem),
                    "file": path.name,
                    "signal_peptide_id": data.get("signal_peptide_id", ""),
                    "leader_id": data.get("leader_id", ""),
                    "mature_id": data.get("mature_id", ""),
                    "saved_at": data.get("saved_at", ""),
                }
            )
        except (json.JSONDecodeError, OSError):
            pass
    return configs


def delete_target_config(name: str, paths: ProjectPaths) -> bool:
    """Delete a saved user configuration by name. Returns True if deleted."""
    safe = _safe_id(name)
    path = _config_dir(paths) / f"{safe}.json"
    if path.exists():
        path.unlink()
        return True
    return False


__all__ = [
    "BuiltinTarget",
    "TargetBuildConfig",
    "delete_target_config",
    "known_leaders",
    "known_mature_proteins",
    "known_signal_peptides",
    "list_saved_configs",
    "load_target_config",
    "save_target_config",
]
