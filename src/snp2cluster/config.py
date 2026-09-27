from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml


@dataclass
class SNP2ClusterConfig:
    raw: dict[str, Any]
    path: Path

    @property
    def project_name(self) -> str:
        return self.raw.get("project", {}).get("name", "snp2cluster_project")

    @property
    def run_id(self) -> str | None:
        return self.raw.get("project", {}).get("run_id")

    @property
    def inputs(self) -> dict[str, str]:
        return self.raw.get("inputs", {})

    @property
    def columns(self) -> dict[str, str]:
        return self.raw.get("columns", {})
    

    @property
    def variables(self) -> dict[str, str]:
        return self.raw.get("variables", {})


    @property
    def analysis(self) -> dict[str, Any]:
        return self.raw.get("analysis", {})

    @property
    def outputs(self) -> dict[str, Any]:
        return self.raw.get("outputs", {})


def load_config(config_path: str | Path) -> SNP2ClusterConfig:
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}

    return SNP2ClusterConfig(raw=raw, path=path)
