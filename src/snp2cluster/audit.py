from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import json
import platform
import subprocess
import sys
from typing import Any

from . import __version__
from .config import SNP2ClusterConfig


def make_run_id(project_name: str) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    safe_project = project_name.lower().replace(" ", "_")
    return f"{safe_project}_{stamp}"


def create_run_dirs(config: SNP2ClusterConfig) -> dict[str, Path]:
    run_id = config.run_id or make_run_id(config.project_name)
    base = Path(config.outputs.get("out_dir", "runs")) / run_id
    dirs = {
        "base": base,
        "audit": base / "audit",
        "tables": base / "tables",
        "figures": base / "figures",
        "report": base / "report",
    }
    for p in dirs.values():
        p.mkdir(parents=True, exist_ok=True)
    return dirs


def dependency_versions() -> dict[str, str]:
    packages = ["pandas", "numpy", "sklearn", "scipy", "networkx", "matplotlib", "seaborn"]
    versions = {}
    for package in packages:
        try:
            mod = __import__(package)
            versions[package] = getattr(mod, "__version__", "unknown")
        except Exception:
            versions[package] = "not_importable"
    return versions


def write_manifest(config: SNP2ClusterConfig, dirs: dict[str, Path], validation_messages: list[str]) -> Path:
    manifest: dict[str, Any] = {
        "tool": "SNP2Cluster-Py",
        "version": __version__,
        "project_name": config.project_name,
        "run_directory": str(dirs["base"]),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "input_files": config.inputs,
        "columns": config.columns,
        "analysis_parameters": config.analysis,
        "validation_messages": validation_messages,
        "dependency_versions": dependency_versions(),
    }
    path = dirs["audit"] / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path


def write_validation_log(dirs: dict[str, Path], messages: list[str]) -> Path:
    path = dirs["audit"] / "validation_log.txt"
    if messages:
        content = "\n".join(messages)
    else:
        content = "Validation passed with no issues."
    path.write_text(content + "\n", encoding="utf-8")
    return path
