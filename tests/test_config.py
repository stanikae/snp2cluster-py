from pathlib import Path
from snp2cluster.config import load_config


def test_load_example_config():
    cfg = load_config(Path("configs/example_config.yaml"))
    assert cfg.project_name == "example_hospital"
    assert "metadata" in cfg.inputs
