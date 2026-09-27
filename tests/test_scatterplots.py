import matplotlib
matplotlib.use("Agg")

from pathlib import Path

import pandas as pd
import pytest

from snp2cluster.visualization.scatterplots import (
    generate_transmission_scatter_from_config,
    plot_transmission_scatter,
)


def _df() -> pd.DataFrame:
    return pd.DataFrame({
        "SampleID": ["S1", "S2", "S3", "S4"],
        "ST": [307, 307, 152, 152],
        "TakenDate": ["2019-10-07", "2019-11-25", "2019-11-25", "2020-01-06"],
        "Clusters": [None, 4, 3, 3],
        "WardType": ["neonatal", "other", "neonatal", "neonatal"],
    })


def _config(cluster_type: str = "Transmission") -> dict:
    return {
        "columns": {"sequence_type": "ST"},
        "variables": {"main_var": "FacilityName", "var_01": "WardType", "var_02": "TakenDate"},
        "analysis": {
            "cluster_type": cluster_type,
            "transmission_level": "Facility",
            "snp_cutoff": 20,
            "days_cutoff": 45,
            "random_seed": 42,
            "date_format": "ymd",
        },
        "outputs": {"generate_figures": True},
    }


def test_config_entry_point_generates_transmission_outputs(tmp_path: Path) -> None:
    paths = generate_transmission_scatter_from_config(_df(), _config(), tmp_path)
    assert paths["png"].exists()
    assert paths["pdf"].exists()


def test_config_entry_point_skips_core_mode(tmp_path: Path) -> None:
    paths = generate_transmission_scatter_from_config(_df(), _config("Core"), tmp_path)
    assert paths == {}
    assert list(tmp_path.iterdir()) == []


def test_config_entry_point_skips_when_figures_disabled(tmp_path: Path) -> None:
    config = _config()
    config["outputs"]["generate_figures"] = False
    assert generate_transmission_scatter_from_config(_df(), config, tmp_path) == {}


def test_plotter_accepts_configured_metadata_names(tmp_path: Path) -> None:
    renamed = _df().rename(columns={"ST": "SequenceType", "WardType": "LocationClass"})
    paths = plot_transmission_scatter(
        renamed,
        tmp_path,
        sequence_type_column="SequenceType",
        collection_date_column="TakenDate",
        cluster_column="Clusters",
        grouping_column="LocationClass",
        date_format="ymd",
    )
    assert all(path.exists() for path in paths.values())


def test_missing_configured_column_is_reported(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Missing columns"):
        generate_transmission_scatter_from_config(
            _df().drop(columns="TakenDate"), _config(), tmp_path
        )
