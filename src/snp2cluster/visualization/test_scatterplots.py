import matplotlib
matplotlib.use("Agg")

from pathlib import Path

import pandas as pd
import pytest

from snp2cluster.visualization.scatterplots import plot_transmission_scatter


def _example_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S3", "S4", "S5"],
            "Epiweek": ["2019.41", "2019.48", "2019.48", "2020.01", "2020.01"],
            "ST": [307, 307, 152, 152, 152],
            "Clusters": [None, 4, 3, 3, 3],
            "WardType": ["neonatal", "neonatal", "other", "neonatal", "neonatal"],
        }
    )


def test_plot_transmission_scatter_writes_png_and_pdf(tmp_path: Path) -> None:
    paths = plot_transmission_scatter(
        _example_df(),
        tmp_path,
        snp_cutoff=20,
        day_window=45,
    )

    assert paths["png"].exists()
    assert paths["pdf"].exists()
    assert paths["png"].name == "02_transmission_scatter.png"
    assert paths["pdf"].name == "02_transmission_scatter.pdf"


def test_plot_transmission_scatter_supports_no_shape_column(tmp_path: Path) -> None:
    paths = plot_transmission_scatter(
        _example_df().drop(columns="WardType"),
        tmp_path,
        shape_column=None,
    )
    assert all(path.exists() for path in paths.values())


def test_plot_transmission_scatter_reports_missing_columns(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Missing columns"):
        plot_transmission_scatter(
            _example_df().drop(columns="ST"),
            tmp_path,
        )
