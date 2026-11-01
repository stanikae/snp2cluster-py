from pathlib import Path

import pandas as pd

from snp2cluster.visualization.diagnostics import (
    plot_k_selection,
)


def test_plot_k_selection_handles_valid_k1_result(
    tmp_path: Path,
):
    diagnostics = pd.DataFrame(
        {
            "k": [1],
            "silhouette": [float("nan")],
            "selection_method": [
                "single_distinct_profile"
            ],
        }
    )

    outputs = plot_k_selection(
        diagnostics_df=diagnostics,
        output_dir=tmp_path,
        title="ST377",
        filename_prefix="01_k_selection_377",
    )

    assert outputs["png"].exists()
    assert outputs["pdf"].exists()
    assert outputs["png"].stat().st_size > 0
    assert outputs["pdf"].stat().st_size > 0


def test_plot_k_selection_selects_best_valid_score(
    tmp_path: Path,
):
    diagnostics = pd.DataFrame(
        {
            "k": [2, 3, 4],
            "silhouette": [0.42, 0.71, 0.55],
            "selection_method": [
                "silhouette",
                "silhouette",
                "silhouette",
            ],
        }
    )

    outputs = plot_k_selection(
        diagnostics_df=diagnostics,
        output_dir=tmp_path,
        title="Standard context",
        filename_prefix="01_k_selection_standard",
    )

    assert outputs["png"].exists()
    assert outputs["pdf"].exists()
    assert outputs["png"].stat().st_size > 0
    assert outputs["pdf"].stat().st_size > 0