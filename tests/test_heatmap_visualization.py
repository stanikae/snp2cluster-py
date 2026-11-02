from __future__ import annotations

from pathlib import Path

import pandas as pd

from snp2cluster.visualization.heatmaps import (
    _legend_grid,
    _plan_legend_layout,
    plot_annotated_snp_heatmap,
)


def test_legend_grid_wraps_high_cardinality_values() -> None:
    assert _legend_grid(
        6,
        max_rows=14,
    ) == (1, 6)

    assert _legend_grid(
        17,
        max_rows=14,
    ) == (2, 9)

    assert _legend_grid(
        43,
        max_rows=14,
    ) == (4, 11)


def test_legend_layout_uses_multiple_side_panels() -> None:
    placements = _plan_legend_layout(
        [6, 2, 43, 28, 5],
        max_rows=14,
    )

    assert len(placements) == 5

    # ST and Core SNP cluster legends must use multiple columns.
    assert int(
        placements[2]["ncol"]
    ) > 1

    assert int(
        placements[3]["ncol"]
    ) > 1

    # At least one legend must move to a new horizontal panel.
    x_positions = {
        float(
            placement["x"]
        )
        for placement in placements
    }

    assert len(x_positions) > 1


def test_legend_layout_does_not_overlap_within_panel() -> None:
    placements = _plan_legend_layout(
        [6, 2, 17, 24, 4],
        max_rows=14,
    )

    grouped: dict[
        float,
        list[dict[str, float | int]],
    ] = {}

    for placement in placements:
        grouped.setdefault(
            float(placement["x"]),
            [],
        ).append(placement)

    for panel in grouped.values():
        for current, following in zip(
            panel,
            panel[1:],
        ):
            current_bottom = (
                float(current["y"])
                - float(current["height"])
            )

            assert (
                float(following["y"])
                < current_bottom
            )


def test_high_cardinality_heatmap_renders(
    tmp_path: Path,
) -> None:
    sample_ids = [
        f"S{index:02d}"
        for index in range(18)
    ]

    matrix = pd.DataFrame(
        [
            [
                abs(row - column)
                for column in range(18)
            ]
            for row in range(18)
        ],
        index=sample_ids,
        columns=sample_ids,
        dtype=float,
    )

    annotations = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "hospital": [
                f"Hospital {index % 3}"
                for index in range(18)
            ],
            "st_type": [
                str(index)
                for index in range(18)
            ],
            "core_snp_cluster": [
                f"cluster_{index}"
                for index in range(18)
            ],
        }
    )

    outputs = plot_annotated_snp_heatmap(
        matrix,
        annotations,
        tmp_path,
        sample_id_column="sample_id",
        annotation_columns=[
            "hospital",
            "st_type",
            "core_snp_cluster",
        ],
        annotation_labels={
            "hospital": "Hospital",
            "st_type": "ST",
            "core_snp_cluster": (
                "Core SNP cluster"
            ),
        },
        filename_prefix=(
            "high_cardinality_heatmap"
        ),
    )

    assert outputs["png"].exists()
    assert outputs["pdf"].exists()
    assert outputs["png"].stat().st_size > 0
    assert outputs["pdf"].stat().st_size > 0
