from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd
import pytest

from snp2cluster.visualization.mst import (
    _legend_shape_svg,
    _shape_map,
    save_interactive_core_mst,
)


def _symbol_key(
    symbol: dict[str, object],
) -> tuple[str, object]:
    border_dashes = symbol[
        "border_dashes"
    ]

    if isinstance(border_dashes, list):
        canonical_border = tuple(
            border_dashes
        )
    else:
        canonical_border = border_dashes

    return (
        str(symbol["shape"]),
        canonical_border,
    )


def test_shape_map_assigns_17_unique_symbols() -> None:
    series = pd.Series(
        [
            "13",
            "14",
            "17",
            "22",
            "25",
            "39",
            "48",
            "152",
            "198",
            "405",
            "464",
            "643",
            "870",
            "985",
            "997",
            "1079",
            "1552",
        ]
    )

    symbols = _shape_map(series)

    assert len(symbols) == 17

    unique_symbols = {
        _symbol_key(symbol)
        for symbol in symbols.values()
    }

    assert len(unique_symbols) == 17


def test_shape_map_is_deterministic() -> None:
    first = _shape_map(
        pd.Series(["152", "13", "39"])
    )

    second = _shape_map(
        pd.Series(["39", "152", "13"])
    )

    assert first == second


def test_shape_map_rejects_more_than_33_categories() -> None:
    series = pd.Series(
        [
            str(index)
            for index in range(34)
        ]
    )

    with pytest.raises(
        ValueError,
        match=(
            "maximum supported number is 33"
        ),
    ):
        _shape_map(series)


def test_legend_svg_supports_all_mst_shapes() -> None:
    shapes = (
        "dot",
        "square",
        "triangle",
        "diamond",
        "star",
        "triangleDown",
        "hexagon",
        "ellipse",
        "circle",
        "box",
        "database",
    )

    rendered = {
        shape: _legend_shape_svg(
            shape,
            [6, 4],
        )
        for shape in shapes
    }

    assert len(set(rendered.values())) == len(shapes)

    for svg in rendered.values():
        assert 'stroke-dasharray="6 4"' in svg


def test_mst_html_contains_shape_border_patterns(
    tmp_path: Path,
) -> None:
    tree = nx.Graph()

    tree.add_edge(
        "A",
        "B",
        weight=0,
    )

    annotations = pd.DataFrame(
        {
            "sample_id": ["A", "B"],
            "core_snp_cluster": [
                "13__1",
                "13__1",
            ],
            "_core_context": [
                "13",
                "14",
            ],
            "st_type": [
                "13",
                "14",
            ],
        }
    )

    output_path = (
        tmp_path
        / "mst.html"
    )

    save_interactive_core_mst(
        tree,
        annotations,
        output_path,
        sample_id_column="sample_id",
        shape_column="st_type",
        split_value="Hospital A",
    )

    document = output_path.read_text(
        encoding="utf-8",
    )

    normalised_document = " ".join(
    document.split()
    )

    assert "shapeProperties" in document
    assert "borderDashes" in document

    assert (
    "Node fill colour represents "
    "Core SNP cluster membership."
    in normalised_document
    )
