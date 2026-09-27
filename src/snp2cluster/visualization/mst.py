from __future__ import annotations

import html
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from matplotlib.colors import to_hex
import networkx as nx
import pandas as pd
import seaborn as sns


_MISSING_LABEL = "Missing"
_MISSING_COLOUR = "#D3D3D3"

# _VIS_SHAPES = (
#     "dot",
#     "square",
#     "triangle",
#     "diamond",
#     "star",
#     "triangleDown",
#     "hexagon",
#     "ellipse",
# )

# TODO: Extend categorical shape support using a secondary style layer,
# such as solid and dashed borders. Combining 12 immediately usable
# VisJS shapes with two border styles would provide 24 distinct category
# symbols without changing the underlying category-mapping logic.

_VIS_SHAPES = (
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
    "text",
)




def _get(
    config: Any,
    *keys: str,
    default: Any = None,
) -> Any:
    """Read nested configuration values."""

    value = config

    for key in keys:

        if isinstance(value, Mapping):

            if key not in value:
                return default

            value = value[key]

        else:

            if not hasattr(value, key):
                return default

            value = getattr(
                value,
                key,
            )

    return value


def _natural_key(
    value: object,
) -> tuple[object, ...]:
    """Return a natural-sort key."""

    return tuple(
        int(part)
        if part.isdigit()
        else part.casefold()
        for part in re.split(
            r"(\d+)",
            str(value),
        )
    )


def _normalise_category(
    value: object,
) -> str:
    """Convert null-like categories to a display label."""

    if pd.isna(value):
        return _MISSING_LABEL

    label = str(value).strip()

    if label.casefold() in {
        "",
        "na",
        "nan",
        "none",
        "null",
    }:
        return _MISSING_LABEL

    return label


def _safe_prefix(
    value: object,
) -> str:
    """Create a filesystem-safe output label."""

    cleaned = str(value).strip()

    for character in (
        " ",
        "/",
        "\\",
        ":",
    ):

        cleaned = cleaned.replace(
            character,
            "_",
        )

    return cleaned or "UNKNOWN"


def _category_colour_map(
    series: pd.Series,
) -> dict[str, str]:
    """Assign deterministic colours to categories."""

    labels = sorted(
        {
            _normalise_category(value)
            for value in series
        },
        key=_natural_key,
    )

    observed = [
        label
        for label in labels
        if label != _MISSING_LABEL
    ]

    palette = sns.color_palette(
        "husl",
        n_colors=max(
            len(observed),
            1,
        ),
    )

    colours = {
        label: to_hex(
            palette[index]
        )
        for index, label
        in enumerate(observed)
    }

    if _MISSING_LABEL in labels:
        colours[
            _MISSING_LABEL
        ] = _MISSING_COLOUR

    return colours


def _shape_map(
    series: pd.Series,
) -> dict[str, str]:
    """Assign deterministic VisJS shapes to categories."""

    labels = sorted(
        {
            _normalise_category(value)
            for value in series
        },
        key=_natural_key,
    )

    return {
        label: _VIS_SHAPES[
            index % len(_VIS_SHAPES)
        ]
        for index, label
        in enumerate(labels)
    }


def _resolve_sample_id_column(
    annotation_df: pd.DataFrame,
    node_ids: set[str],
    candidates: Sequence[str | None],
) -> str:
    """Choose the available sample-ID column with greatest node coverage."""

    available = [
        candidate
        for candidate in candidates
        if (
            candidate
            and candidate
            in annotation_df.columns
        )
    ]

    if not available:
        raise ValueError(
            "No usable sample-ID column "
            "found in MST annotations"
        )

    coverage = {
        column: len(
            set(
                annotation_df[column]
                .dropna()
                .astype(str)
            )
            & node_ids
        )
        for column in available
    }

    return max(
        available,
        key=lambda column: coverage[column],
    )


def prepare_mst_annotations(
    tree: nx.Graph,
    annotation_df: pd.DataFrame,
    *,
    sample_id_column: str,
    core_cluster_column: str,
    core_context_column: str = "_core_context",
    shape_column: str | None = None,
    tooltip_columns: Sequence[str] = (),
) -> pd.DataFrame:
    """Align enriched isolate metadata with MST nodes."""

    if not isinstance(
        tree,
        nx.Graph,
    ):
        raise TypeError(
            "tree must be a NetworkX graph"
        )

    if tree.number_of_nodes() == 0:
        raise ValueError(
            "tree contains no nodes"
        )

    if not isinstance(
        annotation_df,
        pd.DataFrame,
    ):
        raise TypeError(
            "annotation_df must be a pandas DataFrame"
        )

    if annotation_df.empty:
        raise ValueError(
            "annotation_df is empty"
        )

    required = {
        sample_id_column,
        core_cluster_column,
        core_context_column,
    }

    missing = (
        required
        - set(annotation_df.columns)
    )

    if missing:
        raise ValueError(
            "Missing MST annotation columns: "
            f"{missing}"
        )

    selected = [
        sample_id_column,
        core_cluster_column,
        core_context_column,
    ]

    optional_columns = [
        shape_column,
        *tooltip_columns,
    ]

    for column in optional_columns:

        if (
            column
            and column
            in annotation_df.columns
            and column not in selected
        ):
            selected.append(
                column
            )

    annotations = annotation_df.loc[
        :,
        selected,
    ].copy()

    annotations[
        sample_id_column
    ] = (
        annotations[
            sample_id_column
        ]
        .astype(str)
    )

    annotations = (
        annotations.drop_duplicates(
            subset=[
                sample_id_column,
            ],
            keep="first",
        )
    )

    annotations = (
        annotations.set_index(
            sample_id_column
        )
    )

    node_ids = [
        str(node)
        for node in tree.nodes
    ]

    annotations = (
        annotations.reindex(
            node_ids
        )
    )

    annotations.index.name = (
        sample_id_column
    )

    annotations[
        "_core_cluster"
    ] = (
        annotations[
            core_cluster_column
        ]
        .map(
            _normalise_category
        )
    )

    annotations[
        "_core_context_label"
    ] = (
        annotations[
            core_context_column
        ]
        .map(
            _normalise_category
        )
    )

    if (
        shape_column
        and shape_column
        in annotations.columns
    ):

        annotations[
            "_shape_group"
        ] = (
            annotations[
                shape_column
            ]
            .map(
                _normalise_category
            )
        )

    else:

        annotations[
            "_shape_group"
        ] = "Isolate"

    return annotations



def _legend_shape_svg(
    shape: str,
) -> str:
    """Return an inline SVG representation of a VisJS node shape."""

    svg_shapes = {
        "dot": (
            '<circle cx="10" cy="10" r="7" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.5" />'
        ),
        "square": (
            '<rect x="3" y="3" width="14" height="14" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.5" />'
        ),
        "triangle": (
            '<polygon points="10,2 18,17 2,17" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.5" />'
        ),
        "triangleDown": (
            '<polygon points="2,3 18,3 10,18" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.5" />'
        ),
        "diamond": (
            '<polygon points="10,2 18,10 10,18 2,10" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.5" />'
        ),
        "star": (
            '<polygon '
            'points="10,1 12.5,7 19,7.5 14,11.5 '
            '15.5,18 10,14.5 4.5,18 6,11.5 '
            '1,7.5 7.5,7" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.2" />'
        ),
        "hexagon": (
            '<polygon points="5,2 15,2 19,10 '
            '15,18 5,18 1,10" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.5" />'
        ),
        "ellipse": (
            '<ellipse cx="10" cy="10" rx="8" ry="5.5" '
            'fill="#777777" stroke="#444444" '
            'stroke-width="1.5" />'
        ),
    }

    shape_markup = svg_shapes.get(
        shape,
        svg_shapes["dot"],
    )

    return (
        '<svg class="mst-shape-svg" '
        'viewBox="0 0 20 20" '
        'aria-hidden="true">'
        f"{shape_markup}"
        "</svg>"
    )

def _build_html_legend(
    *,
    title: str,
    split_label: str,
    split_value: str,
    core_colours: Mapping[str, str],
    shape_column: str | None,
    shape_values: Mapping[str, str],
) -> str:
    """Build the fixed title and legend panel for the HTML MST."""

    core_items = "\n".join(
        (
            '<div class="mst-legend-item">'
            f'<span class="mst-colour" '
            f'style="background:{html.escape(colour)};"></span>'
            f"<span>{html.escape(label)}</span>"
            "</div>"
        )
        for label, colour
        in core_colours.items()
    )

    if shape_column:

        shape_items = "\n".join(
            (
                '<div class="mst-legend-item">'
                '<span class="mst-shape-symbol">'
                f"{_legend_shape_svg(shape)}"
                "</span>"
                f"<span>{html.escape(label)}</span>"
                "</div>"
            )
            for label, shape
            in shape_values.items()
        )

        shape_section = f"""
        <div class="mst-legend-section">
            <div class="mst-legend-heading">
                {html.escape(shape_column)}
            </div>
            {shape_items}
        </div>
        """

    else:

        shape_section = ""

    return f"""
    <style>
        #mst-info-panel {{
            position: fixed;
            top: 18px;
            left: 18px;
            z-index: 1000;
            width: 280px;
            max-height: calc(100vh - 36px);
            overflow-y: auto;
            background: rgba(255, 255, 255, 0.96);
            border: 1px solid #d0d0d0;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.12);
            padding: 14px 16px;
            font-family: Arial, sans-serif;
            color: #222222;
        }}

        #mst-info-panel h1 {{
            margin: 0 0 5px 0;
            font-size: 18px;
            line-height: 1.25;
        }}

        #mst-info-panel .mst-subtitle {{
            margin-bottom: 14px;
            color: #555555;
            font-size: 13px;
        }}

        .mst-legend-section {{
            margin-top: 12px;
        }}

        .mst-legend-heading {{
            margin-bottom: 7px;
            font-weight: bold;
            font-size: 13px;
        }}

        .mst-legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
            margin: 4px 0;
            font-size: 12px;
        }}

        .mst-colour {{
            display: inline-block;
            width: 13px;
            height: 13px;
            border: 1px solid #444444;
            border-radius: 50%;
            flex: 0 0 auto;
        }}

        .mst-shape-symbol {{
            display: inline-flex;
            width: 20px;
            height: 20px;
            align-items: center;
            justify-content: center;
            flex: 0 0 auto;
        }}

        .mst-shape-svg {{
            display: block;
            width: 18px;
            height: 18px;
            overflow: visible;
        }}

        .mst-edge-note {{
            margin-top: 12px;
            padding-top: 10px;
            border-top: 1px solid #dddddd;
            color: #555555;
            font-size: 11px;
            line-height: 1.35;
        }}

        #mynetwork {{
            height: 100vh !important;
        }}
    </style>

    <div id="mst-info-panel">
        <h1>{html.escape(title)}</h1>

        <div class="mst-subtitle">
            {html.escape(split_label)}:
            <strong>{html.escape(split_value)}</strong>
        </div>

        <div class="mst-legend-section">
            <div class="mst-legend-heading">
                Core SNP cluster
            </div>
            {core_items}
        </div>

        {shape_section}

        <div class="mst-edge-note">
            Lines represent minimum-spanning-tree connections.
            Hover over an edge to view the SNP distance.
            Hover over a node to view isolate metadata.
        </div>
    </div>
    """


def _inject_html_panel(
    output_path: Path,
    panel_html: str,
) -> None:
    """Inject a custom title and legend into PyVis HTML."""

    document = output_path.read_text(
        encoding="utf-8",
    )

    body_marker = "<body>"

    if body_marker not in document:
        raise ValueError(
            "Could not locate the HTML body "
            "in generated PyVis output"
        )

    document = document.replace(
        body_marker,
        body_marker + panel_html,
        1,
    )

    output_path.write_text(
        document,
        encoding="utf-8",
    )


def save_interactive_core_mst(
    tree: nx.Graph,
    annotation_df: pd.DataFrame,
    output_path: str | Path,
    *,
    sample_id_column: str,
    core_cluster_column: str = "core_snp_cluster",
    core_context_column: str = "_core_context",
    shape_column: str | None = None,
    tooltip_columns: Sequence[str] = (),
    title: str = "Core SNP minimum spanning tree",
    split_label: str = "Analysis group",
    split_value: str,
) -> Path:
    """Save an interactive, Core-cluster-focused MST."""

    try:

        from pyvis.network import Network

    except ImportError as exc:

        raise RuntimeError(
            "Interactive MST output requires pyvis"
        ) from exc

    annotations = prepare_mst_annotations(
        tree,
        annotation_df,
        sample_id_column=sample_id_column,
        core_cluster_column=core_cluster_column,
        core_context_column=core_context_column,
        shape_column=shape_column,
        tooltip_columns=tooltip_columns,
    )

    core_colours = (
        _category_colour_map(
            annotations[
                "_core_cluster"
            ]
        )
    )

    shape_values = (
        _shape_map(
            annotations[
                "_shape_group"
            ]
        )
    )

    output_path = Path(
        output_path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    network = Network(
        height="100vh",
        width="100%",
        directed=False,
        bgcolor="#ffffff",
        font_color="#222222",
        heading="",
        cdn_resources="remote",
    )

    node_count = tree.number_of_nodes()

    show_labels = (
        node_count <= 35
    )

    for node in tree.nodes:

        node_id = str(node)

        row = annotations.loc[
            node_id
        ]

        tooltip = [
            f"<b>Sample:</b> "
            f"{html.escape(node_id)}",
            (
                "<b>Core context:</b> "
                f"{html.escape(row['_core_context_label'])}"
            ),
            (
                "<b>Core SNP cluster:</b> "
                f"{html.escape(row['_core_cluster'])}"
            ),
        ]

        for column in tooltip_columns:

            if column in annotations.columns:

                tooltip.append(
                    f"<b>{html.escape(str(column))}:</b> "
                    f"{html.escape(_normalise_category(row[column]))}"
                )

        node_colour = core_colours[
            row["_core_cluster"]
        ]

        node_shape = shape_values[
            row["_shape_group"]
        ]

        network.add_node(
            node_id,
            label=(
                node_id
                if show_labels
                else ""
            ),
            title="<br>".join(
                tooltip
            ),
            shape=node_shape,
            size=18,
            color={
                "background": node_colour,
                "border": "#444444",
                "highlight": {
                    "background": node_colour,
                    "border": "#000000",
                },
            },
            borderWidth=1.5,
        )

    for (
        node_a,
        node_b,
        edge_data,
    ) in tree.edges(
        data=True
    ):

        snp_distance = float(
            edge_data.get(
                "weight",
                0.0,
            )
        )

        visual_length = (
            65
            + min(
                max(
                    snp_distance,
                    0.0,
                ),
                100.0,
            )
            * 2.5
        )

        network.add_edge(
            str(node_a),
            str(node_b),
            title=(
                "SNP distance: "
                f"{snp_distance:g}"
            ),
            label=(
                f"{snp_distance:g}"
                if (
                    snp_distance > 0
                    and node_count <= 30
                )
                else ""
            ),
            length=visual_length,
            color="#888888",
            width=1.4,
        )

    network.set_options(
        """
        {
          "physics": {
            "enabled": true,
            "solver": "forceAtlas2Based",
            "forceAtlas2Based": {
              "gravitationalConstant": -70,
              "centralGravity": 0.008,
              "springLength": 120,
              "springConstant": 0.05,
              "damping": 0.55,
              "avoidOverlap": 0.8
            },
            "stabilization": {
              "enabled": true,
              "iterations": 700,
              "fit": true
            }
          },
          "interaction": {
            "hover": true,
            "navigationButtons": true,
            "keyboard": true,
            "multiselect": true
          },
          "nodes": {
            "font": {
              "size": 12,
              "face": "Arial"
            }
          },
          "edges": {
            "smooth": false,
            "font": {
              "size": 10,
              "align": "middle"
            }
          }
        }
        """
    )

    network.write_html(
        str(output_path),
        open_browser=False,
    )

    panel_html = _build_html_legend(
        title=title,
        split_label=split_label,
        split_value=split_value,
        core_colours=core_colours,
        shape_column=shape_column,
        shape_values=shape_values,
    )

    _inject_html_panel(
        output_path,
        panel_html,
    )

    return output_path



def generate_core_mst_from_config(
    tree: nx.Graph,
    annotation_df: pd.DataFrame,
    config: Any,
    output_path: str | Path,
    split_value: str,
    split_label: str | None = None,
    shape_column: str | None = None,
    core_cluster_column: str = "core_snp_cluster",
    core_context_column: str = "_core_context",
    title: str = "Core SNP minimum spanning tree",
) -> Path:
    """Config-aware wrapper for one interactive Core MST."""

    configured_sample_id = _get(
        config,
        "columns",
        "sample_id",
    )

    sequence_type_column = _get(
        config,
        "columns",
        "sequence_type",
    )

    main_var = _get(
        config,
        "variables",
        "main_var",
    )

    var_01 = _get(
        config,
        "variables",
        "var_01",
    )

    var_02 = _get(
        config,
        "variables",
        "var_02",
    )

    transmission_level = str(
        _get(
            config,
            "analysis",
            "transmission_level",
            default="Facility",
        )
    ).strip().casefold()

    node_ids = {
        str(node)
        for node in tree.nodes
    }

    sample_id_column = _resolve_sample_id_column(
        annotation_df,
        node_ids,
        candidates=(
            configured_sample_id,
            "sample_id",
        ),
    )

    if shape_column is not None:

        resolved_shape_column = (
            str(shape_column)
            if str(shape_column) in annotation_df.columns
            else None
        )

    elif transmission_level == "community":

        resolved_shape_column = (
            str(var_01)
            if (
                var_01
                and str(var_01) in annotation_df.columns
            )
            else None
        )

    else:

        resolved_shape_column = (
            str(sequence_type_column)
            if (
                sequence_type_column
                and str(sequence_type_column)
                in annotation_df.columns
            )
            else None
        )

    tooltip_columns = [
        column
        for column in (
            main_var,
            var_01,
            sequence_type_column,
            var_02,
        )
        if (
            column
            and column in annotation_df.columns
        )
    ]

    resolved_split_label = (
        str(split_label)
        if split_label
        else (
            str(main_var)
            if main_var
            else "Analysis group"
        )
    )

    return save_interactive_core_mst(
        tree,
        annotation_df,
        output_path,
        sample_id_column=sample_id_column,
        core_cluster_column=core_cluster_column,
        core_context_column=core_context_column,
        shape_column=resolved_shape_column,
        tooltip_columns=tooltip_columns,
        title=title,
        split_label=resolved_split_label,
        split_value=str(split_value),
    )