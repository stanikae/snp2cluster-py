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

from ..mst import build_minimum_spanning_tree, mst_edges_table


_MISSING_LABEL = "Missing"
_MISSING_COLOUR = "#D3D3D3"

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
    """Read nested configuration values from mappings or attribute objects."""
    value = config

    for key in keys:
        if isinstance(value, Mapping):
            if key not in value:
                return default
            value = value[key]
        else:
            if not hasattr(value, key):
                return default
            value = getattr(value, key)

    return value


def _natural_key(value: object) -> tuple[object, ...]:
    """Return a natural-sort key for numeric and mixed labels."""
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", str(value))
    )


def _normalise_category(value: object) -> str:
    """Convert null-like values to a stable display label."""
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


def _is_clustered(value: object) -> bool:
    """Return whether a transmission-cluster value is assigned."""
    return _normalise_category(value) != _MISSING_LABEL


def _safe_name(value: object) -> str:
    """Return a filesystem-safe label."""
    cleaned = str(value).strip()

    for character in (
        " ",
        "/",
        "\\",
        ":",
    ):
        cleaned = cleaned.replace(character, "_")

    return cleaned or "UNKNOWN"


def _category_colour_map(series: pd.Series) -> dict[str, str]:
    """Assign deterministic colours to categories."""
    labels = sorted(
        {_normalise_category(value) for value in series},
        key=_natural_key,
    )
    observed = [label for label in labels if label != _MISSING_LABEL]
    palette = sns.color_palette("husl", n_colors=max(len(observed), 1))
    colours = {
        label: to_hex(palette[index])
        for index, label in enumerate(observed)
    }

    if _MISSING_LABEL in labels:
        colours[_MISSING_LABEL] = _MISSING_COLOUR

    return colours


def _shape_map(series: pd.Series) -> dict[str, str]:
    """Assign deterministic VisJS shapes to categories."""
    labels = sorted(
        {_normalise_category(value) for value in series},
        key=_natural_key,
    )

    return {
        label: _VIS_SHAPES[index % len(_VIS_SHAPES)]
        for index, label in enumerate(labels)
    }


def _resolve_sample_id_column(
    transmission_df: pd.DataFrame,
    matrix_ids: set[str],
    candidates: Sequence[str | None],
) -> str:
    """Choose the sample-ID field with greatest SNP-matrix coverage."""
    available = [
        candidate
        for candidate in candidates
        if candidate and candidate in transmission_df.columns
    ]

    if not available:
        raise ValueError(
            "No usable sample-ID column found in transmission-network data"
        )

    coverage = {
        column: len(
            set(transmission_df[column].dropna().astype(str)) & matrix_ids
        )
        for column in available
    }

    return max(available, key=lambda column: coverage[column])


def _parse_collection_dates(
    series: pd.Series,
    *,
    date_format: str | None,
) -> pd.Series:
    """Parse collection dates using the configured date ordering."""
    normalised_format = str(date_format or "ymd").strip().casefold()
    format_map = {
        "ymd": "%Y-%m-%d",
        "ydm": "%Y-%d-%m",
        "dmy": "%d-%m-%Y",
        "mdy": "%m-%d-%Y",
    }

    parsed = pd.to_datetime(
        series,
        errors="coerce",
        format=format_map.get(normalised_format),
    )

    unresolved = parsed.isna() & series.notna()

    if unresolved.any():
        parsed.loc[unresolved] = pd.to_datetime(
            series.loc[unresolved],
            errors="coerce",
            dayfirst=normalised_format in {"dmy", "ydm"},
        )

    return parsed


def prepare_transmission_network_data(
    transmission_df: pd.DataFrame,
    snp_matrix: pd.DataFrame,
    *,
    sample_id_column: str,
    main_var: str,
    collection_date_column: str,
    cluster_column: str = "Clusters",
    shape_column: str | None = None,
    date_format: str | None = "ymd",
) -> pd.DataFrame:
    """Validate and prepare clustered isolates for network generation."""
    if not isinstance(transmission_df, pd.DataFrame):
        raise TypeError("transmission_df must be a pandas DataFrame")
    if transmission_df.empty:
        raise ValueError("transmission_df is empty")
    if not isinstance(snp_matrix, pd.DataFrame):
        raise TypeError("snp_matrix must be a pandas DataFrame")
    if snp_matrix.empty:
        raise ValueError("snp_matrix is empty")

    required = {
        sample_id_column,
        main_var,
        collection_date_column,
        cluster_column,
    }
    missing = required - set(transmission_df.columns)

    if missing:
        raise ValueError(
            "Missing transmission-network columns: "
            f"{sorted(missing)}"
        )

    matrix_ids = set(snp_matrix.index.astype(str)) & set(
        snp_matrix.columns.astype(str)
    )

    selected_columns = list(transmission_df.columns)
    data = transmission_df.loc[:, selected_columns].copy()
    data[sample_id_column] = data[sample_id_column].astype(str)
    data["_transmission_cluster"] = data[cluster_column].map(
        _normalise_category
    )
    data["_network_group"] = data[main_var].map(_normalise_category)
    data["_shape_group"] = (
        data[shape_column].map(_normalise_category)
        if shape_column and shape_column in data.columns
        else "Isolate"
    )
    data["_collection_date"] = _parse_collection_dates(
        data[collection_date_column],
        date_format=date_format,
    )

    duplicate_assignments = (
        data.groupby(sample_id_column, dropna=False)["_transmission_cluster"]
        .nunique(dropna=False)
    )
    conflicting_ids = duplicate_assignments[duplicate_assignments > 1].index.tolist()

    if conflicting_ids:
        raise ValueError(
            "Duplicate sample IDs have conflicting transmission-cluster "
            f"assignments: {conflicting_ids[:10]}"
        )

    data = data.drop_duplicates(subset=[sample_id_column], keep="first")
    data = data[data["_transmission_cluster"] != _MISSING_LABEL].copy()
    data = data[data[sample_id_column].isin(matrix_ids)].copy()

    if data.empty:
        return data

    data = data.set_index(sample_id_column, drop=False)
    return data


def build_transmission_cluster_network(
    group_df: pd.DataFrame,
    snp_matrix: pd.DataFrame,
    *,
    sample_id_column: str,
) -> tuple[nx.Graph, pd.DataFrame]:
    """Build one sparse network by combining cluster-specific SNP MSTs."""
    graph = nx.Graph()
    edge_tables: list[pd.DataFrame] = []

    for cluster_label, cluster_df in group_df.groupby(
        "_transmission_cluster",
        dropna=False,
        sort=False,
    ):
        sample_ids = [
            sample_id
            for sample_id in cluster_df[sample_id_column].astype(str).tolist()
            if sample_id in snp_matrix.index and sample_id in snp_matrix.columns
        ]

        graph.add_nodes_from(
            (sample_id, {"transmission_cluster": cluster_label})
            for sample_id in sample_ids
        )

        if len(sample_ids) < 2:
            continue

        cluster_matrix = snp_matrix.loc[sample_ids, sample_ids].copy()
        cluster_tree = build_minimum_spanning_tree(cluster_matrix)

        for node_a, node_b, edge_data in cluster_tree.edges(data=True):
            date_a = group_df.loc[str(node_a), "_collection_date"]
            date_b = group_df.loc[str(node_b), "_collection_date"]
            day_difference = (
                abs((date_a - date_b).days)
                if pd.notna(date_a) and pd.notna(date_b)
                else pd.NA
            )
            graph.add_edge(
                str(node_a),
                str(node_b),
                transmission_cluster=cluster_label,
                snp_distance=float(edge_data.get("weight", 0.0)),
                day_difference=day_difference,
            )

        cluster_edges = mst_edges_table(cluster_tree).rename(
            columns={"snp_distance": "snp_distance"}
        )
        cluster_edges["Clusters"] = cluster_label
        cluster_edges["day_difference"] = cluster_edges.apply(
            lambda row: _edge_day_difference(
                group_df,
                str(row["sample_id_1"]),
                str(row["sample_id_2"]),
            ),
            axis=1,
        )
        edge_tables.append(cluster_edges)

    if edge_tables:
        edges = pd.concat(edge_tables, ignore_index=True)
    else:
        edges = pd.DataFrame(
            columns=[
                "sample_id_1",
                "sample_id_2",
                "snp_distance",
                "Clusters",
                "day_difference",
            ]
        )

    return graph, edges


def _edge_day_difference(
    group_df: pd.DataFrame,
    sample_a: str,
    sample_b: str,
) -> object:
    """Return the absolute difference between two collection dates."""
    date_a = group_df.loc[sample_a, "_collection_date"]
    date_b = group_df.loc[sample_b, "_collection_date"]

    if pd.isna(date_a) or pd.isna(date_b):
        return pd.NA

    return abs((date_a - date_b).days)


def _legend_shape_svg(shape: str) -> str:
    """Return an inline SVG representation of a VisJS node shape."""
    svg_shapes = {
        "dot": '<circle cx="10" cy="10" r="7" />',
        "circle": '<circle cx="10" cy="10" r="7" />',
        "square": '<rect x="3" y="3" width="14" height="14" />',
        "box": '<rect x="2" y="4" width="16" height="12" rx="2" />',
        "database": (
            '<ellipse cx="10" cy="4" rx="7" ry="3" />'
            '<path d="M3 4v11c0 2 14 2 14 0V4" />'
            '<path d="M3 10c0 2 14 2 14 0" />'
        ),
        "triangle": '<polygon points="10,2 18,17 2,17" />',
        "triangleDown": '<polygon points="2,3 18,3 10,18" />',
        "diamond": '<polygon points="10,2 18,10 10,18 2,10" />',
        "star": (
            '<polygon points="10,1 12.5,7 19,7.5 14,11.5 '
            '15.5,18 10,14.5 4.5,18 6,11.5 1,7.5 7.5,7" />'
        ),
        "hexagon": '<polygon points="5,2 15,2 19,10 15,18 5,18 1,10" />',
        "ellipse": '<ellipse cx="10" cy="10" rx="8" ry="5.5" />',
        "text": '<text x="10" y="14" text-anchor="middle">T</text>',
    }
    shape_markup = svg_shapes.get(shape, svg_shapes["dot"])

    return (
        '<svg class="network-shape-svg" viewBox="0 0 20 20" '
        'aria-hidden="true">'
        '<g fill="#777777" stroke="#444444" stroke-width="1.4">'
        f"{shape_markup}"
        "</g></svg>"
    )

def _build_html_panel(
    *,
    title: str,
    group_label: str,
    group_value: str,
    cluster_colours: Mapping[str, str],
    shape_column: str | None,
    shape_values: Mapping[str, str],
    node_count: int,
    days_cutoff: float,
) -> str:
    """Build a fixed title, summary, and legend panel."""
    cluster_items = "\n".join(
        (
            '<div class="network-legend-item">'
            f'<span class="network-colour" style="background:{html.escape(colour)};">'
            "</span>"
            f"<span>{html.escape(label)}</span>"
            "</div>"
        )
        for label, colour in cluster_colours.items()
    )

    shape_section = ""
    if shape_column:
        shape_items = "\n".join(
            (
                '<div class="network-legend-item">'
                '<span class="network-shape-symbol">'
                f"{_legend_shape_svg(shape)}"
                "</span>"
                f"<span>{html.escape(label)}</span>"
                "</div>"
            )
            for label, shape in shape_values.items()
        )
        shape_section = f"""
        <div class="network-legend-section">
            <div class="network-legend-heading">{html.escape(shape_column)}</div>
            {shape_items}
        </div>
        """

    return f"""
    <style>
        #transmission-network-panel {{
            position: fixed;
            top: 18px;
            left: 18px;
            z-index: 1000;
            width: 300px;
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
        #transmission-network-panel h1 {{
            margin: 0 0 5px 0;
            font-size: 18px;
            line-height: 1.25;
        }}
        .network-subtitle {{
            margin-bottom: 8px;
            color: #555555;
            font-size: 13px;
        }}
        .network-summary {{
            margin-bottom: 12px;
            color: #555555;
            font-size: 12px;
        }}
        .network-legend-section {{ margin-top: 12px; }}
        .network-legend-heading {{
            margin-bottom: 7px;
            font-weight: bold;
            font-size: 13px;
        }}
        .network-legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
            margin: 4px 0;
            font-size: 12px;
        }}
        .network-colour {{
            display: inline-block;
            width: 13px;
            height: 13px;
            border: 1px solid #444444;
            border-radius: 50%;
            flex: 0 0 auto;
        }}
        .network-shape-symbol {{
            display: inline-flex;
            width: 20px;
            height: 20px;
            align-items: center;
            justify-content: center;
            flex: 0 0 auto;
        }}
        .network-shape-svg {{
            display: block;
            width: 18px;
            height: 18px;
            overflow: visible;
        }}
        .network-note {{
            margin-top: 12px;
            padding-top: 10px;
            border-top: 1px solid #dddddd;
            color: #555555;
            font-size: 11px;
            line-height: 1.35;
        }}
        #mynetwork {{ height: 100vh !important; }}
    </style>
    <div id="transmission-network-panel">
        <h1>{html.escape(title)}</h1>
        <div class="network-subtitle">
            {html.escape(group_label)}:
            <strong>{html.escape(group_value)}</strong>
        </div>
        <div class="network-summary">
            {node_count:,} clustered isolates · {len(cluster_colours):,} transmission clusters
        </div>
        <div class="network-legend-section">
            <div class="network-legend-heading">Transmission cluster</div>
            {cluster_items}
        </div>
        {shape_section}
        <div class="network-note">
            <strong>How to interpret this network</strong><br><br>
                Nodes represent isolates assigned to the same inferred transmission
                cluster. Edges form a minimum spanning tree based on SNP distance
                for visualization and do not represent direct transmission or
                necessarily satisfy the configured temporal cutoff of
                {days_cutoff:g} days.<br><br>
                Transmission-cluster membership may reflect a qualifying chain
                through intermediate isolates. Solid edges connect isolates whose
                direct collection-date difference is within the temporal cutoff.
                Dashed edges connect isolates whose direct collection-date
                difference exceeds the cutoff.<br><br>
                Review the corresponding transmission scatterplot for the temporal
                distribution and possible time-based chain of cases. Hover over
                nodes and edges for metadata and distance details.

        </div>
    </div>
    """

def _inject_html_panel(output_path: Path, panel_html: str) -> None:
    """Inject a custom title and legend panel into generated PyVis HTML."""
    document = output_path.read_text(encoding="utf-8")
    body_marker = "<body>"

    if body_marker not in document:
        raise ValueError("Could not locate the HTML body in PyVis output")

    document = document.replace(body_marker, body_marker + panel_html, 1)
    output_path.write_text(document, encoding="utf-8")


def save_interactive_transmission_network(
    graph: nx.Graph,
    group_df: pd.DataFrame,
    output_path: str | Path,
    *,
    sample_id_column: str,
    group_label: str,
    group_value: str,
    cluster_column: str = "Clusters",
    shape_column: str | None = None,
    tooltip_columns: Sequence[str] = (),
    title: str = "Transmission-cluster network",
    days_cutoff: float,
) -> Path:
    """Render one interactive transmission-cluster network."""
    try:
        from pyvis.network import Network
    except ImportError as exc:
        raise RuntimeError(
            "Interactive transmission networks require pyvis"
        ) from exc

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cluster_colours = _category_colour_map(group_df["_transmission_cluster"])
    shape_values = _shape_map(group_df["_shape_group"])
    show_labels = graph.number_of_nodes() <= 35

    network = Network(
        height="100vh",
        width="100%",
        directed=False,
        bgcolor="#ffffff",
        font_color="#222222",
        heading="",
        cdn_resources="remote",
    )

    for node in graph.nodes:
        node_id = str(node)
        row = group_df.loc[node_id]
        cluster_label = row["_transmission_cluster"]
        shape_label = row["_shape_group"]
        tooltip = [
            f"<b>Sample:</b> {html.escape(node_id)}",
            (
                "<b>Transmission cluster:</b> "
                f"{html.escape(cluster_label)}"
            ),
        ]

        for column in tooltip_columns:
            if column in group_df.columns:
                tooltip.append(
                    f"<b>{html.escape(str(column))}:</b> "
                    f"{html.escape(_normalise_category(row[column]))}"
                )

        network.add_node(
            node_id,
            label=node_id if show_labels else "",
            title="<br>".join(tooltip),
            shape=shape_values[shape_label],
            size=18,
            color={
                "background": cluster_colours[cluster_label],
                "border": "#444444",
                "highlight": {
                    "background": cluster_colours[cluster_label],
                    "border": "#000000",
                },
            },
            borderWidth=1.5,
        )

    for node_a, node_b, edge_data in graph.edges(data=True):
        snp_distance = float(edge_data.get("snp_distance", 0.0))
        day_difference = edge_data.get("day_difference", pd.NA)
        edge_title = [f"SNP distance: {snp_distance:g}"]

        if pd.notna(day_difference):
            edge_title.append(f"Collection-date difference: {int(day_difference)} days")

        visual_length = 65 + min(max(snp_distance, 0.0), 100.0) * 2.5
        network.add_edge(
            str(node_a),
            str(node_b),
            title="<br>".join(edge_title),
            label=(
                f"{snp_distance:g}"
                if snp_distance > 0 and graph.number_of_nodes() <= 30
                else ""
            ),
            length=visual_length,
            color="#888888",
            width=1.4,
        )

    for node_a, node_b, edge_data in graph.edges(data=True):

        snp_distance = float(
            edge_data.get(
                "snp_distance",
                0.0,
            )
        )

        day_difference = edge_data.get(
            "day_difference",
            pd.NA,
        )

        exceeds_days_cutoff = (
            pd.notna(day_difference)
            and float(day_difference)
            > float(days_cutoff)
        )

        edge_title = [
            "<b>Genetic visualization edge</b>",
            f"SNP distance: {snp_distance:g}",
        ]

        if pd.notna(day_difference):

            edge_title.extend(
                [
                    (
                        "Collection-date difference: "
                        f"{int(day_difference)} days"
                    ),
                    (
                        "Configured temporal cutoff: "
                        f"{days_cutoff:g} days"
                    ),
                    (
                        "Direct temporal criterion met: "
                        f"{'No' if exceeds_days_cutoff else 'Yes'}"
                    ),
                ]
            )

        else:

            edge_title.extend(
                [
                    (
                        "Collection-date difference: "
                        "Unavailable"
                    ),
                    (
                        "Configured temporal cutoff: "
                        f"{days_cutoff:g} days"
                    ),
                ]
            )

        edge_title.append(
            "This post-clustering SNP minimum-spanning-tree "
            "edge is not evidence of direct transmission."
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
            title="<br>".join(
                edge_title
            ),
            label=(
                f"{snp_distance:g}"
                if (
                    snp_distance > 0
                    and graph.number_of_nodes() <= 30
                )
                else ""
            ),
            length=visual_length,
            color={
                "color": (
                    "#B05A5A"
                    if exceeds_days_cutoff
                    else "#888888"
                ),
                "highlight": (
                    "#8F3F3F"
                    if exceeds_days_cutoff
                    else "#555555"
                ),
                "hover": (
                    "#8F3F3F"
                    if exceeds_days_cutoff
                    else "#555555"
                ),
                "inherit": False,
            },
            dashes=(
                [8, 6]
                if exceeds_days_cutoff
                else False
            ),
            width=(
                2.2
                if exceeds_days_cutoff
                else 1.4
            ),
            smooth=False,
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
            "font": {"size": 12, "face": "Arial"}
          },
          "edges": {
            "smooth": false,
            "font": {"size": 10, "align": "middle"}
          }
        }
        """
    )

    network.write_html(str(output_path), open_browser=False)
    panel_html = _build_html_panel(
        title=title,
        group_label=group_label,
        group_value=group_value,
        cluster_colours=cluster_colours,
        shape_column=shape_column,
        shape_values=shape_values,
        node_count=graph.number_of_nodes(),
        days_cutoff=days_cutoff,
    )
    _inject_html_panel(output_path, panel_html)
    return output_path


def generate_transmission_networks_from_config(
    transmission_df: pd.DataFrame,
    snp_matrix: pd.DataFrame,
    config: Any,
    output_dir: str | Path,
    *,
    cluster_column: str = "Clusters",
    filename_prefix: str = "05_transmission_network",
) -> dict[str, object]:
    """Generate all mode-aware transmission networks and output tables."""
    if not bool(_get(config, "outputs", "generate_figures", default=True)):
        return {
            "paths": {},
            "nodes": pd.DataFrame(),
            "edges": pd.DataFrame(),
        }

    cluster_type = str(
        _get(config, "analysis", "cluster_type", default="")
    ).strip()

    if cluster_type.casefold() != "transmission":
        return {
            "paths": {},
            "nodes": pd.DataFrame(),
            "edges": pd.DataFrame(),
        }

    main_var = _get(config, "variables", "main_var")
    var_01 = _get(config, "variables", "var_01")
    var_02 = _get(config, "variables", "var_02")
    sequence_type_column = _get(config, "columns", "sequence_type")
    configured_sample_id = _get(config, "columns", "sample_id")
    transmission_level = str(
        _get(config, "analysis", "transmission_level", default="Facility")
    ).strip().casefold()

    days_cutoff = float(
        _get(
            config,
            "analysis",
            "days_cutoff",
            default=45,
        )
    )

    if not main_var:
        raise ValueError(
            "Transmission networks require variables.main_var"
        )
    if not var_02:
        raise ValueError(
            "Transmission networks require variables.var_02"
        )

    matrix_ids = set(snp_matrix.index.astype(str)) & set(
        snp_matrix.columns.astype(str)
    )
    sample_id_column = _resolve_sample_id_column(
        transmission_df,
        matrix_ids,
        candidates=(configured_sample_id, "sample_id"),
    )

    shape_column = (
        str(var_01)
        if var_01 and var_01 in transmission_df.columns
        else None
    )

    prepared = prepare_transmission_network_data(
        transmission_df,
        snp_matrix,
        sample_id_column=sample_id_column,
        main_var=str(main_var),
        collection_date_column=str(var_02),
        cluster_column=cluster_column,
        shape_column=shape_column,
        date_format=_get(config, "analysis", "date_format", default="ymd"),
    )

    if prepared.empty:
        return {
            "paths": {},
            "nodes": pd.DataFrame(),
            "edges": pd.DataFrame(),
        }

    tooltip_columns = [
        column
        for column in (
            main_var,
            var_01,
            sequence_type_column,
            var_02,
            "_core_context",
            "core_cluster_label",
            "SNPs",
            "Days",
            "Cluster_Cases_count",
        )
        if column and column in prepared.columns
    ]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    node_tables: list[pd.DataFrame] = []
    edge_tables: list[pd.DataFrame] = []

    for group_value, group_df in prepared.groupby(
        "_network_group",
        dropna=False,
        sort=False,
    ):
        if group_value == _MISSING_LABEL or len(group_df) < 1:
            continue

        # graph, group_edges = build_transmission_cluster_network(
        #     group_df,
        #     snp_matrix,
        #     sample_id_column=sample_id_column,
        # )

        # if graph.number_of_nodes() == 0:
        #     continue

        graph, group_edges = build_transmission_cluster_network(
            group_df,
            snp_matrix,
            sample_id_column=sample_id_column,
        )

        # Record the configured temporal cutoff and identify
        # post-clustering MST edges whose direct collection-date
        # difference exceeds that cutoff.
        group_edges[
            "days_cutoff"
        ] = days_cutoff

        group_edges[
            "exceeds_days_cutoff"
        ] = (
            pd.to_numeric(
                group_edges[
                    "day_difference"
                ],
                errors="coerce",
            )
            > days_cutoff
        )

        if graph.number_of_nodes() == 0:
            continue

        group_label = str(main_var)
        output_group = str(group_value)

        if transmission_level == "community" and (
            sequence_type_column and str(main_var) == str(sequence_type_column)
        ):
            filename_group = f"ST_{_safe_name(group_value)}"
        else:
            filename_group = _safe_name(group_value)

        output_path = output_dir / f"{filename_prefix}_{filename_group}.html"
        paths[output_group] = (
            save_interactive_transmission_network(
                graph,
                group_df,
                output_path,
                sample_id_column=sample_id_column,
                group_label=group_label,
                group_value=output_group,
                cluster_column=cluster_column,
                shape_column=shape_column,
                tooltip_columns=tooltip_columns,
                days_cutoff=days_cutoff,
            )
        )

        node_columns = [
            column
            for column in (
                sample_id_column,
                cluster_column,
                main_var,
                var_01,
                sequence_type_column,
                var_02,
                "_core_context",
                "core_cluster_label",
                "Cluster_Cases_count",
            )
            if column and column in group_df.columns
        ]
        group_nodes = group_df.loc[:, node_columns].copy().reset_index(drop=True)
        group_nodes.insert(0, "network_group", output_group)
        node_tables.append(group_nodes)

        group_edges.insert(0, "network_group", output_group)
        edge_tables.append(group_edges)

    nodes = (
        pd.concat(node_tables, ignore_index=True)
        if node_tables
        else pd.DataFrame()
    )
    edges = (
        pd.concat(edge_tables, ignore_index=True)
        if edge_tables
        else pd.DataFrame()
    )

    return {
        "paths": paths,
        "nodes": nodes,
        "edges": edges,
    }
