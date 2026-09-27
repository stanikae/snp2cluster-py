from __future__ import annotations

import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from ..date_utils import parse_dates

_BASE_CLUSTER_COLOURS = (
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#F0E442", "#332288", "#88CCEE", "#44AA99",
    "#117733", "#999933",
)
_MARKERS = ("o", "^", "s", "D", "P", "X", "v", "<", ">", "*")
_UNCLUSTERED = "Unclustered"
_UNCLUSTERED_COLOUR = "#B8B8B8"


def _get(config: Any, *keys: str, default: Any = None) -> Any:
    """Read nested configuration from dictionaries or attribute-based models."""
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
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", str(value))
    )


def _ordered_labels(series: pd.Series) -> list[str]:
    non_missing = series.dropna().astype(str)
    if isinstance(series.dtype, pd.CategoricalDtype) and series.dtype.ordered:
        observed = set(non_missing)
        return [str(v) for v in series.cat.categories if str(v) in observed]
    return sorted(non_missing.unique().tolist(), key=_natural_key)


def _normalise_cluster(value: object) -> str:
    if pd.isna(value) or str(value).strip().casefold() in {"", "na", "nan", "none"}:
        return _UNCLUSTERED
    return str(value)


def _cluster_order(series: pd.Series) -> list[str]:
    labels = {_normalise_cluster(v) for v in series}
    ordered = sorted((v for v in labels if v != _UNCLUSTERED), key=_natural_key)
    if _UNCLUSTERED in labels:
        ordered.append(_UNCLUSTERED)
    return ordered


def _cluster_palette(labels: list[str]) -> dict[str, object]:
    palette: dict[str, object] = {}
    fallback = plt.get_cmap("tab20")
    for index, label in enumerate(v for v in labels if v != _UNCLUSTERED):
        palette[label] = (
            _BASE_CLUSTER_COLOURS[index]
            if index < len(_BASE_CLUSTER_COLOURS)
            else fallback(index % fallback.N)
        )
    if _UNCLUSTERED in labels:
        palette[_UNCLUSTERED] = _UNCLUSTERED_COLOUR
    return palette


def _safe_prefix(value: str) -> str:
    cleaned = value.strip()
    for character in (" ", "/", "\\", ":"):
        cleaned = cleaned.replace(character, "_")
    return cleaned or "02_transmission_scatter"


def _derive_epiweek(collection_dates: pd.Series, date_format: str | None) -> pd.Series:
    """Convert collection dates to ISO epidemiological year-week labels."""
    format_map = {
        "ymd": "%Y-%m-%d",
        "dmy": "%d-%m-%Y",
        "mdy": "%m-%d-%Y",
    }
    fmt = format_map.get(str(date_format).casefold()) if date_format else None
    # dates = pd.to_datetime(collection_dates, format=fmt, errors="coerce")
    dates = parse_dates(collection_dates,  date_format=date_format,)

    if dates.notna().sum() == 0:
        raise ValueError("The configured collection-date column contains no parseable dates")
    iso = dates.dt.isocalendar()
    return iso["year"].astype("Int64").astype(str) + "." + iso["week"].astype("Int64").astype(str).str.zfill(2)


# def prepare_transmission_plot_data(
#     transmission_df: pd.DataFrame,
#     *,
#     sequence_type_column: str,
#     collection_date_column: str,
#     cluster_column: str,
#     grouping_column: str | None,
#     epiweek_column: str | None = None,
#     date_format: str | None = None,
# ) -> pd.DataFrame:
#     """Normalise configured metadata fields to the internal plotting schema."""
#     required = {sequence_type_column, collection_date_column, cluster_column}
#     if grouping_column is not None:
#         required.add(grouping_column)
#     if epiweek_column is not None:
#         required.add(epiweek_column)
#     missing = required - set(transmission_df.columns)
#     if missing:
#         raise ValueError(f"Missing columns: {missing}")

#     columns = list(required)
#     plot_df = transmission_df.loc[:, columns].copy()
#     plot_df["_st"] = plot_df[sequence_type_column]
#     plot_df["_cluster"] = plot_df[cluster_column]
#     plot_df["_collection_date"] = pd.to_datetime(
#         plot_df[collection_date_column], errors="coerce"
#     )
#     if epiweek_column is None:
#         plot_df["_epiweek"] = _derive_epiweek(
#             plot_df[collection_date_column], date_format
#         )
#     else:
#         plot_df["_epiweek"] = plot_df[epiweek_column]
#     if grouping_column is not None:
#         plot_df["_group"] = plot_df[grouping_column]
#     return plot_df

def prepare_transmission_plot_data(
    transmission_df: pd.DataFrame,
    *,
    sequence_type_column: str,
    collection_date_column: str,
    cluster_column: str,
    grouping_column: str | None,
    y_axis_column: str | None = None,
    epiweek_column: str | None = None,
    date_format: str | None = None,
) -> pd.DataFrame:
    """Normalise configured metadata fields for transmission plotting."""

    if y_axis_column is None:
        y_axis_column = sequence_type_column

    required = {
        sequence_type_column,
        collection_date_column,
        cluster_column,
        y_axis_column,
    }

    if grouping_column is not None:
        required.add(grouping_column)

    if epiweek_column is not None:
        required.add(epiweek_column)

    missing = required - set(
        transmission_df.columns
    )

    if missing:
        raise ValueError(
            f"Missing columns: {missing}"
        )

    plot_df = transmission_df.loc[
        :,
        list(required),
    ].copy()

    plot_df["_st"] = (
        plot_df[sequence_type_column]
    )

    plot_df["_y"] = (
        plot_df[y_axis_column]
    )

    plot_df["_cluster"] = (
        plot_df[cluster_column]
    )

    plot_df["_collection_date"] = parse_dates(
        plot_df[collection_date_column],
        date_format=date_format or "ymd",
    )

    if epiweek_column is None:
        plot_df["_epiweek"] = _derive_epiweek(
            plot_df[collection_date_column],
            date_format,
        )
    else:
        plot_df["_epiweek"] = (
            plot_df[epiweek_column]
        )

    if grouping_column is not None:
        plot_df["_group"] = (
            plot_df[grouping_column]
        )

    return plot_df


def plot_transmission_scatter(
    transmission_df: pd.DataFrame,
    output_dir: str | Path,
    *,
    sequence_type_column: str,
    collection_date_column: str,
    cluster_column: str,
    grouping_column: str | None,
    y_axis_column: str | None = None,
    y_axis_label: str | None = None,
    epiweek_column: str | None = None,
    date_format: str | None = None,
    transmission_level: str | None = None,
    title: str | None = None,
    filename_prefix: str = "02_transmission_scatter",
    snp_cutoff: int | float | None = None,
    day_window: int | None = None,
    random_seed: int = 12,
    figsize: tuple[float, float] = (14, 7),
) -> dict[str, Path]:
    """Generate the open-mode transmission-cluster scatterplot.

    All metadata field names are explicit arguments. The function therefore
    contains no assumptions about names such as WardType, FacilityName or ST.
    """
    if not isinstance(transmission_df, pd.DataFrame):
        raise TypeError("transmission_df must be a pandas DataFrame")
    if transmission_df.empty:
        raise ValueError("transmission_df is empty")

    # plot_df = prepare_transmission_plot_data(
    #     transmission_df,
    #     sequence_type_column=sequence_type_column,
    #     collection_date_column=collection_date_column,
    #     cluster_column=cluster_column,
    #     grouping_column=grouping_column,
    #     epiweek_column=epiweek_column,
    #     date_format=date_format,
    # )
    plot_df = prepare_transmission_plot_data(
        transmission_df,
        sequence_type_column=sequence_type_column,
        collection_date_column=collection_date_column,
        cluster_column=cluster_column,
        grouping_column=grouping_column,
        y_axis_column=y_axis_column,
        epiweek_column=epiweek_column,
        date_format=date_format,
    )
    # plot_df = plot_df.dropna(subset=["_epiweek", "_st"])
    # if plot_df.empty:
    #     raise ValueError("No rows remain after removing missing Epiweek or sequence type values")

    # epiweek_order = _ordered_labels(plot_df["_epiweek"])
    # st_order = _ordered_labels(plot_df["_st"])
    # plot_df["_epiweek_label"] = plot_df["_epiweek"].astype(str)
    # plot_df["_st_label"] = plot_df["_st"].astype(str)

    plot_df = plot_df.dropna(
    subset=[
        "_epiweek",
        "_y",
        ]
    )

    if plot_df.empty:
        raise ValueError(
            "No rows remain after removing missing "
            "Epiweek or y-axis values"
        )

    epiweek_order = _ordered_labels(
        plot_df["_epiweek"]
    )

    y_order = _ordered_labels(
        plot_df["_y"]
    )

    plot_df["_epiweek_label"] = (
        plot_df["_epiweek"].astype(str)
    )

    plot_df["_y_label"] = (
        plot_df["_y"].astype(str)
    )
    plot_df["_cluster_label"] = plot_df["_cluster"].map(_normalise_cluster)
    if grouping_column is None:
        plot_df["_group_label"] = "Isolate"
        grouping_title = None
    else:
        plot_df["_group_label"] = plot_df["_group"].map(
            lambda v: "Unknown" if pd.isna(v) or not str(v).strip() else str(v)
        )
        grouping_title = grouping_column

    # x_positions = {label: i for i, label in enumerate(epiweek_order)}
    # # y_positions = {label: i for i, label in enumerate(st_order)}
    # # plot_df["_x"] = plot_df["_epiweek_label"].map(x_positions).astype(float)
    # # plot_df["_y"] = plot_df["_st_label"].map(y_positions).astype(float)

    # y_positions = {
    #     label: i
    #     for i, label in enumerate(y_order)
    # }

    # plot_df["_y_position"] = (
    #     plot_df["_y_label"]
    #     .map(y_positions)
    #     .astype(float)
    # )




    # rng = np.random.default_rng(random_seed)
    # plot_df["_xj"] = plot_df["_x"] + rng.uniform(-0.15, 0.15, len(plot_df))
    # # plot_df["_yj"] = plot_df["_y"] + rng.uniform(-0.07, 0.07, len(plot_df))
    # plot_df["_yj"] = (
    #     plot_df["_y_position"]
    #     + rng.uniform(
    #         -0.07,
    #         0.07,
    #         len(plot_df),
    #     )
    # )

    x_positions = {
        label: index
        for index, label in enumerate(
            epiweek_order
        )
    }

    y_positions = {
        label: index
        for index, label in enumerate(
            y_order
        )
    }

    plot_df["_x"] = (
        plot_df["_epiweek_label"]
        .map(x_positions)
        .astype(float)
    )

    plot_df["_y_position"] = (
        plot_df["_y_label"]
        .map(y_positions)
        .astype(float)
    )

    rng = np.random.default_rng(
        random_seed
    )

    plot_df["_xj"] = (
        plot_df["_x"]
        + rng.uniform(
            -0.15,
            0.15,
            len(plot_df),
        )
    )

    plot_df["_yj"] = (
        plot_df["_y_position"]
        + rng.uniform(
            -0.07,
            0.07,
            len(plot_df),
        )
    )

    clusters = _cluster_order(plot_df["_cluster"])
    colours = _cluster_palette(clusters)
    groups = sorted(plot_df["_group_label"].unique().tolist(), key=_natural_key)
    markers = {value: _MARKERS[i % len(_MARKERS)] for i, value in enumerate(groups)}

    fig, ax = plt.subplots(figsize=figsize)
    draw_order = sorted(clusters, key=lambda v: (v != _UNCLUSTERED, _natural_key(v)))
    for cluster in draw_order:
        for group in groups:
            subset = plot_df.loc[
                (plot_df["_cluster_label"] == cluster)
                & (plot_df["_group_label"] == group)
            ]
            if subset.empty:
                continue
            ax.scatter(
                subset["_xj"], subset["_yj"],
                c=[colours[cluster]], marker=markers[group], s=72,
                alpha=0.72 if cluster == _UNCLUSTERED else 0.88,
                edgecolors="white", linewidths=0.55,
                zorder=2 if cluster == _UNCLUSTERED else 3,
            )

    ax.set_xticks(range(len(epiweek_order)), epiweek_order, rotation=90)
    # ax.set_yticks(range(len(st_order)), st_order)
    ax.set_xlabel("Epidemiological week")
    # ax.set_ylabel("Sequence type")
    ax.set_xlim(-0.55, len(epiweek_order) - 0.45)
    # ax.set_ylim(-0.55, len(st_order) - 0.45)
    ax.set_yticks(
        range(len(y_order)),
        y_order,
    )

    ax.set_ylabel(
        y_axis_label
        or y_axis_column
        or sequence_type_column
    )

    ax.set_ylim(
        -0.55,
        len(y_order) - 0.45,
    )

    # if title is None:
    #     heading = "Transmission clusters"
    #     if transmission_level:
    #         heading = f"{transmission_level} transmission clusters"
    #     details = []
    #     if snp_cutoff is not None:
    #         details.append(f"SNP threshold <= {snp_cutoff}")
    #     if day_window is not None:
    #         details.append(f"temporal window = {day_window} days")
    #     title = heading + ("\n" + " | ".join(details) if details else "")

    details = []

    if snp_cutoff is not None:
        details.append(f"SNP threshold <= {snp_cutoff}")

    if day_window is not None:
        details.append(f"temporal window = {day_window} days")

    if title is None:

        heading = "Transmission clusters"

        if transmission_level:
            heading = f"{transmission_level} transmission clusters"

    else:

        if transmission_level:
            heading = (
                f"{title}\n"
                f"{transmission_level} transmission clusters"
            )
        else:
            heading = (
                f"{title}\n"
                f"Transmission clusters"
            )

    title = heading + (
        "\n" + " | ".join(details)
        if details else ""
    )

    ax.set_title(title, loc="left", pad=12)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    cluster_handles = [
        Line2D([0], [0], marker="o", linestyle="none",
               markerfacecolor=colours[v], markeredgecolor="white",
               markersize=8, label=v)
        for v in clusters
    ]
    first_legend = ax.legend(
        handles=cluster_handles, title="Transmission cluster",
        loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
        borderaxespad=0,
    )
    ax.add_artist(first_legend)
    # if grouping_title is not None:
    #     group_handles = [
    #         Line2D([0], [0], marker=markers[v], linestyle="none",
    #                markerfacecolor="#4D4D4D", markeredgecolor="white",
    #                markersize=8, label=v)
    #         for v in groups
    #     ]
    #     second_legend = ax.legend(
    #         handles=group_handles, title=grouping_title,
    #         loc="upper left", bbox_to_anchor=(1.01, 0.56), frameon=False,
    #         borderaxespad=0,
    #     )
    #     ax.add_artist(second_legend)

    # fig.tight_layout()
    # output_dir = Path(output_dir)
    # output_dir.mkdir(parents=True, exist_ok=True)
    # prefix = _safe_prefix(filename_prefix)
    # paths = {
    #     "png": output_dir / f"{prefix}.png",
    #     "pdf": output_dir / f"{prefix}.pdf",
    # }
    # fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    # # fig.savefig(paths["pdf"], bbox_inches="tight")
    # fig.savefig(
    #     paths["pdf"],
    #     bbox_inches="tight",
    #     bbox_extra_artists=[
    #         first_legend,
    #         second_legend,
    #     ],
    # )
    # plt.close(fig)
    # return paths

    second_legend = None

    if grouping_title is not None:

        group_handles = [
            Line2D(
                [0],
                [0],
                marker=markers[value],
                linestyle="none",
                markerfacecolor="#4D4D4D",
                markeredgecolor="white",
                markersize=8,
                label=value,
            )
            for value in groups
        ]

        second_legend = ax.legend(
            handles=group_handles,
            title=grouping_title,
            loc="upper left",
            bbox_to_anchor=(1.01, 0.56),
            frameon=False,
            borderaxespad=0,
        )

        ax.add_artist(
            second_legend
        )


    # Include only legends that were created.
    extra_artists = [
        first_legend,
    ]

    if second_legend is not None:

        extra_artists.append(
            second_legend
        )


    fig.tight_layout()

    output_dir = Path(
        output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    prefix = _safe_prefix(
        filename_prefix
    )

    paths = {
        "png": output_dir / f"{prefix}.png",
        "pdf": output_dir / f"{prefix}.pdf",
    }

    fig.savefig(
        paths["png"],
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        paths["pdf"],
        bbox_inches="tight",
        bbox_extra_artists=extra_artists,
    )

    plt.close(
        fig
    )

    return paths



# def generate_transmission_scatter_from_config(
#     transmission_df: pd.DataFrame,
#     config: Any,
#     output_dir: str | Path,
#     *,
#     cluster_column: str = "Clusters",
#     y_axis_column: str | None = None,
#     y_axis_label: str | None = None,
#     grouping_column: str | None = None,
#     epiweek_column: str | None = None,
#     filename_prefix: str = "02_transmission_scatter",
#     title: str | None = None,
# ) -> dict[str, Path]:
#     """Config-aware pipeline entry point.

#     Returns an empty dictionary when figures are disabled or cluster_type is
#     not Transmission. This makes the function safe to call unconditionally
#     from pipeline.py while ensuring temporal plots are never produced in Core
#     mode.
#     """
#     generate_figures = bool(_get(config, "outputs", "generate_figures", default=True))
#     cluster_type = str(_get(config, "analysis", "cluster_type", default="")).strip()
#     if not generate_figures or cluster_type.casefold() != "transmission":
#         return {}

#     sequence_type_column = _get(config, "columns", "sequence_type")
#     collection_date_column = _get(config, "variables", "var_02")
#     grouping_column = _get(config, "variables", "var_01")
#     missing_config = [
#         name for name, value in {
#             "columns.sequence_type": sequence_type_column,
#             "variables.var_02": collection_date_column,
#         }.items() if not value
#     ]
#     if missing_config:
#         raise ValueError(f"Missing configuration values: {set(missing_config)}")

#     return plot_transmission_scatter(
#         transmission_df,
#         output_dir,
#         sequence_type_column=str(sequence_type_column),
#         collection_date_column=str(collection_date_column),
#         cluster_column=cluster_column,
#         grouping_column=str(grouping_column) if grouping_column else None,
#         epiweek_column=epiweek_column,
#         date_format=_get(config, "analysis", "date_format"),
#         transmission_level=_get(config, "analysis", "transmission_level"),
#         filename_prefix=filename_prefix,
#         title=title,
#         snp_cutoff=_get(config, "analysis", "snp_cutoff"),
#         day_window=_get(config, "analysis", "days_cutoff"),
#         random_seed=int(_get(config, "analysis", "random_seed", default=12)),
#     )


def generate_transmission_scatter_from_config(
    transmission_df: pd.DataFrame,
    config: Any,
    output_dir: str | Path,
    *,
    cluster_column: str = "Clusters",
    y_axis_column: str | None = None,
    y_axis_label: str | None = None,
    grouping_column: str | None = None,
    epiweek_column: str | None = None,
    filename_prefix: str = "02_transmission_scatter",
    title: str | None = None,
) -> dict[str, Path]:
    """Generate a transmission scatterplot using configuration values.

    Facility and Area modes:
        - y-axis: configured sequence-type column
        - marker shape: variables.var_01

    Community mode:
        - y-axis: variables.var_01, normally FacilityName
        - marker shape: disabled

    Returns an empty dictionary when figure generation is disabled or
    cluster_type is not Transmission.

    Parameters
    ----------
    transmission_df : pd.DataFrame
        Transmission-context dataset containing clustered and unclustered
        isolates.

    config : Any
        SNP2Cluster configuration object or mapping.

    output_dir : str | Path
        Directory where the PNG and PDF files will be written.

    cluster_column : str
        Column containing transmission-cluster assignments.

    y_axis_column : str | None
        Optional explicit y-axis column. When omitted, the column is resolved
        from transmission_level.

    y_axis_label : str | None
        Optional display label for the y-axis.

    grouping_column : str | None
        Optional explicit marker-shape column for Facility and Area modes.
        Community mode disables marker-shape encoding.

    epiweek_column : str | None
        Optional existing epidemiological-week column. When omitted,
        epidemiological weeks are derived from the configured date column.

    filename_prefix : str
        Prefix used for the PNG and PDF filenames.

    title : str | None
        Optional facility, area, or sequence-type title.

    Returns
    -------
    dict[str, Path]
        Paths to the generated PNG and PDF files.
    """

    generate_figures = bool(
        _get(
            config,
            "outputs",
            "generate_figures",
            default=True,
        )
    )

    cluster_type = str(
        _get(
            config,
            "analysis",
            "cluster_type",
            default="",
        )
    ).strip()

    if (
        not generate_figures
        or cluster_type.casefold() != "transmission"
    ):
        return {}

    transmission_level = str(
        _get(
            config,
            "analysis",
            "transmission_level",
            default="Facility",
        )
    ).strip()

    transmission_level_normalised = (
        transmission_level.casefold()
    )

    supported_levels = {
        "facility",
        "area",
        "community",
    }

    if transmission_level_normalised not in supported_levels:
        raise ValueError(
            "Unsupported transmission_level: "
            f"{transmission_level}. "
            "Expected Facility, Area, or Community."
        )

    sequence_type_column = _get(
        config,
        "columns",
        "sequence_type",
    )

    collection_date_column = _get(
        config,
        "variables",
        "var_02",
    )

    configured_var_01 = _get(
        config,
        "variables",
        "var_01",
    )

    missing_config = [
        name
        for name, value in {
            "columns.sequence_type": sequence_type_column,
            "variables.var_02": collection_date_column,
        }.items()
        if not value
    ]

    if missing_config:
        raise ValueError(
            "Missing configuration values: "
            f"{set(missing_config)}"
        )

    # ---------------------------------------------------------
    # Resolve Community-mode visualization fields
    # ---------------------------------------------------------

    if transmission_level_normalised == "community":

        if not configured_var_01:
            raise ValueError(
                "Community transmission plots require "
                "variables.var_01 to define the y-axis, "
                "for example FacilityName."
            )

        if y_axis_column is None:
            y_axis_column = str(
                configured_var_01
            )

        if y_axis_label is None:
            y_axis_label = str(
                configured_var_01
            )

        # Community mode displays var_01 on the y-axis.
        # Do not encode the same field again using marker shapes.
        resolved_grouping_column = None

    # ---------------------------------------------------------
    # Resolve Facility- and Area-mode visualization fields
    # ---------------------------------------------------------

    else:

        if y_axis_column is None:
            y_axis_column = str(
                sequence_type_column
            )

        if y_axis_label is None:
            y_axis_label = "Sequence type"

        if grouping_column is not None:
            resolved_grouping_column = str(
                grouping_column
            )

        elif configured_var_01:
            resolved_grouping_column = str(
                configured_var_01
            )

        else:
            resolved_grouping_column = None

    return plot_transmission_scatter(
        transmission_df,
        output_dir,
        sequence_type_column=str(
            sequence_type_column
        ),
        collection_date_column=str(
            collection_date_column
        ),
        cluster_column=cluster_column,
        grouping_column=resolved_grouping_column,
        y_axis_column=str(
            y_axis_column
        ),
        y_axis_label=(
            str(y_axis_label)
            if y_axis_label is not None
            else None
        ),
        epiweek_column=epiweek_column,
        date_format=_get(
            config,
            "analysis",
            "date_format",
        ),
        transmission_level=transmission_level,
        filename_prefix=filename_prefix,
        title=title,
        snp_cutoff=_get(
            config,
            "analysis",
            "snp_cutoff",
        ),
        day_window=_get(
            config,
            "analysis",
            "days_cutoff",
        ),
        random_seed=int(
            _get(
                config,
                "analysis",
                "random_seed",
                default=12,
            )
        ),
    )
