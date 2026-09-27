from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


_MISSING_LABEL = "Missing"
_MISSING_COLOUR = "#D3D3D3"
_DEFAULT_REFERENCE_LABELS = ("reference",)


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
    """Return a natural-sort key for numeric and mixed-text category labels."""
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", str(value))
    )


def _normalise_category(value: object) -> str:
    """Convert missing values to a stable display label."""
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


def _safe_prefix(value: str) -> str:
    """Return a filesystem-safe output filename prefix."""
    cleaned = value.strip()

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

    return cleaned or "03_core_heatmap"


def _resolve_sample_id_column(
    annotation_df: pd.DataFrame,
    matrix_ids: set[str],
    candidates: Sequence[str | None],
) -> str:
    """Choose the available sample-ID column with best SNP-matrix coverage."""
    available = []

    for candidate in candidates:
        if (
            candidate
            and candidate in annotation_df.columns
            and candidate not in available
        ):
            available.append(candidate)

    if not available:
        raise ValueError(
            "No usable sample-ID column was found in the annotation data."
        )

    coverage = {
        column: len(
            set(
                annotation_df[column]
                .dropna()
                .astype(str)
            )
            & matrix_ids
        )
        for column in available
    }

    return max(
        available,
        key=lambda column: coverage[column],
    )


def _category_colour_map(
    series: pd.Series,
) -> dict[str, object]:
    """Assign deterministic colours to observed annotation categories."""
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
        n_colors=max(len(observed), 1),
    )

    colour_map: dict[str, object] = {
        label: palette[index]
        for index, label in enumerate(observed)
    }

    if _MISSING_LABEL in labels:
        colour_map[_MISSING_LABEL] = _MISSING_COLOUR

    return colour_map


def prepare_heatmap_data(
    snp_matrix: pd.DataFrame,
    annotation_df: pd.DataFrame,
    *,
    sample_id_column: str,
    annotation_columns: Sequence[str],
    reference_labels: Sequence[str] = _DEFAULT_REFERENCE_LABELS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Align a symmetric SNP matrix with one-row-per-isolate annotations.

    The returned matrix and annotation table contain the same isolates in the
    same order. Reference rows and columns are removed by label.
    """
    if not isinstance(snp_matrix, pd.DataFrame):
        raise TypeError(
            "snp_matrix must be a pandas DataFrame"
        )

    if not isinstance(annotation_df, pd.DataFrame):
        raise TypeError(
            "annotation_df must be a pandas DataFrame"
        )

    if snp_matrix.empty:
        raise ValueError(
            "snp_matrix is empty"
        )

    if annotation_df.empty:
        raise ValueError(
            "annotation_df is empty"
        )

    if snp_matrix.shape[0] != snp_matrix.shape[1]:
        raise ValueError(
            "snp_matrix must be square"
        )

    matrix = snp_matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)

    if set(matrix.index) != set(matrix.columns):
        raise ValueError(
            "SNP matrix row and column identifiers do not match"
        )

    matrix = matrix.loc[
        matrix.index,
        matrix.index,
    ]

    try:
        matrix = matrix.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "snp_matrix contains non-numeric values"
        ) from exc

    reference_set = {
        str(value).casefold()
        for value in reference_labels
    }

    retained_ids = [
        sample_id
        for sample_id in matrix.index
        if sample_id.casefold() not in reference_set
    ]

    matrix = matrix.loc[
        retained_ids,
        retained_ids,
    ]

    if matrix.empty:
        raise ValueError(
            "No isolates remain after removing reference labels"
        )

    if not np.allclose(
        matrix.to_numpy(),
        matrix.to_numpy().T,
        equal_nan=True,
    ):
        raise ValueError(
            "snp_matrix must be symmetric"
        )

    missing_annotations = (
        set(annotation_columns)
        - set(annotation_df.columns)
    )

    if missing_annotations:
        raise ValueError(
            "Missing annotation columns: "
            f"{missing_annotations}"
        )

    if sample_id_column not in annotation_df.columns:
        raise ValueError(
            "Sample-ID column not found in annotations: "
            f"{sample_id_column}"
        )

    annotations = annotation_df.loc[
        :,
        [
            sample_id_column,
            *annotation_columns,
        ],
    ].copy()

    annotations[sample_id_column] = (
        annotations[sample_id_column]
        .astype(str)
    )

    annotations = annotations.drop_duplicates(
        subset=[sample_id_column],
        keep="first",
    )

    annotations = annotations.set_index(
        sample_id_column
    )

    shared_ids = [
        sample_id
        for sample_id in matrix.index
        if sample_id in annotations.index
    ]

    if not shared_ids:
        raise ValueError(
            "No sample identifiers are shared between the SNP matrix "
            "and annotation data"
        )

    matrix = matrix.loc[
        shared_ids,
        shared_ids,
    ]

    annotations = annotations.loc[
        shared_ids,
        list(annotation_columns),
    ]

    return matrix, annotations


def plot_annotated_snp_heatmap(
    snp_matrix: pd.DataFrame,
    annotation_df: pd.DataFrame,
    output_dir: str | Path,
    *,
    sample_id_column: str,
    annotation_columns: Sequence[str],
    annotation_labels: Mapping[str, str] | None = None,
    title: str = "Core SNP distance heatmap",
    filename_prefix: str = "03_core_heatmap",
    reference_labels: Sequence[str] = _DEFAULT_REFERENCE_LABELS,
    cmap: str = "viridis",
    linkage_method: str = "average",
    figsize: tuple[float, float] = (12, 12),
    show_sample_labels: bool | None = None,
) -> dict[str, Path]:
    """Generate a clustered SNP-distance heatmap with metadata tracks."""
    matrix, annotations = prepare_heatmap_data(
        snp_matrix,
        annotation_df,
        sample_id_column=sample_id_column,
        annotation_columns=annotation_columns,
        reference_labels=reference_labels,
    )

    if show_sample_labels is None:
        show_sample_labels = len(matrix) <= 40

    display_labels = {
        column: (
            annotation_labels.get(column, column)
            if annotation_labels
            else column
        )
        for column in annotation_columns
    }

    colour_maps: dict[str, dict[str, object]] = {}
    colour_tracks = pd.DataFrame(
        index=annotations.index
    )

    for column in annotation_columns:
        normalised = annotations[column].map(
            _normalise_category
        )

        colour_map = _category_colour_map(
            normalised
        )

        colour_maps[column] = colour_map
        colour_tracks[display_labels[column]] = normalised.map(
            colour_map
        )

    if len(matrix) > 1:
        condensed = squareform(
            matrix.to_numpy(),
            checks=True,
        )

        shared_linkage = linkage(
            condensed,
            method=linkage_method,
        )
    else:
        shared_linkage = None

    cluster_grid = sns.clustermap(
        matrix,
        row_linkage=shared_linkage,
        col_linkage=shared_linkage,
        row_cluster=shared_linkage is not None,
        col_cluster=shared_linkage is not None,
        row_colors=colour_tracks,
        col_colors=colour_tracks,
        cmap=cmap,
        xticklabels=show_sample_labels,
        yticklabels=show_sample_labels,
        figsize=figsize,
        linewidths=0,
        cbar_kws={
            "label": "Pairwise core SNP distance",
        },
    )

    cluster_grid.ax_heatmap.set_xlabel(
        "Isolates"
    )

    cluster_grid.ax_heatmap.set_ylabel(
        "Isolates"
    )

    if show_sample_labels:
        cluster_grid.ax_heatmap.tick_params(
            axis="x",
            labelrotation=90,
            labelsize=7,
        )

        cluster_grid.ax_heatmap.tick_params(
            axis="y",
            labelsize=7,
        )

    cluster_grid.fig.suptitle(
        title,
        y=1.02,
    )


    # ---------------------------------------------------------
    # Metadata legends
    #
    # Anchor legends to the figure rather than the heatmap axes.
    # bbox_inches="tight" will expand the exported canvas to
    # include the external legend column.
    # ---------------------------------------------------------

    legend_y = 0.94
    legend_artists = []

    for column in annotation_columns:

        handles = [
            Patch(
                facecolor=colour,
                edgecolor="none",
                label=label,
            )
            for label, colour
            in colour_maps[column].items()
        ]

        legend = cluster_grid.fig.legend(
            handles=handles,
            title=display_labels[column],
            loc="upper left",
            bbox_to_anchor=(
                1.02,
                legend_y,
            ),
            bbox_transform=cluster_grid.fig.transFigure,
            frameon=False,
            borderaxespad=0,
            fontsize=8,
            title_fontsize=9,
            handlelength=1.4,
            handletextpad=0.6,
            labelspacing=0.35,
        )

        legend_artists.append(
            legend
        )

        # Reserve vertical space according to the number
        # of categories in the current legend.
        legend_height = min(
            0.055
            + 0.030 * len(handles),
            0.27,
        )

        legend_y -= legend_height

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

    cluster_grid.fig.savefig(
        paths["png"],
        dpi=300,
        bbox_inches="tight",
        bbox_extra_artists=legend_artists,
    )

    cluster_grid.fig.savefig(
        paths["pdf"],
        bbox_inches="tight",
        bbox_extra_artists=legend_artists,
    )

    plt.close(
        cluster_grid.fig
    )

    return paths


def generate_annotated_heatmap_from_config(
    snp_matrix: pd.DataFrame,
    annotation_df: pd.DataFrame,
    config: Any,
    output_dir: str | Path,
    *,
    core_cluster_column: str = "snp_cluster_id",
    transmission_cluster_column: str = "Clusters",
    filename_prefix: str = "03_core_heatmap",
    title: str = "Core SNP distance heatmap",
) -> dict[str, Path]:
    """Config-aware pipeline entry point for annotated SNP heatmaps."""
    generate_figures = bool(
        _get(
            config,
            "outputs",
            "generate_figures",
            default=True,
        )
    )

    if not generate_figures:
        return {}

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

    matrix_ids = {
        str(value)
        for value in snp_matrix.index
    }

    sample_id_column = _resolve_sample_id_column(
        annotation_df,
        matrix_ids,
        candidates=(
            configured_sample_id,
            "sample_id",
        ),
    )

    annotation_candidates = [
        (main_var, str(main_var) if main_var else None),
        (var_01, str(var_01) if var_01 else None),
        (sequence_type_column, "ST"),
        (core_cluster_column, "Core SNP cluster"),
        (
            transmission_cluster_column,
            "Transmission cluster",
        ),
    ]

    annotation_columns = []
    annotation_labels: dict[str, str] = {}

    for column, label in annotation_candidates:
        if (
            column
            and column in annotation_df.columns
            and column not in annotation_columns
        ):
            annotation_columns.append(
                str(column)
            )

            annotation_labels[str(column)] = (
                label or str(column)
            )

    if not annotation_columns:
        raise ValueError(
            "No configured heatmap annotation columns are available"
        )

    return plot_annotated_snp_heatmap(
        snp_matrix,
        annotation_df,
        output_dir,
        sample_id_column=sample_id_column,
        annotation_columns=annotation_columns,
        annotation_labels=annotation_labels,
        title=title,
        filename_prefix=filename_prefix,
    )
