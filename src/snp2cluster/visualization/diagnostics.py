from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd



def plot_k_selection(
    diagnostics_df: pd.DataFrame,
    output_dir: str | Path,
    title: str = "Optimal Number of Clusters",
    filename_prefix: str = "01_k_selection",
) -> dict[str, Path]:
    """
    Plot K-selection diagnostics.

    For ordinary multi-cluster solutions, the figure shows candidate
    K values and their average silhouette widths.

    Where data-driven K selection returns K=1, a numerical silhouette
    coefficient is not defined. In that case, the figure records the
    selected K and the reason for the single-cluster solution without
    assigning an artificial silhouette value.

    Parameters
    ----------
    diagnostics_df
        Required columns are ``k`` and ``silhouette``. The optional
        ``selection_method`` column records how K was selected.

    output_dir
        Directory where PNG and PDF figures will be saved.

    title
        Plot title.

    filename_prefix
        Prefix used for the PNG and PDF filenames.

    Returns
    -------
    dict[str, Path]
        Paths to the generated PNG and PDF figures.
    """
    required = {
        "k",
        "silhouette",
    }

    missing = required - set(diagnostics_df.columns)

    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}"
        )

    if diagnostics_df.empty:
        raise ValueError(
            "K-selection diagnostics are empty."
        )

    diagnostics = diagnostics_df.copy()

    diagnostics["k"] = pd.to_numeric(
        diagnostics["k"],
        errors="coerce",
    )

    diagnostics["silhouette"] = pd.to_numeric(
        diagnostics["silhouette"],
        errors="coerce",
    )

    if diagnostics["k"].isna().all():
        raise ValueError(
            "K-selection diagnostics contain no valid K values."
        )

    diagnostics = (
        diagnostics
        .dropna(subset=["k"])
        .sort_values("k")
        .reset_index(drop=True)
    )

    diagnostics["k"] = diagnostics["k"].astype(int)

    output_dir = Path(output_dir)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    safe_prefix = (
        filename_prefix
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
    )

    png_path = output_dir / f"{safe_prefix}.png"
    pdf_path = output_dir / f"{safe_prefix}.pdf"

    valid_diagnostics = diagnostics.dropna(
        subset=["silhouette"]
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    if valid_diagnostics.empty:
        # A K=1 result is a valid data-driven boundary condition.
        # A silhouette coefficient is not defined because no second
        # cluster exists for comparison.
        selected_k = int(
            diagnostics["k"].iloc[0]
        )

        if "selection_method" in diagnostics.columns:
            raw_method = diagnostics[
                "selection_method"
            ].iloc[0]

            selection_method = (
                str(raw_method)
                if pd.notna(raw_method)
                else "single_cluster_solution"
            )
        else:
            selection_method = (
                "single_cluster_solution"
            )

        method_labels = {
            "single_distinct_profile": (
                "One distinct SNP-distance profile"
            ),
            "no_valid_multicluster_solution": (
                "No valid multi-cluster solution"
            ),
            "single_candidate_cluster": (
                "One supported initial cluster"
            ),
            "single_cluster_solution": (
                "One supported initial cluster"
            ),
            "insufficient_isolates": (
                "Insufficient isolates for clustering"
            ),
        }

        reason = method_labels.get(
            selection_method,
            selection_method.replace("_", " "),
        )

        ax.text(
            0.5,
            0.60,
            f"Selected k = {selected_k}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
        )

        ax.text(
            0.5,
            0.43,
            reason,
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )

        ax.text(
            0.5,
            0.31,
            (
                "No numerical silhouette width is defined "
                "for a one-cluster solution."
            ),
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10,
        )

        ax.set_title(
            f"{title}: selected k = {selected_k}"
        )

        ax.set_axis_off()

    else:
        optimal_idx = (
            valid_diagnostics["silhouette"]
            .idxmax()
        )

        optimal_row = diagnostics.loc[
            optimal_idx
        ]

        optimal_k = int(
            optimal_row["k"]
        )

        optimal_score = float(
            optimal_row["silhouette"]
        )

        if len(valid_diagnostics) == 1:
            ax.scatter(
                valid_diagnostics["k"],
                valid_diagnostics["silhouette"],
                s=120,
                zorder=3,
            )
        else:
            ax.plot(
                valid_diagnostics["k"],
                valid_diagnostics["silhouette"],
                marker="o",
                linewidth=2,
            )

        ax.axvline(
            optimal_k,
            linestyle="--",
            linewidth=1.5,
        )

        y_min = float(
            valid_diagnostics["silhouette"].min()
        )

        y_max = float(
            valid_diagnostics["silhouette"].max()
        )

        y_range = y_max - y_min

        if y_range == 0:
            y_range = 0.05

        if optimal_score > (
            y_min + 0.8 * y_range
        ):
            annotation_offset = (8, -32)
            vertical_alignment = "top"
        else:
            annotation_offset = (8, 8)
            vertical_alignment = "bottom"

        ax.annotate(
            f"k={optimal_k}\n{optimal_score:.2f}",
            xy=(
                optimal_k,
                optimal_score,
            ),
            xytext=annotation_offset,
            textcoords="offset points",
            ha="left",
            va=vertical_alignment,
        )

        if len(valid_diagnostics) == 1:
            ax.set_title(
                f"{title}: optimal k = {optimal_k} "
                "(single valid k evaluated)"
            )
        else:
            ax.set_title(
                f"{title}: optimal k = {optimal_k}"
            )

        ax.set_xlabel(
            "Number of clusters k"
        )

        ax.set_ylabel(
            "Average silhouette width"
        )

        ax.grid(
            alpha=0.3
        )

        padding = y_range * 0.15

        ax.set_ylim(
            y_min - padding,
            y_max + padding,
        )

    fig.tight_layout()

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    plt.close(fig)

    return {
        "png": png_path,
        "pdf": pdf_path,
    }