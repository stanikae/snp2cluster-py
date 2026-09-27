from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# def plot_k_selection(
#     diagnostics_df: pd.DataFrame,
#     output_dir: str | Path,
#     title: str = "Optimal Number of Clusters",
# ) -> dict[str, Path]:
    
def plot_k_selection(
    diagnostics_df: pd.DataFrame,
    output_dir: str | Path,
    title: str = "Optimal Number of Clusters",
    filename_prefix: str = "01_k_selection",
) -> dict[str, Path]:
    """
    Plot silhouette-based K selection diagnostics.

    Parameters
    ----------
    diagnostics_df : pd.DataFrame

        Expected columns:

            k
            silhouette

    output_dir : str | Path

        Directory where figures will be saved.

    title : str

        Plot title.

    Returns
    -------
    dict

        Paths to generated figures.
    """

    required = {
        "k",
        "silhouette",
    }

    missing = required - set(
        diagnostics_df.columns
    )

    if missing:

        raise ValueError(
            f"Missing columns: {missing}"
        )

    output_dir = Path(output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    optimal_row = diagnostics_df.loc[
        diagnostics_df["silhouette"].idxmax()
    ]

    optimal_k = int(
        optimal_row["k"]
    )

    optimal_score = float(
        optimal_row["silhouette"]
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    if len(diagnostics_df) == 1:

        ax.scatter(
            diagnostics_df["k"],
            diagnostics_df["silhouette"],
            s=120,
            zorder=3,
        )

        ax.text(
            diagnostics_df["k"].iloc[0],
            diagnostics_df["silhouette"].iloc[0],
            "single k evaluated",
            ha="left",
            va="bottom",
        )

    else:

        ax.plot(
            diagnostics_df["k"],
            diagnostics_df["silhouette"],
            marker="o",
            linewidth=2,
        )

    # ax.plot(
    #     diagnostics_df["k"],
    #     diagnostics_df["silhouette"],
    #     marker="o",
    #     linewidth=2,
    # )

    ax.axvline(
        optimal_k,
        linestyle="--",
        linewidth=1.5,
    )

    # ax.scatter(
    #     optimal_k,
    #     optimal_score,
    #     s=120,
    #     zorder=3,
    # )

    
    y_min = diagnostics_df["silhouette"].min()
    y_max = diagnostics_df["silhouette"].max()
    y_range = y_max - y_min

    if y_range == 0:
        y_range = 0.05

    # If the optimal point is close to the top, place label below it.
    if optimal_score > (y_min + 0.8 * y_range):

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
    

    # ax.set_title(
    #     f"{title}: {optimal_k}"
    # )

    if len(diagnostics_df) == 1:

        ax.set_title(
            f"{title}: {optimal_k} "
            "(single k evaluated)"
        )

    else:

        # ax.set_title(
        #     f"{title}: {optimal_k}"
        # )
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

    safe_prefix = (
        filename_prefix
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
    )

    png_path = (
        output_dir
        / f"{safe_prefix}.png"
    )

    pdf_path = (
        output_dir
        / f"{safe_prefix}.pdf"
    )


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