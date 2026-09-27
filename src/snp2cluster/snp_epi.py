from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

from .r_compat import (
    CoreSNPResult,
    get_core_snp_clusters_r_compatible,
    matrix_to_mdf,
    rleid,
)

from .date_utils import parse_dates

@dataclass(frozen=True)
class CoreSNPAnalysisResult:
    """Python equivalent of R `run_core_snp_cluster_analysis()` return list."""

    snp_clust: pd.DataFrame
    snp_clust_id: pd.DataFrame
    df2mat: pd.DataFrame | None
    vec_keep_names: list[str]
    mat_optimal_centres: pd.DataFrame | None
    core_result: CoreSNPResult


def _fit_kmeans_centers(m: pd.DataFrame, k: int, random_seed: int = 42) -> pd.DataFrame:
    """Fit KMeans and return centers with R-like row/column labels.

    In R, `res.km$centers` has rows as K-means clusters and columns as original
    matrix features/sample IDs. Because SNP2Cluster treats the SNP distance matrix
    directly as the feature matrix, sklearn's KMeans centers are a valid analogue.
    """
    x = m.astype(float).values
    model = KMeans(n_clusters=k, random_state=random_seed, n_init=25)
    model.fit(x)
    centers = pd.DataFrame(
        model.cluster_centers_,
        index=[str(i + 1) for i in range(k)],
        columns=m.columns.astype(str),
    )
    return centers


def run_core_snp_cluster_analysis_r_compatible(
    mat: pd.DataFrame,
    sigN: int,
    sigNN: int,
    snpco: float,
    datesJoin: pd.DataFrame | None = None,
    random_seed: int = 42,
) -> CoreSNPAnalysisResult | None:
    """Python port of R `run_core_snp_cluster_analysis()`.

    R equivalent:
    ```r
    res_list <- get_core_snp_clusters(m=mat, k=sigN, max=sigNN,
                                      dates=datesJoin, snpco=snpco, orig=T)
    snpClust <- res_list[[1]]
    ...
    core_snp_results_list[[1]] <- snpClust
    core_snp_results_list[[2]] <- snpClustID
    core_snp_results_list[[3]] <- df2mat
    core_snp_results_list[[4]] <- vec_keep_names
    core_snp_results_list[[5]] <- mat_optimal_centres
    ```

    Returns None where the R function would return NA.
    """
    mat = mat.copy()
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)
    mat = mat.loc[mat.index, mat.index].astype(float)

    core = get_core_snp_clusters_r_compatible(
        m=mat,
        k=sigN,
        max_k=sigNN,
        dates=datesJoin,
        snpco=snpco,
        orig=True,
        random_seed=random_seed,
    )

    if core is None:
        return None

    snp_clust = core.snp_clusters.copy()
    if snp_clust.empty:
        return None

    # R: snpClustID <- snpClust %>% select(2,3) %>% distinct()
    # snpClust columns are expected as: name, cluster, km_cluster
    snp_clust["km_cluster"] = snp_clust["km_cluster"].astype(str)
    snp_clust["cluster"] = snp_clust["cluster"].astype(str)
    snp_clust_id = snp_clust[["cluster", "km_cluster"]].drop_duplicates().reset_index(drop=True)

    vec_keep_names = snp_clust["name"].astype(str).drop_duplicates().tolist()

    # R: centers_mat <- res_list[[2]]$centers
    centers_mat = _fit_kmeans_centers(mat, k=sigN, random_seed=random_seed)
    mat_optimal_centres = centers_mat.loc[:, centers_mat.columns.isin(vec_keep_names)].copy()

    # R df2mat: centers rows joined to snpClustID by km_cluster, rownames become SNP cluster IDs.
    tmp = mat_optimal_centres.reset_index(names="rowID")
    tmp = tmp.merge(snp_clust_id, left_on="rowID", right_on="km_cluster", how="inner")

    if tmp.empty:
        df2mat = None
    else:
        df2mat = tmp.drop(columns=["rowID", "km_cluster"]).set_index("cluster")

    return CoreSNPAnalysisResult(
        snp_clust=snp_clust.reset_index(drop=True),
        snp_clust_id=snp_clust_id,
        df2mat=df2mat,
        vec_keep_names=vec_keep_names,
        mat_optimal_centres=mat_optimal_centres,
        core_result=core,
    )


def _make_mdf_bidirectional(mdf: pd.DataFrame) -> pd.DataFrame:
    """Ensure mdf can match both X1-X2 and X2-X1 joins.

    The R function joins exactly by `name = X1` and `ID2 = X2`. In practice, if
    `mdf` was generated from a full matrix it may already include both directions.
    If it was generated from combinations only, this helper prevents false breaks
    due solely to pair orientation.
    """
    required = {"X1", "X2", "X3"}
    missing = required - set(mdf.columns)
    if missing:
        raise ValueError(f"mdf is missing required columns: {sorted(missing)}")

    fwd = mdf[["X1", "X2", "X3"]].copy()
    rev = fwd.rename(columns={"X1": "X2", "X2": "X1"})[["X1", "X2", "X3"]]
    both = pd.concat([fwd, rev], ignore_index=True).drop_duplicates(subset=["X1", "X2"])
    both["X1"] = both["X1"].astype(str)
    both["X2"] = both["X2"].astype(str)
    return both


# def calculate_snp_epi_clusters_r_compatible(
#     snpClust: pd.DataFrame,
#     epiwkDF: pd.DataFrame,
#     mdf: pd.DataFrame,
#     excl_vec: list[str] | None = None,
#     daysco: int = 45,
#     bidirectional_mdf: bool = True,
# ) -> pd.DataFrame | None:
def calculate_snp_epi_clusters_r_compatible(
    snpClust: pd.DataFrame,
    epiwkDF: pd.DataFrame,
    mdf: pd.DataFrame,
    sample_id_col: str,
    date_col: str,
    st_col: str,
    date_format: str,
    excl_vec: list[str] | None = None,
    daysco: int = 45,
    bidirectional_mdf: bool = True,
) -> pd.DataFrame | None:
    """
    Python port of R calculate_SNP_Epi_clusters().

    Parameters
    ----------
    snpClust:
        Must contain:
            name
            cluster

    epiwkDF:
        Metadata dataframe containing:
            sample_id_col
            date_col
            st_col

    mdf:
        Long SNP matrix containing:
            X1
            X2
            X3

    sample_id_col:
        Sample identifier column name from config

    date_col:
        Date column from config (var_02)

    st_col:
        Sequence type column

    date_format:
        YAML date_format value
    """

    excl_vec = excl_vec or []

    required_snp = {"name", "cluster"}

    required_epi = {
        sample_id_col,
        date_col,
        st_col,
    }

    required_mdf = {
        "X1",
        "X2",
        "X3",
    }

    if missing := (required_snp - set(snpClust.columns)):
        raise ValueError(
            f"snpClust missing required columns: {sorted(missing)}"
        )

    if missing := (required_epi - set(epiwkDF.columns)):
        raise ValueError(
            f"epiwkDF missing required columns: {sorted(missing)}"
        )

    if missing := (required_mdf - set(mdf.columns)):
        raise ValueError(
            f"mdf missing required columns: {sorted(missing)}"
        )

    # ------------------------------------------------------------------
    # Prepare inputs
    # ------------------------------------------------------------------

    snp = snpClust.copy()

    epi = epiwkDF.copy()

    snp["name"] = snp["name"].astype(str)

    snp["cluster"] = snp["cluster"].astype(str)

    epi[sample_id_col] = epi[sample_id_col].astype(str)

    epi[date_col] = parse_dates(
        epi[date_col],
        date_format=date_format,
    )

    epi[st_col] = epi[st_col].astype(str)

    # ------------------------------------------------------------------
    # Prepare MDF
    # ------------------------------------------------------------------

    mdf_use = (
        _make_mdf_bidirectional(mdf)
        if bidirectional_mdf
        else mdf.copy()
    )

    mdf_use["X1"] = mdf_use["X1"].astype(str)

    mdf_use["X2"] = mdf_use["X2"].astype(str)

    # ------------------------------------------------------------------
    # Process each SNP cluster
    # ------------------------------------------------------------------

    episnp_frames: list[pd.DataFrame] = []

    vec_clust = (
        snp["cluster"]
        .drop_duplicates()
        .tolist()
    )

    for cluster_id in vec_clust:

        temp = (
            snp[snp["cluster"] == str(cluster_id)]
            .copy()
        )

        temp = temp.merge(
            epi,
            left_on="name",
            right_on=sample_id_col,
            how="inner",
        )

        if temp.empty:
            continue

        # --------------------------------------------------------------
        # Equivalent to:
        # group_by(ST)
        # arrange(TakenDate)
        # lead(name)
        # lead(TakenDate)
        # --------------------------------------------------------------

        temp = (
            temp
            .sort_values(
                [st_col, date_col]
            )
            .reset_index(drop=True)
        )

        temp["ID2"] = (
            temp
            .groupby(st_col)["name"]
            .shift(-1)
        )

        temp["Date2"] = (
            temp
            .groupby(st_col)[date_col]
            .shift(-1)
        )

        temp["Days"] = (
            temp["Date2"]
            - temp[date_col]
        ).dt.days

        temp["epicumsum"] = np.where(
            temp["Days"] <= daysco,
            1,
            0,
        )

        temp.loc[
            temp["Days"].isna(),
            "epicumsum"
        ] = 0

        # --------------------------------------------------------------
        # Join MDF
        # --------------------------------------------------------------

        joined = temp.merge(
            mdf_use,
            left_on=["name", "ID2"],
            right_on=["X1", "X2"],
            how="left",
        )

        joined.loc[
            joined["X3"].isna(),
            "epicumsum"
        ] = 0

        joined = (
            joined
            .sort_values(
                [st_col, date_col]
            )
            .reset_index(drop=True)
        )

        joined["CG"] = rleid(
            joined["epicumsum"].tolist()
        )

        joined = joined[
            joined["epicumsum"] != 0
        ].copy()

        if not joined.empty:
            episnp_frames.append(joined)

    # ------------------------------------------------------------------
    # No transmission clusters
    # ------------------------------------------------------------------

    if not episnp_frames:
        return None

    clusterSet3 = pd.concat(
        episnp_frames,
        ignore_index=True,
    )

    # ------------------------------------------------------------------
    # Create cluster IDs
    # ------------------------------------------------------------------

    clusterSet3["num"] = (
        clusterSet3["cluster"].astype(str)
        + "_ST"
        + clusterSet3[st_col].astype(str)
        + "_"
        + clusterSet3["CG"].astype(str)
    )

    clusterSet3["Clusters"] = rleid(
        clusterSet3["num"].tolist()
    )

    # ------------------------------------------------------------------
    # pivot_longer(name, ID2)
    # ------------------------------------------------------------------

    drop_if_present = [
        sample_id_col,
    ]

    for col in drop_if_present:
        if col in clusterSet3.columns:
            clusterSet3 = clusterSet3.drop(
                columns=[col]
            )

    value_cols = [
        "name",
        "ID2",
    ]

    id_cols = [
        c
        for c in clusterSet3.columns
        if c not in value_cols
    ]

    long_df = (
        clusterSet3
        .melt(
            id_vars=id_cols,
            value_vars=value_cols,
            value_name="sample_id",
        )
        .drop(columns=["variable"])
    )

    long_df = (
        long_df
        .dropna(subset=["sample_id"])
        .copy()
    )

    long_df["sample_id"] = (
        long_df["sample_id"]
        .astype(str)
    )

    # ------------------------------------------------------------------
    # Remove excluded columns
    # ------------------------------------------------------------------

    drop_cols = [
        c
        for c in excl_vec
        if c in long_df.columns
    ]

    if drop_cols:
        long_df = long_df.drop(
            columns=drop_cols
        )

    long_df = long_df.rename(
        columns={
            "X3": "SNPs",
        }
    )

    # ------------------------------------------------------------------
    # R:
    # group_by(Clusters)
    # distinct(sampleID)
    # ------------------------------------------------------------------

    long_df = long_df.drop_duplicates(
        subset=["Clusters", "sample_id"],
        keep="first",
    )

    counts = (
        long_df
        .groupby("Clusters")["sample_id"]
        .transform("count")
    )

    long_df["Cluster_Cases_count"] = counts

    long_df["Clusters"] = (
        long_df["Clusters"]
        .astype(str)
    )

    final = long_df[
        [
            "sample_id",
            "Days",
            "SNPs",
            "Clusters",
            "Cluster_Cases_count",
        ]
    ].copy()

    if final.empty:
        return None

    return final.reset_index(drop=True)
