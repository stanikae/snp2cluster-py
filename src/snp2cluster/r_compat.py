from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score


@dataclass(frozen=True)
class CoreSNPResult:
    """R-compatible core SNP clustering result."""

    snp_clusters: pd.DataFrame
    kmeans_assignments: pd.DataFrame
    silhouette_widths: pd.DataFrame
    cluster_avg_widths: pd.DataFrame
    kept_kmeans_clusters: list[str]
    excluded_kmeans_clusters: list[str]


# def determine_optimal_k_r_style(
#     mat: pd.DataFrame,
# ) -> tuple[int, int]:
#     """
#     Port of R fviz_nbclust silhouette search.
#     """

#     n = len(mat)

#     kmx = round(n / 2)

#     if kmx < 2:
#         return 2, 4

#     best_k = 2
#     best_score = -1

#     for k in range(2, kmx + 1):

#         model = KMeans(
#             n_clusters=k,
#             n_init=25,
#             random_state=42,
#         )

#         labels = model.fit_predict(mat)

#         if len(set(labels)) < 2:
#             continue

#         score = silhouette_score(mat, labels)

#         if score > best_score:
#             best_score = score
#             best_k = k

#     sigN = best_k
#     sigNN = sigN + 2

#     return sigN, sigNN



def determine_optimal_k_r_style(
    mat: pd.DataFrame,
    random_seed: int = 42,
) -> tuple[int, int, pd.DataFrame]:
    """
    Port of R fviz_nbclust silhouette search.

    Returns:
        sigN
            Optimal number of clusters.

        sigNN
            Upper clustering threshold used downstream.

        diagnostics
            DataFrame with k and silhouette values.
    """

    n = len(mat)

    kmx = round(n / 2)

    if kmx < 2:

        diagnostics = pd.DataFrame(
            {
                "k": [1],
                "silhouette": [float("nan")],
            }
        )

        return 2, 4, diagnostics

    best_k = 2
    best_score = -1

    rows = []

    for k in range(2, kmx + 1):

        model = KMeans(
            n_clusters=k,
            n_init=25,
            random_state=random_seed,
        )

        labels = model.fit_predict(mat)

        if len(set(labels)) < 2:

            score = float("nan")

        else:

            try:

                score = silhouette_score(
                    mat,
                    labels,
                )

            except Exception:

                score = float("nan")

        rows.append(
            {
                "k": k,
                "silhouette": score,
            }
        )

        if (
            not pd.isna(score)
            and score > best_score
        ):

            best_score = score
            best_k = k

    diagnostics = (
        pd.DataFrame(rows)
        .sort_values("k")
        .reset_index(drop=True)
    )

    sigN = best_k
    sigNN = sigN + 2

    return sigN, sigNN, diagnostics




def rleid(values) -> list[int]:
    """Small Python equivalent of data.table::rleid."""
    out: list[int] = []
    current = object()
    run_id = 0
    for value in values:
        if value != current:
            run_id += 1
            current = value
        out.append(run_id)
    return out


def cumsum_group_r_port(values, threshold: float) -> list[int]:
    """Python port of SNP2Cluster R helper `cumsum_group`.

    Important: although the helper is named cumsum_group, your R comment shows that
    the current algorithm no longer accumulates distances continuously. It applies
    the SNP cut-off pair-wise to avoid breaking valid cluster chains.
    """
    x = [-555 if pd.isna(v) else v for v in values]
    grp = 0
    group = 1
    result: list[int] = []

    for i, value in enumerate(x):
        current = abs(float(value))
        cumsum_value = current  # R: cumsum <- abs(x[i])

        if abs(cumsum_value) > threshold:
            previous = abs(float(x[i - 1])) if i > 0 else 0.0

            if (current + previous) <= threshold:
                group = group + 1
                cumsum_value = value
            elif current <= threshold:
                group = group + 1
                cumsum_value = value
            else:
                grp = group
                group = 0
                cumsum_value = 0

        result.append(int(group))

        if group == 0 and i != 0:
            group = grp + 1
        elif group == 0 and i == 0:
            group = 1

    return result


def connect_samples_r_port(
    df: pd.DataFrame,
    threshold: float = 20,
    rn_col: str = "rn",
    sample_col: str = "sampleID",
    id2_col: str = "ID2",
    snp_col: str = "X3",
) -> list[list[Any] | None]:
    """Python port of SNP2Cluster R helper `connect_samples`.

    Expected columns after R renaming:
    sampleID, ID2, X3, clust, num, rn
    """
    required = {rn_col, sample_col, id2_col, snp_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns for connect_samples_r_port: {sorted(missing)}")

    x = df.reset_index(drop=True).copy()
    vec_id1 = x[sample_col].astype(str).tolist()
    vec_id2 = x[id2_col].astype(str).tolist()
    vec_x3 = x[snp_col].tolist()
    rn_idx = list(range(len(x)))

    list_checks: list[list[Any] | None] = []
    pos_idx = 0
    grp_cnt = 0

    for row_pos in rn_idx:
        i = pos_idx if pos_idx != 0 else row_pos
        raw = vec_x3[i]
        snp_val = float(raw) if not pd.isna(raw) else np.inf
        vec_id3: list[Any] | None = None

        if int(snp_val) <= threshold:
            grp_cnt += 1
            id1 = vec_id1[i]
            val2 = vec_id2[i]
            vec_in_pair = [id1, val2]

            if i != 0:
                valid_checks = [v for v in list_checks if v is not None and len(v) >= 2]

                duplicate_pair = False
                for checked in valid_checks:
                    id_pair_vec = [checked[0], checked[1]]
                    if all(v in id_pair_vec for v in vec_in_pair):
                        duplicate_pair = True
                        break

                if duplicate_pair:
                    vec_id3 = None
                else:
                    vec_col1_vals = list(dict.fromkeys([v[0] for v in valid_checks]))
                    vec_col2_vals = list(dict.fromkeys([v[1] for v in valid_checks]))

                    if id1 in vec_col1_vals:
                        id_a, id_b = id1, val2
                    elif val2 in vec_col2_vals:
                        id_a, id_b = val2, id1
                    else:
                        id_a, id_b = id1, val2

                    vec_id3 = [id_a, id_b, snp_val, grp_cnt]
                    pos_idx = 0
            else:
                vec_id3 = [vec_in_pair[0], vec_in_pair[1], snp_val, grp_cnt]
                pos_idx = 0
        else:
            grp_cnt += 2
            vec_id3 = None

        list_checks.append(vec_id3)
        grp_cnt = 0

    return list_checks


def list_checks_to_df(list_checks: list[list[Any] | None]) -> pd.DataFrame:
    """Convert R-like list output into V1/V2/V3/V4 dataframe."""
    rows = []
    for item in list_checks:
        if item is None:
            continue
        if len(item) < 4:
            continue
        rows.append({"V1": item[0], "V2": item[1], "V3": float(item[2]), "V4": str(item[3])})
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def matrix_to_mdf(distance_matrix: pd.DataFrame) -> pd.DataFrame:
    """Create R-style long SNP matrix table with columns X1, X2, X3.

    The R algorithm joins all combinations of IDs against a long `mdf` table.
    This function creates one unordered pair per isolate pair using the matrix index order.
    """
    mat = distance_matrix.copy()
    mat.index = mat.index.astype(str)
    mat.columns = mat.columns.astype(str)
    mat = mat.loc[mat.index, mat.index].astype(float)

    rows = []
    ids = list(mat.index)
    for x1, x2 in combinations(ids, 2):
        rows.append({"X1": x1, "X2": x2, "X3": float(mat.loc[x1, x2])})
    return pd.DataFrame(rows)


def _fit_kmeans_and_silhouette(m: pd.DataFrame, k: int, random_seed: int = 42):
    x = m.astype(float).values
    labels = KMeans(n_clusters=k, random_state=random_seed, n_init=25).fit_predict(x)
    clusters = labels + 1

    if len(set(labels)) > 1 and len(labels) > k:
        widths = silhouette_samples(x, labels)
    else:
        widths = np.zeros(len(labels))

    widths_df = pd.DataFrame(
        {
            "ID": m.index.astype(str),
            "cluster": clusters.astype(str),
            "sil_width": widths,
        }
    )
    avg_df = (
        widths_df.groupby("cluster", as_index=False)["sil_width"]
        .mean()
        .rename(columns={"sil_width": "cluster_avg_width"})
    )
    return labels, widths_df, avg_df


def _prepare_keep_df1(compare_df: pd.DataFrame, mdf: pd.DataFrame, snpco: float) -> pd.DataFrame:
    check_df1 = compare_df.merge(mdf, on=["X1", "X2"], how="inner")
    check_df1 = check_df1[(check_df1["X1"] != check_df1["X2"]) & (check_df1["X3"] <= snpco)]
    if check_df1.empty:
        return pd.DataFrame()

    keep_df1 = check_df1.sort_values("X3").copy()
    keep_df1["snpcumsum"] = cumsum_group_r_port(keep_df1["X3"].tolist(), snpco)
    keep_df1 = keep_df1[keep_df1["snpcumsum"] != 0].copy()
    if keep_df1.empty:
        return keep_df1
    keep_df1["num"] = (keep_df1["X3"].astype(str) + keep_df1["snpcumsum"].astype(str)).astype(float)
    return keep_df1


def _select_non_overlapping_groups(df_filty: pd.DataFrame) -> pd.DataFrame:
    """Port of the list_dfs/nm_vec section that avoids assigning samples twice."""
    grps_vec = list(pd.unique(df_filty["snpcumsum3"]))
    list_dfs: dict[int, pd.DataFrame] = {}
    nm_vec: set[str] = set()

    for z_idx, z in enumerate(grps_vec, start=1):
        if z_idx == 1:
            k_values = [z]
            grp_df = df_filty[df_filty["snpcumsum3"].isin(k_values)].copy()
            list_dfs[int(z)] = grp_df
            nm_vec = set(pd.unique(pd.concat([grp_df["V1"], grp_df["V2"]]).astype(str)))

        if z_idx >= len(grps_vec):
            break

        j = grps_vec[z_idx]
        keep = df_filty[df_filty["snpcumsum3"] == j].copy()
        keep = keep[~keep["V1"].astype(str).isin(nm_vec)]
        keep = keep[~keep["V2"].astype(str).isin(nm_vec)]

        if not keep.empty:
            list_dfs[int(j)] = keep

        if list_dfs:
            bnd_df = pd.concat(list_dfs.values(), ignore_index=True)
            nm_vec = set(pd.unique(pd.concat([bnd_df["V1"], bnd_df["V2"]]).astype(str)))

    if not list_dfs:
        return pd.DataFrame()
    return pd.concat(list_dfs.values(), ignore_index=True)


def _clusters_from_pairs(df_pairs: pd.DataFrame, max_v: int = 0) -> pd.DataFrame:
    if df_pairs.empty:
        return pd.DataFrame(columns=["sampleID", "SNPs", "clst", "Clusters"])

    df_pairs = df_pairs.copy()
    df_pairs["Clusters"] = [x + max_v for x in rleid(df_pairs["snpcumsum3"].tolist())]
    df_pairs = df_pairs.rename(columns={"V1": "sampleID", "V2": "ID2", "V3": "SNPs", "snpcumsum3": "clst"})
    df_pairs["sampleID"] = df_pairs["sampleID"].astype(str)
    df_pairs["ID2"] = df_pairs["ID2"].astype(str)
    df_pairs["SNPs"] = pd.to_numeric(df_pairs["SNPs"])
    df_pairs["clst"] = pd.to_numeric(df_pairs["clst"])

    long_df = df_pairs.melt(
    id_vars=["SNPs","clst","Clusters"],
    value_vars=["sampleID","ID2"],
    value_name="sampleID_value",
    )

    long_df = long_df.rename(columns={"sampleID_value": "sampleID"})
    long_df = long_df.drop(columns=["variable"]).drop_duplicates(subset=["Clusters", "sampleID"])
    return long_df.reset_index(drop=True)


def get_core_snp_clusters_r_compatible(
    m: pd.DataFrame,
    k: int,
    max_k: int,
    snpco: float,
    dates: pd.DataFrame | None = None,
    orig: bool = True,
    random_seed: int = 42,
    silhouette_keep_threshold: float = 0.5,
) -> CoreSNPResult | None:
    """Python port of your R `get_core_snp_clusters` function, core/orig branch.

    This implements the main R logic you shared:
    - K-means with nstart=25 equivalent
    - calculate sample and average cluster silhouette widths
    - keep only K-means clusters with abs(avg silhouette) >= 0.5
    - for each kept K-means cluster, build all pairwise comparisons
    - filter pairs to SNP <= snpco
    - apply `cumsum_group`
    - apply `connect_samples`
    - remove overlapping sample assignments
    - return isolate-level SNP clusters joined to K-means cluster IDs

    The non-orig/date-ordered branch is intentionally not implemented here yet;
    share the remaining SNP-Epi functions and we can port that as the next layer.
    """
    if not orig:
        raise NotImplementedError("The date-ordered orig=FALSE branch will be added with SNP-Epi porting.")

    m = m.copy()
    m.index = m.index.astype(str)
    m.columns = m.columns.astype(str)
    # m = m.loc[m.index, m.index].astype(float)

    m = m.loc[m.index, m.index].astype(float)

    # scaled_values = StandardScaler().fit_transform(m)

    # m = pd.DataFrame(
    #     scaled_values,
    #     index=m.index,
    #     columns=m.columns,
    # )

    if m.shape[0] <= max_k:
        max_k = m.shape[0] - 1
    if k >= m.shape[0]:
        k = max(1, m.shape[0] - 1)
    if k < 2:
        return None

    labels, widths_df, avg_df = _fit_kmeans_and_silhouette(m, k=k, random_seed=random_seed)
    kmeans_assignments = widths_df[["ID", "cluster"]].rename(columns={"ID": "name", "cluster": "km_cluster"})

    avg_df["keep"] = avg_df["cluster_avg_width"].abs() >= silhouette_keep_threshold
    vec_keep = avg_df.loc[avg_df["keep"], "cluster"].astype(str).tolist()
    vec_excl = avg_df.loc[~avg_df["keep"], "cluster"].astype(str).tolist()

    if not vec_keep:
        return None

    mdf = matrix_to_mdf(m)
    list_clusters: list[pd.DataFrame] = []
    max_v_global = 0

    for grp in vec_keep:
        names_vec = widths_df.loc[widths_df["cluster"].astype(str) == str(grp), "ID"].astype(str).tolist()
        if len(names_vec) <= 2:
            continue

        clster_ids_df = widths_df.loc[widths_df["cluster"].astype(str) == str(grp), ["ID", "cluster"]].copy()
        clster_ids_df = clster_ids_df.rename(columns={"cluster": "km_cluster"})
        clster_ids_df["km_cluster"] = clster_ids_df["km_cluster"].astype(str)

        compare_df = pd.DataFrame(list(combinations(names_vec, 2)), columns=["X1", "X2"])
        keep_df1 = _prepare_keep_df1(compare_df, mdf, snpco)
        if keep_df1.empty:
            continue

    
        # if str(grp) == "3":
        #     print("\n===== KMEANS CLUSTER 3 =====")
        #     print(sorted(names_vec))

        #     print("\n===== VALID SNP PAIRS <= SNP CUTOFF =====")
        #     print(
        #         keep_df1[
        #             ["X1", "X2", "X3"]
        #         ]
        #         .sort_values("X3")
        #         .to_string(index=False)
        #     )

        keep_df2 = keep_df1.rename(columns={"snpcumsum": "clst", "X3": "SNPs"}).copy()
        df_filtx = keep_df2.sort_values("SNPs").reset_index(drop=True).copy()
        df_filtx["rn"] = np.arange(1, len(df_filtx) + 1)
        df_filtx = df_filtx.rename(columns={"X1": "sampleID", "X2": "ID2", "SNPs": "X3"})
        df_filtx = df_filtx[["sampleID", "ID2", "X3", "clst", "num", "rn"]]

        list_checks = connect_samples_r_port(df_filtx, threshold=snpco)
        df_filty = list_checks_to_df(list_checks)
        if df_filty.empty:
            continue

        df_filty["snpcumsum8"] = [str(x) for x in cumsum_group_r_port(df_filty["V3"].astype(float).tolist(), snpco)]
        df_filty["clstGRP"] = df_filty["V4"].astype(str) + df_filty["snpcumsum8"].astype(str)
        df_filty["snpcumsum3"] = rleid(df_filty["clstGRP"].tolist())

        selected_pairs = _select_non_overlapping_groups(df_filty)
        if selected_pairs.empty:
            continue

        keep_df3 = _clusters_from_pairs(selected_pairs, max_v=max_v_global)
        if keep_df3.empty:
            continue


        # if str(grp) == "3":
        #     print("\n===== FINAL SNP CLUSTER MEMBERS =====")
        #     print(
        #         keep_df3.to_string(index=False)
        #     )
        
        if "Clusters" in keep_df3.columns:
            keep_df = keep_df3
        else:
            # Defensive fallback; normally not reached after _clusters_from_pairs.
            keep_df = keep_df3.copy()
            keep_df["Clusters"] = rleid(keep_df["clst"].tolist())
            keep_df = keep_df.groupby("Clusters").filter(lambda x: len(x) > 1)

        if keep_df.empty:
            continue

        max_v_global = int(pd.to_numeric(keep_df["Clusters"]).max())

        out = (
            keep_df[["sampleID", "Clusters"]]
            .rename(columns={"sampleID": "name", "Clusters": "cluster"})
            .drop_duplicates()
            .merge(clster_ids_df, left_on="name", right_on="ID", how="inner")
            .drop(columns=["ID"])
        )
        out["cluster"] = out["cluster"].astype(str)
        out["km_cluster"] = out["km_cluster"].astype(str)
        list_clusters.append(out)

    if not list_clusters:
        return None

    snp_clust = pd.concat(list_clusters, ignore_index=True).drop_duplicates()
    if snp_clust.empty:
        return None

    return CoreSNPResult(
        snp_clusters=snp_clust,
        kmeans_assignments=kmeans_assignments,
        silhouette_widths=widths_df,
        cluster_avg_widths=avg_df,
        kept_kmeans_clusters=vec_keep,
        excluded_kmeans_clusters=vec_excl,
    )
