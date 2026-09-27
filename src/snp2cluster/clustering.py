from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def choose_optimal_k(
    distance_matrix: pd.DataFrame,
    max_k: int | None = None,
    random_seed: int = 42,
) -> tuple[int, pd.DataFrame]:
    """
    R-parity version of K selection.

    Mirrors:

        kmx <- round(nrow(mat)/2)

        fviz_nbclust(
            mat,
            kmeans,
            method = "silhouette"
        )

        sigN <- cluster count with highest silhouette

    Returns:
        best_k
        diagnostics dataframe
    """

    x = distance_matrix.astype(float).values

    n = x.shape[0]

    if n < 3:
        return 1, pd.DataFrame(
            {
                "k": [1],
                "silhouette": [np.nan],
            }
        )

    #
    # R:
    # kmx <- round(nrow(mat)/2)
    #
    r_kmax = max(2, round(n / 2))

    if max_k is not None:
        r_kmax = min(r_kmax, max_k)

    rows = []

    best_k = 2
    best_score = -1.0

    for k in range(2, r_kmax + 1):

        model = KMeans(
            n_clusters=k,
            random_state=random_seed,
            n_init=25,  # R parity
        )

        labels = model.fit_predict(x)

        if len(set(labels)) < 2:

            score = np.nan

        else:

            try:
                score = silhouette_score(
                    x,
                    labels,
                    metric="euclidean",
                )
            except Exception:
                score = np.nan

        rows.append(
            {
                "k": k,
                "silhouette": score,
            }
        )

        if not np.isnan(score) and score > best_score:

            best_score = score
            best_k = k

    diagnostics = (
        pd.DataFrame(rows)
        .sort_values("k")
        .reset_index(drop=True)
    )

    return best_k, diagnostics


def assign_kmeans_clusters(
    distance_matrix: pd.DataFrame,
    k: int,
    random_seed: int = 42,
) -> pd.DataFrame:
    x = distance_matrix.astype(float).values
    if k <= 1:
        labels = np.zeros(x.shape[0], dtype=int)
    else:
        model = KMeans(n_clusters=k, random_state=random_seed, n_init="auto")
        labels = model.fit_predict(x)

    return pd.DataFrame(
        {
            "sample_id": distance_matrix.index.astype(str),
            "kmeans_cluster": [f"K{label + 1}" for label in labels],
        }
    )


def generate_snp_chains(
    distance_matrix: pd.DataFrame,
    assignments: pd.DataFrame,
    snp_cutoff: int,
) -> pd.DataFrame:
    """Generate preliminary SNP cluster chains.

    NOTE: This is a conservative placeholder.
    The exact R cumulative SNP-chain algorithm should be ported after reviewing the R function body.
    Current behaviour: within each KMeans group, connected components are formed where pairwise SNP distance
    is <= snp_cutoff.
    """
    import networkx as nx

    df = assignments.copy()
    cluster_rows = []

    for k_group, sub in df.groupby("kmeans_cluster"):
        ids = sub["sample_id"].astype(str).tolist()
        graph = nx.Graph()
        graph.add_nodes_from(ids)
        for i, a in enumerate(ids):
            for b in ids[i + 1 :]:
                d = float(distance_matrix.loc[a, b])
                if d <= snp_cutoff:
                    graph.add_edge(a, b, snp_distance=d)

        for comp_i, component in enumerate(nx.connected_components(graph), start=1):
            snp_cluster = f"{k_group}_SNP{comp_i}"
            for sample_id in sorted(component):
                cluster_rows.append(
                    {
                        "sample_id": sample_id,
                        "kmeans_cluster": k_group,
                        "snp_cluster_id": snp_cluster,
                        "snp_cutoff": snp_cutoff,
                    }
                )

    return pd.DataFrame(cluster_rows)
