from __future__ import annotations

from collections.abc import Sequence

import networkx as nx
import numpy as np
import pandas as pd


_DEFAULT_REFERENCE_LABELS = ("reference",)


def prepare_mst_distance_matrix(
    distance_matrix: pd.DataFrame,
    *,
    reference_labels: Sequence[str] = _DEFAULT_REFERENCE_LABELS,
) -> pd.DataFrame:
    """Validate and prepare a pairwise SNP-distance matrix for an MST."""
    if not isinstance(distance_matrix, pd.DataFrame):
        raise TypeError("distance_matrix must be a pandas DataFrame")
    if distance_matrix.empty:
        raise ValueError("distance_matrix is empty")
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise ValueError("distance_matrix must be square")

    matrix = distance_matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)

    if set(matrix.index) != set(matrix.columns):
        raise ValueError("Distance-matrix row and column identifiers do not match")

    matrix = matrix.loc[matrix.index, matrix.index]
    try:
        matrix = matrix.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("distance_matrix contains non-numeric values") from exc

    if not np.allclose(matrix.to_numpy(), matrix.to_numpy().T, equal_nan=True):
        raise ValueError("distance_matrix must be symmetric")
    if not np.allclose(np.diag(matrix.to_numpy()), 0.0):
        raise ValueError("distance_matrix diagonal must contain zeros")

    reference_set = {str(value).casefold() for value in reference_labels}
    retained_ids = [
        sample_id
        for sample_id in matrix.index
        if sample_id.casefold() not in reference_set
    ]
    matrix = matrix.loc[retained_ids, retained_ids]
    if matrix.empty:
        raise ValueError("No isolates remain after removing reference labels")
    return matrix


def build_minimum_spanning_tree(
    distance_matrix: pd.DataFrame,
    *,
    reference_labels: Sequence[str] = _DEFAULT_REFERENCE_LABELS,
    algorithm: str = "kruskal",
) -> nx.Graph:
    """Build a weighted minimum spanning tree from pairwise SNP distances."""
    matrix = prepare_mst_distance_matrix(
        distance_matrix,
        reference_labels=reference_labels,
    )
    graph = nx.Graph()
    sample_ids = list(matrix.index)
    graph.add_nodes_from(sample_ids)

    for index, sample_a in enumerate(sample_ids):
        for sample_b in sample_ids[index + 1 :]:
            snp_distance = float(matrix.loc[sample_a, sample_b])
            if not np.isfinite(snp_distance):
                raise ValueError(
                    f"Non-finite SNP distance between {sample_a} and {sample_b}"
                )
            graph.add_edge(sample_a, sample_b, weight=snp_distance)

    return nx.minimum_spanning_tree(
        graph,
        weight="weight",
        algorithm=algorithm,
    )




def mst_edges_table(tree: nx.Graph) -> pd.DataFrame:
    """Convert weighted MST edges to a tabular representation."""
    rows = [
        {
            "sample_id_1": sample_a,
            "sample_id_2": sample_b,
            "snp_distance": edge_data.get("weight"),
        }
        for sample_a, sample_b, edge_data in tree.edges(data=True)
    ]
    return pd.DataFrame(
        rows,
        columns=["sample_id_1", "sample_id_2", "snp_distance"],
    )
