import pandas as pd

from snp2cluster.clustering import assign_kmeans_clusters, generate_snp_chains


def test_assign_kmeans_clusters_small_matrix():
    mat = pd.DataFrame(
        [[0, 1, 10], [1, 0, 11], [10, 11, 0]],
        index=["A", "B", "C"],
        columns=["A", "B", "C"],
    )
    out = assign_kmeans_clusters(mat, k=2, random_seed=42)
    assert set(out["sample_id"]) == {"A", "B", "C"}
    assert "kmeans_cluster" in out.columns


def test_generate_snp_chains_connected_component():
    mat = pd.DataFrame(
        [[0, 1, 10], [1, 0, 11], [10, 11, 0]],
        index=["A", "B", "C"],
        columns=["A", "B", "C"],
    )
    assignments = pd.DataFrame({"sample_id": ["A", "B", "C"], "kmeans_cluster": ["K1", "K1", "K1"]})
    out = generate_snp_chains(mat, assignments, snp_cutoff=2)
    assert out["snp_cluster_id"].nunique() == 2
