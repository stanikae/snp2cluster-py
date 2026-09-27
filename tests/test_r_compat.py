import pandas as pd

from snp2cluster.r_compat import (
    cumsum_group_r_port,
    connect_samples_r_port,
    get_core_snp_clusters_r_compatible,
    matrix_to_mdf,
)


def test_cumsum_group_r_port_matches_pairwise_breaking():
    assert cumsum_group_r_port([3, 5, 25, 2, 4], 20) == [1, 1, 0, 2, 2]


def test_connect_samples_r_port_basic_threshold():
    df = pd.DataFrame(
        {
            "sampleID": ["A", "A", "B"],
            "ID2": ["B", "C", "C"],
            "X3": [5, 25, 7],
            "clst": [1, 0, 1],
            "num": [51, 250, 71],
            "rn": [1, 2, 3],
        }
    )
    out = connect_samples_r_port(df, threshold=20)
    retained = [x for x in out if x is not None]
    assert len(retained) == 2
    assert retained[0][0:2] == ["A", "B"]


def test_matrix_to_mdf_has_expected_pairs():
    mat = pd.DataFrame(
        [[0, 1, 10], [1, 0, 11], [10, 11, 0]],
        index=["A", "B", "C"],
        columns=["A", "B", "C"],
    )
    mdf = matrix_to_mdf(mat)
    assert list(mdf.columns) == ["X1", "X2", "X3"]
    assert len(mdf) == 3


def test_get_core_snp_clusters_r_compatible_runs():
    # Six samples in two obvious SNP groups.
    mat = pd.DataFrame(
        [
            [0, 2, 3, 50, 51, 52],
            [2, 0, 4, 49, 50, 51],
            [3, 4, 0, 48, 49, 50],
            [50, 49, 48, 0, 2, 3],
            [51, 50, 49, 2, 0, 4],
            [52, 51, 50, 3, 4, 0],
        ],
        index=["A", "B", "C", "D", "E", "F"],
        columns=["A", "B", "C", "D", "E", "F"],
    )
    result = get_core_snp_clusters_r_compatible(mat, k=2, max_k=4, snpco=20)
    assert result is not None
    assert not result.snp_clusters.empty
    assert {"name", "cluster", "km_cluster"}.issubset(result.snp_clusters.columns)
