import pandas as pd

from snp2cluster.r_compat import (
    cumsum_group_r_port,
    connect_samples_r_port,
    determine_optimal_k_r_style,
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


def test_identical_profiles_select_k1():
    """
    Multiple eligible isolates with identical SNP-distance profiles
    must produce a data-driven initial K=1 result.
    """
    ids = ["A", "B", "C", "D", "E"]

    mat = pd.DataFrame(
        0.0,
        index=ids,
        columns=ids,
    )

    sig_n, sig_nn, diagnostics = (
        determine_optimal_k_r_style(
            mat,
            random_seed=42,
        )
    )

    assert sig_n == 1
    assert sig_nn == 3

    assert diagnostics["k"].tolist() == [1]
    assert diagnostics["silhouette"].isna().all()

    assert (
        diagnostics["selection_method"].tolist()
        == ["single_distinct_profile"]
    )


def test_identical_five_isolates_form_one_core_cluster():
    """
    A valid K=1 solution must continue through SNP-cutoff refinement.

    Five isolates with all pairwise distances equal to zero must
    remain one final Core SNP cluster containing all five isolates.
    """
    ids = ["A", "B", "C", "D", "E"]

    mat = pd.DataFrame(
        0.0,
        index=ids,
        columns=ids,
    )

    result = get_core_snp_clusters_r_compatible(
        m=mat,
        k=1,
        max_k=3,
        snpco=20,
        random_seed=42,
    )

    assert result is not None

    assert set(
        result.snp_clusters["name"]
    ) == set(ids)

    assert len(
        result.snp_clusters
    ) == 5

    assert (
        result.snp_clusters["cluster"].nunique()
        == 1
    )

    assert (
        result.snp_clusters["km_cluster"].nunique()
        == 1
    )

    assert (
        result.kmeans_assignments[
            "km_cluster"
        ].nunique()
        == 1
    )

    assert (
        result.silhouette_widths[
            "sil_width"
        ].isna().all()
    )

    assert (
        result.cluster_avg_widths[
            "selection_method"
        ].tolist()
        == ["single_cluster_solution"]
    )

    assert (
        result.cluster_avg_widths[
            "cluster_avg_width"
        ].isna().all()
    )

    assert (
        result.cluster_avg_widths[
            "keep"
        ].tolist()
        == [True]
    )

    assert (
        result.cluster_avg_widths[
            "selection_method"
        ].tolist()
        == ["single_cluster_solution"]
    )

    assert result.kept_kmeans_clusters == ["1"]
    assert result.excluded_kmeans_clusters == []


def test_two_isolates_within_cutoff_form_valid_cluster():
    """
    A two-isolate context is eligible for Core clustering.

    If its pairwise SNP distance is within the configured cutoff,
    both isolates must be retained in one final Core SNP cluster.
    """
    mat = pd.DataFrame(
        [
            [0, 5],
            [5, 0],
        ],
        index=["A", "B"],
        columns=["A", "B"],
    )

    result = get_core_snp_clusters_r_compatible(
        m=mat,
        k=1,
        max_k=1,
        snpco=20,
        random_seed=42,
    )

    assert result is not None

    assert set(
        result.snp_clusters["name"]
    ) == {"A", "B"}

    assert len(
        result.snp_clusters
    ) == 2

    assert (
        result.snp_clusters["cluster"].nunique()
        == 1
    )

    cluster_sizes = (
        result.snp_clusters
        .groupby("cluster")["name"]
        .nunique()
    )

    assert cluster_sizes.min() == 2


def test_two_isolates_above_cutoff_form_no_valid_cluster():
    """
    Initial K=1 must not automatically make a final Core SNP cluster.

    If the only pairwise distance exceeds the configured SNP cutoff,
    the pair must not be reported as a valid Core SNP cluster.
    """
    mat = pd.DataFrame(
        [
            [0, 25],
            [25, 0],
        ],
        index=["A", "B"],
        columns=["A", "B"],
    )

    result = get_core_snp_clusters_r_compatible(
        m=mat,
        k=1,
        max_k=1,
        snpco=20,
        random_seed=42,
    )

    assert result is None

def test_single_isolate_context_returns_no_cluster():
    """
    A Core-clustering context containing fewer than two isolates is
    not eligible to produce a valid Core SNP cluster.
    """
    mat = pd.DataFrame(
        [[0]],
        index=["A"],
        columns=["A"],
    )

    result = get_core_snp_clusters_r_compatible(
        m=mat,
        k=1,
        max_k=1,
        snpco=20,
        random_seed=42,
    )

    assert result is None
