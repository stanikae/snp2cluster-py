import pandas as pd

from snp2cluster.r_compat import (
    get_core_snp_clusters_r_compatible,
    matrix_to_mdf,
)
from snp2cluster.snp_epi import calculate_snp_epi_clusters_r_compatible
from snp2cluster.date_utils import parse_dates


def test_calculate_snp_epi_clusters_r_compatible_basic_chain():
    snp_clust = pd.DataFrame(
        {
            "name": ["A", "B", "C", "D"],
            "cluster": ["1", "1", "1", "1"],
            "km_cluster": ["1", "1", "1", "1"],
        }
    )
    epi = pd.DataFrame(
        {
            "sampleID": ["A", "B", "C", "D"],
            "TakenDate": ["2024-01-01", "2024-01-05", "2024-03-01", "2024-03-03"],
            "ST": ["15", "15", "20", "20"],
        }
    )
    mat = pd.DataFrame(
        [
            [0, 3, 50, 51],
            [3, 0, 49, 50],
            [50, 49, 0, 4],
            [51, 50, 4, 0],
        ],
        index=["A", "B", "C", "D"],
        columns=["A", "B", "C", "D"],
    )
    mdf = matrix_to_mdf(mat)

    out = calculate_snp_epi_clusters_r_compatible(
        snpClust=snp_clust,
        epiwkDF=epi,
        mdf=mdf,
        sample_id_col="sampleID",
        date_col="TakenDate",
        st_col="ST",
        date_format="ymd",
        excl_vec=[],
        daysco=45,
    )

    assert out is not None
    # assert set(out["sampleID"]) == {"A", "B", "C", "D"}
    assert set(out["sample_id"]) == {"A", "B", "C", "D"}
    assert out["Clusters"].nunique() == 2
    assert set(out["Cluster_Cases_count"]) == {2}


def test_calculate_snp_epi_clusters_breaks_by_days_cutoff():
    snp_clust = pd.DataFrame(
        {
            "name": ["A", "B"],
            "cluster": ["1", "1"],
            "km_cluster": ["1", "1"],
        }
    )
    epi = pd.DataFrame(
        {
            "sampleID": ["A", "B"],
            "TakenDate": ["2024-01-01", "2024-04-01"],
            "ST": ["15", "15"],
        }
    )
    mdf = pd.DataFrame({"X1": ["A"], "X2": ["B"], "X3": [3]})


    out = calculate_snp_epi_clusters_r_compatible(
        snpClust=snp_clust,
        epiwkDF=epi,
        mdf=mdf,
        sample_id_col="sampleID",
        date_col="TakenDate",
        st_col="ST",
        date_format="ymd",
        excl_vec=[],
        daysco=45,
    )

    assert out is None


def test_parse_dates_ymd_slash():

    s = pd.Series(
        ["2019/12/17"]
    )

    out = parse_dates(
        s,
        date_format="ymd",
    )

    assert str(out.iloc[0].date()) == "2019-12-17"


def test_parse_dates_dmy_slash():

    s = pd.Series(
        ["17/12/2019"]
    )

    out = parse_dates(
        s,
        date_format="dmy",
    )

    assert str(out.iloc[0].date()) == "2019-12-17"

def test_st307_core_cluster_produces_four_isolate_transmission_chain():
    """
    Six Hospital C ST307 isolates form one Core SNP cluster at an
    11-SNP cutoff, but only the four temporally connected isolates form
    a transmission cluster at a 14-day cutoff.

    BG32 is separated from BG315 by 48 days. BG317 is separated from
    BG254 by 42 days. The middle four isolates have consecutive gaps
    of 7, 1 and 5 days.
    """
    ids = [
        "BG258",
        "BG254",
        "BG32",
        "BG317",
        "BG315",
        "BG259",
    ]

    matrix = pd.DataFrame(
        [
            [0, 5, 4, 5, 2, 2],
            [5, 0, 1, 6, 3, 3],
            [4, 1, 0, 5, 2, 2],
            [5, 6, 5, 0, 7, 3],
            [2, 3, 2, 7, 0, 4],
            [2, 3, 2, 3, 4, 0],
        ],
        index=ids,
        columns=ids,
        dtype=float,
    )

    metadata = pd.DataFrame(
        {
            "sampleID": ids,
            "TakenDate": [
                "12/5/19",
                "12/10/19",
                "10/10/19",
                "1/21/20",
                "11/27/19",
                "12/4/19",
            ],
            "ST": ["307"] * 6,
        }
    )

    core_result = get_core_snp_clusters_r_compatible(
        m=matrix,
        k=2,
        max_k=4,
        snpco=11,
        random_seed=42,
    )

    assert core_result is not None

    assert set(
        core_result.snp_clusters["name"]
    ) == set(ids)

    assert (
        core_result.snp_clusters[
            "cluster"
        ].nunique()
        == 1
    )

    transmission_result = (
        calculate_snp_epi_clusters_r_compatible(
            snpClust=core_result.snp_clusters,
            epiwkDF=metadata,
            mdf=matrix_to_mdf(matrix),
            sample_id_col="sampleID",
            date_col="TakenDate",
            st_col="ST",
            date_format="mdy",
            excl_vec=[],
            daysco=14,
        )
    )

    assert transmission_result is not None

    expected_transmission_ids = {
        "BG254",
        "BG258",
        "BG259",
        "BG315",
    }

    observed_transmission_ids = set(
        transmission_result[
            "sample_id"
        ].astype(str)
    )

    assert (
        observed_transmission_ids
        == expected_transmission_ids
    )

    assert (
        transmission_result[
            "Clusters"
        ].nunique()
        == 1
    )

    assert (
        transmission_result[
            "Cluster_Cases_count"
        ]
        .astype(int)
        .eq(4)
        .all()
    )

    assert "BG32" not in observed_transmission_ids
    assert "BG317" not in observed_transmission_ids


def test_st307_transmission_membership_is_input_order_invariant():
    """
    Reordering metadata and SNP-matrix rows must not change ST307 Core
    or transmission-cluster membership.
    """
    original_ids = [
        "BG258",
        "BG254",
        "BG32",
        "BG317",
        "BG315",
        "BG259",
    ]

    matrix = pd.DataFrame(
        [
            [0, 5, 4, 5, 2, 2],
            [5, 0, 1, 6, 3, 3],
            [4, 1, 0, 5, 2, 2],
            [5, 6, 5, 0, 7, 3],
            [2, 3, 2, 7, 0, 4],
            [2, 3, 2, 3, 4, 0],
        ],
        index=original_ids,
        columns=original_ids,
        dtype=float,
    )

    dates_by_sample = {
        "BG32": "10/10/19",
        "BG315": "11/27/19",
        "BG259": "12/4/19",
        "BG258": "12/5/19",
        "BG254": "12/10/19",
        "BG317": "1/21/20",
    }

    expected_core_ids = set(original_ids)

    expected_transmission_ids = {
        "BG254",
        "BG258",
        "BG259",
        "BG315",
    }

    orders = [
        original_ids,
        [
            "BG317",
            "BG32",
            "BG259",
            "BG254",
            "BG315",
            "BG258",
        ],
    ]

    for sample_order in orders:
        ordered_matrix = matrix.loc[
            sample_order,
            sample_order,
        ]

        metadata = pd.DataFrame(
            {
                "sampleID": sample_order,
                "TakenDate": [
                    dates_by_sample[sample_id]
                    for sample_id in sample_order
                ],
                "ST": ["307"] * len(sample_order),
            }
        )

        core_result = (
            get_core_snp_clusters_r_compatible(
                m=ordered_matrix,
                k=2,
                max_k=4,
                snpco=11,
                random_seed=42,
            )
        )

        assert core_result is not None

        assert set(
            core_result.snp_clusters["name"]
        ) == expected_core_ids

        assert (
            core_result.snp_clusters[
                "cluster"
            ].nunique()
            == 1
        )

        transmission_result = (
            calculate_snp_epi_clusters_r_compatible(
                snpClust=core_result.snp_clusters,
                epiwkDF=metadata,
                mdf=matrix_to_mdf(
                    ordered_matrix
                ),
                sample_id_col="sampleID",
                date_col="TakenDate",
                st_col="ST",
                date_format="mdy",
                excl_vec=[],
                daysco=14,
            )
        )

        assert transmission_result is not None

        assert set(
            transmission_result[
                "sample_id"
            ].astype(str)
        ) == expected_transmission_ids

        assert (
            transmission_result[
                "Clusters"
            ].nunique()
            == 1
        )

        assert (
            transmission_result[
                "Cluster_Cases_count"
            ]
            .astype(int)
            .eq(4)
            .all()
        )
