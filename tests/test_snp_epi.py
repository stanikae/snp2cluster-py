import pandas as pd

from snp2cluster.r_compat import matrix_to_mdf
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