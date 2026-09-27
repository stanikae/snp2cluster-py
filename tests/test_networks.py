from __future__ import annotations

import networkx as nx
import pandas as pd
import pytest

from snp2cluster.visualization.networks import (
    build_transmission_cluster_network,
    generate_transmission_networks_from_config,
    prepare_transmission_network_data,
)


def _matrix() -> pd.DataFrame:
    ids = ["A", "B", "C", "D"]
    return pd.DataFrame(
        [
            [0, 1, 8, 9],
            [1, 0, 7, 8],
            [8, 7, 0, 2],
            [9, 8, 2, 0],
        ],
        index=ids,
        columns=ids,
    )


def _transmission_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "SampleID": ["A", "B", "C", "D"],
            "FacilityName": ["H1", "H1", "H1", "H1"],
            "WardType": ["NICU", "NICU", "Other", "Other"],
            "TakenDate": [
                "2026-01-01",
                "2026-01-03",
                "2026-02-01",
                "2026-02-04",
            ],
            "ST": [152, 152, 39, 39],
            "Clusters": [1, 1, 2, 2],
            "_core_context": ["152", "152", "39", "39"],
            "core_cluster_label": ["152__1", "152__1", "39__1", "39__1"],
        }
    )


def test_prepare_excludes_unclustered_isolates() -> None:
    data = _transmission_data()
    data.loc[len(data)] = [
        "E",
        "H1",
        "Other",
        "2026-03-01",
        25,
        pd.NA,
        "25",
        "25__1",
    ]
    matrix = _matrix()

    prepared = prepare_transmission_network_data(
        data,
        matrix,
        sample_id_column="SampleID",
        main_var="FacilityName",
        collection_date_column="TakenDate",
        shape_column="WardType",
    )

    assert set(prepared.index) == {"A", "B", "C", "D"}


def test_cluster_specific_msts_have_expected_edges() -> None:
    prepared = prepare_transmission_network_data(
        _transmission_data(),
        _matrix(),
        sample_id_column="SampleID",
        main_var="FacilityName",
        collection_date_column="TakenDate",
        shape_column="WardType",
    )

    graph, edges = build_transmission_cluster_network(
        prepared,
        _matrix(),
        sample_id_column="SampleID",
    )

    assert isinstance(graph, nx.Graph)
    assert graph.number_of_nodes() == 4
    assert graph.number_of_edges() == 2
    assert len(edges) == 2
    assert set(edges["Clusters"].astype(str)) == {"1", "2"}


def test_edge_day_difference_is_calculated() -> None:
    prepared = prepare_transmission_network_data(
        _transmission_data(),
        _matrix(),
        sample_id_column="SampleID",
        main_var="FacilityName",
        collection_date_column="TakenDate",
        shape_column="WardType",
    )

    _, edges = build_transmission_cluster_network(
        prepared,
        _matrix(),
        sample_id_column="SampleID",
    )

    assert sorted(edges["day_difference"].tolist()) == [2, 3]


def test_conflicting_duplicate_cluster_assignments_fail() -> None:
    data = _transmission_data()
    duplicate = data.iloc[[0]].copy()
    duplicate["Clusters"] = 9
    data = pd.concat([data, duplicate], ignore_index=True)

    with pytest.raises(ValueError, match="conflicting"):
        prepare_transmission_network_data(
            data,
            _matrix(),
            sample_id_column="SampleID",
            main_var="FacilityName",
            collection_date_column="TakenDate",
            shape_column="WardType",
        )


def test_core_mode_returns_no_network_outputs(tmp_path) -> None:
    config = {
        "analysis": {"cluster_type": "Core"},
        "outputs": {"generate_figures": True},
    }

    outputs = generate_transmission_networks_from_config(
        _transmission_data(),
        _matrix(),
        config,
        tmp_path,
    )

    assert outputs["paths"] == {}
    assert outputs["nodes"].empty
    assert outputs["edges"].empty
