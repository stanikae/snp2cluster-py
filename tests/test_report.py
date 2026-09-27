from __future__ import annotations

from pathlib import Path

import pandas as pd

from snp2cluster.report import write_html_report


def test_write_html_report_discovers_outputs_and_escapes_html(tmp_path: Path) -> None:
    figures = tmp_path / "figures"
    tables = tmp_path / "tables"
    report_dir = tmp_path / "report"
    figures.mkdir()
    tables.mkdir()

    (figures / "03_core_heatmap.png").write_bytes(b"png")
    (figures / "04_mst_Hospital_A.html").write_text("<html></html>", encoding="utf-8")
    (figures / "05_transmission_network_Hospital_A.html").write_text(
        "<html></html>",
        encoding="utf-8",
    )
    (tables / "core_snp_clusters.csv").write_text("sample_id\nA\n", encoding="utf-8")

    assignments = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C"],
            "core_cluster_label": ["ST152__1", "ST152__1", pd.NA],
        }
    )
    transmission = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C"],
            "Clusters": [1, 1, pd.NA],
        }
    )
    config = {
        "analysis": {
            "cluster_type": "Transmission",
            "transmission_level": "Facility",
            "core_cluster_context": ["ST"],
            "snp_cutoff": 20,
            "days_cutoff": 45,
        },
        "variables": {
            "main_var": "FacilityName",
            "var_01": "WardType",
            "var_02": "TakenDate",
        },
        "columns": {
            "sample_id": "sample_id",
            "sequence_type": "ST",
        },
    }

    output = write_html_report(
        report_dir / "snp2cluster_report.html",
        project_name="Demo <Project>",
        assignments=assignments,
        transmission=transmission,
        config=config,
        figures_dir=figures,
        tables_dir=tables,
        software_version="1.2.3",
    )

    document = output.read_text(encoding="utf-8")

    assert "Demo &lt;Project&gt;" in document
    assert "Core SNP clusters" in document
    assert "Transmission clusters" in document
    assert "../figures/03_core_heatmap.png" in document
    assert "../figures/04_mst_Hospital_A.html" in document
    assert "../figures/05_transmission_network_Hospital_A.html" in document
    assert "../tables/core_snp_clusters.csv" in document
    assert "1.2.3" in document


def test_write_html_report_supports_core_only_mode(tmp_path: Path) -> None:
    assignments = pd.DataFrame(
        {
            "SampleID": ["A", "B"],
            "snp_cluster_id": [1, 1],
        }
    )
    config = {
        "analysis": {"cluster_type": "Core"},
        "columns": {"sample_id": "SampleID"},
    }

    output = write_html_report(
        tmp_path / "report" / "snp2cluster_report.html",
        project_name="Core run",
        assignments=assignments,
        transmission=None,
        config=config,
        figures_dir=tmp_path / "figures",
        tables_dir=tmp_path / "tables",
    )

    document = output.read_text(encoding="utf-8")

    assert "Core run" in document
    assert "No interactive outputs were generated for this run." in document
    assert "No output tables were found." in document
