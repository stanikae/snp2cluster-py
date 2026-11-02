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

    run_manifest = {
        "tool": "SNP2Cluster-Py",
        "version": "1.2.3",
        "project_name": "Demo Project",
        "run_directory": (
            "runs/demo_project_20261101T210835Z"
        ),
        "timestamp_utc": (
            "2026-11-01T21:08:35+00:00"
        ),
        "python_version": "3.13.9",
        "platform": "macOS-arm64",
        "input_files": {
            "metadata": "data/metadata.csv",
            "snp_matrix": "data/snp_matrix.csv",
            "mlst_profile": "data/mlst.xlsx",
        },
        "columns": {
            "sample_id": "sample_id",
            "sequence_type": "ST",
        },
        "analysis_parameters": {
            "cluster_type": "Transmission",
            "core_cluster_context": ["ST"],
            "snp_cutoff": 20,
            "days_cutoff": 45,
            "random_seed": 42,
            "engine": "r_compatible",
        },
        "validation_messages": [],
        "dependency_versions": {
            "pandas": "3.0.6",
            "numpy": "2.5.3",
            "sklearn": "1.9.1",
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
        run_manifest=run_manifest,
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
    assert "No run-manifest metadata was supplied." not in document
    assert "Run Directory" in document
    assert "runs/demo_project_20261101T210835Z" in document
    assert "Python Version" in document
    assert "3.13.9" in document
    assert "Platform" in document
    assert "macOS-arm64" in document
    assert "Run identity" in document
    assert "Input files" in document
    assert "Analysis parameters" in document
    assert "Configured columns" in document
    assert "Validation" in document
    assert "Software dependencies" in document
    assert "Software Version" in document
    assert "data/metadata.csv" in document
    assert "data/snp_matrix.csv" in document
    assert "data/mlst.xlsx" in document
    assert "Core Cluster Context" in document
    assert "ST" in document
    assert "Snp Cutoff" in document
    assert "20" in document
    assert "Days Cutoff" in document
    assert "45" in document
    assert "Random Seed" in document
    assert "42" in document
    assert "r_compatible" in document
    assert (
        "Validation passed with no reported issues."
        in document
    )
    assert "Pandas" in document
    assert "3.0.6" in document

    assert "Numpy" in document
    assert "2.5.3" in document

    assert "Sklearn" in document
    assert "1.9.1" in document

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
    assert "No run-manifest metadata was supplied." in document


def test_write_html_report_renders_validation_messages(tmp_path: Path,) -> None:
    assignments = pd.DataFrame(
        {
            "sample_id": ["A", "B"],
            "core_cluster_label": [
                "ST1__1",
                "ST1__1",
            ],
        }
    )

    manifest = {
        "tool": "SNP2Cluster-Py",
        "version": "1.2.3",
        "validation_messages": [
            "Example validation notice.",
        ],
    }

    output = write_html_report(
        (
            tmp_path
            / "report"
            / "snp2cluster_report.html"
        ),
        project_name="Validation run",
        assignments=assignments,
        transmission=None,
        config={
            "analysis": {
                "cluster_type": "Core",
            },
            "columns": {
                "sample_id": "sample_id",
            },
        },
        figures_dir=tmp_path / "figures",
        tables_dir=tmp_path / "tables",
        software_version="1.2.3",
        run_manifest=manifest,
    )

    document = output.read_text(
        encoding="utf-8"
    )

    assert "Example validation notice." in document

    assert (
        "Validation passed with no reported issues."
        not in document
    )
