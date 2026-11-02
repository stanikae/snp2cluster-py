from __future__ import annotations
import base64
import html
import os
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


_MISSING_VALUES = {"", "na", "nan", "none", "null", "missing"}


def _get(config: Any, *keys: str, default: Any = None) -> Any:
    """Read nested values from mappings or attribute-based config objects."""
    value = config

    for key in keys:
        if value is None:
            return default

        if isinstance(value, Mapping):
            if key not in value:
                return default
            value = value[key]
        else:
            if not hasattr(value, key):
                return default
            value = getattr(value, key)

    return value


def _normalise(value: object, *, default: str = "Not available") -> str:
    """Return a safe human-readable scalar value."""
    if value is None or pd.isna(value):
        return default

    text = str(value).strip()
    return default if text.casefold() in _MISSING_VALUES else text


def _assigned_mask(series: pd.Series) -> pd.Series:
    """Return a Boolean mask for non-missing assignments."""
    values = series.astype("string").str.strip()
    return values.notna() & ~values.str.casefold().isin(_MISSING_VALUES)


def _count_unique_assigned(df: pd.DataFrame | None, candidates: Sequence[str]) -> int | None:
    """Count non-missing unique values in the first available candidate column."""
    if df is None or df.empty:
        return None

    for column in candidates:
        if column in df.columns:
            mask = _assigned_mask(df[column])
            return int(df.loc[mask, column].astype(str).nunique())

    return None


def _count_assigned_rows(df: pd.DataFrame | None, candidates: Sequence[str]) -> int | None:
    """Count rows with an assignment in the first available candidate column."""
    if df is None or df.empty:
        return None

    for column in candidates:
        if column in df.columns:
            return int(_assigned_mask(df[column]).sum())

    return None


def _resolve_sample_id_column(
    assignments: pd.DataFrame,
    transmission: pd.DataFrame | None,
    config: Any,
) -> str | None:
    """Resolve the configured or conventional sample-ID field."""
    configured = _get(config, "columns", "sample_id")
    candidates = [configured, "sample_id", "SampleID", "isolate_id"]

    for column in candidates:
        if not column:
            continue
        if column in assignments.columns:
            return str(column)
        if transmission is not None and column in transmission.columns:
            return str(column)

    return None


def _relative_url(target: Path, report_dir: Path) -> str:
    """Return a browser-compatible relative URL."""
    return Path(os.path.relpath(target, report_dir)).as_posix()


def _discover_files(directory: Path | None, patterns: Sequence[str]) -> list[Path]:
    """Discover unique files matching ordered glob patterns."""
    if directory is None or not directory.exists():
        return []

    found: list[Path] = []
    seen: set[Path] = set()

    for pattern in patterns:
        for path in sorted(directory.glob(pattern)):
            if path.is_file() and path not in seen:
                found.append(path)
                seen.add(path)

    return found


def _display_label(path: Path) -> str:
    """Create a concise, readable label from an output filename."""
    stem = path.stem

    if stem.startswith("03_core_heatmap"):
        return "Annotated Core SNP heatmap"

    if stem.startswith("02_transmission_scatter_"):
        context = stem.removeprefix(
            "02_transmission_scatter_"
        )

        context_label = (
            context
            .replace("__", " / ")
            .replace("_", " ")
            .strip()
        )

        return (
            f"Transmission scatterplot: {context_label}"
            if context_label
            else "Transmission scatterplot"
        )

    if stem.startswith("02_transmission_scatter"):
        return "Transmission scatterplot"

    if stem.startswith("01_k_selection_"):
        context = stem.removeprefix(
            "01_k_selection_"
        )

        context_label = (
            context
            .replace("__", " / ")
            .replace("_", " ")
            .strip()
        )

        return (
            f"K-selection diagnostic: {context_label}"
            if context_label
            else "K-selection diagnostic"
        )

    if stem.startswith("01_k_selection"):
        return "K-selection diagnostic"

    if stem.startswith("04_mst_"):
        context = stem.removeprefix("04_mst_")

        context_label = (
            context
            .replace("__", " / ")
            .replace("_", " ")
            .strip()
        )

        return (
            f"Core SNP MST: {context_label}"
            if context_label
            else "Core SNP MST"
        )

    if stem.startswith("05_transmission_network_"):
        context = stem.removeprefix(
            "05_transmission_network_"
        )

        context_label = (
            context
            .replace("__", " / ")
            .replace("_", " ")
            .strip()
        )

        return (
            f"Transmission network: {context_label}"
            if context_label
            else "Transmission network"
        )

    return (
        stem
        .replace("__", " / ")
        .replace("_", " ")
        .strip()
        or path.stem
    )


def _metric_card(
    label: str,
    value: object,
    detail: str | None = None,
) -> str:
    """
    Render a summary metric card with optional supporting detail.
    """
    detail_html = (
        (
            '<span class="metric-detail">'
            f"{html.escape(detail)}"
            "</span>"
        )
        if detail
        else ""
    )

    return (
        '<div class="metric">'
        f'<span class="metric-label">{html.escape(label)}</span>'
        f'<strong class="metric-value">{html.escape(str(value))}</strong>'
        f"{detail_html}"
        "</div>"
    )


def _config_rows(config: Any) -> list[tuple[str, str]]:
    """Build the open scientific configuration summary."""
    core_context = _get(config, "analysis", "core_cluster_context")

    if isinstance(core_context, (list, tuple)):
        core_context_display = " → ".join(str(value) for value in core_context)
    else:
        core_context_display = _normalise(core_context)

    values = [
        ("Analysis mode", _get(config, "analysis", "cluster_type")),
        ("Transmission level", _get(config, "analysis", "transmission_level")),
        ("Core cluster context", core_context_display),
        ("SNP cutoff", _get(config, "analysis", "snp_cutoff")),
        ("Temporal cutoff (days)", _get(config, "analysis", "days_cutoff")),
        ("Bootstraps", _get(config, "analysis", "n_bootstraps")),
        ("Random seed", _get(config, "analysis", "random_seed")),
        ("Date format", _get(config, "analysis", "date_format")),
        ("Clustering engine", _get(config, "analysis", "engine")),
        ("Maximum k", _get(config, "analysis", "max_k")),
        ("Main variable", _get(config, "variables", "main_var")),
        ("Secondary variable", _get(config, "variables", "var_01")),
        ("Collection-date variable", _get(config, "variables", "var_02")),
        ("Sample-ID column", _get(config, "columns", "sample_id")),
        ("Sequence-type column", _get(config, "columns", "sequence_type")),
    ]

    return [
        (label, _normalise(value))
        for label, value in values
    ]


def _render_table(rows: Sequence[tuple[str, str]]) -> str:
    body = "".join(
        "<tr>"
        f"<th>{html.escape(label)}</th>"
        f"<td>{html.escape(value)}</td>"
        "</tr>"
        for label, value in rows
    )
    return f'<div class="table-wrap"><table>{body}</table></div>'


def _render_static_gallery(files: Sequence[Path], report_dir: Path) -> str:
    if not files:
        return '<p class="empty-state">No static figures were generated for this run.</p>'

    figures = []

    for path in files:
        url = html.escape(_relative_url(path, report_dir), quote=True)
        label = html.escape(_display_label(path))
        figures.append(
            '<figure class="figure-card">'
            f'<a href="{url}" target="_blank" rel="noopener">'
            f'<img src="{url}" alt="{label}">'
            "</a>"
            f"<figcaption>{label}</figcaption>"
            "</figure>"
        )

    return '<div class="figure-grid">' + "".join(figures) + "</div>"


def _render_link_cards(
    files: Sequence[Path],
    report_dir: Path,
    *,
    button_text: str,
) -> str:
    if not files:
        return '<p class="empty-state">No interactive outputs were generated for this run.</p>'

    cards = []

    for path in files:
        url = html.escape(_relative_url(path, report_dir), quote=True)
        label = html.escape(_display_label(path))
        cards.append(
            '<article class="link-card">'
            f"<h3>{label}</h3>"
            f'<a class="button" href="{url}" target="_blank" rel="noopener">'
            f"{html.escape(button_text)}</a>"
            "</article>"
        )

    return '<div class="link-grid">' + "".join(cards) + "</div>"


def _output_file_description(
    path: Path,
) -> str:
    """
    Return a concise description for a generated output table.
    """
    descriptions = {
        "core_snp_clusters.csv": (
            "Retained context-aware Core SNP clusters and their "
            "member isolates."
        ),
        "isolate_cluster_assignments.csv": (
            "Isolate-level Core SNP and transmission-cluster "
            "assignments for the run."
        ),
        "k_selection_diagnostics.csv": (
            "K-means cluster-level diagnostic results, including "
            "average silhouette widths, retention status and "
            "selection method."
        ),
        "kmeans_cluster_assignments.csv": (
            "Initial K-means cluster assignments generated before "
            "SNP-cutoff refinement."
        ),
        "minimum_spanning_tree_edges.csv": (
            "Core SNP minimum-spanning-tree edges used in the "
            "interactive MST outputs."
        ),
        "transmission_clusters.csv": (
            "Inferred transmission-cluster assignments produced "
            "after genomic, epidemiological and temporal analysis."
        ),
        "transmission_context.csv": (
            "Isolate-level epidemiological and temporal context used "
            "during transmission-cluster analysis."
        ),
        "transmission_network_edges.csv": (
            "Edges used to visualise genetic relationships within "
            "inferred transmission clusters."
        ),
        "transmission_network_nodes.csv": (
            "Node attributes used in the interactive transmission "
            "network outputs."
        ),
    }

    return descriptions.get(
        path.name,
        "Machine-readable output generated for this analysis run.",
    )


def _render_downloads(
    paths: Sequence[Path],
    report_dir: Path,
) -> str:
    """
    Render downloadable output files with concise descriptions.
    """
    if not paths:
        return (
            '<p class="empty-state">'
            "No output tables were found."
            "</p>"
        )

    cards: list[str] = []

    for path in paths:
        relative_path = os.path.relpath(
            path,
            start=report_dir,
        )

        href = html.escape(
            relative_path,
            quote=True,
        )

        label = html.escape(
            _display_label(path)
        )

        description = html.escape(
            _output_file_description(path)
        )

        cards.append(
            '<article class="download-card">'
            '<div class="download-card-content">'
            f"<h3>{label}</h3>"
            f'<p class="download-description">{description}</p>'
            "</div>"
            f'<a class="button" href="{href}" download>'
            "Download"
            "</a>"
            "</article>"
        )

    return (
        '<div class="download-list">'
        + "".join(cards)
        + "</div>"
    )


def _humanize_key(value: object) -> str:
    """
    Convert a machine-readable manifest key into a report label.
    """
    return (
        str(value)
        .strip()
        .replace("_", " ")
        .title()
    )


def _format_manifest_value(value: object) -> str:
    """
    Format and HTML-escape a manifest value.
    """
    if value is None:
        return "Not specified"

    if isinstance(value, bool):
        return "Yes" if value else "No"

    if isinstance(value, Sequence) and not isinstance(
        value,
        (str, bytes),
    ):
        if not value:
            return "None"

        return ", ".join(
            html.escape(str(item))
            for item in value
        )

    text = str(value).strip()

    if not text:
        return "Not specified"

    return html.escape(text)


def _render_key_value_table(
    values: Mapping[str, object] | None,
    *,
    empty_message: str,
) -> str:
    """
    Render a mapping as a two-column HTML table.
    """
    if not values:
        return (
            '<p class="empty-state">'
            f"{html.escape(empty_message)}"
            "</p>"
        )

    rows: list[str] = []

    for key, value in values.items():
        rows.append(
            "<tr>"
            f"<th>{html.escape(_humanize_key(key))}</th>"
            f"<td>{_format_manifest_value(value)}</td>"
            "</tr>"
        )

    return (
        '<div class="table-wrap">'
        "<table>"
        + "".join(rows)
        + "</table>"
        "</div>"
    )

def _render_manifest(
    run_manifest: Mapping[str, object] | None,
    software_version: str | None,
) -> str:
    """
    Render structured run provenance from the run manifest.
    """
    if not run_manifest:
        return (
            '<p class="empty-state">'
            "No run-manifest metadata was supplied."
            "</p>"
        )

    identity_keys = (
        "tool",
        "project_name",
        "run_directory",
        "timestamp_utc",
        "python_version",
        "platform",
    )

    identity: dict[str, object] = {}

    manifest_version = run_manifest.get("version")

    if software_version:
        identity["software_version"] = software_version
    elif manifest_version:
        identity["software_version"] = manifest_version

    for key in identity_keys:
        if key in run_manifest:
            identity[key] = run_manifest[key]

    input_files = run_manifest.get("input_files")
    configured_columns = run_manifest.get("columns")
    analysis_parameters = run_manifest.get(
        "analysis_parameters"
    )
    dependency_versions = run_manifest.get(
        "dependency_versions"
    )
    validation_messages = run_manifest.get(
        "validation_messages"
    )

    if (
        isinstance(validation_messages, Sequence)
        and not isinstance(
            validation_messages,
            (str, bytes),
        )
    ):
        validation_items = list(validation_messages)
    elif validation_messages:
        validation_items = [validation_messages]
    else:
        validation_items = []

    if validation_items:
        validation_html = (
            '<ul class="validation-list">'
            + "".join(
                f"<li>{html.escape(str(item))}</li>"
                for item in validation_items
            )
            + "</ul>"
        )
    else:
        validation_html = (
            '<p class="validation-success">'
            "Validation passed with no reported issues."
            "</p>"
        )

    return f"""
    <div class="provenance-block">
      <h3>Run identity</h3>
      {_render_key_value_table(
          identity,
          empty_message=(
              "Run identity metadata is unavailable."
          ),
      )}

      <details class="provenance-details" open>
        <summary>Input files</summary>
        {_render_key_value_table(
            input_files
            if isinstance(input_files, Mapping)
            else None,
            empty_message=(
                "Input-file metadata is unavailable."
            ),
        )}
      </details>

      <details class="provenance-details">
        <summary>Analysis parameters</summary>
        {_render_key_value_table(
            analysis_parameters
            if isinstance(
                analysis_parameters,
                Mapping,
            )
            else None,
            empty_message=(
                "Analysis-parameter metadata is unavailable."
            ),
        )}
      </details>

      <details class="provenance-details">
        <summary>Configured columns</summary>
        {_render_key_value_table(
            configured_columns
            if isinstance(
                configured_columns,
                Mapping,
            )
            else None,
            empty_message=(
                "Configured-column metadata is unavailable."
            ),
        )}
      </details>

      <details class="provenance-details">
        <summary>Validation</summary>
        {validation_html}
      </details>

      <details class="provenance-details">
        <summary>Software dependencies</summary>
        {_render_key_value_table(
            dependency_versions
            if isinstance(
                dependency_versions,
                Mapping,
            )
            else None,
            empty_message=(
                "Dependency-version metadata is unavailable."
            ),
        )}
      </details>
    </div>
    """


def _load_logo_data_uri(
    logo_path: str | Path | None = None,
) -> str | None:
    """
    Load a local logo file and return an embedded Base64 data URI.

    If no explicit path is supplied, use the packaged SNP2Cluster-Py
    logo. Returning None allows the report to use a text fallback.
    """
    if logo_path is None:
        path = (
            Path(__file__).resolve().parent
            / "assets"
            / "snp2cluster-py.png"
        )
    else:
        path = Path(logo_path)

    if not path.is_file():
        return None

    mime_types = {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".svg": "image/svg+xml",
        ".webp": "image/webp",
    }

    mime_type = mime_types.get(
        path.suffix.lower()
    )

    if mime_type is None:
        return None

    try:
        encoded = base64.b64encode(
            path.read_bytes()
        ).decode("ascii")
    except OSError:
        return None

    return f"data:{mime_type};base64,{encoded}"


def _format_percentage(
    numerator: int,
    denominator: int,
) -> str:
    """
    Format a count as a percentage of the supplied denominator.
    """
    if denominator <= 0:
        return "0.0%"

    percentage = (
        float(numerator)
        / float(denominator)
        * 100.0
    )

    return f"{percentage:.1f}%"


def write_html_report(
    output_path: str | Path,
    *,
    project_name: str,
    assignments: pd.DataFrame,
    transmission: pd.DataFrame | None = None,
    config: Any = None,
    figures_dir: str | Path | None = None,
    tables_dir: str | Path | None = None,
    run_manifest: Mapping[str, object] | None = None,
    software_version: str | None = None,
    logo_path: str | Path | None = None,
) -> Path:
    """Write the open scientific SNP2Cluster-Py HTML run report."""
    logo_data_uri = _load_logo_data_uri(
        logo_path
    )

    if logo_data_uri:
        logo_html = (
            '<div class="hero-logo-panel">'
            '<img class="hero-logo" '
            f'src="{logo_data_uri}" '
            'alt="SNP2Cluster-Py logo">'
            "</div>"
        )
    else:
        logo_html = (
            '<div class="hero-logo-fallback">'
            "SNP2Cluster-Py"
            "</div>"
        )
    if not isinstance(assignments, pd.DataFrame):
        raise TypeError("assignments must be a pandas DataFrame")
    if transmission is not None and not isinstance(transmission, pd.DataFrame):
        raise TypeError("transmission must be a pandas DataFrame or None")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir = output_path.parent
    figures_path = Path(figures_dir) if figures_dir is not None else None
    tables_path = Path(tables_dir) if tables_dir is not None else None

    sample_id_column = _resolve_sample_id_column(assignments, transmission, config)
    sample_source = assignments

    if sample_id_column and sample_id_column not in assignments.columns:
        sample_source = transmission if transmission is not None else assignments

    if sample_id_column and sample_id_column in sample_source.columns:
        n_samples = int(sample_source[sample_id_column].astype(str).nunique())
    else:
        n_samples = int(len(assignments))

    n_core_clusters = _count_unique_assigned(
        assignments,
        ("core_cluster_label", "snp_cluster_id"),
    )
    n_core_clustered = _count_assigned_rows(
        assignments,
        ("core_cluster_label", "snp_cluster_id"),
    )
    n_transmission_clusters = _count_unique_assigned(
        transmission,
        ("Clusters", "transmission_cluster_id"),
    )
    n_transmission_clustered = _count_assigned_rows(
        transmission,
        ("Clusters", "transmission_cluster_id"),
    )

    core_clustered_percentage = _format_percentage(
        n_core_clustered,
        n_samples,
    )

    transmission_clustered_percentage = (
        _format_percentage(
            n_transmission_clustered,
            n_samples,
        )
        if transmission is not None
        else None
    )

    k_selection_figures = _discover_files(
        figures_path,
        ("01_k_selection*.png",),
    )

    transmission_scatterplots = _discover_files(
        figures_path,
        ("02_transmission_scatter*.png",),
    )

    core_heatmaps = _discover_files(
        figures_path,
        ("03_core_heatmap*.png",),
    )

    core_msts = _discover_files(figures_path, ("04_mst_*.html",))
    transmission_networks = _discover_files(
        figures_path,
        ("05_transmission_network_*.html",),
    )
    table_files = _discover_files(tables_path, ("*.csv", "*.tsv"))

    core_heatmap_html = _render_static_gallery(
        core_heatmaps,
        report_dir,
    )

    transmission_scatterplot_html = (
        _render_static_gallery(
            transmission_scatterplots,
            report_dir,
        )
        if transmission is not None
        else ""
    )

    k_selection_html = _render_static_gallery(
        k_selection_figures,
        report_dir,
    )


    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    cluster_type = _normalise(_get(config, "analysis", "cluster_type"))
    transmission_level = _normalise(
        _get(config, "analysis", "transmission_level")
    )

    metric_cards = [
        _metric_card(
            "Samples analysed",
            n_samples,
        ),
        _metric_card(
            "Core SNP clusters",
            n_core_clusters,
        ),
        _metric_card(
            "Core-clustered isolates",
            n_core_clustered,
            (
                f"{core_clustered_percentage} "
                "of analysed isolates"
            ),
        ),
    ]

    if transmission is not None:
        metric_cards.extend(
            [
                _metric_card(
                    "Transmission clusters",
                    n_transmission_clusters,
                ),
                _metric_card(
                    "Transmission-clustered isolates",
                    n_transmission_clustered,
                    (
                        f"{transmission_clustered_percentage} "
                        "of analysed isolates"
                    ),
                ),
            ]
        )

    metrics = "".join(metric_cards)

    navigation_items = [
        ("summary", "Summary"),
        ("configuration", "Configuration"),
        ("core-results", "Core SNP results"),
    ]

    if transmission is not None:
        navigation_items.append(
            (
                "transmission-results",
                "Transmission results",
            )
        )

    navigation_items.extend(
        [
            (
                "k-selection",
                "K-selection diagnostics",
            ),
            (
                "interactive-outputs",
                "Interactive outputs",
            ),
            ("downloads", "Downloads"),
            ("methods", "Methods"),
            ("provenance", "Provenance"),
        ]
    )

    navigation_html = "".join(
        (
            f'<a href="#{html.escape(section_id)}">'
            f"{html.escape(label)}"
            "</a>"
        )
        for section_id, label in navigation_items
    )


    if transmission is not None:
        transmission_results_html = f"""
  <section
    class="section report-section"
    id="transmission-results"
    aria-labelledby="transmission-results-heading"
  >
    <h2 id="transmission-results-heading">
      Transmission results
    </h2>

    <p>
      Transmission scatterplots show temporal and SNP-distance
      patterns within the configured epidemiological analysis groups.
    </p>

    <div class="transmission-figure-gallery">
      {transmission_scatterplot_html}
    </div>
  </section>
        """
    else:
        transmission_results_html = ""

    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SNP2Cluster-Py Report · {html.escape(project_name)}</title>
  <style>
    :root {{
      --ink: #1f2937;
      --muted: #5f6b7a;
      --line: #d9e0e7;
      --surface: #ffffff;
      --surface-alt: #f4f7fa;
      --brand: #245c6a;
      --brand-dark: #173f49;
      --accent: #dcecef;
      --warning: #fff7df;
      --warning-line: #d8a827;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background: var(--surface-alt);
      font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    a {{ color: var(--brand); }}
    .page {{ max-width: 1240px; margin: 0 auto; padding: 28px 22px 60px; }}
    .hero {{
      padding: 30px;
      color: white;
      background: linear-gradient(135deg, var(--brand-dark), var(--brand));
      border-radius: 16px;
      box-shadow: 0 8px 28px rgba(23, 63, 73, 0.18);
    }}
    .hero-layout {{
      display: grid;
      grid-template-columns: minmax(190px, 280px) minmax(0, 1fr);
      gap: 28px;
      align-items: center;
    }}

    .hero-logo-panel {{
      display: flex;
      align-items: center;
      justify-content: center;
      min-height: 150px;
      padding: 14px;
      background: #ffffff;
      border: 1px solid rgba(255, 255, 255, 0.5);
      border-radius: 14px;
      box-shadow: 0 8px 24px rgba(0, 0, 0, 0.14);
    }}

    .hero-logo {{
      display: block;
      width: 100%;
      max-width: 250px;
      max-height: 180px;
      height: auto;
      object-fit: contain;
    }}

    .hero-logo-fallback {{
      display: flex;
      align-items: center;
      justify-content: center;
      min-height: 120px;
      padding: 20px;
      color: var(--brand-dark);
      background: #ffffff;
      border-radius: 14px;
      font-size: 1.2rem;
      font-weight: 800;
      text-align: center;
    }}

    .hero-content {{
      min-width: 0;
    }}

    .hero-content h1 {{
      margin-top: 0;
    }}
    .hero h1 {{ margin: 0 0 8px; font-size: clamp(28px, 4vw, 44px); }}
    .hero p {{ margin: 4px 0; color: #e8f3f5; }}
    .section {{
      margin-top: 24px;
      padding: 24px;
      background: var(--surface);
      border: 1px solid var(--line);
      border-radius: 14px;
    }}
    .section h2 {{ margin: 0 0 16px; color: var(--brand-dark); }}
    .section h3 {{ color: var(--brand-dark); }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
      gap: 14px;
    }}
    .metric-card {{
      min-height: 112px;
      padding: 18px;
      background: var(--surface);
      border: 1px solid var(--line);
      border-radius: 12px;
      box-shadow: 0 4px 14px rgba(31, 41, 55, 0.06);
    }}
    .metric-value {{ font-size: 30px; font-weight: 750; color: var(--brand); }}
    .metric-label {{ margin-top: 5px; color: var(--muted); }}
    .table-wrap {{ overflow-x: auto; }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ padding: 10px 12px; border-bottom: 1px solid var(--line); text-align: left; }}
    th {{ width: 34%; color: var(--brand-dark); background: #f8fafc; }}
    .figure-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 18px;
    }}
    .figure-card {{ margin: 0; border: 1px solid var(--line); border-radius: 12px; overflow: hidden; }}
    .figure-card img {{ display: block; width: 100%; max-height: 560px; object-fit: contain; background: white; }}
    figcaption {{ padding: 10px 12px; color: var(--muted); background: #f8fafc; }}
    .link-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 14px;
    }}
    .link-card {{ padding: 18px; border: 1px solid var(--line); border-radius: 12px; }}
    .link-card h3 {{ margin: 0 0 14px; font-size: 16px; }}
    .button {{
      display: inline-block;
      padding: 9px 13px;
      color: white;
      background: var(--brand);
      border-radius: 8px;
      text-decoration: none;
      font-weight: 650;
    }}
    .button:hover {{ background: var(--brand-dark); }}
    .callout {{
      padding: 16px 18px;
      background: var(--warning);
      border-left: 4px solid var(--warning-line);
      border-radius: 8px;
    }}
    .download-list {{ columns: 2 300px; padding-left: 20px; }}
    .download-list li {{ margin-bottom: 7px; break-inside: avoid; }}
    .empty-state {{ color: var(--muted); font-style: italic; }}
    code {{ padding: 2px 5px; background: var(--surface-alt); border-radius: 4px; }}
    footer {{ margin-top: 24px; color: var(--muted); text-align: center; font-size: 13px; }}
    @media print {{
      body {{ background: white; }}
      .page {{ max-width: none; padding: 0; }}
      .hero, .section, .metric-card {{ box-shadow: none; break-inside: avoid; }}
      .button {{ color: var(--brand-dark); background: transparent; border: 1px solid var(--brand); }}
    }}
    .provenance-block h3 {{
      margin: 22px 0 10px;
    }}

    .provenance-block h3:first-child {{
      margin-top: 0;
    }}

    .provenance-details {{
      margin-top: 14px;
      border: 1px solid var(--line);
      border-radius: 10px;
      background: #fbfcfd;
      overflow: hidden;
    }}

    .provenance-details summary {{
      padding: 13px 16px;
      color: var(--brand-dark);
      background: #f3f7f8;
      cursor: pointer;
      font-weight: 700;
    }}

    .provenance-details summary:focus-visible {{
      outline: 3px solid rgba(20, 94, 128, 0.25);
      outline-offset: -3px;
    }}

    .provenance-details[open] summary {{
      border-bottom: 1px solid var(--line);
    }}

    .provenance-details .table-wrap,
    .provenance-details .empty-state,
    .provenance-details .validation-success,
    .provenance-details .validation-list {{
      margin: 14px 16px;
    }}

    .provenance-details td {{
      overflow-wrap: anywhere;
      word-break: break-word;
    }}

    .validation-success {{
      padding: 12px 14px;
      color: #155724;
      background: #edf8f0;
      border-left: 4px solid #2f855a;
      border-radius: 7px;
    }}

    .validation-list {{
      padding-left: 22px;
    }}

    .validation-list li + li {{
      margin-top: 7px;
    }}

    html {{
      scroll-behavior: smooth;
      scroll-padding-top: 84px;
    }}

    .report-nav {{
      position: sticky;
      top: 0;
      z-index: 100;
      margin: 18px 0 24px;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 5px 18px rgba(20, 46, 58, 0.08);
      backdrop-filter: blur(8px);
    }}

    .report-nav-inner {{
      display: flex;
      gap: 6px;
      align-items: center;
      overflow-x: auto;
      padding: 9px;
      scrollbar-width: thin;
    }}

    .report-nav a {{
      flex: 0 0 auto;
      padding: 8px 11px;
      color: var(--brand-dark);
      border-radius: 8px;
      font-size: 0.9rem;
      font-weight: 700;
      line-height: 1.2;
      text-decoration: none;
      white-space: nowrap;
    }}

    .report-nav a:hover {{
      color: #ffffff;
      background: var(--brand-dark);
    }}

    .report-nav a:focus-visible {{
      color: #ffffff;
      background: var(--brand-dark);
      outline: 3px solid rgba(20, 94, 128, 0.25);
      outline-offset: 2px;
    }}

    .report-section {{
      scroll-margin-top: 84px;
    }}

    .metric-detail {{
      display: block;
      margin-top: 5px;
      color: #5f6f76;
      font-size: 0.83rem;
      font-weight: 600;
      line-height: 1.25;
    }}

    .primary-figure .gallery {{
      grid-template-columns: 1fr;
    }}

    .primary-figure .figure-card {{
      max-width: none;
    }}

    .primary-figure img {{
      width: 100%;
      max-height: none;
      object-fit: contain;
    }}

    .transmission-figure-gallery .gallery {{
      grid-template-columns: repeat(
        auto-fit,
        minmax(340px, 1fr)
      );
    }}

    .figure-disclosure {{
      margin-top: 18px;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: #fbfcfd;
      overflow: hidden;
    }}

    .figure-disclosure summary {{
      padding: 14px 16px;
      color: var(--brand-dark);
      background: #f3f7f8;
      cursor: pointer;
      font-weight: 700;
    }}

    .figure-disclosure summary:hover {{
      background: #eaf2f4;
    }}

    .figure-disclosure summary:focus-visible {{
      outline: 3px solid rgba(20, 94, 128, 0.25);
      outline-offset: -3px;
    }}

    .figure-disclosure[open] summary {{
      border-bottom: 1px solid var(--line);
    }}

    .figure-disclosure-content {{
      padding: 16px;
    }}

    .download-list {{
      display: grid;
      gap: 12px;
      margin-top: 18px;
    }}

    .download-card {{
      display: flex;
      gap: 18px;
      align-items: center;
      justify-content: space-between;
      padding: 16px;
      border: 1px solid var(--line);
      border-radius: 10px;
      background: #fbfcfd;
    }}

    .download-card-content {{
      min-width: 0;
    }}

    .download-card h3 {{
      margin: 0 0 5px;
      color: var(--brand-dark);
      font-size: 1rem;
      overflow-wrap: anywhere;
    }}

    .download-description {{
      margin: 0;
      color: #5f6f76;
      font-size: 0.9rem;
      line-height: 1.45;
    }}

    @media (max-width: 760px) {{
      .hero-layout {{
        grid-template-columns: 1fr;
      }}

      .hero-logo-panel {{
        width: 100%;
        max-width: 320px;
        min-height: auto;
      }}

      .hero-logo {{
        max-height: 140px;
      }}

      .transmission-figure-gallery .gallery {{
        grid-template-columns: 1fr;
      }}

      .figure-disclosure-content {{
        padding: 10px;
      }}

      .download-card {{
        align-items: flex-start;
        flex-direction: column;
      }}

    }}

    @media print {{
      .provenance-details {{
        break-inside: avoid;
      }}

      .provenance-details summary {{
        background: #ffffff;
      }}

      .hero-layout {{
        grid-template-columns: 180px 1fr;
      }}

      .hero-logo-panel {{
        box-shadow: none;
        border: 1px solid #cccccc;
      }}

      .report-nav {{
        display: none;
      }}

      .report-section {{
        scroll-margin-top: 0;
      }}

      .figure-disclosure {{
        break-inside: avoid;
      }}

      .figure-disclosure summary {{
        background: #ffffff;
      }}
    }}
  </style>
</head>
<body>
<main class="page">
  <header class="hero">
    <div class="hero-layout">
      {logo_html}

      <div class="hero-content">
        <h1>SNP2Cluster-Py Analysis Report</h1>

        <p>
          <strong>Project:</strong>
          {html.escape(project_name)}
        </p>

        <p>
          <strong>Analysis:</strong>
          {html.escape(cluster_type)}
          ·
          <strong>Transmission level:</strong>
          {html.escape(transmission_level)}
        </p>

        <p>
          <strong>Report generated:</strong>
          {generated_at}
        </p>
      </div>
    </div>
  </header>

  <nav
    class="report-nav"
    aria-label="Report sections"
  >
    <div class="report-nav-inner">
      {navigation_html}
    </div>
  </nav>

    <section
    class="section report-section"
    id="summary"
    aria-labelledby="summary-heading"
  >
    <h2 id="summary-heading">Analysis summary</h2>
    <div class="metrics">{metrics}</div>
  </section>

  <section
    class="section report-section"
    id="configuration"
    aria-labelledby="configuration-heading"
  >
    <h2 id="configuration-heading">
      Analysis configuration
    </h2>
    {_render_table(_config_rows(config))}
  </section>

  <section
    class="section report-section"
    id="core-results"
    aria-labelledby="core-results-heading"
  >
    <h2 id="core-results-heading">
      Core SNP results
    </h2>

    <p>
      The annotated Core SNP heatmap summarises pairwise SNP-distance
      relationships and the retained context-aware Core SNP cluster
      assignments.
    </p>

    <div class="primary-figure">
      {core_heatmap_html}
    </div>
  </section>

  {transmission_results_html}

  <section
    class="section report-section"
    id="k-selection"
    aria-labelledby="k-selection-heading"
  >
    <h2 id="k-selection-heading">
      K-selection diagnostics
    </h2>

    <p>
      These figures document the initial data-driven cluster-number
      selection undertaken within eligible Core-clustering contexts.
      Where candidate multi-cluster solutions are supported,
      silhouette diagnostics are used to compare candidate K values.
      An initial K of 1 is valid where the supplied SNP-distance data
      support one initial group. That group still proceeds through the
      configured SNP-cutoff refinement before final Core SNP clusters
      are reported.
    </p>

    <details class="figure-disclosure">
      <summary>
        View K-selection diagnostic figures
        ({len(k_selection_figures)})
      </summary>

      <div class="figure-disclosure-content">
        {k_selection_html}
      </div>
    </details>
  </section>

  <section
    class="section report-section"
    id="interactive-outputs"
    aria-labelledby="mst-heading"
  >
    <h2 id="mst-heading">
      Interactive Core SNP minimum spanning trees
    </h2>
    <p>
      Core MSTs summarize SNP-distance relationships and
      context-aware Core SNP clusters. Open each analysis
      group in a separate browser tab.
    </p>
    {_render_link_cards(
        core_msts,
        report_dir,
        button_text="Open Core SNP MST",
    )}
  </section>

  <section
    class="section report-section"
    aria-labelledby="network-heading"
  >
    <h2 id="network-heading">
      Interactive transmission-cluster networks
    </h2>

    <div class="callout">
      Nodes represent isolates assigned to inferred transmission
      clusters. Network edges are post-clustering SNP
      minimum-spanning-tree connections used to visualize genetic
      relationships. They do not represent direct transmission and
      do not necessarily satisfy the configured temporal cutoff.
      Transmission-cluster membership may reflect qualifying chains
      through intermediate isolates. Review the transmission
      scatterplots for temporal patterns.
    </div>

    <div style="margin-top:16px">
      {_render_link_cards(
          transmission_networks,
          report_dir,
          button_text="Open transmission network",
      )}
    </div>
  </section>

  <section
    class="section report-section"
    id="downloads"
    aria-labelledby="tables-heading"
  >
    <h2 id="tables-heading">Output tables</h2>
    <p>
      Download the machine-readable analytical outputs generated for
      this run. The descriptions below distinguish initial K-means
      assignments, retained Core SNP clusters, transmission-cluster
      outputs and network-visualisation data.
    </p>
    {_render_downloads(
        table_files,
        report_dir,
    )}
  </section>

  <section
    class="section report-section"
    id="methods"
    aria-labelledby="methods-heading"
  >
    <h2 id="methods-heading">
      Methods and interpretation guidance
    </h2>

    <h3>Core-clustering contexts</h3>

    <p>
      Core SNP clustering is performed within the configured ordered
      Core-clustering context. Depending on the configuration, the
      context may represent sequence type, the main epidemiological
      variable, a secondary variable, or an ordered combination of
      variables. A context is eligible for Core clustering when it
      contains at least two isolates.
    </p>

    <h3>Initial data-driven K clustering</h3>

    <p>
      For each eligible context, SNP2Cluster-Py uses the supplied Core
      SNP-distance matrix to perform an initial data-driven clustering
      step. Candidate multi-cluster solutions are assessed using
      silhouette-based diagnostics where the data support comparison
      between two or more initial clusters.
    </p>

    <p>
      An initial K of 1 is a valid data-driven result when the
      SNP-distance data support one initial group. This includes
      contexts containing multiple isolates with one distinct
      SNP-distance profile. A numerical silhouette width is not
      defined for a one-cluster solution, but the K=1 group remains
      eligible for the subsequent SNP-cutoff refinement step.
    </p>

    <h3>SNP-cutoff refinement</h3>

    <p>
      After initial K assignment, each initial group is evaluated
      using the unscaled pairwise Core SNP distances and the configured
      SNP cutoff. Pairwise comparisons above the cutoff provide points
      for splitting the initial group according to the implemented
      chain-refinement logic. Only final Core SNP clusters containing
      at least two isolates are retained and reported.
    </p>

    <h3>Context-aware cluster labels</h3>

    <p>
      Final Core SNP cluster labels retain the configured context so
      that cluster identifiers remain interpretable across analytical
      groups. Core SNP cluster membership therefore reflects both the
      selected analytical context and the post-K-selection SNP-cutoff
      refinement.
    </p>

    <h3>Transmission clustering</h3>

    <p>
      In Transmission mode, SNP2Cluster-Py subsequently integrates
      Core SNP cluster membership with the configured epidemiological
      grouping, specimen-collection dates, SNP threshold and temporal
      threshold. Transmission-cluster assignments should be
      interpreted alongside the transmission scatterplots and the
      available epidemiological information.
    </p>

    <h3>Interactive networks</h3>

    <p>
      Minimum-spanning-tree and transmission-network edges visualise
      genetic relationships between isolates. Network connections are
      analytical visualisation aids and must not be interpreted as
      confirmed direct transmission events.
    </p>

    <h3>Important interpretation limitation</h3>

    <p>
      This report summarises analytical outputs produced by
      SNP2Cluster-Py. Results should be interpreted alongside the
      effective configuration, laboratory information,
      epidemiological metadata, sampling strategy and relevant
      public-health context.
    </p>
  </section>

  <section
    class="section report-section"
    id="provenance"
    aria-labelledby="provenance-heading"
  >
    <h2 id="provenance-heading">
      Reproducibility and provenance
    </h2>
    {_render_manifest(
        run_manifest,
        software_version,
    )}
  </section>

  <footer>
    Generated by SNP2Cluster-Py · Open scientific analysis report
  </footer>
</main>
</body>
</html>
"""

    output_path.write_text(document, encoding="utf-8")
    return output_path
