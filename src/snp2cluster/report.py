from __future__ import annotations

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
    """Create a concise label from a generated output filename."""
    stem = path.stem

    for prefix in (
        "01_k_selection_",
        "01_k_selection",
        "02_transmission_scatter_",
        "02_transmission_scatter",
        "03_core_heatmap",
        "04_mst_",
        "05_transmission_network_",
    ):
        if stem.startswith(prefix):
            stem = stem[len(prefix) :]
            break

    return stem.replace("__", " / ").replace("_", " ").strip() or path.stem


def _metric_card(label: str, value: object) -> str:
    display = "NA" if value is None else str(value)
    return (
        '<div class="metric-card">'
        f'<div class="metric-value">{html.escape(display)}</div>'
        f'<div class="metric-label">{html.escape(label)}</div>'
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


def _render_downloads(files: Sequence[Path], report_dir: Path) -> str:
    if not files:
        return '<p class="empty-state">No output tables were found.</p>'

    items = []

    for path in files:
        url = html.escape(_relative_url(path, report_dir), quote=True)
        items.append(
            "<li>"
            f'<a href="{url}" download>{html.escape(path.name)}</a>'
            "</li>"
        )

    return '<ul class="download-list">' + "".join(items) + "</ul>"


def _render_manifest(
    run_manifest: Mapping[str, object] | None,
    software_version: str | None,
) -> str:
    rows: list[tuple[str, str]] = []

    if software_version:
        rows.append(("SNP2Cluster-Py version", str(software_version)))

    if run_manifest:
        for key, value in run_manifest.items():
            if isinstance(value, (str, int, float, bool)) or value is None:
                rows.append((str(key).replace("_", " ").title(), _normalise(value)))

    if not rows:
        return '<p class="empty-state">No run-manifest metadata was supplied.</p>'

    return _render_table(rows)


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
) -> Path:
    """Write the open scientific SNP2Cluster-Py HTML run report."""
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

    static_figures = _discover_files(
        figures_path,
        (
            "01_k_selection*.png",
            "02_transmission_scatter*.png",
            "03_core_heatmap*.png",
        ),
    )
    core_msts = _discover_files(figures_path, ("04_mst_*.html",))
    transmission_networks = _discover_files(
        figures_path,
        ("05_transmission_network_*.html",),
    )
    table_files = _discover_files(tables_path, ("*.csv", "*.tsv"))

    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    cluster_type = _normalise(_get(config, "analysis", "cluster_type"))
    transmission_level = _normalise(
        _get(config, "analysis", "transmission_level")
    )

    metrics = "".join(
        (
            _metric_card("Samples analysed", n_samples),
            _metric_card("Core SNP clusters", n_core_clusters),
            _metric_card("Core-clustered isolates", n_core_clustered),
            _metric_card("Transmission clusters", n_transmission_clusters),
            _metric_card("Transmission-clustered isolates", n_transmission_clustered),
        )
    )

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
  </style>
</head>
<body>
<main class="page">
  <header class="hero">
    <h1>SNP2Cluster-Py Analysis Report</h1>
    <p><strong>Project:</strong> {html.escape(project_name)}</p>
    <p><strong>Analysis:</strong> {html.escape(cluster_type)} · <strong>Transmission level:</strong> {html.escape(transmission_level)}</p>
    <p><strong>Report generated:</strong> {generated_at}</p>
  </header>

  <section class="section" aria-labelledby="summary-heading">
    <h2 id="summary-heading">Analysis summary</h2>
    <div class="metrics">{metrics}</div>
  </section>

  <section class="section" aria-labelledby="configuration-heading">
    <h2 id="configuration-heading">Analysis configuration</h2>
    {_render_table(_config_rows(config))}
  </section>

  <section class="section" aria-labelledby="figures-heading">
    <h2 id="figures-heading">Static scientific figures</h2>
    {_render_static_gallery(static_figures, report_dir)}
  </section>

  <section class="section" aria-labelledby="mst-heading">
    <h2 id="mst-heading">Interactive Core SNP minimum spanning trees</h2>
    <p>Core MSTs summarize SNP-distance relationships and context-aware Core SNP clusters. Open each analysis group in a separate browser tab.</p>
    {_render_link_cards(core_msts, report_dir, button_text="Open Core SNP MST")}
  </section>

  <section class="section" aria-labelledby="network-heading">
    <h2 id="network-heading">Interactive transmission-cluster networks</h2>
    <div class="callout">
      Nodes represent isolates assigned to inferred transmission clusters. Network edges are post-clustering SNP minimum-spanning-tree connections used to visualize genetic relationships. They do not represent direct transmission and do not necessarily satisfy the configured temporal cutoff. Transmission-cluster membership may reflect qualifying chains through intermediate isolates. Review the transmission scatterplots for temporal patterns.
    </div>
    <div style="margin-top:16px">
      {_render_link_cards(transmission_networks, report_dir, button_text="Open transmission network")}
    </div>
  </section>

  <section class="section" aria-labelledby="tables-heading">
    <h2 id="tables-heading">Output tables</h2>
    <p>Download the machine-readable results generated for this run.</p>
    {_render_downloads(table_files, report_dir)}
  </section>

  <section class="section" aria-labelledby="methods-heading">
    <h2 id="methods-heading">Methods and interpretation guidance</h2>
    <h3>Core SNP clustering</h3>
    <p>Core SNP clustering groups isolates within the configured ordered Core context. The context determines which isolates are eligible to be clustered together, while the context-aware cluster label keeps cluster identifiers interpretable across analytical groups.</p>
    <h3>Transmission clustering</h3>
    <p>Transmission analysis uses the configured epidemiological grouping, collection dates, Core SNP clusters, SNP threshold, and temporal threshold. Cluster membership should be interpreted alongside the temporal scatterplots and available epidemiological information.</p>
    <h3>Important limitation</h3>
    <p>This report summarizes analytical outputs produced by SNP2Cluster-Py. Results should be interpreted alongside the configured parameters, laboratory information, epidemiological metadata, and public-health context. Network connections are visualization aids and should not be interpreted as confirmed direct transmission events.</p>
  </section>

  <section class="section" aria-labelledby="provenance-heading">
    <h2 id="provenance-heading">Reproducibility and provenance</h2>
    {_render_manifest(run_manifest, software_version)}
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
