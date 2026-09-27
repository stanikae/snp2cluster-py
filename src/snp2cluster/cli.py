from __future__ import annotations

from pathlib import Path
import typer
from rich.console import Console
from rich.table import Table

from .config import load_config
from .pipeline import run_pipeline
from .validation import validate_all, ValidationError

app = typer.Typer(help="SNP2Cluster-Py: genomic transmission cluster detection CLI")
console = Console()


@app.command()
def validate(config: Path = typer.Option(..., "--config", "-c", help="Path to YAML config")):
    """Validate configuration and input files."""
    cfg = load_config(config)
    messages = validate_all(cfg, strict_files=True)

    if messages:
        console.print("[bold red]Validation failed[/bold red]")
        for msg in messages:
            console.print(f"- {msg}")
        raise typer.Exit(code=1)

    console.print("[bold green]Validation passed[/bold green]")


@app.command()
def run(config: Path = typer.Option(..., "--config", "-c", help="Path to YAML config")):
    """Run SNP2Cluster analysis."""
    cfg = load_config(config)
    try:
        dirs = run_pipeline(cfg)
    except ValidationError as exc:
        console.print(f"[bold red]{exc}[/bold red]")
        raise typer.Exit(code=1)

    table = Table(title="SNP2Cluster-Py run completed")
    table.add_column("Output", style="cyan")
    table.add_column("Path")
    for key, path in dirs.items():
        table.add_row(key, str(path))
    console.print(table)


@app.command()
def version():
    """Show package version."""
    from . import __version__
    console.print(__version__)

if __name__ == "__main__":
    app()
