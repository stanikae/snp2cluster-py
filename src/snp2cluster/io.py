from __future__ import annotations

from pathlib import Path
import pandas as pd


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path)
    if suffix in {".tsv"}:
        return pd.read_csv(path, sep="\t")
    if suffix in {".xlsx", ".xlsm"}:
        return pd.read_excel(path, engine="openpyxl")
    if suffix == ".xls":
        return pd.read_excel(path, engine="xlrd")
    if suffix == ".parquet":
        return pd.read_parquet(path)

    raise ValueError(f"Unsupported input format: {path}")


def read_snp_matrix(path: str | Path) -> pd.DataFrame:
    """Read pairwise SNP distance matrix.

    Expected common format:
    - first column contains sample IDs, or
    - row index already contains sample IDs.
    """
    df = read_table(path)

    # If first column looks like sample IDs, set it as index.
    first_col = df.columns[0]
    if not pd.api.types.is_numeric_dtype(df[first_col]):
        df = df.set_index(first_col)

    # Make column labels strings for stable matching.
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def write_csv(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
