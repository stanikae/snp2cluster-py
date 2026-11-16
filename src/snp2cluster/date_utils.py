from __future__ import annotations

import pandas as pd


DATE_FORMAT_PATTERNS = {
    "ymd": (
        "%Y/%m/%d",
        "%y/%m/%d",
    ),
    "ydm": (
        "%Y/%d/%m",
        "%y/%d/%m",
    ),
    "dmy": (
        "%d/%m/%Y",
        "%d/%m/%y",
    ),
    "mdy": (
        "%m/%d/%Y",
        "%m/%d/%y",
    ),
}


def parse_dates(
    series: pd.Series,
    date_format: str = "ymd",
) -> pd.Series:
    """Parse dates using the configured component order.

    Supported orders are ``ymd``, ``ydm``, ``dmy`` and ``mdy``.
    Dash, slash and dot separators are accepted, with two-digit or
    four-digit years.

    Missing, empty and invalid values are returned as ``NaT``.
    """
    patterns = DATE_FORMAT_PATTERNS.get(date_format)

    if patterns is None:
        raise ValueError(
            f"Unsupported date format: {date_format}"
        )

    normalized = (
        series.astype("string")
        .str.strip()
        .str.replace(
            r"[.-]",
            "/",
            regex=True,
        )
        .replace("", pd.NA)
    )

    parsed = pd.Series(
        pd.NaT,
        index=series.index,
        dtype="datetime64[ns]",
        name=series.name,
    )

    for pattern in patterns:
        unresolved = (
            parsed.isna()
            & normalized.notna()
        )

        if not unresolved.any():
            break

        candidate = pd.to_datetime(
            normalized.loc[unresolved],
            format=pattern,
            errors="coerce",
        )

        successful_indices = candidate.index[
            candidate.notna()
        ]

        if len(successful_indices) > 0:
            parsed.loc[successful_indices] = (
                candidate.loc[successful_indices]
            )

    return parsed
