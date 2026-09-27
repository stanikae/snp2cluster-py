from __future__ import annotations

import pandas as pd

DATE_FORMAT_MAP = {
    "ymd": {
        "yearfirst": True,
        "dayfirst": False,
    },
    "ydm": {
        "yearfirst": True,
        "dayfirst": True,
    },
    "dmy": {
        "yearfirst": False,
        "dayfirst": True,
    },
    "mdy": {
        "yearfirst": False,
        "dayfirst": False,
    },
}


def parse_dates(
    series: pd.Series,
    date_format: str = "ymd",
) -> pd.Series:
    """
    Lubridate-style date parsing.

    Accepts:
        2019-12-17
        2019/12/17
        2019.12.17

    Supported:
        ymd
        ydm
        dmy
        mdy
    """

    settings = DATE_FORMAT_MAP.get(date_format)

    if settings is None:
        raise ValueError(
            f"Unsupported date format: {date_format}"
        )

    return pd.to_datetime(
        series,
        yearfirst=settings["yearfirst"],
        dayfirst=settings["dayfirst"],
        errors="coerce",
    )
