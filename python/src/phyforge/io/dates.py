"""Tip dates for time trees.

physher wants a decimal year per taxon.  Dates come either from the taxon
names (the common case for sequence data) or from a CSV file.
"""

from __future__ import annotations

import csv
import re
from typing import Dict, Iterable, Optional

DAYS_IN_MONTH = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


class DateParseError(ValueError):
    """Raised when a date cannot be recovered for one or more taxa."""


def decimal_year(day: int, month: int, year: int) -> float:
    """Convert a calendar date to a decimal year."""
    days = list(DAYS_IN_MONTH)
    if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        days[1] = 29
    offset = sum(days[: month - 1]) + day - 1
    return year + offset / sum(days)


def dates_from_names(
    taxa: Iterable[str],
    regex: str = r"_(\d+\.?\d*)$",
    date_format: Optional[str] = None,
) -> Dict[str, float]:
    """Extract a date from each taxon name using ``regex``.

    With no ``date_format`` the first capture group is read as a decimal year.
    With one (e.g. ``yyyy-MM-dd``) the pattern is expected to capture the day,
    month and year in groups named by that format's field order.
    """
    pattern = re.compile(regex)
    fields = _format_fields(date_format)
    dates: Dict[str, float] = {}
    missing = []
    for taxon in taxa:
        match = pattern.search(taxon)
        if match is None:
            missing.append(taxon)
            continue
        if fields is None or len(match.groups()) < 3:
            dates[taxon] = float(match.group(1))
        else:
            day, month, year = fields
            dates[taxon] = decimal_year(
                int(match.group(day)), int(match.group(month)), int(match.group(year))
            )
    if missing:
        raise DateParseError(
            f"could not extract a date from {len(missing)} taxon name(s) using "
            f"{regex!r}; first is {missing[0]!r}"
        )
    return dates


def dates_from_csv(
    path: str,
    date_format: Optional[str] = None,
    name_column: str = "strain",
    date_column: str = "date",
) -> Dict[str, float]:
    """Read dates from a CSV file with a header naming the two columns."""
    dates: Dict[str, float] = {}
    with open(path, newline="") as fp:
        reader = csv.DictReader(fp, skipinitialspace=True)
        if reader.fieldnames is None:
            raise DateParseError(f"{path} has no header row")
        for column in (name_column, date_column):
            if column not in reader.fieldnames:
                raise DateParseError(
                    f"{path} has no {column!r} column (found "
                    f"{', '.join(reader.fieldnames)})"
                )
        for row in reader:
            dates[row[name_column]] = _parse_date(row[date_column], date_format)
    if not dates:
        raise DateParseError(f"{path} contains no rows")
    return dates


def _parse_date(value: str, date_format: Optional[str]) -> float:
    value = value.strip()
    if date_format is None:
        return float(value)
    day, month, year = _format_fields(date_format)
    parts = re.split(r"[/-]", value)
    return decimal_year(
        int(parts[day - 1]), int(parts[month - 1]), int(parts[year - 1])
    )


def _format_fields(date_format: Optional[str]):
    """Map a format such as ``yyyy-MM-dd`` to 1-based (day, month, year) slots."""
    if date_format is None:
        return None
    parts = re.split(r"[/-]", date_format)
    try:
        return parts.index("dd") + 1, parts.index("MM") + 1, parts.index("yyyy") + 1
    except ValueError:
        raise DateParseError(
            f"date format {date_format!r} must contain dd, MM and yyyy"
        ) from None
