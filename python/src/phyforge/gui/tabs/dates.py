"""Tab 2 -- Dates."""

from __future__ import annotations

from typing import Dict, Optional

import streamlit as st

from ...io.dates import DateParseError, dates_from_csv, dates_from_names
from ...io.fasta import read_fasta
from .. import state, uploads
from ..state import Analysis
from ..widgets import bind, bind_session, optional_text, show_problems

SOURCES = {
    "regex": "From taxon names (regex)",
    "csv": "From a CSV file",
    "zero": "All contemporaneous",
}


def taxa() -> Optional[list]:
    entry = uploads.info("alignment")
    if entry is None:
        return None
    try:
        return [sequence.taxon for sequence in read_fasta(entry["path"])]
    except (OSError, ValueError):
        return None


def tip_dates(analysis: Analysis) -> Optional[Dict[str, float]]:
    """The dates the build would read, or None while they cannot be read."""
    names = taxa()
    if names is None or not state.is_time_tree():
        return None
    dates = analysis.spec.clock.dates
    try:
        if dates.source == "zero":
            return {name: 0.0 for name in names}
        if dates.source == "csv":
            if dates.path is None:
                return None
            table = dates_from_csv(dates.path, dates.date_format)
            return {name: table[name] for name in names if name in table}
        return dates_from_names(names, dates.regex, dates.date_format)
    except (DateParseError, OSError, ValueError, KeyError, IndexError):
        return None


def render(analysis: Analysis) -> None:
    st.header("Dates")
    names = taxa()
    if names is None:
        st.info("Load an alignment first (Tab 1).")
        return

    bind_session(
        st.checkbox,
        state.TIME_TREE,
        "Use tip dates (for time trees)",
        help="A time tree has a molecular clock (Tab 5) and node heights instead "
        "of branch lengths; without dates the tree is unrooted.",
    )
    if not state.is_time_tree():
        st.info("Dates disabled – the tree is treated as unrooted.")
        return

    source = bind(
        st.radio,
        "clock.dates.source",
        "Date source",
        analysis,
        options=list(SOURCES),
        format_func=SOURCES.get,
        horizontal=True,
    )
    if source == "regex":
        bind(
            st.text_input,
            "clock.dates.regex",
            "Regex (first capture group = decimal year)",
            analysis,
        )
    elif source == "csv":
        uploads.upload_widget("dates", "CSV with 'strain' and 'date' columns")
        show_problems(analysis, "clock.dates.path")
    if source != "zero":
        optional_text(
            "clock.dates.date_format",
            "Date format",
            analysis,
            placeholder="decimal year; or e.g. yyyy-MM-dd",
            help="Leave blank when the dates are decimal years. With a format "
            "such as yyyy-MM-dd, a regex has to capture day, month and year.",
        )
    show_problems(analysis, "clock.dates.source")

    dates = tip_dates(analysis)
    if dates is None:
        error = _date_error(analysis, names)
        if error:
            st.error(error, icon="⚠️")
        return
    if source != "zero":
        missing = len(names) - len(dates)
        low, high = min(dates.values()), max(dates.values())
        message = (
            f"Dates loaded for **{len(dates)}/{len(names)}** taxa · range "
            f"**{low:.3f} – {high:.3f}**"
        )
        (st.warning if missing else st.success)(message)

    import pandas as pd

    with st.expander("Date table", expanded=True):
        newest = max(dates.values()) if dates else 0.0
        frame = pd.DataFrame(
            {
                "Taxon": list(dates),
                "Date": list(dates.values()),
                "Height": [newest - value for value in dates.values()],
            }
        )
        st.dataframe(frame, width="stretch", hide_index=True)

    if source != "zero":
        with st.expander("Date histogram", expanded=False):
            import plotly.graph_objects as go

            count = len(dates)
            bins = st.slider(
                "Number of bins",
                min_value=1,
                max_value=min(100, max(2, count)),
                value=min(30, max(5, count // 3 or 1)),
                key="dates_hist_nbins",
            )
            figure = go.Figure(
                go.Histogram(
                    x=list(dates.values()),
                    nbinsx=bins,
                    marker=dict(color="#1f77b4", line=dict(color="#0d4c80", width=1)),
                )
            )
            figure.update_layout(
                xaxis_title="Date",
                yaxis_title="Count",
                bargap=0.02,
                height=320,
                margin=dict(l=40, r=20, t=20, b=40),
            )
            st.plotly_chart(figure, width="stretch")


def _date_error(analysis: Analysis, names) -> Optional[str]:
    """Why the dates cannot be read, in the parser's own words."""
    dates = analysis.spec.clock.dates
    try:
        if dates.source == "csv":
            if dates.path is None:
                return None
            table = dates_from_csv(dates.path, dates.date_format)
            missing = [name for name in names if name not in table]
            if missing:
                return (
                    f"{len(missing)} taxa are not in the file; first is {missing[0]!r}"
                )
            return None
        dates_from_names(names, dates.regex, dates.date_format)
    except (DateParseError, OSError, ValueError, IndexError) as error:
        return str(error)
    return None
