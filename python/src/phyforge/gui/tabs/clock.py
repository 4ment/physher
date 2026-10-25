"""Tab 5 -- Clock Model."""

from __future__ import annotations

import streamlit as st

from ...build.branchmodel import ROLE_BRANCH_RATES, ROLE_CLOCK_SCALE, ROLE_RATE
from ...spec import CLOCK_MODELS
from .. import state
from ..state import Analysis
from ..widgets import bind, init_table, show_problems

CLOCKS = {"strict": "Strict", "ncln": "Uncorrelated lognormal (non-centred)"}


def render(analysis: Analysis) -> None:
    st.header("Clock Model")
    if not state.is_time_tree():
        st.info(
            "The tree is unrooted, so there is no clock. Turn on tip dates (Tab 2) "
            "for a time tree."
        )
        return

    bind(
        st.radio,
        "clock.kind",
        "Clock",
        analysis,
        options=list(CLOCK_MODELS),
        format_func=CLOCKS.get,
        horizontal=True,
    )

    key = "w:@fix_rate"
    if key not in st.session_state:
        st.session_state[key] = state.values()["clock.rate"] is not None

    def changed() -> None:
        values = state.values()
        values["clock.rate"] = (
            values["clock.rate_init"] if st.session_state[key] else None
        )
        st.session_state.pop("w:clock.rate", None)

    fixed = st.checkbox("Fix the substitution rate", key=key, on_change=changed)
    if fixed:
        bind(
            st.number_input,
            "clock.rate",
            "Substitution rate",
            analysis,
            min_value=0.0,
            format="%g",
        )
    show_problems(analysis, "clock.rate")
    show_problems(analysis, "clock.kind")

    init_table(
        analysis, (ROLE_RATE, ROLE_BRANCH_RATES, ROLE_CLOCK_SCALE), "Clock parameters"
    )
