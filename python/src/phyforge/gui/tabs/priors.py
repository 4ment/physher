"""Tab 6 -- Priors: the tree prior, and the prior on every parameter."""

from __future__ import annotations

import streamlit as st

from ...build.coalescent import (
    PIECEWISE_MODELS,
    ROLE_GROWTH,
    ROLE_PRECISION,
    ROLE_THETA,
)
from ...spec import COALESCENT_MODELS
from .. import state
from ..state import Analysis
from ..widgets import bind, init_table, prior_table, show_problems


def render(analysis: Analysis) -> None:
    st.header("Priors")
    if not analysis.spec.uses_priors:
        st.info(
            "Maximum likelihood fits the likelihood alone, so it has no priors. "
            "Pick map, mcmc or advi in the sidebar to set them."
        )
        return

    if state.is_time_tree():
        st.subheader("Tree prior")
        kind = bind(
            st.selectbox,
            "coalescent.kind",
            "Coalescent",
            analysis,
            options=[None, *COALESCENT_MODELS],
            format_func=lambda k: "None (uniform node heights)" if k is None else k,
        )
        if kind in ("skygrid", "skyglide"):
            columns = st.columns(2)
            with columns[0]:
                bind(
                    st.number_input,
                    "coalescent.grid",
                    "Grid points",
                    analysis,
                    value=None,
                    min_value=2,
                    step=1,
                )
            with columns[1]:
                bind(
                    st.number_input,
                    "coalescent.cutoff",
                    "Cutoff (oldest grid point)",
                    analysis,
                    value=None,
                    min_value=0.0,
                    format="%g",
                )
        if kind in PIECEWISE_MODELS:
            bind(
                st.checkbox,
                "coalescent.gmrf",
                "Smooth the population sizes with a GMRF",
                analysis,
            )
        show_problems(analysis, "coalescent.kind")
        init_table(
            analysis, (ROLE_THETA, ROLE_GROWTH, ROLE_PRECISION), "Coalescent parameters"
        )

    st.subheader("Parameter priors")
    st.caption(
        "Each estimated parameter starts with the prior phyforge would write; pick "
        "another distribution to replace it, or 'none' to drop it. A vector "
        "hyperparameter takes one value per element."
    )
    prior_table(analysis)
