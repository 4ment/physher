"""Tab 4 -- Site Model: substitution model and rate variation across sites."""

from __future__ import annotations

import streamlit as st

from ...build.sitemodel import (
    ROLE_PROPORTION_INVARIANT,
    ROLE_SHAPE,
    ROLE_SITE_PROPORTIONS,
    ROLE_SITE_RATES,
    ROLE_TOP_RATE,
)
from ...build.substitution import (
    EQUAL_FREQUENCY_MODELS,
    ROLE_FREQUENCIES,
    ROLE_KAPPA,
    ROLE_RATES,
)
from ...spec import CAT_ASSIGNMENTS, RATE_PARAMETERIZATIONS, SUBSTITUTION_MODELS
from .. import state
from ..state import Analysis
from ..widgets import bind, init_table, optional_text, select, show_problems

#: How rates vary across sites, as one choice; "none" is a single category.
RATE_MODELS = {
    "none": "None (one rate for every site)",
    "gamma": "Discretized gamma (+G)",
    "weibull": "Discretized Weibull",
    "lognormal": "Discretized lognormal",
    "discrete": "Free rates (+R)",
    "cat": "Empirical CAT",
}

QUADRATURES = [None, "median", "mean", "gausslaguerre"]


def _rate_model(analysis: Analysis) -> str:
    values = state.values()
    key = "w:@rate_model"
    if key not in st.session_state:
        none = values["site.categories"] <= 1 and values["site.distribution"] != "cat"
        st.session_state[key] = "none" if none else values["site.distribution"]

    def changed() -> None:
        choice = st.session_state[key]
        if choice == "none":
            values["site.categories"] = 1
            values["site.distribution"] = "gamma"
        else:
            values["site.distribution"] = choice
            if values["site.categories"] < 2:
                values["site.categories"] = 20 if choice == "cat" else 4
        st.session_state.pop("w:site.categories", None)

    choice = st.selectbox(
        "Rate variation across sites",
        list(RATE_MODELS),
        key=key,
        on_change=changed,
        format_func=RATE_MODELS.get,
    )
    show_problems(analysis, "site.distribution")
    return choice


def render(analysis: Analysis) -> None:
    st.header("Site Model")

    st.subheader("Substitution model")
    datatype = st.session_state[state.DATATYPE]
    models = [
        name
        for name, kind in SUBSTITUTION_MODELS.items()
        if kind == datatype  # codon models are refused by validate anyway
    ]
    left, right = st.columns(2)
    with left:
        model = select("substitution.kind", "Model", models, analysis)
    with right:
        if model.upper() not in EQUAL_FREQUENCY_MODELS:
            bind(
                st.radio,
                "substitution.frequencies",
                "Equilibrium frequencies",
                analysis,
                options=["estimate", "equal"],
                horizontal=True,
            )

    st.subheader("Rate heterogeneity")
    choice = _rate_model(analysis)
    if choice != "none":
        columns = st.columns(3)
        with columns[0]:
            bind(
                st.number_input,
                "site.categories",
                "Categories",
                analysis,
                min_value=2,
                step=1,
            )
        if choice != "cat":
            with columns[1]:
                bind(st.checkbox, "site.invariant", "+ invariant sites (+I)", analysis)
        if choice in ("gamma", "weibull", "lognormal"):
            with columns[2]:
                bind(
                    st.selectbox,
                    "site.quadrature",
                    "Quadrature",
                    analysis,
                    options=QUADRATURES,
                    format_func=lambda q: "default (median)" if q is None else q,
                )
        if choice == "discrete":
            bind(
                st.selectbox,
                "site.rate_parameterization",
                "Rate parameterization",
                analysis,
                options=list(RATE_PARAMETERIZATIONS),
                help="rate_increments and rate_ratios order the categories; "
                "rate_shape and mean_contribution leave them exchangeable",
            )
            bind(
                st.checkbox,
                "site.free_rate_em",
                "Fit the rates and weights by EM",
                analysis,
            )
        if choice == "cat":
            _cat(analysis)
    else:
        bind(
            st.checkbox,
            "site.invariant",
            "Proportion of invariant sites (+I)",
            analysis,
        )
    show_problems(analysis, "site.categories")

    init_table(
        analysis,
        (
            ROLE_FREQUENCIES,
            ROLE_RATES,
            ROLE_KAPPA,
            ROLE_SHAPE,
            ROLE_PROPORTION_INVARIANT,
            ROLE_SITE_RATES,
            ROLE_SITE_PROPORTIONS,
            ROLE_TOP_RATE,
        ),
    )


def _cat(analysis: Analysis) -> None:
    columns = st.columns(2)
    with columns[0]:
        bind(
            st.selectbox,
            "site.cat.assignment",
            "Assignment",
            analysis,
            options=list(CAT_ASSIGNMENTS),
            help="argmax is FastTree's rule; posterior_mean is more accurate on "
            "branch lengths",
        )
        bind(
            st.checkbox,
            "site.cat.mixture",
            "Report the mixture log-likelihood",
            analysis,
        )
    with columns[1]:
        optional_text(
            "site.cat.prior",
            "Prior over the rate grid",
            analysis,
            placeholder="physher's default",
            help="'npmle' to estimate it from the alignment, or the shape of a "
            "mean-one Gamma",
        )
        bind(
            st.number_input,
            "site.cat.probe",
            "Probe grid size",
            analysis,
            value=None,
            min_value=1,
            step=1,
            placeholder="physher's default",
        )
