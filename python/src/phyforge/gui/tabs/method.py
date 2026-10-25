"""Tab 7 -- Algorithm: how the model is fitted."""

from __future__ import annotations

import streamlit as st

from ...spec import GRADIENT_ALGORITHMS, LADDERS, MARGINAL_POSTERIOR, MARGINAL_TEMPERED
from ..state import Analysis
from ..widgets import (
    bind,
    optional_text,
    show_problems,
    variational_table,
    weight_table,
)

ALGORITHMS = {
    "meta": "meta (cycle over per-parameter optimizers)",
    "bfgs": "BFGS (gradient)",
    "cg": "conjugate gradient",
    "sg": "stochastic gradient",
}

MARGINAL_LABELS = {
    "ss": "stepping stone",
    "ps": "path sampling",
    "ps2": "modified path sampling",
    "hm": "harmonic mean",
    "shm": "stabilized harmonic mean",
    "is": "importance sampling",
}


def render(analysis: Analysis) -> None:
    st.header("Algorithm")
    method = analysis.spec.method
    st.caption(f"Method **{method}** — change it in the sidebar.")
    if method in ("ml", "map"):
        _optimizer(analysis)
    elif method == "mcmc":
        _mcmc(analysis)
    else:
        _advi(analysis)


def _optimizer(analysis: Analysis) -> None:
    columns = st.columns(3)
    with columns[0]:
        bind(
            st.selectbox,
            "optimizer.algorithm",
            "Optimizer",
            analysis,
            options=list(ALGORITHMS),
            format_func=ALGORITHMS.get,
        )
    with columns[1]:
        bind(
            st.number_input, "optimizer.iterations", "Iterations", analysis, min_value=1
        )
    with columns[2]:
        bind(
            st.number_input,
            "optimizer.precision",
            "Precision",
            analysis,
            min_value=0.0,
            format="%g",
        )
    if analysis.spec.optimizer.algorithm in GRADIENT_ALGORITHMS:
        st.caption(
            "Gradient optimizers move every parameter on the unconstrained scale."
        )

    with st.expander("Bootstrap", expanded=analysis.spec.bootstrap.enabled):
        columns = st.columns(2)
        with columns[0]:
            bind(
                st.number_input,
                "bootstrap.replicates",
                "Replicates",
                analysis,
                min_value=0,
            )
            bind(st.number_input, "bootstrap.threads", "Threads", analysis, min_value=1)
        with columns[1]:
            bind(
                st.checkbox,
                "bootstrap.include_observed",
                "Prepend a fit on the observed data",
                analysis,
            )
            bind(
                st.checkbox,
                "bootstrap.reset",
                "Restart every replicate from the estimate",
                analysis,
            )
            bind(
                st.checkbox,
                "bootstrap.compact",
                "Drop the patterns a replicate did not draw",
                analysis,
            )
        show_problems(analysis, "bootstrap.replicates")
        show_problems(analysis, "bootstrap.threads")


def _mcmc(analysis: Analysis) -> None:
    columns = st.columns(2)
    with columns[0]:
        bind(
            st.number_input,
            "mcmc.length",
            "Chain length",
            analysis,
            min_value=1,
            step=1000,
        )
    with columns[1]:
        bind(
            st.number_input, "mcmc.every", "Log every", analysis, min_value=1, step=100
        )

    with st.expander("Predictive scores (pattern trace)"):
        columns = st.columns(2)
        with columns[0]:
            bind(st.checkbox, "mcmc.cpo", "Conditional predictive ordinate", analysis)
            bind(
                st.number_input,
                "mcmc.cpo_burnin",
                "CPO burn-in (logged rows)",
                analysis,
                value=None,
                min_value=0,
                placeholder="a tenth",
            )
        with columns[1]:
            bind(st.checkbox, "mcmc.lppd", "lppd and WAIC", analysis)
            bind(
                st.number_input,
                "mcmc.lppd_burnin",
                "lppd burn-in (logged rows)",
                analysis,
                value=None,
                min_value=0,
                placeholder="a tenth",
            )

    with st.expander("Marginal likelihood"):
        options = list(MARGINAL_TEMPERED + MARGINAL_POSTERIOR)
        _marginal(analysis, options)
        tempered = set(analysis.spec.marginal.algorithms) & set(MARGINAL_TEMPERED)
        if tempered:
            columns = st.columns(4)
            with columns[0]:
                bind(st.number_input, "marginal.steps", "Steps", analysis, min_value=0)
            with columns[1]:
                bind(
                    st.selectbox,
                    "marginal.ladder",
                    "Ladder",
                    analysis,
                    options=list(LADDERS),
                )
            if analysis.spec.marginal.ladder == "beta":
                with columns[2]:
                    bind(
                        st.number_input,
                        "marginal.ladder_alpha",
                        "α",
                        analysis,
                        format="%g",
                    )
                with columns[3]:
                    bind(
                        st.number_input,
                        "marginal.ladder_beta",
                        "β",
                        analysis,
                        format="%g",
                    )
            show_problems(analysis, "marginal.steps")
        if analysis.spec.marginal.algorithms:
            bind(
                st.number_input,
                "marginal.burnin",
                "Burn-in (logged rows)",
                analysis,
                value=None,
                min_value=0,
                placeholder="a tenth",
            )

    st.subheader("Operator weights")
    weight_table(analysis)


def _advi(analysis: Analysis) -> None:
    columns = st.columns(3)
    with columns[0]:
        bind(
            st.number_input,
            "variational.iterations",
            "Iterations",
            analysis,
            min_value=1,
        )
        bind(st.number_input, "variational.eta", "Learning rate", analysis, format="%g")
    with columns[1]:
        bind(
            st.number_input,
            "variational.elbo_samples",
            "ELBO samples",
            analysis,
            min_value=1,
        )
        bind(
            st.number_input,
            "variational.grad_samples",
            "Gradient samples",
            analysis,
            min_value=1,
        )
    with columns[2]:
        bind(st.number_input, "variational.tol", "Tolerance", analysis, format="%g")
        bind(
            st.number_input,
            "variational.samples",
            "Samples drawn at the end",
            analysis,
            min_value=0,
        )
    with st.expander("Advanced"):
        columns = st.columns(3)
        with columns[0]:
            bind(st.text_input, "variational.update", "Update rule", analysis)
            bind(st.checkbox, "variational.entropy", "Use the entropy", analysis)
        with columns[1]:
            bind(
                st.number_input,
                "variational.sigma_init",
                "Initial σ of every factor",
                analysis,
                min_value=0.0,
                format="%g",
            )
        with columns[2]:
            optional_text("variational.checkpoint", "Checkpoint file", analysis)

    with st.expander("Marginal likelihood"):
        _marginal(analysis, ["is"])
        if "is" in analysis.spec.marginal.algorithms:
            bind(st.number_input, "marginal.samples", "Draws", analysis, min_value=1)

    st.subheader("Variational approximation")
    variational_table(analysis)


def _marginal(analysis: Analysis, options) -> None:
    bind(
        st.multiselect,
        "marginal.algorithms",
        "Estimators",
        analysis,
        options=options,
        format_func=lambda a: f"{a} — {MARGINAL_LABELS[a]}",
        encode=lambda chosen: [a for a in chosen if a in options],
        decode=tuple,
    )
