"""Tab 3 -- Tree."""

from __future__ import annotations

import streamlit as st
import streamlit.components.v1 as components

from ...build.tree import ROLE_BRANCH_LENGTHS, ROLE_PROPORTIONS, ROLE_ROOT_HEIGHT
from ...io.trees import read_newick
from .. import state, uploads
from ..root_to_tip import CRITERIA, root_to_tip
from ..state import Analysis
from ..viewers import phylocanvas_html
from ..widgets import bind, init_table, show_problems
from .dates import tip_dates

#: Where an unrooted tree's starting branch lengths come from, and the two
#: spec fields each choice sets.
STARTS = {
    "parsimony": ("A parsimony fit of the topology", None, False),
    "flat": ("0.1 on every branch", False, False),
    "file": ("The lengths in the tree file", None, True),
}


@st.dialog("Tree viewer", width="large")
def _view(newick: str) -> None:
    height = 600
    components.html(phylocanvas_html(newick, height=height), height=height + 50)


def _start() -> str:
    values = state.values()
    if values["data.keep_branch_lengths"]:
        return "file"
    return "flat" if values["data.init_parsimony"] is False else "parsimony"


def _starting_lengths(analysis: Analysis) -> None:
    key = "w:@start"
    if key not in st.session_state:
        st.session_state[key] = _start()

    def changed() -> None:
        _, parsimony, keep = STARTS[st.session_state[key]]
        state.values()["data.init_parsimony"] = parsimony
        state.values()["data.keep_branch_lengths"] = keep

    choice = st.radio(
        "Starting branch lengths",
        list(STARTS),
        key=key,
        on_change=changed,
        format_func=lambda name: STARTS[name][0],
        horizontal=True,
    )
    show_problems(analysis, "data.init_parsimony")
    if choice == "parsimony":
        bind(
            st.number_input,
            "data.init_parsimony_min",
            "Parsimony floor",
            analysis,
            value=None,
            min_value=0.0,
            format="%g",
            placeholder="physher's default",
            help="Applied to the parsimony estimates, which are otherwise zero on "
            "a branch with no inferred substitution",
        )


def _root_to_tip(newick: str, dates) -> None:
    with st.expander("Root-to-tip regression", expanded=False):
        left, right = st.columns([3, 1])
        criterion = left.selectbox(
            "Optimisation criterion",
            list(CRITERIA),
            format_func=CRITERIA.get,
            key="r2t_criterion",
        )
        right.write("")
        if right.button("Compute", key="r2t_compute", width="stretch"):
            try:
                with st.spinner("Searching for the optimal root…"):
                    st.session_state["r2t_result"] = root_to_tip(
                        newick, dates, criterion
                    )
            except ValueError as error:
                st.session_state.pop("r2t_result", None)
                st.error(f"Root-to-tip failed: {error}")
        result = st.session_state.get("r2t_result")
        if not result:
            st.caption(
                "Hit **Compute** to search every branch for the rooting that "
                "best fits a molecular clock (linear date vs root-to-tip)."
            )
            return

        tmrca = f"{result['tmrca']:.4g}" if result["tmrca"] is not None else "—"
        metrics = st.columns(4)
        metrics[0].metric("Rate (slope)", f"{result['slope']:.4g}")
        metrics[1].metric("TMRCA (x-intercept)", tmrca)
        metrics[2].metric("r²", f"{result['r2']:.4f}")
        metrics[3].metric("n tips", str(result["n"]))

        import plotly.graph_objects as go

        xs, ys = result["dates"], result["root_to_tip"]
        low, high = min(xs), max(xs)
        figure = go.Figure()
        figure.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="markers",
                text=result["leaf_names"],
                hovertemplate="%{text}<br>date %{x}<br>r2t %{y:.5g}<extra></extra>",
                marker=dict(color="#1f77b4", size=8),
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[low, high],
                y=[result["slope"] * x + result["intercept"] for x in (low, high)],
                mode="lines",
                line=dict(color="#d62728"),
            )
        )
        figure.update_layout(
            xaxis_title="Date",
            yaxis_title="Root-to-tip distance",
            showlegend=False,
            height=360,
            margin=dict(l=40, r=20, t=20, b=40),
        )
        st.plotly_chart(figure, width="stretch")


def render(analysis: Analysis) -> None:
    st.header("Tree")
    if uploads.info("alignment") is None:
        st.info("Load an alignment first (Tab 1).")
        return

    uploads.upload_widget("tree", "Upload tree file (Newick or NEXUS)")
    show_problems(analysis, "data.tree")
    entry = uploads.info("tree")
    if entry is None:
        st.info("No tree loaded yet.")
        return
    try:
        newick = read_newick(entry["path"])
    except (OSError, ValueError) as error:
        st.error(f"Could not read {entry['name']}: {error}")
        return

    message, view = st.columns([4, 1])
    message.success(f"Tree loaded from **{entry['name']}**")
    if view.button("👁 View tree", key="view_tree_btn"):
        _view(newick)

    if state.is_time_tree():
        st.caption(
            "A time tree: the topology is kept, and its node heights are "
            "estimated as a root height plus one proportion per internal node."
        )
        dates = tip_dates(analysis)
        if dates:
            _root_to_tip(newick, dates)
    else:
        _starting_lengths(analysis)

    init_table(
        analysis,
        (ROLE_BRANCH_LENGTHS, ROLE_ROOT_HEIGHT, ROLE_PROPORTIONS),
        "Tree parameters",
    )
