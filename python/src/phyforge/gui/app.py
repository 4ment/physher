"""phyforge -- BEAUti-style configuration builder for physher.

Run with ``phyforge gui``, or ``streamlit run`` on this file.
"""

from __future__ import annotations

import os
import sys

# Streamlit runs this file as __main__, outside the package.
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import streamlit as st  # noqa: E402

from phyforge.gui import session, state  # noqa: E402
from phyforge.gui.tabs import (  # noqa: E402
    alignment,
    clock,
    dates,
    method,
    output,
    priors,
    site_model,
    tree,
)
from phyforge.gui.widgets import bind  # noqa: E402
from phyforge.spec import METHODS  # noqa: E402

METHOD_LABELS = {
    "ml": "ml — maximum likelihood",
    "map": "map — maximum a posteriori",
    "mcmc": "mcmc — Metropolis-Hastings",
    "advi": "advi — variational inference",
}

st.set_page_config(
    page_title="phyforge",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

state.init_state()
analysis = state.analyse()

with st.sidebar:
    st.header("Analysis")
    bind(
        st.radio,
        "method",
        "Method",
        analysis,
        options=list(METHODS),
        format_func=METHOD_LABELS.get,
    )
    bind(
        st.text_input,
        "stem",
        "Output stem",
        analysis,
        help="stem of the files physher writes, and of the configuration",
    )
    st.divider()
    st.subheader("Session")
    st.download_button(
        "Save session",
        data=session.dump(),
        file_name="phyforge-session.json",
        mime="application/json",
        width="stretch",
    )
    st.session_state.setdefault("phyforge.session_counter", 0)
    loaded = st.file_uploader(
        "Load session",
        type=["json"],
        key=f"session:{st.session_state['phyforge.session_counter']}",
    )
    if loaded is not None:
        try:
            session.load(loaded.getvalue().decode())
        except (ValueError, KeyError, TypeError) as error:
            st.error(f"Could not load the session: {error}")
        else:
            st.session_state["phyforge.session_counter"] += 1
            st.rerun()

st.title("phyforge")
st.caption(
    "Graphical configuration builder for physher — generates a valid JSON "
    "config, and the phyforge command that writes the same one."
)

tabs = st.tabs(
    [
        "1 · Alignment",
        "2 · Dates",
        "3 · Tree",
        "4 · Site Model",
        "5 · Clock Model",
        "6 · Priors",
        "7 · Algorithm",
        "8 · Output",
    ]
)
for tab, module in zip(
    tabs, (alignment, dates, tree, site_model, clock, priors, method, output)
):
    with tab:
        module.render(analysis)
