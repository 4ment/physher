"""Tab 1 -- Alignment."""

from __future__ import annotations

from typing import List

import streamlit as st
import streamlit.components.v1 as components

from ...io.fasta import Sequence, read_fasta
from .. import state, uploads
from ..state import Analysis
from ..viewers import msaviewer_html
from ..widgets import bind_session, show_problems

DATATYPES = {"nucleotide": "Nucleotide", "aa": "Amino acid"}

_NUCLEOTIDE_CHARACTERS = set("ACGTUN-?RYKMSWBDHV.")


def guess_datatype(sequences: List[Sequence]) -> str:
    """Nucleotide when nearly every character is a nucleotide or ambiguity code."""
    text = "".join(sequence.sequence for sequence in sequences[:50]).upper()
    if not text:
        return "nucleotide"
    nucleotide = sum(char in _NUCLEOTIDE_CHARACTERS for char in text)
    return "nucleotide" if nucleotide / len(text) > 0.9 else "aa"


@st.dialog("Alignment viewer", width="large")
def _view(sequences: List[Sequence], name: str, protein: bool) -> None:
    sites = max(len(s.sequence) for s in sequences)
    st.caption(f"`{name}` — {len(sequences)} taxa × {sites} sites")
    height = 600
    components.html(
        msaviewer_html(
            [(s.taxon, s.sequence) for s in sequences], protein=protein, height=height
        ),
        height=height,
        scrolling=True,
    )


def render(analysis: Analysis) -> None:
    st.header("Alignment")
    uploads.upload_widget(
        "alignment", "Upload FASTA alignment", help="Sequences in FASTA format"
    )
    show_problems(analysis, "data.alignment")

    entry = uploads.info("alignment")
    if entry is None:
        st.info("No alignment loaded yet.")
        return
    try:
        sequences = read_fasta(entry["path"])
    except (OSError, ValueError) as error:
        st.error(f"Could not read {entry['name']}: {error}")
        return

    guessed = st.session_state.setdefault(
        f"phyforge.guessed:{entry['path']}", guess_datatype(sequences)
    )
    if st.session_state.get("phyforge.datatype_for") != entry["path"]:
        # A new alignment: start from what it looks like.
        st.session_state[state.DATATYPE] = guessed
        st.session_state.pop(f"w:@{state.DATATYPE}", None)
        st.session_state["phyforge.datatype_for"] = entry["path"]

    sites = {len(sequence.sequence) for sequence in sequences}
    header = st.columns([3, 3, 1, 1, 2, 1, 1])
    for column, title in zip(header, ["Name", "File", "Taxa", "Sites", "Type"]):
        column.markdown(f"**{title}**")
    row = st.columns([3, 3, 1, 1, 2, 1, 1])
    row[0].write(entry["name"].rsplit("/", 1)[-1].rsplit(".", 1)[0])
    row[1].write(entry["name"])
    row[2].write(len(sequences))
    row[3].write(max(sites))
    bind_session(
        st.selectbox,
        state.DATATYPE,
        "Data type",
        container=row[4],
        options=list(DATATYPES),
        format_func=DATATYPES.get,
        label_visibility="collapsed",
    )
    if row[5].button("👁", key="view_alignment", help="View alignment"):
        _view(sequences, entry["name"], st.session_state[state.DATATYPE] == "aa")
    if row[6].button("🗑️", key="delete_alignment", help="Remove alignment"):
        uploads.forget("alignment")
        st.rerun()

    if len(sites) > 1:
        st.warning(
            f"The sequences are not all the same length ({min(sites)}–{max(sites)} "
            "sites); physher expects an aligned FASTA."
        )
