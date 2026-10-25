"""Tab 8 -- Output: the JSON configuration and the command that writes it."""

from __future__ import annotations

import dataclasses
import json
from typing import Any, Dict

import streamlit as st

from ...build import build
from ...cli.argv import to_command
from ...validate import SpecError
from .. import state, uploads
from ..state import Analysis
from ..widgets import bind


def displayed(spec, analysis: Analysis):
    """The spec as the user's files are called outside this session."""
    data = dataclasses.replace(
        spec.data,
        alignment=uploads.display_name("alignment") or spec.data.alignment,
        tree=uploads.display_name("tree") or spec.data.tree,
    )
    dates = dataclasses.replace(
        spec.clock.dates, path=uploads.display_name("dates") or spec.clock.dates.path
    )
    clock = dataclasses.replace(spec.clock, dates=dates)
    return dataclasses.replace(spec, data=data, clock=clock)


def rename_files(document: Any, names: Dict[str, str]) -> Any:
    """Replace every temporary upload path in ``document`` by the file's name."""
    if isinstance(document, dict):
        return {key: rename_files(value, names) for key, value in document.items()}
    if isinstance(document, list):
        return [rename_files(value, names) for value in document]
    if isinstance(document, str):
        return names.get(document, document)
    return document


def render(analysis: Analysis) -> None:
    st.header("JSON Configuration")
    if uploads.info("alignment") is None or uploads.info("tree") is None:
        st.info("Load an alignment (Tab 1) and a tree (Tab 3) first.")
        return

    columns = st.columns(2)
    with columns[0]:
        bind(
            st.checkbox,
            "data.embed_alignment",
            "Embed the alignment in the JSON",
            analysis,
        )
    with columns[1]:
        bind(st.checkbox, "data.embed_tree", "Embed the tree in the JSON", analysis)
    names = {
        entry["path"]: entry["name"]
        for entry in state.files().values()
        if entry["upload"]
    }
    if names and not (
        analysis.spec.data.embed_alignment and analysis.spec.data.embed_tree
    ):
        st.caption(
            "Files that are not embedded are referenced by the name they were "
            "uploaded under: keep them next to the configuration when running physher."
        )

    _inactive(analysis)

    shown = [p for p in analysis.problems if (p.field or "") not in analysis.shown]
    if analysis.problems:
        st.error(
            "The configuration cannot be built yet:\n\n"
            + "\n".join(f"- {problem}" for problem in analysis.problems),
            icon="⚠️",
        )
        if shown:
            st.caption("Problems not shown next to a setting are listed above too.")
        return

    try:
        document = rename_files(build(analysis.spec), names)
    except (SpecError, ValueError, KeyError, OSError) as error:
        st.error(f"Could not build configuration: {error}")
        return

    stem = analysis.spec.stem or "physher"
    filename = f"{stem}.json"
    text = json.dumps(document, indent=2)
    st.success("Configuration ready.")
    buttons = st.columns(2)
    buttons[0].download_button(
        f"Download {filename}",
        data=text + "\n",
        file_name=filename,
        mime="application/json",
        type="primary",
    )
    if analysis.spec.parameters:
        buttons[1].download_button(
            "Download parameters.json",
            data=json.dumps(analysis.spec.to_dict()["parameters"], indent=2) + "\n",
            file_name="parameters.json",
            mime="application/json",
            help="The per-parameter overrides, for phyforge --params",
        )

    st.subheader("Command line")
    st.caption("The same configuration, written by the phyforge CLI:")
    st.code(
        to_command(displayed(analysis.spec, analysis), output=filename), language="bash"
    )
    st.code(f"physher {filename}", language="bash")

    with st.expander("JSON preview", expanded=False):
        st.json(document, expanded=2)

    with st.expander("Component summary", expanded=False):
        import pandas as pd

        rows = [
            {"id": key, "type": node.get("type", "—")}
            for key, node in document.items()
            if isinstance(node, dict)
        ]
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)


def _inactive(analysis: Analysis) -> None:
    if not analysis.inactive:
        return
    with st.expander(
        f"{len(analysis.inactive)} kept override(s) not used by this model"
    ):
        st.caption(
            "Settings for parameters the current model does not have. They come "
            "back into effect if you switch back to a model that has them."
        )
        for id_ in analysis.inactive:
            columns = st.columns([3, 1])
            columns[0].markdown(f"`{id_}`")
            if columns[1].button("Discard", key=f"discard:{id_}"):
                state.overrides().pop(id_, None)
                st.rerun()
