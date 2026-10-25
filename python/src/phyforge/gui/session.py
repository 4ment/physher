"""Saving and loading a GUI session.

A session file holds the settings, the overrides and the uploaded files
themselves, so it reopens on another machine; files referenced by path are
stored as the path.
"""

from __future__ import annotations

import dataclasses
import json
from typing import Any, Dict

import streamlit as st

from ..spec import ParameterOverride
from . import state, uploads

FORMAT = "phyforge-session"
VERSION = 1


def dump() -> str:
    files: Dict[str, Any] = {}
    for role, entry in state.files().items():
        if entry["upload"]:
            files[role] = {"name": entry["name"], "content": uploads.read_text(role)}
        else:
            files[role] = {"path": entry["path"]}
    return json.dumps(
        {
            "format": FORMAT,
            "version": VERSION,
            "values": state.values(),
            "overrides": {
                id_: dataclasses.asdict(o) for id_, o in state.overrides().items()
            },
            "time_tree": state.is_time_tree(),
            "datatype": st.session_state[state.DATATYPE],
            "files": files,
        },
        indent=2,
    )


def load(text: str) -> None:
    data = json.loads(text)
    if data.get("format") != FORMAT:
        raise ValueError("not a phyforge session file")
    values = state.values()
    for path, value in data["values"].items():
        if path in values:
            values[path] = tuple(value) if path == "marginal.algorithms" else value
    st.session_state[state.OVERRIDES] = {
        id_: ParameterOverride.from_dict(o) for id_, o in data["overrides"].items()
    }
    st.session_state[state.TIME_TREE] = bool(data.get("time_tree"))
    st.session_state[state.DATATYPE] = data.get("datatype", "nucleotide")
    st.session_state[state.FILES] = {}
    for role in uploads.PATHS:
        values[uploads.PATHS[role]] = None
    for role, entry in data.get("files", {}).items():
        if "content" in entry and entry["content"] is not None:
            uploads.store(role, entry["name"], entry["content"].encode())
        elif "path" in entry:
            uploads.use_local(role, entry["path"])
    state.reset_widgets()
