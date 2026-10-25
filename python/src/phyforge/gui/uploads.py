"""Input files: uploads written to disk, or paths on the machine running the GUI.

The spec, and physher after it, refers to files by path, while an upload is a
stream of bytes; each upload is written to a temporary directory that belongs
to the session, under its own name.  The document is then either made
self-contained (the alignment and tree embedded) or pointed at the file name
the user uploaded, which is where it will sit next to the configuration --
never at the temporary path, which means nothing outside this session.
"""

from __future__ import annotations

import os
import tempfile
from typing import Any, Dict, Optional

import streamlit as st

from . import state

DIRECTORY = "phyforge.directory"

#: The spec field each kind of input file fills in.
PATHS = {
    "alignment": "data.alignment",
    "tree": "data.tree",
    "dates": "clock.dates.path",
}


def directory() -> str:
    if DIRECTORY not in st.session_state:
        st.session_state[DIRECTORY] = tempfile.mkdtemp(prefix="phyforge-")
    return st.session_state[DIRECTORY]


def store(role: str, name: str, data: bytes) -> str:
    """Write an upload for ``role`` and point the spec at it."""
    folder = os.path.join(directory(), role)
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, os.path.basename(name))
    with open(path, "wb") as fp:
        fp.write(data)
    state.files()[role] = {"name": os.path.basename(name), "path": path, "upload": True}
    state.values()[PATHS[role]] = path
    return path


def use_local(role: str, path: str) -> None:
    """Point the spec at a file already on disk, referenced as given."""
    state.files()[role] = {"name": path, "path": path, "upload": False}
    state.values()[PATHS[role]] = path


def forget(role: str) -> None:
    state.files().pop(role, None)
    state.values()[PATHS[role]] = None


def info(role: str) -> Optional[Dict[str, Any]]:
    return state.files().get(role)


def read_text(role: str) -> Optional[str]:
    entry = info(role)
    if entry is None or not os.path.exists(entry["path"]):
        return None
    with open(entry["path"], "r", errors="replace") as fp:
        return fp.read()


def display_name(role: str) -> Optional[str]:
    """What the configuration and the command call the file."""
    entry = info(role)
    return None if entry is None else entry["name"]


def upload_widget(role: str, label: str, help: str = "") -> None:
    """An uploader plus a path box; whichever is used last wins."""
    counter = f"phyforge.uploader.{role}"
    st.session_state.setdefault(counter, 0)
    columns = st.columns([3, 2])
    uploaded = columns[0].file_uploader(
        label, key=f"upload:{role}:{st.session_state[counter]}", help=help
    )
    if uploaded is not None:
        store(role, uploaded.name, uploaded.getvalue())
        # A fresh key drops the uploader's own copy of the file, which would
        # otherwise be stored again on every rerun; a stale path would win.
        st.session_state[counter] += 1
        st.session_state[f"local:{role}"] = ""
        st.rerun()
    path = columns[1].text_input(
        "…or a path on this machine",
        key=f"local:{role}",
        placeholder="/path/to/file",
    )
    if path and (info(role) is None or info(role)["path"] != path):
        if os.path.exists(path):
            use_local(role, path)
            st.rerun()
        else:
            columns[1].error(f"{path} does not exist", icon="⚠️")
