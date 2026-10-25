"""``phyforge gui``: start the Streamlit app."""

from __future__ import annotations

import os
import sys
from typing import Sequence


def launch(argv: Sequence[str] = ()) -> int:
    """Run the app, passing ``argv`` on to ``streamlit run``."""
    try:
        from streamlit.web import cli as streamlit_cli
    except ImportError:
        print(
            "phyforge: the GUI needs the optional dependencies; install them "
            "with  pip install 'phyforge[gui]'  (or  uv sync --extra gui)",
            file=sys.stderr,
        )
        return 1
    app = os.path.join(os.path.dirname(os.path.abspath(__file__)), "app.py")
    argv = list(argv)
    if not any(arg.startswith("--server.address") for arg in argv):
        # Streamlit listens on every interface by default; the data loaded here
        # is often unpublished, so stay on this machine unless asked otherwise.
        argv = ["--server.address", "localhost", *argv]
    sys.argv = ["streamlit", "run", app, *argv]
    return streamlit_cli.main()
