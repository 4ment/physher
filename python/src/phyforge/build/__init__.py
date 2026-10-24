"""Turning an :class:`~phyforge.spec.AnalysisSpec` into physher JSON.

Every module here is a pure function of the spec plus a
:class:`~phyforge.build.context.BuildContext`: no argparse, no I/O beyond
reading the alignment and tree the spec names.  That is what lets a GUI reuse
the whole layer.
"""

from .analysis import build, make_context
from .context import BuildContext, FreeParameter

__all__ = ["build", "make_context", "BuildContext", "FreeParameter"]
