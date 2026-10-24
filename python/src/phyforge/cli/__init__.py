"""Command line front end."""

from .args import make_parser, parse, to_spec
from .main import main

__all__ = ["main", "make_parser", "parse", "to_spec"]
