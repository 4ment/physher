"""Build physher JSON configuration files.

The package is layered so that front ends other than the CLI can reuse it::

    phyforge.spec       what the user wants, as plain dataclasses
    phyforge.validate   checks over a spec, reported per field
    phyforge.io         reading alignments, trees and dates
    phyforge.build      spec -> physher JSON
    phyforge.cli        argv -> spec, and printing the result

A front end only needs the first four::

    from phyforge import AnalysisSpec, build
    document = build(AnalysisSpec(method="ml", ...))

``describe(spec)`` lists the parameters of the analysis -- their starting
values, priors and the rest -- which ``AnalysisSpec.parameters`` overrides by
id.
"""

from .build import build, describe, make_context
from .spec import AnalysisSpec, ParameterOverride, PriorChoice
from .validate import SpecError, check, validate

__version__ = "0.1.0"

__all__ = [
    "AnalysisSpec",
    "ParameterOverride",
    "PriorChoice",
    "SpecError",
    "build",
    "check",
    "describe",
    "make_context",
    "validate",
    "__version__",
]
