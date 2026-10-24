"""Creating free parameters: the JSON node and its registry entry together.

Whether a parameter is emitted constrained or backed by an unconstrained one
is a build-wide decision (gradient-based methods need transforms everywhere),
so it is made here rather than repeated in every model builder.
"""

from __future__ import annotations

from typing import Optional, Tuple

from ..nodes import parameter, simplex, transformed_parameter
from .context import (
    POSITIVE,
    REAL,
    SIMPLEX,
    UNIT,
    BuildContext,
    FreeParameter,
    to_unconstrained,
)

__all__ = ["free_parameter", "free_simplex", "fixed_parameter"]

_BOUNDS = {
    REAL: (None, None),
    POSITIVE: (0, None),
    UNIT: (0, 1),
}


def free_parameter(
    context: BuildContext,
    id_: str,
    value: float,
    *,
    role: str,
    kind: str = POSITIVE,
    dimension: Optional[int] = None,
    transformed: Optional[bool] = None,
    log: bool = True,
    bounds: Optional[Tuple[Optional[float], Optional[float]]] = None,
) -> dict:
    """Build a parameter node and register it as free.

    ``transformed`` overrides the spec-wide choice; a time tree's height
    parameters use it because they are always reparameterized.  ``bounds``
    narrows the box the ``kind`` implies, for parameters a line search should
    not be allowed to push all the way onto their boundary.
    """
    if transformed is None:
        transformed = context.spec.uses_unconstrained and kind != REAL
    lower, upper = bounds if bounds is not None else _BOUNDS[kind]

    if not transformed:
        node = parameter(id_, value, dimension=dimension, lower=lower, upper=upper)
        context.register(FreeParameter(id_, role, kind, dimension or 1, log=log))
        return node

    unconstrained_id = context.ids.unconstrained(id_)
    unconstrained_value = to_unconstrained(value, kind)
    node = transformed_parameter(
        id_,
        unconstrained_id,
        unconstrained_value,
        dimension=dimension,
        lower=lower,
        upper=upper,
    )
    context.register(
        FreeParameter(
            id_,
            role,
            kind,
            dimension or 1,
            unconstrained_id,
            unconstrained_value,
            log=log,
        )
    )
    return node


def free_simplex(
    context: BuildContext,
    id_: str,
    dimension: int,
    *,
    role: str,
    transformed: Optional[bool] = None,
    log: bool = True,
) -> dict:
    """Build a free simplex of ``dimension`` elements.

    Transformed, it carries an unconstrained parameter of ``dimension - 1``
    starting at zero, which the stick-breaking map sends to equal proportions.
    """
    if transformed is None:
        # MCMC moves a simplex directly with a Dirichlet operator; everything
        # else needs the unconstrained chart.
        transformed = context.spec.method != "mcmc"

    if not transformed:
        node = simplex(id_, dimension=dimension)
        context.register(FreeParameter(id_, role, SIMPLEX, dimension, log=log))
        return node

    unconstrained_id = context.ids.unconstrained(id_)
    node = simplex(
        id_,
        parameter(unconstrained_id, 0.0, dimension=dimension - 1),
    )
    context.register(
        FreeParameter(id_, role, SIMPLEX, dimension, unconstrained_id, 0.0, log=log)
    )
    return node


def fixed_parameter(id_: str, value, *, dimension: Optional[int] = None) -> dict:
    """A parameter that is never registered, and so is never estimated."""
    return parameter(id_, value, dimension=dimension)
