"""Creating free parameters: the JSON node and its registry entry together.

Whether a parameter is emitted constrained or backed by an unconstrained one
is a build-wide decision (gradient-based methods need transforms everywhere),
so it is made here rather than repeated in every model builder.  So is applying
the user's per-parameter overrides: a model builder passes its default starting
value and gets back the node the override asks for.

A parameter an override fixes keeps the node shape it would have had -- the
models reading it do not care whether anything moves it -- but is registered as
fixed, so no prior, optimizer, operator or variational factor is built for it.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

from ..nodes import parameter, simplex, transformed_parameter
from .context import (
    POSITIVE,
    REAL,
    SIMPLEX,
    UNIT,
    BuildContext,
    FreeParameter,
    Values,
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
    override = context.override(id_)
    default = value
    value = _resolve_init(
        context, id_, value, override.init, dimension or 1, lower, upper
    )
    known = dict(init=value, default_init=default, estimated=override.estimate)

    if not transformed:
        node = parameter(id_, value, dimension=dimension, lower=lower, upper=upper)
        context.register(
            FreeParameter(
                id_,
                role,
                kind,
                dimension or 1,
                # Untransformed, a real parameter is its own unconstrained value.
                unconstrained_init=value if kind == REAL else 0.0,
                log=log,
                **known,
            )
        )
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
            **known,
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
    override = context.override(id_)
    default = [1.0 / dimension] * dimension
    # None leaves physher's own start, equal proportions, which the node
    # spells as no value at all.
    init = _resolve_simplex_init(context, id_, override.init, dimension)
    known = dict(
        init=init or default, default_init=default, estimated=override.estimate
    )

    if not transformed:
        node = simplex(id_, init, dimension=dimension)
        context.register(FreeParameter(id_, role, SIMPLEX, dimension, log=log, **known))
        return node

    unconstrained_id = context.ids.unconstrained(id_)
    unconstrained_value = 0.0 if init is None else to_unconstrained(init, SIMPLEX)
    node = simplex(
        id_,
        parameter(unconstrained_id, unconstrained_value, dimension=dimension - 1),
    )
    context.register(
        FreeParameter(
            id_,
            role,
            SIMPLEX,
            dimension,
            unconstrained_id,
            unconstrained_value,
            log=log,
            **known,
        )
    )
    return node


def _resolve_init(
    context: BuildContext,
    id_: str,
    default: Values,
    init: Optional[Values],
    dimension: int,
    lower: Optional[float],
    upper: Optional[float],
) -> Values:
    """The override's starting value, or ``default`` if it has none or a bad one.

    A bad one is recorded as a problem rather than raised, so a caller hears
    about every override at once.
    """
    if init is None:
        return default
    values = _as_list(init)
    if len(values) not in (1, dimension):
        context.problem(
            f"{id_} has {dimension} element(s), so it takes 1 or {dimension} "
            f"starting value(s), not {len(values)}",
            id_,
            "init",
        )
        return default
    for value in values:
        if (lower is not None and value <= lower) or (
            upper is not None and value >= upper
        ):
            context.problem(
                f"{id_} must start inside ({_bound(lower, '-inf')}, "
                f"{_bound(upper, 'inf')}), not at {value}",
                id_,
                "init",
            )
            return default
    return values[0] if len(values) == 1 else values


def _resolve_simplex_init(
    context: BuildContext, id_: str, init: Optional[Values], dimension: int
) -> Optional[List[float]]:
    """The override's starting proportions, normalized, or None."""
    if init is None:
        return None
    values = _as_list(init)
    if len(values) != dimension:
        context.problem(
            f"{id_} is a simplex of {dimension} proportions, so it takes "
            f"{dimension} starting values, not {len(values)}",
            id_,
            "init",
        )
        return None
    total = sum(values)
    if any(value <= 0 for value in values) or abs(total - 1.0) > 1e-6:
        context.problem(
            f"{id_} must start at positive proportions summing to 1",
            id_,
            "init",
        )
        return None
    return [value / total for value in values]


def _as_list(values: Values) -> List[float]:
    return list(values) if isinstance(values, (list, tuple)) else [values]


def _bound(value: Optional[float], infinite: str) -> str:
    return infinite if value is None else f"{value:g}"


def fixed_parameter(id_: str, value, *, dimension: Optional[int] = None) -> dict:
    """A parameter that is never registered, and so is never estimated."""
    return parameter(id_, value, dimension=dimension)
