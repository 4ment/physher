"""MCMC operators, one set per free parameter.

Operators attach to the *constrained* parameter even when it is transformed:
physher pushes the proposed value back through the transform (this is what the
time-tree examples do with ``&root.height``).
"""

from __future__ import annotations

from typing import List

from ..nodes import param_ref
from .context import SIMPLEX, UNIT, BuildContext, FreeParameter


def build(context: BuildContext) -> List[dict]:
    operators: List[dict] = []
    for free in context.free_parameters:
        operators.extend(_operators_for(free))
    return operators


def _operators_for(free: FreeParameter) -> List[dict]:
    x = param_ref(free.id)
    # An operator is picked in proportion to its weight, and each proposal
    # moves one parameter, so a scalar and a vector of branch lengths sharing a
    # weight would give every branch length 1/N of the attention the scalar
    # gets. Weighting by the dimension the operator moves keeps the per-element
    # proposal rate the same across the whole model.
    weight = free.dimension
    if free.kind == SIMPLEX:
        return [_operator("dirichlet", free.id, x, weight)]
    if free.kind == UNIT:
        return [_operator("beta", free.id, x, weight)]
    # A scaler moves the magnitude, a slider the location; positive and real
    # parameters both benefit from having each available.
    return [
        _operator("scaler", free.id, x, weight),
        _operator("slider", free.id, x, weight),
    ]


def _operator(algorithm: str, id_: str, x: str, weight: int = 1) -> dict:
    return {
        "id": f"op.{algorithm}.{id_}",
        "type": "operator",
        "algorithm": algorithm,
        "x": x,
        "weight": weight,
    }
