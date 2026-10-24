"""Primitive physher JSON nodes.

Every helper here returns a plain ``dict``, so the builder layer above can be
exercised without a physher binary and reused from a GUI.

physher refers to previously declared objects by id with a sigil prefix:
``@id`` for a Model and ``&id`` for a Parameter.  A parameter's value lives
under ``"x"`` -- a number, a list of numbers, or a nested parameter node, the
last of which makes it a *transformed* parameter whose free value is the nested
one (see ``transforms.c``).
"""

from __future__ import annotations

from typing import Any, Optional, Sequence, Union

Number = Union[int, float]
Value = Union[Number, Sequence[Number], dict]

__all__ = [
    "model_ref",
    "param_ref",
    "parameter",
    "transformed_parameter",
    "simplex",
    "compound",
    "distribution",
]


def model_ref(id_: str) -> str:
    """Reference to a Model declared elsewhere in the document."""
    return f"@{id_}"


def param_ref(id_: str) -> str:
    """Reference to a Parameter declared elsewhere in the document."""
    return f"&{id_}"


def parameter(
    id_: str,
    x: Optional[Value] = None,
    *,
    dimension: Optional[int] = None,
    lower: Optional[Union[Number, str]] = None,
    upper: Optional[Union[Number, str]] = None,
) -> dict:
    """A ``"type": "parameter"`` node.

    ``x`` may be a scalar, a list, or another parameter node.  When ``x`` is a
    scalar and ``dimension`` is given, physher broadcasts the scalar across the
    whole vector.
    """
    node: dict = {"id": id_, "type": "parameter"}
    if dimension is not None:
        node["dimension"] = dimension
    if lower is not None:
        node["lower"] = lower
    if upper is not None:
        node["upper"] = upper
    if x is not None:
        node["x"] = x
    return node


def transformed_parameter(
    id_: str,
    unconstrained_id: str,
    unconstrained_x: Value,
    *,
    dimension: Optional[int] = None,
    lower: Optional[Union[Number, str]] = None,
    upper: Optional[Union[Number, str]] = None,
) -> dict:
    """A constrained parameter backed by an unconstrained one.

    The bounds on the outer node select the transform: ``lower`` alone gives
    the log/exp map, ``lower`` with ``upper`` the logit map.  Optimizers and
    variational distributions attach to the *inner* id, which is the leaf of
    the autodiff graph.
    """
    return parameter(
        id_,
        parameter(unconstrained_id, unconstrained_x, dimension=dimension),
        lower=lower,
        upper=upper,
    )


def simplex(
    id_: str, x: Optional[Value] = None, *, dimension: Optional[int] = None
) -> dict:
    """A ``"type": "Simplex"`` node.

    Pass ``x`` as a list to fix the simplex at those values, or as a nested
    parameter node of dimension ``K - 1`` to estimate it.  With neither, a
    ``dimension`` of ``K`` starts from equal values.
    """
    node: dict = {"id": id_, "type": "Simplex"}
    if dimension is not None:
        node["dimension"] = dimension
    if x is not None:
        node["x"] = x
    return node


def compound(id_: str, distributions: Sequence[Any]) -> dict:
    """A ``"type": "compound"`` node summing the log densities of its members."""
    return {"id": id_, "type": "compound", "distributions": list(distributions)}


def distribution(
    id_: str,
    name: str,
    x: Any,
    *,
    parameters: Optional[dict] = None,
    **extra: Any,
) -> dict:
    """A ``"type": "distribution"`` node evaluating ``name`` at ``x``."""
    node: dict = {"id": id_, "type": "distribution", "distribution": name, "x": x}
    if parameters is not None:
        node["parameters"] = parameters
    node.update(extra)
    return node
