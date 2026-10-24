"""The mean-field variational approximation used by ``advi``.

One independent normal per free parameter, placed on the *unconstrained*
value.  A simplex contributes one fewer factor than it has elements, which
``FreeParameter.free_dimension`` already accounts for.
"""

from __future__ import annotations

from typing import List, Tuple

from ..nodes import compound, distribution, model_ref, param_ref, parameter
from .context import BuildContext

__all__ = ["build", "build_jacobian", "build_bounds"]


def build(context: BuildContext) -> Tuple[dict, List[str]]:
    """Return the variational compound and the ids of its own parameters."""
    options = context.spec.variational
    ids = context.ids
    distributions = []
    variational_parameters: List[str] = []

    for free in context.free_parameters:
        dimension = free.free_dimension
        vector = dimension > 1
        mu = parameter(
            f"mu.{free.id}",
            free.unconstrained_init,
            dimension=dimension if vector else None,
        )
        sigma = parameter(
            f"sigma.{free.id}",
            options.sigma_init,
            dimension=dimension if vector else None,
            lower=0,
        )
        distributions.append(
            distribution(
                f"{ids.variational}.{free.id}",
                "normal",
                param_ref(free.free_id),
                parameters={"mu": mu, "sigma": sigma},
            )
        )
        variational_parameters.extend((param_ref(mu["id"]), param_ref(sigma["id"])))

    return compound(ids.variational, distributions), variational_parameters


def build_jacobian(context: BuildContext) -> dict:
    """The log Jacobian determinant of every transformed parameter.

    Needed once, at the top of the joint, because the approximation lives in
    the unconstrained space while the model's densities are defined in the
    constrained one.
    """
    parameters = [
        param_ref(free.id) for free in context.free_parameters if free.is_transformed
    ]
    return {"id": context.ids.jacobian, "type": "jacobian", "parameters": parameters}


def build_bounds(context: BuildContext) -> List[dict]:
    """The two ELBO nodes: many samples for reporting, few for the gradient."""
    options = context.spec.variational
    ids = context.ids
    return [
        {
            "id": ids.elbo,
            "type": "bound",
            "bound": "elbo",
            "joint": model_ref(ids.joint),
            "variational": model_ref(ids.variational),
            "samples": options.elbo_samples,
            "entropy": options.entropy,
        },
        {
            "id": ids.grad_elbo,
            "type": "bound",
            "bound": "elbo",
            "joint": model_ref(ids.joint),
            "variational": model_ref(ids.variational),
            "samples": options.grad_samples,
            "entropy": options.entropy,
        },
    ]
