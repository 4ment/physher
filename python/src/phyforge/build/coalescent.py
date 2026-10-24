"""Coalescent tree priors for time trees."""

from __future__ import annotations

from ..nodes import distribution, model_ref, param_ref
from .context import POSITIVE, REAL, BuildContext
from .parameters import free_parameter

ROLE_THETA = "theta"
ROLE_GROWTH = "growth"
ROLE_PRECISION = "precision"

#: Models whose theta is a vector smoothed by a GMRF prior.
PIECEWISE_MODELS = ("skyride", "skygrid", "skyglide")


def theta_dimension(context: BuildContext) -> int:
    """Number of population size intervals for the spec's coalescent."""
    spec = context.spec.coalescent
    if spec.kind in ("skygrid", "skyglide"):
        if spec.grid is None:
            raise ValueError(f"{spec.kind} requires a number of grid points")
        return spec.grid
    if spec.kind == "skyride":
        # One interval per coalescent event.
        return context.taxa_count - 1
    return 1


def build(context: BuildContext) -> dict:
    """Build the ``coalescent`` node."""
    spec = context.spec.coalescent
    ids = context.ids
    dimension = theta_dimension(context)
    piecewise = spec.kind in PIECEWISE_MODELS

    node = {
        "id": "coalescent",
        "type": "coalescent",
        "kind": spec.kind,
        # A GMRF smooths the *log* population sizes, so a piecewise model
        # always keeps its unconstrained backing for gmrf to point at.
        "theta": free_parameter(
            context,
            ids.theta,
            spec.theta_init,
            role=ROLE_THETA,
            kind=POSITIVE,
            dimension=dimension if piecewise else None,
            transformed=True if piecewise and spec.gmrf else None,
        ),
        "tree": model_ref(ids.tree),
    }

    if spec.kind in ("skygrid", "skyglide"):
        if spec.cutoff is None:
            raise ValueError(f"{spec.kind} requires a cutoff")
        node["cutoff"] = spec.cutoff
    elif spec.kind == "exponential":
        node["growth"] = free_parameter(
            context,
            ids.growth,
            spec.growth_init,
            role=ROLE_GROWTH,
            kind=REAL,
        )
    return node


def build_priors(context: BuildContext) -> list:
    """Hyperpriors on the coalescent's own parameters."""
    spec = context.spec.coalescent
    ids = context.ids
    theta = context.by_role(ROLE_THETA)
    priors = []

    if spec.kind in PIECEWISE_MODELS:
        if spec.gmrf:
            precision = free_parameter(
                context,
                ids.precision,
                spec.precision_init,
                role=ROLE_PRECISION,
                kind=POSITIVE,
            )
            priors.append(
                distribution(
                    "prior.gmrf",
                    "gmrf",
                    param_ref(theta.free_id),
                    parameters={"precision": precision},
                )
            )
            priors.append(
                distribution(
                    "prior.gmrf.precision",
                    "gamma",
                    param_ref(ids.precision),
                    parameters={"shape": 0.001, "rate": 0.001},
                )
            )
    else:
        # The scale-invariant 1/x reference prior on a population size.
        priors.append(distribution("prior.theta", "oneonx", param_ref(theta.id)))
    return priors
