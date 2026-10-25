"""Coalescent tree priors for time trees."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..nodes import distribution, model_ref, param_ref
from ..spec import PriorChoice
from .context import POSITIVE, REAL, BuildContext
from .parameters import free_parameter

if TYPE_CHECKING:
    from .overrides import PriorWriter

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


def build_priors(context: BuildContext, writer: "PriorWriter") -> list:
    """Hyperpriors on the coalescent's own parameters."""
    spec = context.spec.coalescent
    ids = context.ids
    theta = context.by_role(ROLE_THETA)
    priors = []

    if theta is None:
        # Held fixed by an override: there is nothing left to smooth or to
        # put a reference prior on.
        return priors
    if spec.kind in PIECEWISE_MODELS:
        if spec.gmrf:
            precision = free_parameter(
                context,
                ids.precision,
                spec.precision_init,
                role=ROLE_PRECISION,
                kind=POSITIVE,
            )
            writer.structural(theta, "the GMRF smoothing the log population sizes")
            priors.append(
                distribution(
                    "prior.gmrf",
                    "gmrf",
                    param_ref(theta.free_id),
                    parameters={"precision": precision},
                )
            )
            priors.extend(
                writer.write(
                    context.by_role(ROLE_PRECISION),
                    "prior.gmrf.precision",
                    PriorChoice("gamma", {"shape": 0.001, "rate": 0.001}),
                )
            )
    else:
        # The scale-invariant 1/x reference prior on a population size.
        priors.extend(writer.write(theta, "prior.theta", PriorChoice("oneonx")))
    return priors
