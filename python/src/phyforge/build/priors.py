"""Priors over the free parameters of the evolutionary model."""

from __future__ import annotations

from typing import List

from ..nodes import compound, distribution, model_ref, param_ref
from . import coalescent as coalescent_builder
from .branchmodel import (
    NCLN_REFERENCE_MU,
    NCLN_REFERENCE_SIGMA,
    ROLE_BRANCH_RATES,
    ROLE_CLOCK_SCALE,
    ROLE_RATE,
)
from .context import POSITIVE, SIMPLEX, BuildContext
from .sitemodel import (
    ROLE_SHAPE,
    ROLE_SITE_PROPORTIONS,
    ROLE_SITE_RATES,
    ROLE_TOP_RATE,
)
from .substitution import ROLE_FREQUENCIES, ROLE_RATES
from .tree import ROLE_BRANCH_LENGTHS


def build(context: BuildContext) -> dict:
    """Build the ``prior`` compound for the spec.

    Called after the likelihood so that every parameter it covers has already
    been registered.
    """
    spec = context.spec
    distributions: List[dict] = []

    if spec.clock.is_time_tree:
        if spec.coalescent.kind is not None:
            distributions.append(coalescent_builder.build(context))
            distributions.extend(coalescent_builder.build_priors(context))
        if spec.clock.estimates_rate:
            distributions.extend(_clock_prior(context))
        distributions.extend(_relaxed_clock_priors(context))
    else:
        distributions.extend(_branch_length_prior(context))

    distributions.extend(_substitution_priors(context))
    distributions.extend(_site_model_priors(context))

    return compound(context.ids.prior, distributions)


def _branch_length_prior(context: BuildContext) -> List[dict]:
    spec = context.spec.priors
    branch_lengths = context.by_role(ROLE_BRANCH_LENGTHS)
    if branch_lengths is None or spec.branch_lengths == "none":
        return []
    if spec.branch_lengths != "exponential":
        raise ValueError(f"unsupported branch length prior {spec.branch_lengths!r}")
    return [
        distribution(
            "prior.branch_lengths",
            "exponential",
            param_ref(branch_lengths.id),
            parameters={"lambda": 1.0 / spec.branch_lengths_mean},
        )
    ]


def _clock_prior(context: BuildContext) -> List[dict]:
    spec = context.spec.priors
    rate = context.by_role(ROLE_RATE)
    if rate is None:
        return []
    if spec.clock == "ctmcscale":
        # The CTMC scale reference prior, which needs the tree to scale itself.
        return [
            distribution(
                "prior.rate",
                "ctmcscale",
                param_ref(rate.id),
                tree=model_ref(context.ids.tree),
            )
        ]
    if spec.clock == "exponential":
        return [
            distribution(
                "prior.rate",
                "exponential",
                param_ref(rate.id),
                parameters={"lambda": spec.clock_mean},
            )
        ]
    raise ValueError(f"unsupported clock prior {spec.clock!r}")


def _relaxed_clock_priors(context: BuildContext) -> List[dict]:
    """The reference lognormal on an ``ncln`` clock's raw branch rates.

    Its mu and sigma are the constants the branch model's transform is written
    against, so they are pinned rather than estimated: moving them would mean
    the rates are no longer the standardized draws the transform expects.
    """
    rates = context.by_role(ROLE_BRANCH_RATES)
    if rates is None:
        return []
    priors = [
        distribution(
            "prior.branch_rates",
            "lognormal",
            param_ref(rates.id),
            parameters={
                "mu": NCLN_REFERENCE_MU,
                "sigma": NCLN_REFERENCE_SIGMA,
            },
        )
    ]
    scale = context.by_role(ROLE_CLOCK_SCALE)
    if scale is not None:
        priors.append(
            distribution(
                "prior.clock_scale",
                "exponential",
                param_ref(scale.id),
                parameters={"lambda": 1.0},
            )
        )
    return priors


def _substitution_priors(context: BuildContext) -> List[dict]:
    """A flat Dirichlet on each estimated simplex."""
    priors = []
    for role, id_ in (
        (ROLE_FREQUENCIES, "prior.frequencies"),
        (ROLE_RATES, "prior.rates"),
    ):
        parameter = context.by_role(role)
        if parameter is not None:
            priors.append(distribution(id_, "dirichlet", param_ref(parameter.id)))
    return priors


def _site_model_priors(context: BuildContext) -> List[dict]:
    priors = []
    shape = context.by_role(ROLE_SHAPE)
    if shape is not None:
        priors.append(
            distribution(
                "prior.shape",
                "exponential",
                param_ref(shape.id),
                parameters={"lambda": 1.0 / context.spec.priors.shape_mean},
            )
        )
    priors.extend(_free_rate_priors(context))
    return priors


def _free_rate_priors(context: BuildContext) -> List[dict]:
    """Priors on the free category rates and their weights.

    Flat on the simplices, as everywhere else.  The positive parameters get an
    exponential rather than nothing because the ordered charts carry a
    redundant degree of freedom -- the likelihood is *exactly* flat along
    ``theta -> c theta``, since the rates are normalized afterwards -- so an
    improper prior would leave an improper posterior.  A parameter bounded
    into ``(0, 1)`` (the ratios) is already proper, and is left alone.
    """
    priors = []
    for role, id_ in (
        (ROLE_SITE_RATES, "prior.site_rates"),
        (ROLE_SITE_PROPORTIONS, "prior.site_proportions"),
        (ROLE_TOP_RATE, "prior.top_rate"),
    ):
        parameter = context.by_role(role)
        if parameter is None:
            continue
        if parameter.kind == SIMPLEX:
            priors.append(distribution(id_, "dirichlet", param_ref(parameter.id)))
        elif parameter.kind == POSITIVE:
            priors.append(
                distribution(
                    id_,
                    "exponential",
                    param_ref(parameter.id),
                    parameters={"lambda": 1.0},
                )
            )
    return priors
