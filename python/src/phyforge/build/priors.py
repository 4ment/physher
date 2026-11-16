"""Priors over the free parameters of the evolutionary model.

Each builder below names the *default* prior on the parameters it covers; the
:class:`~phyforge.build.overrides.PriorWriter` swaps in the user's choice where
an override gives one, and gives the parameters with no default a prior of
their own when asked.
"""

from __future__ import annotations

from typing import List

from ..nodes import compound, distribution, param_ref
from ..spec import PriorChoice
from . import coalescent as coalescent_builder
from .branchmodel import (
    NCLN_REFERENCE_MU,
    NCLN_REFERENCE_SIGMA,
    ROLE_BRANCH_RATES,
    ROLE_CLOCK_SCALE,
    ROLE_RATE,
)
from .context import POSITIVE, SIMPLEX, BuildContext
from .overrides import PriorWriter
from .sitemodel import (
    ROLE_SHAPE,
    ROLE_SITE_PROPORTIONS,
    ROLE_SITE_RATES,
    ROLE_TOP_RATE,
)
from .substitution import ROLE_FREQUENCIES, ROLE_RATES
from .tree import ROLE_BRANCH_LENGTHS, ROLE_PROPORTIONS

#: The flat prior every estimated simplex gets.
FLAT_DIRICHLET = PriorChoice("dirichlet")


def build(context: BuildContext) -> dict:
    """Build the ``prior`` compound for the spec.

    Called after the likelihood so that every parameter it covers has already
    been registered.
    """
    spec = context.spec
    writer = PriorWriter(context)
    distributions: List[dict] = []

    if spec.clock.is_time_tree:
        if spec.coalescent.kind is not None:
            distributions.append(coalescent_builder.build(context))
            # The coalescent is the density of the node heights, so it is what
            # the height proportions are distributed as.
            writer.structural(
                context.by_role(ROLE_PROPORTIONS), "the coalescent tree prior"
            )
            distributions.extend(coalescent_builder.build_priors(context, writer))
        if spec.clock.estimates_rate:
            distributions.extend(_clock_prior(context, writer))
        distributions.extend(_relaxed_clock_priors(context, writer))
    else:
        distributions.extend(_branch_length_prior(context, writer))

    distributions.extend(_substitution_priors(context, writer))
    distributions.extend(_site_model_priors(context, writer))
    distributions.extend(writer.remaining())

    return compound(context.ids.prior, distributions)


def _branch_length_prior(context: BuildContext, writer: PriorWriter) -> List[dict]:
    spec = context.spec.priors
    if spec.branch_lengths == "none":
        default = None
    elif spec.branch_lengths == "exponential":
        default = PriorChoice("exponential", {"lambda": 1.0 / spec.branch_lengths_mean})
    else:
        raise ValueError(f"unsupported branch length prior {spec.branch_lengths!r}")
    return writer.write(
        context.by_role(ROLE_BRANCH_LENGTHS), "prior.branch_lengths", default
    )


def _clock_prior(context: BuildContext, writer: PriorWriter) -> List[dict]:
    spec = context.spec.priors
    if spec.clock == "ctmcscale":
        # The CTMC scale reference prior, which needs the tree to scale itself.
        default = PriorChoice("ctmcscale")
    elif spec.clock == "exponential":
        default = PriorChoice("exponential", {"lambda": 1.0 / spec.clock_mean})
    else:
        raise ValueError(f"unsupported clock prior {spec.clock!r}")
    return writer.write(context.by_role(ROLE_RATE), "prior.rate", default)


def _relaxed_clock_priors(context: BuildContext, writer: PriorWriter) -> List[dict]:
    """The reference lognormal on an ``ncln`` clock's raw branch rates.

    Its mu and sigma are the constants the branch model's transform is written
    against, so they are pinned rather than estimated: moving them would mean
    the rates are no longer the standardized draws the transform expects.
    """
    rates = context.by_role(ROLE_BRANCH_RATES)
    if rates is None:
        return []
    reference = PriorChoice(
        "lognormal", {"mu": NCLN_REFERENCE_MU, "sigma": NCLN_REFERENCE_SIGMA}
    )
    # Written directly rather than through the writer: an override must not
    # replace it, but the parameter should still report what it is.
    rates.prior = rates.default_prior = reference
    writer.structural(rates, "the reference lognormal of the ncln clock")
    priors = [
        _node("prior.branch_rates", reference, rates.id),
    ]
    priors.extend(
        writer.write(
            context.by_role(ROLE_CLOCK_SCALE),
            "prior.clock_scale",
            PriorChoice("exponential", {"lambda": 1.0}),
        )
    )
    return priors


def _substitution_priors(context: BuildContext, writer: PriorWriter) -> List[dict]:
    """A flat Dirichlet on each estimated simplex."""
    priors = []
    for role, id_ in (
        (ROLE_FREQUENCIES, "prior.frequencies"),
        (ROLE_RATES, "prior.rates"),
    ):
        priors.extend(writer.write(context.by_role(role), id_, FLAT_DIRICHLET))
    return priors


def _site_model_priors(context: BuildContext, writer: PriorWriter) -> List[dict]:
    priors = writer.write(
        context.by_role(ROLE_SHAPE),
        "prior.shape",
        PriorChoice("exponential", {"lambda": 1.0 / context.spec.priors.shape_mean}),
    )
    priors.extend(_free_rate_priors(context, writer))
    return priors


def _free_rate_priors(context: BuildContext, writer: PriorWriter) -> List[dict]:
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
            default = FLAT_DIRICHLET
        elif parameter.kind == POSITIVE:
            default = PriorChoice("exponential", {"lambda": 1.0})
        else:
            default = None
        priors.extend(writer.write(parameter, id_, default))
    return priors


def _node(id_: str, choice: PriorChoice, x: str) -> dict:
    return distribution(
        id_, choice.distribution, param_ref(x), parameters=dict(choice.parameters)
    )
