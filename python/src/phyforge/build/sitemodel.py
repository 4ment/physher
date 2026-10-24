"""The site model node -- rate heterogeneity across sites."""

from __future__ import annotations

from typing import List, Union

from ..nodes import model_ref, param_ref, parameter
from .context import POSITIVE, UNIT, BuildContext
from .parameters import free_parameter, free_simplex

ROLE_SHAPE = "shape"
ROLE_PROPORTION_INVARIANT = "proportion_invariant"
ROLE_SITE_RATES = "site_rates"
ROLE_SITE_PROPORTIONS = "site_proportions"
ROLE_TOP_RATE = "top_rate"

#: The free parameters EM estimates for itself.  It maximizes over the rates
#: and their weights jointly, so a coordinate sweep over the same parameters
#: is redundant and the schedule leaves them out.
EM_ROLES = (ROLE_SITE_RATES, ROLE_SITE_PROPORTIONS)

#: Each rate distribution names its own parameter: gamma and Weibull take a
#: "shape", the lognormal a "scale" (see new_SiteModel_from_json).
DISTRIBUTION_PARAMETER = {
    "gamma": "shape",
    "weibull": "shape",
    "lognormal": "scale",
}

#: Bounds on a free-rate parameter, per parameterization.  A line search on the
#: constrained value has to be kept off the boundary: an increment of exactly 0
#: collapses two categories onto one rate, and a ratio of exactly 0 or 1 does
#: the same at either end of the chain.
RATE_INCREMENT_BOUNDS = (0.001, 100.0)
RATE_RATIO_BOUNDS = (1e-08, 0.99)
TOP_RATE_BOUNDS = (0.001, 100.0)


def build(context: BuildContext) -> dict:
    """Build the ``sitemodel`` node.

    ``categories`` counts the *variable* categories, so +I adds its class on
    top of them (see ``_check_proportions_dimension`` in ``sitemodel.c``).
    """
    spec = context.spec.site
    ids = context.ids
    node = {"id": ids.site_model, "type": "sitemodel"}

    if spec.is_cat:
        return _cat(context, node)
    if spec.is_free_rates:
        return _free_rates(context, node)

    if spec.categories > 1:
        node["distribution"] = spec.distribution
        node["categories"] = spec.categories
        key = DISTRIBUTION_PARAMETER[spec.distribution]
        node[key] = free_parameter(
            context,
            ids.shape if key == "shape" else ids.scale,
            spec.shape_init,
            role=ROLE_SHAPE,
            kind=POSITIVE,
        )
        if spec.quadrature is not None:
            node["quadrature"] = spec.quadrature
    elif spec.invariant:
        # +I on its own: a point mass at zero plus a single variable rate.
        node["distribution"] = "discrete"
        node["categories"] = 1

    if spec.invariant:
        node["proportion_invariant"] = free_parameter(
            context,
            ids.proportion_invariant,
            spec.proportion_invariant_init,
            role=ROLE_PROPORTION_INVARIANT,
            kind=UNIT,
        )
    return node


# --------------------------------------------------------------------------
# free rates (+R): "distribution": "discrete" with every category rate free
# --------------------------------------------------------------------------


def _free_rates(context: BuildContext, node: dict) -> dict:
    """A free-rate mixture, in whichever chart the spec asks for.

    The unit-mean constraint leaves ``categories - 1`` degrees of freedom
    whatever the chart, and the *key* names the chart: physher checks the shape
    of the parameter against it rather than guessing, so a simplex under
    ``rate_increments`` is a parse error rather than a different model.
    """
    spec = context.spec.site
    ids = context.ids
    categories = spec.categories

    node["distribution"] = "discrete"
    node["categories"] = categories
    if spec.invariant:
        # Not "proportion_invariant", which can only build a 2-element simplex:
        # here category 0 is one more element of the weights simplex, pinned at
        # rate 0 and left out of the rates.
        node["invariant"] = True

    node[spec.rate_parameterization] = _rate_parameter(context, categories)
    if spec.rate_parameterization == "rate_ratios":
        # The ratios alone do not pin the sequence; it is built downwards from
        # the rate of the last category.
        node["top_rate"] = free_parameter(
            context,
            ids.top_rate,
            1.0,
            role=ROLE_TOP_RATE,
            kind=POSITIVE,
            bounds=TOP_RATE_BOUNDS,
        )

    node["proportions"] = free_simplex(
        context,
        ids.site_proportions,
        categories + (1 if spec.invariant else 0),
        role=ROLE_SITE_PROPORTIONS,
    )
    return node


def _rate_parameter(context: BuildContext, categories: int) -> dict:
    """The free parameter of a ``discrete`` site model.

    Its dimension counts the *variable* categories: an invariant class sits at
    rate 0 and contributes nothing to the weighted mean, so it lengthens the
    proportions simplex but not the rates.
    """
    spec = context.spec.site
    ids = context.ids
    parameterization = spec.rate_parameterization

    if parameterization in ("rate_shape", "mean_contribution"):
        return free_simplex(context, ids.site_rates, categories, role=ROLE_SITE_RATES)
    if parameterization == "rate_increments":
        # Equal increments start the rates at 1, 2, ... K, which normalization
        # then centres on 1 -- spread out enough to be told apart, close enough
        # that no category starts on a rate the data cannot support.
        return free_parameter(
            context,
            ids.site_rates,
            1.0,
            role=ROLE_SITE_RATES,
            kind=POSITIVE,
            dimension=categories,
            bounds=RATE_INCREMENT_BOUNDS,
        )
    if parameterization == "rate_ratios":
        return free_parameter(
            context,
            ids.site_rates,
            0.5,
            role=ROLE_SITE_RATES,
            kind=UNIT,
            dimension=categories - 1,
            bounds=RATE_RATIO_BOUNDS,
        )
    raise ValueError(f"unknown rate parameterization {parameterization!r}")


def em_optimizer(context: BuildContext, target: str, prefix: str = "opt") -> dict:
    """One EM step over the free-rate mixture, as a schedule entry.

    The E-step reads the posterior class membership of every pattern off the
    current fit and the M-step re-solves each category rate against it, so
    both halves need the branch lengths and follow the sweep that fits them.
    Left to itself the step takes one iteration per category, which is what
    IQ-TREE does; the enclosing schedule decides how many times it is called.

    EM also rescales the branch lengths, by the factor the site model divides
    out of the raw rates -- the same likelihood, written differently -- which
    is why it is an unrooted-tree step only.
    """
    return {
        "id": f"{prefix}.em",
        "type": "optimizer",
        "algorithm": "em",
        "target": target,
        # Not the objective: the mixture EM splits hangs off the tree
        # likelihood, and EM maximizes the likelihood alone even when the
        # target is a posterior that also carries the free-rate priors.
        "model": model_ref(context.ids.tree_likelihood),
    }


# --------------------------------------------------------------------------
# empirical CAT
# --------------------------------------------------------------------------


def _cat(context: BuildContext, node: dict) -> dict:
    """FastTree's CAT model: one rate per site pattern, read off the alignment.

    The rates are estimated by the ``cat`` action rather than by an optimizer,
    so they are *not* registered as a free parameter: nothing downstream --
    optimizer schedule, MCMC operators, variational factors -- should move
    them.  Their starting values are overwritten by the first assignment.
    """
    spec = context.spec.site
    ids = context.ids
    node["sitepattern"] = param_ref(ids.site_pattern)
    node["categories"] = spec.categories
    node["parameters"] = parameter(
        ids.cat_rates, 1.0, dimension=spec.categories, lower=0.0
    )
    return node


def cat_optimizer(context: BuildContext, target: str) -> dict:
    """A CAT reassignment pass, as an entry in the meta optimizer's schedule.

    It rewrites the category rates and the per-pattern assignment from the
    current branch lengths, so it is only worth anything after a first
    branch-length fit; inside the schedule it alternates with that fit under
    one convergence criterion instead of being unrolled by hand.
    """
    spec = context.spec.site.cat
    node = {
        "id": "opt.cat",
        "type": "optimizer",
        "algorithm": "cat",
        "target": target,
        # Not the objective: the site model the pass rewrites lives here, and
        # the objective may well be a posterior wrapping it.
        "model": model_ref(context.ids.tree_likelihood),
        "assignment": spec.assignment,
    }
    if spec.prior is not None:
        node["prior"] = cat_prior(spec.prior)
    if spec.probe is not None:
        node["probe"] = spec.probe
    if spec.verbosity:
        node["verbosity"] = spec.verbosity
    return node


def cat_report(context: BuildContext) -> dict:
    """Sum the category label out and report the mixture log-likelihood.

    What the tree likelihood prints under CAT is a *score*, not a likelihood:
    it sits above one by the mutual information between a column and the
    category it was handed, and must not be compared with a gamma or free-rate
    fit.  This node changes nothing -- it leaves the assignment, the rates and
    the partials as it found them -- and reports the number that can be.
    """
    node = {
        "id": "catmixture",
        "type": "cat",
        "model": model_ref(context.ids.tree_likelihood),
        "assign": False,
    }
    verbosity = context.spec.site.cat.verbosity
    if verbosity:
        node["verbosity"] = verbosity
    return node


def cat_prior(prior: str) -> Union[str, float]:
    """physher takes the CAT prior as a Gamma shape or the string ``npmle``."""
    if prior == "npmle":
        return prior
    return float(prior)


def log_columns(context: BuildContext) -> List[dict]:
    """Logger columns for site-model quantities that are not free parameters.

    Only CAT has any: its rates are written by the estimator, so they are not
    in the free-parameter registry the trace logger maps over, and the site
    model has to be asked for them by name.
    """
    if not context.spec.site.is_cat:
        return []
    return [{"ref": model_ref(context.ids.site_model), "quantity": "rates"}]
