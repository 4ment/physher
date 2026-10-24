"""The ``physher`` run list: what to do with the model once it is built."""

from __future__ import annotations

from typing import List, Optional

from ..nodes import model_ref, param_ref
from ..spec import GRADIENT_ALGORITHMS
from . import loggers, sitemodel
from .context import BuildContext
from .operators import build as build_operators
from .tree import ROLE_BRANCH_LENGTHS

__all__ = ["build", "target", "GRADIENT_ALGORITHMS"]

#: Where the per-pattern log likelihoods go, and where the estimators that
#: read them back (cpo, lppd) look for them; the three have to agree.
_PATTERN_FILE = "{stem}.patterns.txt"


def build(
    context: BuildContext, variational_parameters: Optional[List[str]] = None
) -> List[dict]:
    """Build the run list for the spec's method.

    ``variational_parameters`` are the ids the stochastic optimizer moves; the
    analysis builder supplies them because they only exist once the variational
    approximation has been built.
    """
    method = context.spec.method
    if method in ("ml", "map"):
        return _optimize(context)
    if method == "mcmc":
        return _mcmc(context)
    if method == "advi":
        return _advi(context, variational_parameters or [])
    raise ValueError(f"unknown method {method!r}")


def target(context: BuildContext) -> str:
    """The density being optimized or sampled."""
    ids = context.ids
    return model_ref(ids.joint if context.spec.uses_priors else ids.tree_likelihood)


# --------------------------------------------------------------------------
# maximum likelihood / maximum a posteriori
# --------------------------------------------------------------------------


def _optimize(context: BuildContext) -> List[dict]:
    stem = context.spec.stem

    actions = [
        _optimizer(context, "opt"),
        # tree_length: what the run just fitted for an unrooted tree, summed
        # into one number. The branch lengths themselves are a nuisance vector
        # that only the tree logger writes out.
        loggers.trace(
            "log", context, targets=loggers.target_refs(context), tree_length=True
        ),
        loggers.tree("treelog", context, file=f"{stem}.tree"),
    ]
    if context.spec.site.is_cat and context.spec.site.cat.mixture:
        # Last, because the number worth printing is the one at the tree that
        # was written out.
        actions.append(sitemodel.cat_report(context))
    actions.extend(_bootstrap(context))
    return actions


def _optimizer(context: BuildContext, id_: str, *, prefix: str = "opt") -> dict:
    """The optimizer that fits the point estimate.

    Built twice: once for the run itself, once for each bootstrap replicate,
    which has to be refitted the same way to mean anything. ``prefix`` names
    the sub-optimizers of a ``meta`` schedule, whose ids have to stay distinct
    between the two copies.
    """
    options = context.spec.optimizer
    if options.algorithm in GRADIENT_ALGORITHMS:
        return {
            "id": id_,
            "type": "optimizer",
            "algorithm": options.algorithm,
            "target": target(context),
            "parameters": [param_ref(free.free_id) for free in context.free_parameters],
            "iterations": options.iterations,
            "precision": options.precision,
        }
    optimizer = {
        "id": id_,
        "type": "optimizer",
        "algorithm": options.algorithm,
        "target": target(context),
        "precision": options.precision,
        "iterations": options.iterations,
        "list": _optimizer_list(context, prefix),
    }
    if context.spec.site.is_cat:
        # The first sweep fits the branch lengths the assignment is read
        # from, so a schedule that stops there has never reassigned
        # anything: force a second one.
        optimizer["min"] = 2
    return optimizer


def _optimizer_list(context: BuildContext, prefix: str = "opt") -> List[dict]:
    """Sub-optimizers for the ``meta`` algorithm, which cycles over them.

    Branch lengths come first. The meta optimizer stops once a whole sweep
    stops improving, so a schedule that tunes a site-model parameter against
    the arbitrary starting branch lengths can settle into a local optimum it
    never leaves -- fitting the gamma shape to branches still at their initial
    value pins the shape, and every later sweep is flat.
    """
    optimizers = []
    uses_em = context.spec.site.uses_em
    ordered = sorted(
        context.free_parameters, key=lambda free: free.role != ROLE_BRANCH_LENGTHS
    )
    for free in ordered:
        if uses_em and free.role in sitemodel.EM_ROLES:
            # EM maximizes over the category rates and their weights jointly,
            # in the step appended below; sweeping over the same parameters
            # one coordinate at a time would only pay for the answer twice.
            continue
        if free.role == ROLE_BRANCH_LENGTHS:
            # Branch lengths have a dedicated routine that exploits the
            # likelihood's structure rather than treating them as N scalars.
            optimizers.append(
                {
                    "id": f"{prefix}.branch_lengths",
                    "type": "optimizer",
                    "algorithm": "serial",
                    "target": target(context),
                    "treelikelihood": model_ref(context.ids.tree_likelihood),
                }
            )
        else:
            optimizers.append(
                {
                    "id": f"{prefix}.{free.id}",
                    "type": "optimizer",
                    "algorithm": "serial",
                    "target": target(context),
                    "parameters": param_ref(free.free_id),
                }
            )
    if context.spec.site.is_cat:
        # Not a parameter sweep: it rewrites the category rates and the
        # per-pattern assignment from the branch lengths the sweep above just
        # fitted, and so has to follow them.
        optimizers.append(sitemodel.cat_optimizer(context, target(context)))
    if uses_em:
        # Both halves of the step are read off the current fit, so like the
        # CAT pass it follows the branch lengths rather than leading them.
        optimizers.append(sitemodel.em_optimizer(context, target(context), prefix))
    return optimizers


# --------------------------------------------------------------------------
# MCMC
# --------------------------------------------------------------------------


def _mcmc(context: BuildContext) -> List[dict]:
    options = context.spec.mcmc
    stem = context.spec.stem
    targets = loggers.target_refs(context)

    log = [
        loggers.trace("screenlog", context, every=options.every * 10, targets=targets),
        loggers.trace(
            "log", context, file=f"{stem}.log", every=options.every, targets=targets
        ),
        loggers.tree("treelog", context, file=f"{stem}.trees", every=options.every),
    ]
    if options.logs_patterns():
        log.append(
            loggers.patterns(
                "patternlog",
                context,
                file=_PATTERN_FILE.format(stem=stem),
                every=options.every,
            )
        )

    mcmc = {
        "id": "mcmc",
        "type": "mcmc",
        "target": target(context),
        "length": options.length,
        "log": log,
        "operators": build_operators(context),
    }

    marginal = context.spec.marginal
    if not marginal.uses_ladder:
        return [mcmc] + _pattern_estimators(context) + _marginal(context)

    # A ladder of temperatures, the same chain run at each rung. Every logger
    # inside is redirected per step, so `<stem>.log` becomes "0<stem>.log",
    # "1<stem>.log" and so on -- which is what the estimators below read.
    ladder = {"id": "mmcmc", "type": "mmcmc", **_ladder(marginal), "mcmc": mcmc}
    return [ladder] + _marginal(context)


def _ladder(marginal) -> dict:
    """The temperature ladder, spelled the same way for mmcmc and the estimator.

    Both rebuild the rungs from these keys rather than sharing them, so a
    mismatch would have the estimator weighting the traces by temperatures no
    chain was run at.
    """
    ladder = {"steps": marginal.steps, "distribution": marginal.ladder}
    if marginal.ladder == "beta":
        ladder["alpha"] = marginal.ladder_alpha
        ladder["beta"] = marginal.ladder_beta
    return ladder


def _pattern_estimators(context: BuildContext) -> List[dict]:
    """What reads the pattern trace back once the chain is done.

    Both run after the chain because each is a pass over the file the chain
    wrote, not something accumulated as the chain goes. They read the same
    file: cpo takes the harmonic mean of each pattern's likelihood, lppd the
    arithmetic one plus its variance.
    """
    options = context.spec.mcmc
    filename = _PATTERN_FILE.format(stem=context.spec.stem)
    model = model_ref(context.ids.tree_likelihood)

    nodes: List[dict] = []
    if options.cpo:
        nodes.append(
            {
                "id": "cpo",
                "type": "cpo",
                "filename": filename,
                "model": model,
                "burnin": options.resolved_cpo_burnin(),
            }
        )
    if options.lppd:
        nodes.append(
            {
                "id": "lppd",
                "type": "predictive",
                "filename": filename,
                "model": model,
                "burnin": options.resolved_lppd_burnin(),
            }
        )
    return nodes


# --------------------------------------------------------------------------
# marginal likelihood
# --------------------------------------------------------------------------


def _marginal(context: BuildContext) -> List[dict]:
    """The estimators that read the chain's trace back.

    One node per algorithm rather than one node computing several: physher
    reads the trace files itself, so the cost is a re-read, and keeping them
    apart means the output names which estimator produced which number.
    """
    marginal = context.spec.marginal
    stem = context.spec.stem
    ladder = _ladder(marginal) if marginal.uses_ladder else {}
    # A ladder logs `samples` rows per rung; a plain chain logs them once.
    burnin = marginal.resolved_burnin(context.spec.mcmc.logged_samples())

    nodes: List[dict] = []
    for algorithm in marginal.tempered + marginal.posterior:
        nodes.append(
            {
                "id": f"marginal.{algorithm}",
                "type": "marginallikelihood",
                "file": f"{stem}.log",
                "treelikelihood": context.ids.tree_likelihood,
                "burnin": burnin,
                "algorithm": algorithm,
                **ladder,
            }
        )
    return nodes


def _importance_sampling(context: BuildContext) -> List[dict]:
    """Marginal likelihood by reweighting draws from the fitted approximation.

    ``normalize`` is off because the joint carries a log Jacobian term the
    self-normalized form would drop: it splits the compound into a likelihood
    and a prior and ignores everything after them.
    """
    marginal = context.spec.marginal
    ids = context.ids
    if not marginal.importance:
        return []
    return [
        {
            "id": "marginal.is",
            "type": "vbis",
            "model": model_ref(ids.joint),
            "distribution": model_ref(ids.variational),
            "samples": marginal.samples,
            "normalize": False,
            # Snapshotted before the draws and put back afterwards: the
            # sampling writes into the unconstrained leaves the approximation
            # was fitted on.
            "parameters": [param_ref(free.free_id) for free in context.free_parameters],
        }
    ]


# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------


def _bootstrap(context: BuildContext) -> List[dict]:
    """Refit the model on resampled alignments, from the point estimate.

    It runs after the loggers above, so what they report is the fit on the
    observed data; the replicates are written to files of their own. The
    optimizer is a second copy of the one that produced the estimate --
    replicates fitted a different way would not be comparable with it.
    """
    options = context.spec.bootstrap
    stem = context.spec.stem
    if not options.enabled:
        return []
    bootstrap = [
        {
            "id": "bootstrap",
            "type": "bootstrap",
            "replicates": options.replicates,
            "compact": options.compact,
            "reset": options.reset,
            "include_observed": options.include_observed,
            "treelikelihood": model_ref(context.ids.tree_likelihood),
            # What is snapshotted at the estimate and restored before each
            # replicate. Given explicitly because a meta schedule's branch
            # length optimizer declares no parameters of its own.
            "parameters": [param_ref(free.free_id) for free in context.free_parameters],
            "optimizer": _optimizer(context, "bootstrap.opt", prefix="bootstrap.opt"),
            "loggers": [
                loggers.trace(
                    "bootstrap.log",
                    context,
                    file=f"{stem}.bootstrap.csv",
                    every=1,
                    targets=loggers.target_refs(context),
                    tree_length=True,
                ),
                loggers.tree(
                    "bootstrap.treelog",
                    context,
                    file=f"{stem}.bootstrap.trees",
                    every=1,
                ),
            ],
        }
    ]
    if options.threads > 1:
        bootstrap[0]["threads"] = options.threads
    return bootstrap


# --------------------------------------------------------------------------
# variational inference
# --------------------------------------------------------------------------


def _advi(context: BuildContext, variational_parameters: List[str]) -> List[dict]:
    options = context.spec.variational
    ids = context.ids
    stem = context.spec.stem

    optimizer = {
        "id": "sg",
        "type": "optimizer",
        "algorithm": "sg",
        "update": options.update,
        "target": model_ref(ids.grad_elbo),
        "eta": options.eta,
        "tol": options.tol,
        "iterations": options.iterations,
        "parameters": variational_parameters,
        "logger": {
            "id": "elbolog",
            "type": "logger",
            "columns": [model_ref(ids.elbo)],
        },
    }
    if options.checkpoint is not None:
        optimizer["checkpoint"] = options.checkpoint

    actions = [optimizer]
    if options.samples > 0:
        actions.append(
            {
                "id": "sampler",
                "type": "sampler",
                "model": model_ref(ids.variational),
                "samples": options.samples,
                "loggers": [
                    loggers.trace(
                        "log",
                        context,
                        file=f"{stem}.csv",
                        targets=[model_ref(ids.joint), model_ref(ids.variational)],
                    ),
                    loggers.tree("treelog", context, file=f"{stem}.trees"),
                ],
            }
        )
    actions.extend(_importance_sampling(context))
    return actions
