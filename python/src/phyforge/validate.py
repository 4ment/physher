"""Spec-level checks, shared by every front end.

The CLI turns these into error messages; a GUI can attach them to the field
named by :attr:`SpecError.field`.  Checking here rather than inside the
builders means a caller learns about every problem before any JSON is written.
"""

from __future__ import annotations

import os
from typing import List, Optional

from .spec import (
    CAT_ASSIGNMENTS,
    CLOCK_MODELS,
    COALESCENT_MODELS,
    GRADIENT_ALGORITHMS,
    LADDERS,
    MARGINAL_ALGORITHMS,
    MARGINAL_IMPORTANCE,
    MARGINAL_TEMPERED,
    METHODS,
    NO_PRIOR,
    PRIOR_DISTRIBUTIONS,
    RATE_PARAMETERIZATIONS,
    RATE_PARAMETERIZATIONS_WITH_INVARIANT,
    REAL_HYPERPARAMETERS,
    SITE_MODELS,
    SUBSTITUTION_MODELS,
    AnalysisSpec,
    PriorChoice,
)

__all__ = ["SpecError", "check", "check_prior", "validate"]


class SpecError(ValueError):
    """A problem with an analysis spec, tied to the field that caused it.

    ``field`` is the dotted path of the spec attribute at fault, e.g.
    ``site.categories`` or ``parameters.kappa.prior``, so a front end can put
    the message next to the widget or flag that sets it.
    """

    def __init__(self, message: str, field: Optional[str] = None):
        super().__init__(message)
        self.field = field

    def __str__(self) -> str:
        message = super().__str__()
        return f"{self.field}: {message}" if self.field else message


def validate(spec: AnalysisSpec) -> None:
    """Raise the first :class:`SpecError` found in ``spec``."""
    problems = check(spec)
    if problems:
        raise problems[0]


def check(spec: AnalysisSpec) -> List[SpecError]:
    """Return every problem with ``spec``, in the order fields are declared."""
    problems: List[SpecError] = []

    if spec.method not in METHODS:
        problems.append(
            SpecError(
                f"unknown method {spec.method!r} (choose from {', '.join(METHODS)})",
                "method",
            )
        )

    problems.extend(_check_data(spec))
    problems.extend(_check_substitution(spec))
    problems.extend(_check_site(spec))
    problems.extend(_check_clock(spec))
    problems.extend(_check_priors(spec))
    problems.extend(_check_coalescent(spec))
    problems.extend(_check_mcmc(spec))
    problems.extend(_check_marginal(spec))
    problems.extend(_check_bootstrap(spec))
    problems.extend(_check_parameters(spec))
    return problems


def _check_priors(spec: AnalysisSpec) -> List[SpecError]:
    # the exponential priors are written as lambda = 1 / mean
    means = (
        ("branch_lengths_mean", "the mean of the branch length prior"),
        ("clock_mean", "the mean of the clock rate prior"),
        ("shape_mean", "the mean of the shape prior"),
    )
    return [
        SpecError(f"{label} must be positive", f"priors.{name}")
        for name, label in means
        if not getattr(spec.priors, name) > 0
    ]


def _check_mcmc(spec: AnalysisSpec) -> List[SpecError]:
    options = spec.mcmc
    problems = _check_pattern_estimator(
        spec,
        "cpo",
        "the conditional predictive ordinate",
        options.cpo,
        options.cpo_burnin,
    )
    problems.extend(
        _check_pattern_estimator(
            spec,
            "lppd",
            "the log pointwise predictive density",
            options.lppd,
            options.lppd_burnin,
        )
    )
    return problems


def _check_pattern_estimator(
    spec: AnalysisSpec,
    name: str,
    label: str,
    enabled: bool,
    burnin: Optional[int],
) -> List[SpecError]:
    """Check one of the estimators that read the pattern trace back.

    ``cpo`` and ``lppd`` differ only in what they do with the trace, so they
    are subject to the same conditions: a posterior sample to read, a chain
    that actually writes one, and a burnin that leaves something behind.
    """
    problems: List[SpecError] = []
    options = spec.mcmc
    if not enabled:
        return problems

    if spec.method != "mcmc":
        problems.append(
            SpecError(
                f"{label} is estimated from a "
                f"posterior sample, which {spec.method} does not draw",
                f"mcmc.{name}",
            )
        )
        return problems

    if spec.marginal.uses_ladder:
        # Under stepping stone every logger is redirected per step, to
        # "0<file>", "1<file>" and so on, so the file the estimator names is
        # never written -- and a tempered chain is not a posterior sample
        # anyway.
        problems.append(
            SpecError(
                f"{name} cannot be combined with stepping stone (--steps)",
                f"mcmc.{name}",
            )
        )

    samples = options.pattern_samples()
    if samples == 0:
        problems.append(
            SpecError(
                f"a chain of {options.length} sampled every {options.every} "
                f"writes no pattern log likelihoods for {name} to read",
                f"mcmc.{name}",
            )
        )
    elif burnin is not None:
        if burnin < 0:
            problems.append(
                SpecError(f"{name} burnin cannot be negative", f"mcmc.{name}_burnin")
            )
        elif burnin >= samples:
            problems.append(
                SpecError(
                    f"a burnin of {burnin} leaves nothing of the "
                    f"{samples} logged samples",
                    f"mcmc.{name}_burnin",
                )
            )
    return problems


def _check_marginal(spec: AnalysisSpec) -> List[SpecError]:
    """Marginal likelihood estimators, each against what it needs to read."""
    marginal = spec.marginal
    problems: List[SpecError] = []

    unknown = [a for a in marginal.algorithms if a not in MARGINAL_ALGORITHMS]
    for algorithm in unknown:
        problems.append(
            SpecError(
                f"unknown marginal likelihood estimator {algorithm!r} (choose "
                f"from {', '.join(MARGINAL_ALGORITHMS)})",
                "marginal.algorithms",
            )
        )
    if unknown:
        return problems

    # A chain-based estimator needs a chain, and importance sampling needs the
    # variational approximation it reweights.
    chain = tuple(a for a in marginal.algorithms if a not in MARGINAL_IMPORTANCE)
    if chain and spec.method != "mcmc":
        problems.append(
            SpecError(
                f"{', '.join(chain)} {'estimate' if len(chain) > 1 else 'estimates'} "
                f"the marginal likelihood from a chain, which {spec.method} does "
                "not run",
                "marginal.algorithms",
            )
        )
    if marginal.importance and spec.method != "advi":
        problems.append(
            SpecError(
                "importance sampling reweights draws from a variational "
                f"approximation, which {spec.method} does not fit",
                "marginal.algorithms",
            )
        )

    tempered = marginal.tempered
    if tempered and marginal.steps < 2:
        problems.append(
            SpecError(
                f"{', '.join(tempered)} integrates over a ladder of tempered "
                "chains, so it needs at least 2 steps",
                "marginal.steps",
            )
        )
    if marginal.steps == 1:
        problems.append(
            SpecError(
                "a ladder of 1 step has nothing to integrate over", "marginal.steps"
            )
        )
    elif marginal.steps < 0:
        problems.append(
            SpecError("the number of steps cannot be negative", "marginal.steps")
        )
    elif marginal.steps > 1 and not tempered:
        # Running the ladder is the expensive part; without an estimator that
        # reads it, every chain but the last is thrown away.
        problems.append(
            SpecError(
                "a ladder of tempered chains is only read by "
                f"{', '.join(MARGINAL_TEMPERED)}; none was asked for",
                "marginal.steps",
            )
        )

    if marginal.ladder not in LADDERS:
        problems.append(
            SpecError(
                f"unknown temperature ladder {marginal.ladder!r} (choose from "
                f"{', '.join(LADDERS)})",
                "marginal.ladder",
            )
        )
    if marginal.ladder == "beta":
        if marginal.ladder_alpha <= 0:
            problems.append(
                SpecError(
                    "the beta ladder's alpha must be positive", "marginal.ladder_alpha"
                )
            )
        if marginal.ladder_beta <= 0:
            problems.append(
                SpecError(
                    "the beta ladder's beta must be positive", "marginal.ladder_beta"
                )
            )

    if marginal.importance and marginal.samples < 1:
        problems.append(
            SpecError(
                "importance sampling needs at least one draw to reweight",
                "marginal.samples",
            )
        )

    if chain and spec.method == "mcmc":
        samples = spec.mcmc.logged_samples()
        if samples == 0:
            problems.append(
                SpecError(
                    f"a chain of {spec.mcmc.length} sampled every "
                    f"{spec.mcmc.every} writes nothing for the marginal "
                    "likelihood estimators to read",
                    "marginal.algorithms",
                )
            )
        elif marginal.burnin is not None:
            if marginal.burnin < 0:
                problems.append(
                    SpecError(
                        "the marginal burnin cannot be negative", "marginal.burnin"
                    )
                )
            elif marginal.burnin >= samples:
                problems.append(
                    SpecError(
                        f"a burnin of {marginal.burnin} leaves nothing of the "
                        f"{samples} logged samples",
                        "marginal.burnin",
                    )
                )
    return problems


def _check_bootstrap(spec: AnalysisSpec) -> List[SpecError]:
    """The bootstrap resamples the data and refits, so it needs a fit to redo."""
    options = spec.bootstrap
    problems: List[SpecError] = []
    if options.replicates < 0:
        problems.append(
            SpecError(
                "the number of replicates cannot be negative", "bootstrap.replicates"
            )
        )
    if options.threads < 1:
        problems.append(
            SpecError("the number of threads must be at least 1", "bootstrap.threads")
        )
    if not options.enabled:
        return problems
    # physher shares the replicates out only when each one is independent of the
    # others and of the pattern set the threads all read
    if options.threads > 1 and not options.reset:
        problems.append(
            SpecError(
                "threads need every replicate to start from the point estimate, "
                "which --bootstrap-continue turns off",
                "bootstrap.threads",
            )
        )
    if options.threads > 1 and not options.compact:
        problems.append(
            SpecError(
                "threads need the compacted pattern sets: without them a replicate "
                "rewrites the weights every thread reads",
                "bootstrap.threads",
            )
        )
    if spec.method not in ("ml", "map"):
        problems.append(
            SpecError(
                "the bootstrap refits a point estimate on resampled alignments, "
                f"and {spec.method} does not produce one",
                "bootstrap.replicates",
            )
        )
    return problems


def _check_data(spec: AnalysisSpec) -> List[SpecError]:
    problems = []
    if spec.data.alignment is None:
        problems.append(SpecError("an alignment is required", "data.alignment"))
    elif not os.path.exists(spec.data.alignment):
        problems.append(
            SpecError(f"{spec.data.alignment} does not exist", "data.alignment")
        )

    if spec.data.tree is None:
        problems.append(SpecError("a tree is required", "data.tree"))
    elif not os.path.exists(spec.data.tree):
        problems.append(SpecError(f"{spec.data.tree} does not exist", "data.tree"))

    problems.extend(_check_parsimony(spec))
    return problems


def _check_parsimony(spec: AnalysisSpec) -> List[SpecError]:
    problems: List[SpecError] = []
    # Only an explicit --init-parsimony conflicts with anything: left unset it
    # is the default, which stands aside wherever it does not apply.
    if spec.data.init_parsimony:
        if spec.clock.is_time_tree:
            # The initializer writes branch lengths, and a time tree has none:
            # it carries node heights, reparameterized as a root height and
            # proportions.  physher errors out rather than ignoring it.
            problems.append(
                SpecError(
                    "parsimony initializes branch lengths, which a time tree "
                    "does not have",
                    "data.init_parsimony",
                )
            )

        if spec.data.keep_branch_lengths:
            problems.append(
                SpecError(
                    "--init-parsimony and --keep-branch-lengths both set the "
                    "starting branch lengths; parsimony would overwrite the ones "
                    "read from the tree file",
                    "data.init_parsimony",
                )
            )

    if spec.data.init_parsimony_min is not None:
        if not spec.uses_parsimony:
            problems.append(
                SpecError(
                    "a parsimony branch length floor needs the parsimony "
                    "initializer, which is off for this analysis",
                    "data.init_parsimony_min",
                )
            )
        elif spec.data.init_parsimony_min < 0:
            problems.append(
                SpecError(
                    "the parsimony branch length floor cannot be negative",
                    "data.init_parsimony_min",
                )
            )
    return problems


def _check_substitution(spec: AnalysisSpec) -> List[SpecError]:
    kind = spec.substitution.kind.upper()
    if kind not in SUBSTITUTION_MODELS:
        return [
            SpecError(
                f"unknown substitution model {spec.substitution.kind!r} (choose "
                f"from {', '.join(sorted(SUBSTITUTION_MODELS))})",
                "substitution.kind",
            )
        ]
    if SUBSTITUTION_MODELS[kind] == "codon":
        # physher's codon data type needs a genetic code alongside it, which
        # its JSON reader cannot express as a bare string; refuse rather than
        # emit a config that will not load.
        return [
            SpecError(
                f"{kind} is a codon model, which phyforge cannot configure yet",
                "substitution.kind",
            )
        ]
    if spec.substitution.frequencies not in ("equal", "estimate"):
        return [
            SpecError(
                f"frequencies must be 'equal' or 'estimate', not "
                f"{spec.substitution.frequencies!r}",
                "substitution.frequencies",
            )
        ]
    return []


def _check_site(spec: AnalysisSpec) -> List[SpecError]:
    site = spec.site
    problems = []
    if site.categories < 1:
        problems.append(
            SpecError(
                "the number of rate categories must be at least 1", "site.categories"
            )
        )
    if site.distribution not in SITE_MODELS:
        return problems + [
            SpecError(
                f"unknown site rate model {site.distribution!r} "
                f"(choose from {', '.join(SITE_MODELS)})",
                "site.distribution",
            )
        ]
    if site.distribution in ("discrete", "cat") and site.categories < 2:
        problems.append(
            SpecError(
                f"{site.distribution!r} is a mixture over rate categories, so it "
                "needs at least 2 of them",
                "site.categories",
            )
        )
    if site.is_cat:
        problems.extend(_check_cat(spec))
    elif site.is_free_rates:
        problems.extend(_check_free_rates(spec))
    elif site.free_rate_em:
        problems.append(
            SpecError(
                "EM splits a free-rate mixture, and there is none here (it "
                "needs 'discrete' with at least 2 categories)",
                "site.free_rate_em",
            )
        )
    return problems


def _check_free_rates(spec: AnalysisSpec) -> List[SpecError]:
    """Free rates (+R): every category rate estimated, and no gradient."""
    site = spec.site
    problems = []
    if site.rate_parameterization not in RATE_PARAMETERIZATIONS:
        problems.append(
            SpecError(
                f"unknown rate parameterization {site.rate_parameterization!r} "
                f"(choose from {', '.join(RATE_PARAMETERIZATIONS)})",
                "site.rate_parameterization",
            )
        )
    elif (
        site.invariant
        and site.rate_parameterization not in RATE_PARAMETERIZATIONS_WITH_INVARIANT
    ):
        # The other two write their first rate into category 0, which is where
        # the invariant class would sit.
        problems.append(
            SpecError(
                f"{site.rate_parameterization} has no invariant form (use "
                f"{' or '.join(RATE_PARAMETERIZATIONS_WITH_INVARIANT)})",
                "site.invariant",
            )
        )
    if site.quadrature is not None:
        problems.append(
            SpecError(
                "a quadrature discretizes a continuous distribution; free rates "
                "have none to discretize",
                "site.quadrature",
            )
        )
    if site.free_rate_em:
        problems.extend(_check_free_rate_em(spec))
    problems.extend(
        _check_derivative_free(spec, "the free-rate parameterizations have")
    )
    return problems


def _check_free_rate_em(spec: AnalysisSpec) -> List[SpecError]:
    """What physher's free-rate EM (``em.c``) refuses to run on.

    Every one of these is a guard the C code checks and then declines the
    whole step over, leaving the mixture at its starting values -- which a
    schedule reports as a converged sweep rather than as a failure, so it has
    to be caught here.
    """
    site = spec.site
    problems = []
    if site.rate_parameterization != "rate_increments":
        # The M-step returns the rates sorted, which is what makes the
        # increments positive; no other chart is written from that.
        problems.append(
            SpecError(
                "EM is implemented for the rate_increments chart only, not "
                f"{site.rate_parameterization}",
                "site.free_rate_em",
            )
        )
    if spec.method not in ("ml", "map"):
        problems.append(
            SpecError(
                "EM is an entry in the optimizer schedule, and only ml and map "
                "run an optimizer",
                "site.free_rate_em",
            )
        )
    if spec.clock.is_time_tree:
        # The M-step holds the branch lengths fixed and the site model then
        # renormalizes the rates, so the tree has to absorb the same factor.
        # A time tree's lengths are heights x clock rate and own no such
        # vector to scale.
        problems.append(
            SpecError(
                "EM re-imposes the unit-mean constraint by scaling branch "
                "lengths, which a time tree does not own",
                "site.free_rate_em",
            )
        )
    if (site.categories + 1) * site.categories * 1e-4 >= 1.0:
        # The weight floor EM applies is paid for out of the largest weight,
        # which only works while that weight stays above the floor itself.
        problems.append(
            SpecError(
                f"{site.categories} categories cannot all clear the weight "
                "floor EM applies",
                "site.categories",
            )
        )
    return problems


def _check_cat(spec: AnalysisSpec) -> List[SpecError]:
    """Empirical CAT: rates and assignment read off the alignment."""
    cat = spec.site.cat
    problems = []

    datatype = SUBSTITUTION_MODELS.get(spec.substitution.kind.upper())
    if datatype is not None and datatype != "nucleotide":
        problems.append(
            SpecError(
                "the CAT likelihood kernels are implemented for nucleotides only, "
                f"and {spec.substitution.kind} is a {datatype} model",
                "site.distribution",
            )
        )
    if spec.method != "ml":
        # The score CAT reports is not a likelihood -- it sits above one by the
        # mutual information between a column and the category it was handed --
        # so there is no posterior to be had from it, and physher estimates the
        # categories from an optimizer action rather than from a chain.
        problems.append(
            SpecError(
                "the CAT site model is only available for ml: its score is not a "
                "likelihood, and its categories are estimated by an optimizer "
                "action rather than sampled",
                "method",
            )
        )
    if spec.site.invariant:
        problems.append(
            SpecError(
                "CAT reads its categories off the alignment and takes no "
                "invariant class",
                "site.invariant",
            )
        )
    if spec.site.quadrature is not None:
        problems.append(
            SpecError(
                "CAT reads its categories off the alignment rather than "
                "discretizing a distribution",
                "site.quadrature",
            )
        )
    if cat.assignment not in CAT_ASSIGNMENTS:
        problems.append(
            SpecError(
                f"unknown CAT assignment {cat.assignment!r} (choose from "
                f"{', '.join(CAT_ASSIGNMENTS)})",
                "site.cat.assignment",
            )
        )
    if cat.prior is not None and cat.prior != "npmle":
        try:
            shape = float(cat.prior)
        except ValueError:
            problems.append(
                SpecError(
                    f"the CAT prior is 'npmle' or a Gamma shape, not {cat.prior!r}",
                    "site.cat.prior",
                )
            )
        else:
            if shape < 0:
                problems.append(
                    SpecError("a Gamma shape cannot be negative", "site.cat.prior")
                )
    if cat.probe is not None:
        if cat.assignment == "argmax":
            problems.append(
                SpecError(
                    "the probe grid is what the posterior mean is taken over; "
                    "argmax reads the categories themselves",
                    "site.cat.probe",
                )
            )
        elif cat.probe < 1:
            problems.append(
                SpecError("the probe grid must be positive", "site.cat.probe")
            )
    if spec.optimizer.algorithm in GRADIENT_ALGORITHMS:
        # physher does not guard its branch-length gradient against CAT: a
        # gradient-based run reads past the end of the partials buffer rather
        # than refusing, so this has to be caught here.
        problems.append(
            SpecError(
                "the empirical CAT rates are not differentiable, so the site "
                "model has no gradient; use the meta optimizer",
                "optimizer.algorithm",
            )
        )
    return problems


def _check_derivative_free(spec: AnalysisSpec, reason: str) -> List[SpecError]:
    """Refuse the methods that would differentiate a site model that cannot be.

    physher does not guard the gradient against these models -- the free-rate
    parameters simply never accumulate one, which reaches an optimizer as a
    zero, indistinguishable from a converged parameter -- so this has to be
    caught here.
    """
    problems = []
    if spec.method == "advi":
        problems.append(SpecError(f"{reason} no gradient, which advi needs", "method"))
    elif spec.optimizer.algorithm in GRADIENT_ALGORITHMS:
        problems.append(
            SpecError(
                f"{reason} no gradient; use the meta optimizer",
                "optimizer.algorithm",
            )
        )
    return problems


def _check_clock(spec: AnalysisSpec) -> List[SpecError]:
    clock = spec.clock
    if clock.kind is None:
        return []
    problems = []
    if clock.kind not in CLOCK_MODELS:
        problems.append(
            SpecError(
                f"unknown clock {clock.kind!r} (choose from {', '.join(CLOCK_MODELS)})",
                "clock.kind",
            )
        )
    if clock.rate is not None and clock.rate <= 0:
        problems.append(SpecError("a fixed rate must be positive", "clock.rate"))
    if clock.dates.source not in ("regex", "csv", "zero"):
        problems.append(
            SpecError(
                f"unknown date source {clock.dates.source!r}", "clock.dates.source"
            )
        )
    if clock.dates.source == "csv":
        if clock.dates.path is None:
            problems.append(SpecError("a date file is required", "clock.dates.path"))
        elif not os.path.exists(clock.dates.path):
            problems.append(
                SpecError(f"{clock.dates.path} does not exist", "clock.dates.path")
            )
    return problems


def _check_coalescent(spec: AnalysisSpec) -> List[SpecError]:
    coalescent = spec.coalescent
    if coalescent.kind is None:
        return []
    problems = []
    if not spec.clock.is_time_tree:
        problems.append(
            SpecError(
                "a coalescent is a prior on a time tree, so it needs a clock",
                "coalescent.kind",
            )
        )
    if coalescent.kind not in COALESCENT_MODELS:
        problems.append(
            SpecError(
                f"unknown coalescent {coalescent.kind!r} (choose from "
                f"{', '.join(COALESCENT_MODELS)})",
                "coalescent.kind",
            )
        )
    if coalescent.kind in ("skygrid", "skyglide"):
        if coalescent.grid is None:
            problems.append(
                SpecError(
                    f"{coalescent.kind} needs a number of grid points",
                    "coalescent.grid",
                )
            )
        if coalescent.cutoff is None:
            problems.append(
                SpecError(f"{coalescent.kind} needs a cutoff", "coalescent.cutoff")
            )
    if not spec.uses_priors and coalescent.kind is not None:
        problems.append(
            SpecError(
                "a coalescent prior has no effect on a maximum likelihood run",
                "coalescent.kind",
            )
        )
    return problems


def _check_parameters(spec: AnalysisSpec) -> List[SpecError]:
    """Per-parameter overrides, as far as they can be judged on their own.

    Whether an id exists, and whether a prior or a starting value suits the
    parameter it names, depends on the model; the builder checks those once
    it has declared every parameter (see ``build.overrides``).
    """
    problems: List[SpecError] = []
    for id_, override in spec.parameters.items():
        where = f"parameters.{id_}"
        if not override.estimate:
            given = [
                name
                for name in ("prior", "mu", "sigma", "weight")
                if getattr(override, name) is not None
            ]
            if given:
                problems.append(
                    SpecError(
                        f"{id_} is fixed, so it has no {' or '.join(given)}",
                        f"{where}.{given[0]}",
                    )
                )
        if override.prior is not None:
            if not spec.uses_priors:
                problems.append(
                    SpecError(
                        f"a prior on {id_} has no effect on a maximum likelihood run",
                        f"{where}.prior",
                    )
                )
            problems.extend(check_prior(override.prior, f"{where}.prior"))
        for name in ("mu", "sigma"):
            if getattr(override, name) is not None and spec.method != "advi":
                problems.append(
                    SpecError(
                        f"{name} places the variational approximation of {id_}, "
                        f"and only advi fits one",
                        f"{where}.{name}",
                    )
                )
        if override.sigma is not None and any(
            value <= 0 for value in _as_list(override.sigma)
        ):
            problems.append(
                SpecError(f"the scale of {id_} must be positive", f"{where}.sigma")
            )
        if override.weight is not None:
            if spec.method != "mcmc":
                problems.append(
                    SpecError(
                        f"weight sets how often an operator moves {id_}, and only "
                        "mcmc has operators",
                        f"{where}.weight",
                    )
                )
            elif override.weight <= 0:
                problems.append(
                    SpecError(
                        f"the operator weight of {id_} must be positive",
                        f"{where}.weight",
                    )
                )
    return problems


def check_prior(prior: PriorChoice, field: str) -> List[SpecError]:
    """Check a prior's distribution and hyperparameters against physher's."""
    if prior.distribution == NO_PRIOR:
        if prior.parameters:
            return [SpecError("'none' takes no hyperparameters", field)]
        return []
    family = PRIOR_DISTRIBUTIONS.get(prior.distribution)
    if family is None:
        return [
            SpecError(
                f"unknown prior {prior.distribution!r} (choose from "
                f"{', '.join(list(PRIOR_DISTRIBUTIONS) + [NO_PRIOR])})",
                field,
            )
        ]
    if set(prior.parameters) not in [set(keys) for keys in family.hyperparameters]:
        accepted = " or ".join(
            "(" + ", ".join(keys) + ")" if keys else "nothing"
            for keys in family.hyperparameters
        )
        return [
            SpecError(
                f"{prior.distribution} takes {accepted}, not "
                f"({', '.join(prior.parameters)})",
                field,
            )
        ]
    problems = []
    for key, value in prior.parameters.items():
        if key not in REAL_HYPERPARAMETERS and any(v <= 0 for v in _as_list(value)):
            problems.append(
                SpecError(f"the {prior.distribution} {key} must be positive", field)
            )
    return problems


def _as_list(values) -> List[float]:
    return list(values) if isinstance(values, (list, tuple)) else [values]
