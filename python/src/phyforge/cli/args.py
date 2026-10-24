"""Mapping command line arguments onto an :class:`AnalysisSpec`.

This module is the *only* place that knows about argparse.  Everything it
produces is a plain spec, so the same analysis can be described by a GUI
without going through argv.
"""

from __future__ import annotations

import argparse
from typing import Optional, Sequence

from ..spec import (
    CAT_ASSIGNMENTS,
    CLOCK_MODELS,
    COALESCENT_MODELS,
    LADDERS,
    MARGINAL_ALGORITHMS,
    METHODS,
    RATE_PARAMETERIZATIONS,
    SITE_MODELS,
    SUBSTITUTION_MODELS,
    AnalysisSpec,
    BootstrapOptions,
    CATSpec,
    ClockSpec,
    CoalescentSpec,
    DataSpec,
    DatesSpec,
    MarginalOptions,
    MCMCOptions,
    OptimizerOptions,
    PriorSpec,
    SiteModelSpec,
    SubstitutionSpec,
    VariationalOptions,
)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phyforge",
        description="Create JSON configuration files for physher.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "method",
        choices=METHODS,
        help="ml: maximum likelihood, map: maximum a posteriori, "
        "mcmc: Metropolis-Hastings, advi: variational inference",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="write the configuration here instead of standard output",
    )
    parser.add_argument(
        "--stem",
        default="physher",
        help="stem for the files physher will write",
    )

    _add_data(parser)
    _add_substitution(parser)
    _add_site(parser)
    _add_clock(parser)
    _add_coalescent(parser)
    _add_priors(parser)
    _add_optimizer(parser)
    _add_mcmc(parser)
    _add_variational(parser)
    _add_marginal(parser)
    _add_bootstrap(parser)
    return parser


def _add_data(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("data")
    group.add_argument("-i", "--input", required=True, help="alignment file (FASTA)")
    group.add_argument(
        "-t", "--tree", required=True, help="tree file (newick or nexus)"
    )
    group.add_argument(
        "--embed-alignment",
        action="store_true",
        help="embed the sequences in the JSON instead of referencing the file",
    )
    group.add_argument(
        "--embed-tree",
        action="store_true",
        help="embed the newick string instead of referencing the file",
    )
    group.add_argument(
        "--keep-branch-lengths",
        action="store_true",
        help="start from the branch lengths in the input tree",
    )
    group.add_argument(
        "--init-parsimony",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="start from a parsimony fit of the input topology; left unset it "
        "is on wherever it applies, which an unrooted tree does, while "
        "--no-init-parsimony starts every branch at 0.1",
    )
    group.add_argument(
        "--init-parsimony-min",
        type=float,
        help="floor applied to the parsimony branch lengths",
    )
    group.add_argument(
        "--tip-states",
        action="store_true",
        help="use tip states instead of tip partials",
    )
    group.add_argument(
        "--no-sse", action="store_true", help="disable the vectorized likelihood cores"
    )


def _add_substitution(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("substitution model")
    group.add_argument(
        "-m",
        "--model",
        default="JC69",
        choices=sorted(SUBSTITUTION_MODELS),
        help="substitution model",
    )
    group.add_argument(
        "-f",
        "--frequencies",
        default="estimate",
        choices=("estimate", "equal"),
        help="equilibrium frequencies",
    )
    group.add_argument(
        "--kappa", type=float, default=1.0, help="initial transition/transversion ratio"
    )


def _add_site(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("site model")
    group.add_argument(
        "-C",
        "--categories",
        type=int,
        default=1,
        metavar="C",
        help="number of variable rate categories",
    )
    group.add_argument(
        "-D",
        "--site-distribution",
        "--site-model",
        dest="site_distribution",
        default="gamma",
        choices=SITE_MODELS,
        help="how rates vary across sites: a distribution to discretize, "
        "'discrete' to estimate every category rate (+R), or 'cat' for the "
        "empirical CAT model",
    )
    group.add_argument(
        "-I",
        "--invariant",
        action="store_true",
        help="add a proportion of invariant sites",
    )
    group.add_argument(
        "--quadrature",
        help="discretization of the rate distribution (e.g. gausslaguerre)",
    )
    group.add_argument(
        "--shape", type=float, default=0.1, help="initial rate distribution shape"
    )
    group.add_argument(
        "--rate-parameterization",
        default="rate_increments",
        choices=RATE_PARAMETERIZATIONS,
        help="coordinate system for the free rates of a 'discrete' site model; "
        "the first two order the categories, the last two leave them "
        "exchangeable",
    )
    group.add_argument(
        "--cat-em",
        action="store_true",
        help="fit the rate categories of a 'discrete' site model -- their "
        "rates and their weights -- by EM instead of sweeping over them one "
        "coordinate at a time; rate_increments and an unrooted tree only",
    )

    group = parser.add_argument_group("empirical CAT")
    group.add_argument(
        "--cat-assignment",
        default="posterior_mean",
        choices=CAT_ASSIGNMENTS,
        help="how a per-pattern rate profile becomes a category; argmax is "
        "FastTree's rule, posterior_mean is more accurate on branch lengths",
    )
    group.add_argument(
        "--cat-prior",
        help="prior over the rate grid: 'npmle' to estimate it from the "
        "alignment, or the shape of a mean-one Gamma",
    )
    group.add_argument(
        "--cat-probe",
        type=int,
        help="size of the rate grid the posterior mean is taken over",
    )
    group.add_argument(
        "--no-cat-mixture",
        action="store_true",
        help="skip the mixture log-likelihood report; what the CAT likelihood "
        "itself prints is a score, and is not comparable with another model's",
    )
    group.add_argument(
        "--cat-verbosity",
        type=int,
        default=1,
        help="what the CAT estimator prints about each assignment",
    )


def _add_clock(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("clock")
    group.add_argument(
        "-c",
        "--clock",
        choices=CLOCK_MODELS,
        help="molecular clock; without one the tree is unrooted",
    )
    group.add_argument("--rate", type=float, help="fix the substitution rate")
    group.add_argument(
        "--rate-init", type=float, default=1e-3, help="initial substitution rate"
    )
    group.add_argument(
        "--root-height",
        type=float,
        default=5.0,
        help="initial root height above the oldest tip",
    )
    group.add_argument(
        "--dates",
        choices=("regex", "csv", "zero"),
        default="regex",
        help="where tip dates come from",
    )
    group.add_argument(
        "--date-file", help="CSV file of dates, with strain/date columns"
    )
    group.add_argument(
        "--date-regex",
        default=r"_(\d+\.?\d*)$",
        help="pattern capturing the date in a taxon name",
    )
    group.add_argument(
        "--date-format",
        help="date layout, e.g. yyyy-MM-dd, when it is not a decimal year",
    )


def _add_coalescent(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("tree prior")
    group.add_argument(
        "--coalescent",
        choices=COALESCENT_MODELS,
        help="coalescent prior on a time tree",
    )
    group.add_argument(
        "--theta", type=float, default=3.0, help="initial population size"
    )
    group.add_argument(
        "--grid", type=int, help="number of grid points for skygrid/skyglide"
    )
    group.add_argument("--cutoff", type=float, help="oldest grid point")
    group.add_argument(
        "--no-gmrf",
        action="store_true",
        help="drop the GMRF smoothing prior on a piecewise coalescent",
    )


def _add_priors(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("priors")
    group.add_argument(
        "--branch-length-prior",
        default="exponential",
        choices=("exponential", "none"),
        help="prior on the branch lengths of an unrooted tree",
    )
    group.add_argument("--branch-length-mean", type=float, default=0.1, help="its mean")
    group.add_argument(
        "--clock-prior",
        default="ctmcscale",
        choices=("ctmcscale", "exponential"),
        help="prior on the clock rate",
    )
    group.add_argument(
        "--clock-prior-mean", type=float, default=1000.0, help="its mean"
    )


def _add_optimizer(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("ml/map")
    group.add_argument(
        "--algorithm",
        default="meta",
        choices=("meta", "bfgs", "cg", "sg"),
        help="optimization algorithm; the last three are gradient based",
    )
    group.add_argument(
        "--iterations", type=int, default=10000, help="maximum optimizer iterations"
    )
    group.add_argument(
        "--precision", type=float, default=1e-3, help="convergence tolerance"
    )


def _add_mcmc(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("mcmc")
    group.add_argument("--length", type=int, default=1000000, help="chain length")
    group.add_argument("--every", type=int, default=1000, help="logging frequency")
    group.add_argument(
        "--cpo",
        action="store_true",
        help="log the per-pattern log likelihoods and report the conditional "
        "predictive ordinate at the end of the chain",
    )
    group.add_argument(
        "--cpo-burnin",
        type=int,
        help="samples of the pattern trace to discard, counted in logged rows "
        "rather than iterations; defaults to a tenth of them",
    )
    group.add_argument(
        "--lppd",
        action="store_true",
        help="log the per-pattern log likelihoods and report the log pointwise "
        "predictive density and p_WAIC at the end of the chain",
    )
    group.add_argument(
        "--lppd-burnin",
        type=int,
        help="samples of the pattern trace the lppd estimator discards; "
        "defaults to a tenth of them",
    )


def _add_marginal(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("marginal likelihood")
    group.add_argument(
        "--marginal",
        nargs="+",
        default=(),
        metavar="ESTIMATOR",
        choices=MARGINAL_ALGORITHMS,
        help="estimate the marginal likelihood: ss (stepping stone), ps and "
        "ps2 (path sampling, and its modified form) integrate over the "
        "ladder --steps builds; hm and shm are the harmonic mean and its "
        "stabilized form, taken over a plain chain; is reweights draws from "
        "an advi approximation",
    )
    group.add_argument(
        "--steps",
        type=int,
        default=0,
        help="rungs of the tempered ladder ss/ps/ps2 integrate over; physher "
        "runs one chain per rung, so this multiplies the run time",
    )
    group.add_argument(
        "--marginal-burnin",
        type=int,
        help="samples of each trace to discard, counted in logged rows rather "
        "than iterations; defaults to a tenth of them",
    )
    group.add_argument(
        "--ladder",
        default="beta",
        choices=LADDERS,
        help="spacing of the inverse temperatures: beta concentrates the "
        "rungs near the prior, where the integrand moves fastest",
    )
    group.add_argument(
        "--ladder-alpha",
        type=float,
        default=0.3,
        help="shape of the beta ladder; smaller crowds the rungs harder "
        "towards the prior",
    )
    group.add_argument(
        "--ladder-beta", type=float, default=1.0, help="its second shape"
    )
    group.add_argument(
        "--is-samples",
        type=int,
        default=100000,
        help="draws taken from the variational approximation by --marginal is",
    )


def _add_bootstrap(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("bootstrap")
    group.add_argument(
        "-b",
        "--bootstrap",
        type=int,
        default=0,
        metavar="N",
        help="refit the model on N alignments resampled from the site "
        "patterns, starting from the estimate the run just reached",
    )
    group.add_argument(
        "--bootstrap-continue",
        action="store_true",
        help="start each replicate from the previous one rather than from the "
        "point estimate",
    )
    group.add_argument(
        "--bootstrap-observed",
        action="store_true",
        help="prepend a replicate fitted on the observed weights, as a control",
    )
    group.add_argument(
        "--no-bootstrap-compact",
        action="store_true",
        help="overwrite the pattern weights in place instead of swapping in a "
        "pattern set with the undrawn patterns dropped",
    )
    group.add_argument(
        "--bootstrap-threads",
        type=int,
        default=1,
        metavar="N",
        help="fit N replicates concurrently, each thread on its own copy of the "
        "model; the output does not depend on N",
    )


def _add_variational(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("advi")
    group.add_argument("--eta", type=float, default=0.1, help="learning rate")
    group.add_argument(
        "--vi-iterations", type=int, default=100000, help="maximum ADVI iterations"
    )
    group.add_argument("--tol", type=float, default=1e-5, help="convergence tolerance")
    group.add_argument(
        "--elbo-samples", type=int, default=100, help="samples for the reported ELBO"
    )
    group.add_argument(
        "--grad-samples", type=int, default=1, help="samples for each gradient estimate"
    )
    group.add_argument(
        "--samples",
        type=int,
        default=1000,
        help="samples drawn from the fitted approximation",
    )
    group.add_argument(
        "--entropy", action="store_true", help="use the entropy rather than log q"
    )
    group.add_argument("--checkpoint", help="write optimizer checkpoints here")


def to_spec(arguments: argparse.Namespace) -> AnalysisSpec:
    """Build an :class:`AnalysisSpec` from parsed arguments."""
    return AnalysisSpec(
        method=arguments.method,
        stem=arguments.stem,
        data=DataSpec(
            alignment=arguments.input,
            tree=arguments.tree,
            embed_alignment=arguments.embed_alignment,
            embed_tree=arguments.embed_tree,
            keep_branch_lengths=arguments.keep_branch_lengths,
            init_parsimony=arguments.init_parsimony,
            init_parsimony_min=arguments.init_parsimony_min,
            tip_states=arguments.tip_states,
            sse=not arguments.no_sse,
        ),
        substitution=SubstitutionSpec(
            kind=arguments.model,
            frequencies=arguments.frequencies,
            kappa_init=arguments.kappa,
        ),
        site=SiteModelSpec(
            categories=arguments.categories,
            distribution=arguments.site_distribution,
            invariant=arguments.invariant,
            quadrature=arguments.quadrature,
            shape_init=arguments.shape,
            rate_parameterization=arguments.rate_parameterization,
            free_rate_em=arguments.cat_em,
            cat=CATSpec(
                assignment=arguments.cat_assignment,
                prior=arguments.cat_prior,
                probe=arguments.cat_probe,
                mixture=not arguments.no_cat_mixture,
                verbosity=arguments.cat_verbosity,
            ),
        ),
        clock=ClockSpec(
            kind=arguments.clock,
            rate=arguments.rate,
            rate_init=arguments.rate_init,
            root_height_init=arguments.root_height,
            dates=DatesSpec(
                source=arguments.dates,
                path=arguments.date_file,
                regex=arguments.date_regex,
                date_format=arguments.date_format,
            ),
        ),
        coalescent=CoalescentSpec(
            kind=arguments.coalescent,
            theta_init=arguments.theta,
            cutoff=arguments.cutoff,
            grid=arguments.grid,
            gmrf=not arguments.no_gmrf,
        ),
        priors=PriorSpec(
            branch_lengths=arguments.branch_length_prior,
            branch_lengths_mean=arguments.branch_length_mean,
            clock=arguments.clock_prior,
            clock_mean=arguments.clock_prior_mean,
        ),
        optimizer=OptimizerOptions(
            iterations=arguments.iterations,
            precision=arguments.precision,
            algorithm=arguments.algorithm,
        ),
        mcmc=MCMCOptions(
            length=arguments.length,
            every=arguments.every,
            cpo=arguments.cpo,
            cpo_burnin=arguments.cpo_burnin,
            lppd=arguments.lppd,
            lppd_burnin=arguments.lppd_burnin,
        ),
        marginal=MarginalOptions(
            # --steps on its own means stepping stone, which is what it meant
            # when it was the only way to ask for a marginal likelihood.
            algorithms=tuple(arguments.marginal)
            or (("ss",) if arguments.steps > 0 else ()),
            steps=arguments.steps,
            burnin=arguments.marginal_burnin,
            ladder=arguments.ladder,
            ladder_alpha=arguments.ladder_alpha,
            ladder_beta=arguments.ladder_beta,
            samples=arguments.is_samples,
        ),
        bootstrap=BootstrapOptions(
            replicates=arguments.bootstrap,
            reset=not arguments.bootstrap_continue,
            include_observed=arguments.bootstrap_observed,
            compact=not arguments.no_bootstrap_compact,
            threads=arguments.bootstrap_threads,
        ),
        variational=VariationalOptions(
            iterations=arguments.vi_iterations,
            eta=arguments.eta,
            tol=arguments.tol,
            elbo_samples=arguments.elbo_samples,
            grad_samples=arguments.grad_samples,
            samples=arguments.samples,
            entropy=arguments.entropy,
            checkpoint=arguments.checkpoint,
        ),
    )


def parse(argv: Optional[Sequence[str]] = None) -> AnalysisSpec:
    """Parse ``argv`` straight into a spec."""
    return to_spec(make_parser().parse_args(argv))
