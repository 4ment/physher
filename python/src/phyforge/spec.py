"""The description of an analysis, independent of how it was collected.

``AnalysisSpec`` is the boundary between the front end and the JSON builders:
the CLI turns argv into one of these, a GUI turns widgets into one, and
``phyforge.build.analysis.build`` turns one into a physher configuration.  The
dataclasses hold plain values only, so ``dataclasses.asdict`` round-trips a
spec through JSON for session save/load.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from typing import Any, Dict, Optional, Tuple

__all__ = [
    "DataSpec",
    "SubstitutionSpec",
    "SiteModelSpec",
    "CATSpec",
    "ClockSpec",
    "DatesSpec",
    "CoalescentSpec",
    "PriorSpec",
    "OptimizerOptions",
    "MCMCOptions",
    "MarginalOptions",
    "BootstrapOptions",
    "VariationalOptions",
    "AnalysisSpec",
    "METHODS",
    "GRADIENT_ALGORITHMS",
    "SUBSTITUTION_MODELS",
    "SITE_DISTRIBUTIONS",
    "SITE_MODELS",
    "RATE_PARAMETERIZATIONS",
    "CAT_ASSIGNMENTS",
    "CLOCK_MODELS",
    "COALESCENT_MODELS",
    "MARGINAL_ALGORITHMS",
    "MARGINAL_TEMPERED",
    "MARGINAL_POSTERIOR",
    "MARGINAL_IMPORTANCE",
    "LADDERS",
]

METHODS = ("ml", "map", "mcmc", "advi")

#: Optimizers that move every parameter at once using gradients.  A model with
#: no analytic gradient has to stay away from them, and from advi.
GRADIENT_ALGORITHMS = ("bfgs", "cg", "sg")

# Substitution model kind -> the data type it applies to.  physher matches
# "kind" case-insensitively (see substmodel.c).
SUBSTITUTION_MODELS: Dict[str, str] = {
    "JC69": "nucleotide",
    "K80": "nucleotide",
    "F81": "nucleotide",
    "HKY": "nucleotide",
    "GTR": "nucleotide",
    "SYM": "nucleotide",
    "LG": "aa",
    "WAG": "aa",
    "DAYHOFF": "aa",
    "GY94": "codon",
    "MG94": "codon",
}

#: Parametric rate distributions, discretized into ``categories`` bins.
SITE_DISTRIBUTIONS = ("gamma", "weibull", "lognormal")

#: Every rate-heterogeneity model: the parametric distributions above, the free
#: rate mixture ("discrete", +R) whose category rates are estimated directly,
#: and the empirical CAT model, which reads its categories off the alignment.
SITE_MODELS = SITE_DISTRIBUTIONS + ("discrete", "cat")

#: Coordinate systems for the category rates of a free-rate model.  They reach
#: the same constraint surface (sum_k p_k r_k = 1) and differ in conditioning
#: and in the prior they induce; the JSON key names the one that is meant.
RATE_PARAMETERIZATIONS = (
    "rate_increments",
    "rate_ratios",
    "rate_shape",
    "mean_contribution",
)

#: Parameterizations that leave category 0 pinned at rate 0, and so can carry
#: an invariant class.  The other two write their first rate there.
RATE_PARAMETERIZATIONS_WITH_INVARIANT = ("rate_increments", "mean_contribution")

#: How the CAT estimator turns a per-pattern rate profile into an assignment.
CAT_ASSIGNMENTS = ("argmax", "posterior_mean")

#: Marginal likelihood estimators, by what they read.  ``ss``/``ps``/``ps2``
#: need a ladder of tempered chains; ``hm``/``shm`` are a pass over a single
#: posterior sample; ``is`` reweights draws from a variational approximation.
MARGINAL_TEMPERED = ("ss", "ps", "ps2")
MARGINAL_POSTERIOR = ("hm", "shm")
MARGINAL_IMPORTANCE = ("is",)
MARGINAL_ALGORITHMS = MARGINAL_TEMPERED + MARGINAL_POSTERIOR + MARGINAL_IMPORTANCE

#: How the inverse temperatures of a stepping-stone ladder are spaced.
LADDERS = ("beta", "uniform")

CLOCK_MODELS = ("strict", "ncln")
COALESCENT_MODELS = ("constant", "exponential", "skyride", "skygrid", "skyglide")


@dataclass
class DataSpec:
    """Where the alignment and tree come from, and how to embed them."""

    alignment: Optional[str] = None
    tree: Optional[str] = None
    datatype: str = "nucleotide"
    genetic_code: Optional[int] = None
    #: embed sequences in the JSON instead of pointing at the alignment file
    embed_alignment: bool = False
    #: embed the newick string instead of pointing at the tree file
    embed_tree: bool = False
    #: start from the branch lengths in the input tree rather than a default
    keep_branch_lengths: bool = False
    #: start from a parsimony fit of the input topology rather than a default;
    #: ``None`` leaves the choice to :attr:`AnalysisSpec.uses_parsimony`, which
    #: turns it on wherever it applies
    init_parsimony: Optional[bool] = None
    #: floor applied to the parsimony estimates; physher's own if left None
    init_parsimony_min: Optional[float] = None
    #: pass tip states rather than tip partials to the likelihood
    tip_states: bool = False
    #: use the SSE/NEON likelihood cores
    sse: bool = True


@dataclass
class SubstitutionSpec:
    kind: str = "JC69"
    #: "equal" pins the frequencies, "estimate" makes them a free simplex
    frequencies: str = "estimate"
    kappa_init: float = 1.0

    @property
    def normalized_kind(self) -> str:
        return self.kind.lower()

    @property
    def has_free_frequencies(self) -> bool:
        return self.frequencies == "estimate" and self.kind.upper() not in (
            "JC69",
            "K80",
            "SYM",
        )

    @property
    def has_free_rates(self) -> bool:
        return self.kind.upper() in ("GTR", "SYM")

    @property
    def has_kappa(self) -> bool:
        return self.kind.upper() in ("K80", "HKY")


@dataclass
class CATSpec:
    """Options of the empirical CAT estimator, physher's ``cat`` action.

    The rates and the per-pattern assignment are read off the alignment rather
    than optimized, so none of this appears in the site model node itself.
    """

    assignment: str = "posterior_mean"
    #: prior over the rate grid the assignment reads: ``"npmle"`` estimates it
    #: from the alignment, a number is the shape of a mean-one Gamma, and
    #: ``None`` leaves physher's own default (3.0 for ``argmax``, flat for
    #: ``posterior_mean``)
    prior: Optional[str] = None
    #: size of the rate grid the posterior mean is taken over
    probe: Optional[int] = None
    #: report the mixture log-likelihood once the fit is done -- the number
    #: that may be compared with another site model's, unlike the CAT score
    mixture: bool = True
    #: how much the estimator prints about each assignment
    verbosity: int = 1


@dataclass
class SiteModelSpec:
    """Rate heterogeneity across sites."""

    #: number of *variable* rate categories; 1 means no rate distribution
    categories: int = 1
    #: a parametric distribution, "discrete" for free rates, or "cat"
    distribution: str = "gamma"
    invariant: bool = False
    quadrature: Optional[str] = None
    shape_init: float = 0.1
    proportion_invariant_init: float = 0.1
    #: which coordinate system the free rates of a "discrete" model live in
    rate_parameterization: str = "rate_increments"
    #: fit the free rates and their weights by EM rather than by a coordinate
    #: sweep over the increments
    free_rate_em: bool = False
    cat: CATSpec = field(default_factory=CATSpec)

    @property
    def is_cat(self) -> bool:
        """The empirical CAT model: one rate per site, read off the alignment."""
        return self.distribution == "cat"

    @property
    def is_free_rates(self) -> bool:
        """A free-rate mixture (+R): every category rate is estimated."""
        return self.distribution == "discrete" and self.categories > 1

    @property
    def uses_em(self) -> bool:
        """Whether the free-rate mixture is fitted by physher's EM step.

        The flag only means anything for a free-rate model; ``validate``
        refuses it elsewhere rather than letting it be silently dropped here.
        """
        return self.free_rate_em and self.is_free_rates

    @property
    def is_parametric(self) -> bool:
        """A discretized continuous distribution, which has a shape/scale."""
        return self.categories > 1 and self.distribution in SITE_DISTRIBUTIONS

    @property
    def has_shape(self) -> bool:
        return self.is_parametric

    @property
    def has_gradient(self) -> bool:
        """Whether physher can differentiate this site model.

        Neither the free-rate parameterizations nor the empirical rates of CAT
        have an analytic gradient, so both are for the derivative-free
        optimizers and the MCMC operators only.
        """
        return not (self.is_cat or self.is_free_rates)


@dataclass
class DatesSpec:
    """How to obtain tip dates for a time tree."""

    #: "regex" parses the taxon name, "csv" reads a file, "zero" makes all
    #: taxa contemporaneous
    source: str = "regex"
    path: Optional[str] = None
    regex: str = r"_(\d+\.?\d*)$"
    date_format: Optional[str] = None


@dataclass
class ClockSpec:
    """Molecular clock.  ``kind is None`` means an unrooted tree."""

    kind: Optional[str] = None
    #: fix the rate at this value instead of estimating it
    rate: Optional[float] = None
    rate_init: float = 1e-3
    #: spread of the per-branch rates of an ``ncln`` clock
    ncln_scale_init: float = 0.5
    root_height_init: float = 1.0
    dates: DatesSpec = field(default_factory=DatesSpec)

    @property
    def is_time_tree(self) -> bool:
        return self.kind is not None

    @property
    def estimates_rate(self) -> bool:
        return self.kind is not None and self.rate is None

    @property
    def is_relaxed(self) -> bool:
        return self.kind == "ncln"


@dataclass
class CoalescentSpec:
    """Tree prior for a time tree.  ``kind is None`` means no coalescent."""

    kind: Optional[str] = None
    theta_init: float = 3.0
    cutoff: Optional[float] = None
    #: number of grid points for skygrid/skyglide
    grid: Optional[int] = None
    growth_init: float = 0.0
    #: put a GMRF smoothing prior on the log thetas of skyride/skygrid
    gmrf: bool = True
    precision_init: float = 0.1


@dataclass
class PriorSpec:
    """Priors used by ``map``, ``mcmc`` and ``advi``."""

    #: prior on the branch lengths of an unrooted tree
    branch_lengths: str = "exponential"
    branch_lengths_mean: float = 0.1
    #: prior on the strict clock rate: "ctmcscale" or "exponential"
    clock: str = "ctmcscale"
    clock_mean: float = 1000.0
    #: prior on the site model shape parameter
    shape_mean: float = 1.0


@dataclass
class OptimizerOptions:
    iterations: int = 10000
    precision: float = 1e-5
    #: "meta" cycles over per-parameter optimizers, "bfgs"/"cg" are gradient
    #: based and run on the unconstrained parameters
    algorithm: str = "meta"


@dataclass
class MCMCOptions:
    length: int = 1000000
    every: int = 1000
    #: log the per-pattern log likelihoods and estimate the conditional
    #: predictive ordinate from them once the chain is done
    cpo: bool = False
    #: rows of the pattern trace to discard, counted in logged samples rather
    #: than iterations; None drops the first tenth of them
    cpo_burnin: Optional[int] = None
    #: estimate the log pointwise predictive density (and the effective number
    #: of parameters) from that same trace; independent of ``cpo``, and the two
    #: share the one logger when both are on
    lppd: bool = False
    #: rows of the pattern trace the lppd estimator discards; None drops the
    #: first tenth of them
    lppd_burnin: Optional[int] = None

    def logs_patterns(self) -> bool:
        """Whether anything needs the per-pattern log likelihood trace."""
        return self.cpo or self.lppd

    def logged_samples(self) -> int:
        """Number of rows each of the chain's file loggers writes."""
        return self.length // self.every

    def pattern_samples(self) -> int:
        """Number of rows the pattern logger writes."""
        return self.logged_samples()

    def resolved_cpo_burnin(self) -> int:
        """Rows of that trace the cpo estimator skips."""
        return self._resolved_burnin(self.cpo_burnin)

    def resolved_lppd_burnin(self) -> int:
        """Rows of that trace the lppd estimator skips."""
        return self._resolved_burnin(self.lppd_burnin)

    def _resolved_burnin(self, burnin: Optional[int]) -> int:
        if burnin is not None:
            return burnin
        return self.pattern_samples() // 10


@dataclass
class MarginalOptions:
    """Marginal likelihood estimation.

    ``algorithms`` is empty when none is wanted, which is the default: every
    estimator here costs either a ladder of extra chains or a second pass over
    a trace, and none of them is free.  The three families need different
    things from the run, which is what :data:`MARGINAL_TEMPERED`,
    :data:`MARGINAL_POSTERIOR` and :data:`MARGINAL_IMPORTANCE` separate.
    """

    algorithms: Tuple[str, ...] = ()
    #: rungs of the tempered ladder ``ss``/``ps``/``ps2`` are estimated from;
    #: physher runs one chain per rung, so this multiplies the run time
    steps: int = 0
    #: rows of each trace to discard, counted in logged samples rather than
    #: iterations; None drops the first tenth of them
    burnin: Optional[int] = None
    #: spacing of the inverse temperatures: "beta" concentrates the rungs near
    #: the prior, where the integrand moves fastest, "uniform" spreads them
    ladder: str = "beta"
    ladder_alpha: float = 0.3
    ladder_beta: float = 1.0
    #: draws taken from the variational approximation by ``is``
    samples: int = 100000

    @property
    def enabled(self) -> bool:
        return len(self.algorithms) > 0

    @property
    def tempered(self) -> Tuple[str, ...]:
        """The requested estimators that read a ladder of tempered chains."""
        return tuple(a for a in self.algorithms if a in MARGINAL_TEMPERED)

    @property
    def posterior(self) -> Tuple[str, ...]:
        """The requested estimators that read one posterior sample."""
        return tuple(a for a in self.algorithms if a in MARGINAL_POSTERIOR)

    @property
    def importance(self) -> Tuple[str, ...]:
        """The requested estimators that reweight variational draws."""
        return tuple(a for a in self.algorithms if a in MARGINAL_IMPORTANCE)

    @property
    def uses_ladder(self) -> bool:
        """Whether the chain has to be run at a ladder of temperatures."""
        return self.steps > 0

    def resolved_burnin(self, samples: int) -> int:
        """Rows of a trace of ``samples`` rows the estimators skip."""
        if self.burnin is not None:
            return self.burnin
        return samples // 10


@dataclass
class BootstrapOptions:
    """Nonparametric bootstrap over alignment sites, run after the fit.

    Each replicate resamples the site patterns with replacement and refits,
    from the point estimate the main optimizer just reached.
    """

    #: resampled replicates; 0 runs no bootstrap
    replicates: int = 0
    #: restart every replicate from the point estimate rather than from the
    #: previous replicate
    reset: bool = True
    #: prepend a round fitted on the observed weights, as a control
    include_observed: bool = False
    #: swap in a compacted pattern set per replicate instead of overwriting the
    #: weights in place; drops the patterns a replicate did not draw
    compact: bool = True
    #: replicates fitted concurrently, each worker on its own copy of the model
    threads: int = 1

    @property
    def enabled(self) -> bool:
        return self.replicates > 0

    @property
    def rounds(self) -> int:
        """Fits the bootstrap performs, the observed control included."""
        return self.replicates + (1 if self.include_observed else 0)


@dataclass
class VariationalOptions:
    iterations: int = 100000
    eta: float = 0.1
    tol: float = 1e-5
    elbo_samples: int = 100
    grad_samples: int = 1
    #: samples drawn from the fitted approximation at the end
    samples: int = 1000
    entropy: bool = False
    update: str = "adam"
    checkpoint: Optional[str] = None
    mu_init: float = 0.0
    sigma_init: float = 0.308


@dataclass
class AnalysisSpec:
    """A complete description of one physher run."""

    method: str = "ml"
    stem: str = "physher"
    data: DataSpec = field(default_factory=DataSpec)
    substitution: SubstitutionSpec = field(default_factory=SubstitutionSpec)
    site: SiteModelSpec = field(default_factory=SiteModelSpec)
    clock: ClockSpec = field(default_factory=ClockSpec)
    coalescent: CoalescentSpec = field(default_factory=CoalescentSpec)
    priors: PriorSpec = field(default_factory=PriorSpec)
    optimizer: OptimizerOptions = field(default_factory=OptimizerOptions)
    mcmc: MCMCOptions = field(default_factory=MCMCOptions)
    marginal: MarginalOptions = field(default_factory=MarginalOptions)
    bootstrap: BootstrapOptions = field(default_factory=BootstrapOptions)
    variational: VariationalOptions = field(default_factory=VariationalOptions)

    @property
    def uses_priors(self) -> bool:
        """Whether the posterior, rather than the likelihood alone, is the target."""
        return self.method in ("map", "mcmc", "advi")

    @property
    def parsimony_applies(self) -> bool:
        """Whether a parsimony initializer means anything for this spec.

        It writes branch lengths, so it has nothing to say about a time tree,
        and it would overwrite the ones ``keep_branch_lengths`` asks to keep.
        """
        return not self.clock.is_time_tree and not self.data.keep_branch_lengths

    @property
    def uses_parsimony(self) -> bool:
        """Whether the tree node carries a parsimony initializer.

        An unset :attr:`DataSpec.init_parsimony` means the default: parsimony
        wherever it applies, since a warm start costs nothing and saves
        optimizer iterations.
        """
        if self.data.init_parsimony is None:
            return self.parsimony_applies
        return self.data.init_parsimony

    @property
    def uses_unconstrained(self) -> bool:
        """Whether free parameters should be backed by unconstrained ones.

        Gradient-based methods walk the unconstrained space, so every free
        parameter needs a transform.  A time tree always reparameterizes its
        node heights, whatever the method.
        """
        return self.method == "advi" or (
            self.method in ("ml", "map")
            and self.optimizer.algorithm in GRADIENT_ALGORITHMS
        )

    def replace(self, **changes: Any) -> "AnalysisSpec":
        return replace(self, **changes)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AnalysisSpec":
        """Rebuild a spec from ``to_dict`` output (for GUI session files)."""
        nested = {
            "data": DataSpec,
            "substitution": SubstitutionSpec,
            "coalescent": CoalescentSpec,
            "priors": PriorSpec,
            "optimizer": OptimizerOptions,
            "mcmc": MCMCOptions,
            "bootstrap": BootstrapOptions,
            "variational": VariationalOptions,
        }
        kwargs: Dict[str, Any] = dict(data)
        for key, klass in nested.items():
            if key in kwargs:
                kwargs[key] = klass(**kwargs[key])
        if "marginal" in kwargs:
            marginal = dict(kwargs["marginal"])
            # JSON has no tuple, so a round-tripped spec brings it back as a
            # list; the rest of the code treats algorithms as a tuple.
            marginal["algorithms"] = tuple(marginal.get("algorithms", ()))
            kwargs["marginal"] = MarginalOptions(**marginal)
        if "site" in kwargs:
            site = dict(kwargs["site"])
            if "cat" in site:
                site["cat"] = CATSpec(**site["cat"])
            kwargs["site"] = SiteModelSpec(**site)
        if "clock" in kwargs:
            clock = dict(kwargs["clock"])
            if "dates" in clock:
                clock["dates"] = DatesSpec(**clock["dates"])
            kwargs["clock"] = ClockSpec(**clock)
        return cls(**kwargs)
