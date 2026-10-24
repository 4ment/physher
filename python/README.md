# phyforge

Build [physher](https://github.com/4ment/physher) JSON configuration files.

```bash
uv sync

uv run phyforge ml -i alignment.fa -t tree.newick -m GTR -C 4 -o gtr-g4.json
physher gtr-g4.json
```

`uv sync` creates `.venv` from `uv.lock` and installs `phyforge` into it in
editable mode, so `uv run phyforge` picks up edits to `src/` immediately. To
use it outside this checkout, `uv tool install .` puts the `phyforge` command
on your PATH, and `uv build` writes a wheel and sdist to `dist/`.

## Methods

| Command | What it writes |
| --- | --- |
| `phyforge ml` | a `meta` optimizer maximizing the tree likelihood |
| `phyforge map` | the same, targeting likelihood × priors |
| `phyforge mcmc` | a chain with an operator per free parameter, plus the pattern-trace estimators `--cpo` / `--lppd` |
| `phyforge advi` | a mean-field approximation, its ELBO, and a stochastic optimizer |

Two things are available on top of a method rather than instead of one:
`--marginal` estimates the marginal likelihood, and `-b/--bootstrap` resamples
the alignment and refits.

## Starting branch lengths

An unrooted analysis starts from a parsimony fit of the input topology unless
told otherwise:

| Flag | Where the initial branch lengths come from |
| --- | --- |
| *(default)* | a parsimony fit of the input topology |
| `--no-init-parsimony` | 0.1 everywhere |
| `--keep-branch-lengths` | the lengths already in the input tree file |

The default emits an `initializer` block on the tree node, which makes
physher count the substitutions parsimony infers on each branch and divide by
the number of sites. On the fluA example that moves the starting log
likelihood from -15782.9 to -4337.7; the optimum is unchanged, so this buys
optimizer iterations rather than a better fit. `--init-parsimony-min` sets the
floor applied to those estimates, which otherwise fall to zero on any branch
with no inferred substitution.

Parsimony writes branch lengths, so it stands aside where there are none to
write: a time tree (`-c/--clock`) carries node heights instead, and
`--keep-branch-lengths` asks for the lengths in the file. Neither needs
`--no-init-parsimony` — the default simply does not apply. Spelling
`--init-parsimony` out next to either of them is an error rather than a silent
preference.

One model is sensitive to the change: empirical CAT (`-D cat`) reads its
per-pattern rate assignment off the starting branch lengths, so the start
picks the optimum rather than just the path to it. On fluA, JC69-CAT20 reaches
-4257.3 from flat lengths and -4265.2 from parsimony; use
`--no-init-parsimony` there.

## Rate heterogeneity

`-D/--site-model` picks how rates vary across sites, and `-C` how many
categories there are:

| `-D` | Model |
| --- | --- |
| `gamma`, `weibull`, `lognormal` | a continuous distribution, discretized into `-C` equal-probability categories (`--quadrature` changes the rule) |
| `discrete` | free rates (+R): every category rate is estimated, under the unit-mean constraint |
| `cat` | empirical CAT: each site pattern is assigned to a single rate category, read off the alignment |

`-I` adds a proportion of invariant sites to any of them except CAT (and
except the two free-rate charts that write a rate into category 0).

**Free rates** have four equivalent parameterizations, chosen with
`--rate-parameterization`; they reach the same constraint surface and differ in
conditioning and in the prior they induce. `rate_increments` (the default) and
`rate_ratios` build an increasing sequence, so the categories are identifiable;
`rate_shape` and `mean_contribution` leave them exchangeable. All four are
derivative-free, so `advi` and the gradient optimizers are refused.

`--cat-em` fits the free rates and their weights by EM (Wang, Li, Susko & Roger
2008) instead of sweeping over them one coordinate at a time: the schedule
keeps its branch-length sweep and gets an `em` entry in place of the rate and
weight sweeps. It needs the default `rate_increments` chart and an unrooted
tree -- EM renormalizes the rates and rescales the branch lengths to match,
which a time tree has no distance vector to absorb -- and it maximizes the
likelihood, so under `map` the free-rate priors are still reported but no
longer fitted.

**CAT** estimates its rates and its assignment with physher's `cat` action,
scheduled after the branch-length sweep and repeated until the fit converges.
`--cat-assignment posterior_mean` (the default) is more accurate on branch
lengths; `argmax` reproduces FastTree. `--cat-prior npmle` estimates the prior
over the rate grid from the alignment instead of assuming a Gamma. What the
tree likelihood prints under CAT is a *score*, not a likelihood, so a final
`cat` action reports the mixture log-likelihood — the number that may be
compared with a gamma or free-rate fit — unless `--no-cat-mixture` is given.
CAT is nucleotide-only and maximum likelihood only.

## MCMC

Every free parameter gets its own operators -- a scaler and a slider, a
Dirichlet for a simplex, a beta for a unit-interval parameter -- each weighted
by the **dimension** of the parameter it moves. Weight is how often physher
picks an operator, and a proposal moves one parameter at a time, so a shared
weight would give each of the 135 branch lengths of a 69-taxon tree the same
attention as `kappa`; weighting by dimension keeps the per-element proposal
rate even across the model.

`--cpo` and `--lppd` are the two estimators that read a **pattern trace**
back: a logger writing one column per site pattern to `<stem>.patterns.txt`,
plus an action that passes over that file once the chain is done.

- `--cpo` adds the conditional predictive ordinate and prints `logCPO`.
- `--lppd` adds the log pointwise predictive density and prints `lppd` and
  `pwaic`, the two halves of the WAIC.

Either flag turns the logger on and both share it, so asking for both writes
one file and reads it twice. `--cpo-burnin` and `--lppd-burnin` discard
leading rows of that file -- rows, not iterations, so they count logged
samples; each defaults to a tenth of them. Neither estimator can be combined
with `--steps`, which redirects every file logger per step, so the file they
name is never written.

```bash
phyforge mcmc -i fluA.fa -t fluA.tree -m HKY --length 2000000 --cpo --lppd
```

## Marginal likelihood

`--marginal` names one or more estimators; without it none is written, since
each costs either a ladder of extra chains or a second pass over a trace.

| Estimator | Method | What it reads |
| --- | --- | --- |
| `ss` | `mcmc` | stepping stone over the `--steps` ladder |
| `ps` | `mcmc` | path sampling over the same ladder |
| `ps2` | `mcmc` | its modified form |
| `hm` | `mcmc` | the harmonic mean of one plain chain |
| `shm` | `mcmc` | its stabilized form |
| `is` | `advi` | importance sampling from the fitted approximation |

```bash
# stepping stone and path sampling over the same 32 chains
phyforge mcmc -i fluA.fa -t fluA.tree -m HKY --marginal ss ps --steps 32

# the harmonic mean, which needs no ladder -- one chain, read back once
phyforge mcmc -i fluA.fa -t fluA.tree -m HKY --marginal shm

# reweighting variational draws
phyforge advi -i fluA.fa -t fluA.tree -m HKY --marginal is --is-samples 100000
```

**The ladder.** `--steps` is the number of rungs `ss`/`ps`/`ps2` integrate
over, and physher runs a full chain at each one, so it multiplies the run
time. Both the chains and the estimator rebuild the rungs from `--ladder`,
`--ladder-alpha` and `--ladder-beta`, so phyforge writes those keys into both
nodes rather than letting them drift apart. The default `beta` ladder crowds
the rungs towards the prior, where the integrand moves fastest; `uniform`
spreads them evenly. `--steps` on its own still means stepping stone, which is
what it meant when it was the only way to ask for a marginal likelihood.

Under a ladder every file logger is redirected per step -- `<stem>.log` becomes
`0<stem>.log`, `1<stem>.log` and so on -- which is why `--cpo` and `--lppd`
cannot be combined with it, and why `hm`/`shm` are refused a ladder they would
only read the last rung of.

`--marginal-burnin` discards leading rows of each trace; like the pattern-trace
burnins it counts logged samples rather than iterations, and defaults to a
tenth of them.

`is` reweights `--is-samples` draws from the fitted approximation by the joint
density, which includes the log Jacobian term -- so the estimate is of the
marginal likelihood itself, and should come out above the ELBO.

## Bootstrap

`-b/--bootstrap N` appends a nonparametric bootstrap to an `ml` or `map` run:
`N` replicates, each resampling the site patterns with replacement and
refitting from the point estimate the run just reached.

```bash
phyforge ml -i fluA.fa -t fluA.tree -m HKY -b 100 --bootstrap-observed
```

It is scheduled after the loggers, so what the run reports on screen is still
the fit on the observed data; the replicates go to `<stem>.bootstrap.csv` and
`<stem>.bootstrap.trees`, one row and one tree each. The replicate optimizer is
a second copy of the one that produced the estimate -- replicates fitted a
different way would not be comparable with it -- under ids of its own, since
physher registers every optimizer it reads.

| Flag | Effect |
| --- | --- |
| `--bootstrap-observed` | prepend a round fitted on the observed weights; it should reproduce the point estimate exactly, which makes it a cheap check that the setup is sound |
| `--bootstrap-continue` | start each replicate from the previous one rather than from the estimate |
| `--no-bootstrap-compact` | overwrite the pattern weights in place instead of swapping in a pattern set with the undrawn patterns dropped |
| `--bootstrap-threads N` | fit N replicates concurrently, each thread on its own copy of the model; the output does not depend on N. Not compatible with `--bootstrap-continue` or `--no-bootstrap-compact` |

## Examples

```bash
# unrooted maximum likelihood, GTR + 4 gamma categories + invariant sites
phyforge ml -i fluA.fa -t fluA.tree -m GTR -C 4 -I

# flat starting branch lengths instead of the parsimony default
phyforge ml -i fluA.fa -t fluA.tree -m GTR --no-init-parsimony

# free rates: four categories whose rates and weights are all estimated
phyforge ml -i fluA.fa -t fluA.tree -D discrete -C 4

# the same model, with the mixture fitted by EM rather than by a coordinate sweep
phyforge ml -i fluA.fa -t fluA.tree -D discrete -C 4 --cat-em

# empirical CAT with 20 categories and a prior estimated from the alignment
phyforge ml -i fluA.fa -t fluA.tree -D cat -C 20 --cat-prior npmle

# time tree with a strict clock, dates parsed from the taxon names
phyforge map -i fluA.fa -t fluA.tree -m HKY -c strict --coalescent constant

# dates from a CSV file with strain/date columns
phyforge mcmc -i fluA.fa -t fluA.tree -c strict --dates csv \
    --date-file dates.csv --date-format yyyy-MM-dd

# variational inference under a skygrid
phyforge advi -i fluA.fa -t fluA.tree -c strict \
    --coalescent skygrid --grid 20 --cutoff 30

# 100 bootstrap replicates of a GTR+G4 fit
phyforge ml -i fluA.fa -t fluA.tree -m GTR -C 4 -b 100

# stepping stone marginal likelihood over 32 tempered chains
phyforge mcmc -i fluA.fa -t fluA.tree -m HKY --marginal ss --steps 32
```

`phyforge --help` lists every option, grouped by the part of the model it
configures.

## Using it without the CLI

The CLI is a thin front end. It parses argv into an `AnalysisSpec` and hands
that to the builders; nothing below `phyforge.cli` knows argparse exists, so a
GUI can drive the same code:

```python
from phyforge import AnalysisSpec, build, check
from phyforge.spec import ClockSpec, DataSpec, SubstitutionSpec

spec = AnalysisSpec(
    method="mcmc",
    data=DataSpec(alignment="fluA.fa", tree="fluA.tree"),
    substitution=SubstitutionSpec(kind="HKY"),
    clock=ClockSpec(kind="strict"),
)

problems = check(spec)  # [] when the spec is usable
document = build(spec)  # the physher JSON, as a dict
```

`check` returns every problem it finds, each carrying the `field` that caused
it, so the same validation drives a CLI error message and a red form field.
`AnalysisSpec.to_dict()` / `from_dict()` round-trip through JSON, which is
enough for a GUI to save and reload a session.

### Layout

```
spec.py       what the user wants, as plain dataclasses
validate.py   checks over a spec, reported per field
ids.py        the ids given to generated objects
nodes.py      primitive physher JSON nodes
io/           reading alignments, trees and dates
build/        spec -> physher JSON
cli/          argv -> spec, and printing the result
```

The piece that ties it together is the **free-parameter registry**
(`build/context.py`). As each model builder creates a parameter it declares
whether that parameter is estimated, its dimension, and whether it is
transformed. Everything downstream — MCMC operators, the variational factors,
the log Jacobian term, which columns to log, which parameters an optimizer
moves — is a map over that list rather than a second guess at the model
structure.

## Adding to it

A new substitution model usually needs only an entry in
`SUBSTITUTION_MODELS` and a branch in `build/substitution.py`. A new rate
distribution needs an entry in `SITE_DISTRIBUTIONS` and in
`DISTRIBUTION_PARAMETER` in `build/sitemodel.py`, since physher names each
family's parameter differently (`shape` for gamma and Weibull, `scale` for the
lognormal). A new inference method needs a branch in `build/actions.py` and,
if it introduces its own top-level models, one in `build/analysis.py`. A new
marginal likelihood estimator needs an entry in one of `MARGINAL_TEMPERED`,
`MARGINAL_POSTERIOR` or `MARGINAL_IMPORTANCE`, which are what `validate.py`
checks against and what `build/actions.py` dispatches on; the three families
differ in what the estimator needs the run to have produced.

A model whose parameters are *not* estimated by an optimizer — CAT, whose
rates the `cat` action writes — must keep them out of the free-parameter
registry, or the operators, variational factors and optimizer schedule will
all try to move them. It then owes the trace logger an explicit column, which
is what `sitemodel.log_columns` is for.

## Tests

```bash
uv run pytest
```

`tests/test_physher.py` runs generated configurations through a real physher
binary (`../build/physher`, or `PHYSHER_BIN`) and checks they reach the
expected optimum on the `examples/fluA` data; it is skipped when no binary is
found.

## Known physher issues

Two things phyforge works around or cannot express:

- **Codon models are refused.** physher's codon data type needs a genetic code
  alongside it, which its JSON reader cannot express as a bare `"codon"`
  string, so `-m GY94`/`-m MG94` report an error rather than emit a config
  that will not load.

- **`vbis` had to be repaired to be usable.** Its backup buffer was sized by
  the *number* of parameters rather than their summed dimension, and its
  destructor `free()`d models the hashtable still owned; both crashed. Fixed
  in `src/phyc/is.c`, along with a null dereference when its optional
  `parameters` key was left out.

- **The temperature ladder's shape was unreachable.** `mmcmc` and
  `marginallikelihood` both read `alpha` and `beta` for a beta-spaced ladder,
  but neither declared them, so the JSON validator rejected any config that
  set them. Added to both schemas.

- **Gradients are not guarded against CAT.** A `bfgs`, `cg`, `sg` or `advi`
  run on a CAT tree likelihood reads past the end of the partials buffer
  rather than refusing, so `check` rejects the combination before physher
  sees it. The free-rate parameterizations are safer — they simply have no
  analytic gradient, which reaches an optimizer as a zero — but are refused
  for the same reason.

