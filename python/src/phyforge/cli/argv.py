"""The inverse of :mod:`phyforge.cli.args`: a spec written back as a command.

A GUI shows this command next to the JSON it builds, so the same analysis can
be regenerated from a script.  Only options that differ from what the parser
would choose on its own are written, which keeps the command short; the round
trip ``parse(to_argv(spec)) == spec`` is what the tests hold it to.
"""

from __future__ import annotations

import shlex
from typing import Any, List, Sequence, Tuple

from ..spec import AnalysisSpec, format_values
from .args import parse

__all__ = ["to_argv", "to_command", "UNEXPRESSED"]

#: How an option is written: with its value, as a bare switch when the field
#: is true ("on") or false ("off"), as a pair of switches for a tri-state
#: (``--x`` / ``--no-x``), or with every element of a list.
VALUE, ON, OFF, TRISTATE, LIST = "value", "on", "off", "tristate", "list"

#: (spec path, flag, how), in the order the parser declares them.
OPTIONS: Sequence[Tuple[str, str, str]] = (
    ("stem", "--stem", VALUE),
    ("data.embed_alignment", "--embed-alignment", ON),
    ("data.embed_tree", "--embed-tree", ON),
    ("data.keep_branch_lengths", "--keep-branch-lengths", ON),
    ("data.init_parsimony", "--init-parsimony", TRISTATE),
    ("data.init_parsimony_min", "--init-parsimony-min", VALUE),
    ("data.tip_states", "--tip-states", ON),
    ("data.sse", "--no-sse", OFF),
    ("substitution.kind", "--model", VALUE),
    ("substitution.frequencies", "--frequencies", VALUE),
    ("substitution.kappa_init", "--kappa", VALUE),
    ("site.categories", "--categories", VALUE),
    ("site.distribution", "--site-model", VALUE),
    ("site.invariant", "--invariant", ON),
    ("site.quadrature", "--quadrature", VALUE),
    ("site.shape_init", "--shape", VALUE),
    ("site.proportion_invariant_init", "--pinv", VALUE),
    ("site.rate_parameterization", "--rate-parameterization", VALUE),
    ("site.free_rate_em", "--cat-em", ON),
    ("site.cat.assignment", "--cat-assignment", VALUE),
    ("site.cat.prior", "--cat-prior", VALUE),
    ("site.cat.probe", "--cat-probe", VALUE),
    ("site.cat.mixture", "--no-cat-mixture", OFF),
    ("site.cat.verbosity", "--cat-verbosity", VALUE),
    ("clock.kind", "--clock", VALUE),
    ("clock.rate", "--rate", VALUE),
    ("clock.rate_init", "--rate-init", VALUE),
    ("clock.ncln_scale_init", "--ncln-scale", VALUE),
    ("clock.root_height_init", "--root-height", VALUE),
    ("clock.dates.source", "--dates", VALUE),
    ("clock.dates.path", "--date-file", VALUE),
    ("clock.dates.regex", "--date-regex", VALUE),
    ("clock.dates.date_format", "--date-format", VALUE),
    ("coalescent.kind", "--coalescent", VALUE),
    ("coalescent.theta_init", "--theta", VALUE),
    ("coalescent.growth_init", "--growth", VALUE),
    ("coalescent.grid", "--grid", VALUE),
    ("coalescent.cutoff", "--cutoff", VALUE),
    ("coalescent.gmrf", "--no-gmrf", OFF),
    ("coalescent.precision_init", "--gmrf-precision", VALUE),
    ("priors.branch_lengths", "--branch-length-prior", VALUE),
    ("priors.branch_lengths_mean", "--branch-length-mean", VALUE),
    ("priors.clock", "--clock-prior", VALUE),
    ("priors.clock_mean", "--clock-prior-mean", VALUE),
    ("priors.shape_mean", "--shape-prior-mean", VALUE),
    ("optimizer.algorithm", "--algorithm", VALUE),
    ("optimizer.iterations", "--iterations", VALUE),
    ("optimizer.precision", "--precision", VALUE),
    ("mcmc.length", "--length", VALUE),
    ("mcmc.every", "--every", VALUE),
    ("mcmc.cpo", "--cpo", ON),
    ("mcmc.cpo_burnin", "--cpo-burnin", VALUE),
    ("mcmc.lppd", "--lppd", ON),
    ("mcmc.lppd_burnin", "--lppd-burnin", VALUE),
    ("variational.eta", "--eta", VALUE),
    ("variational.iterations", "--vi-iterations", VALUE),
    ("variational.tol", "--tol", VALUE),
    ("variational.elbo_samples", "--elbo-samples", VALUE),
    ("variational.grad_samples", "--grad-samples", VALUE),
    ("variational.samples", "--samples", VALUE),
    ("variational.entropy", "--entropy", ON),
    ("variational.checkpoint", "--checkpoint", VALUE),
    ("variational.update", "--vi-update", VALUE),
    ("variational.sigma_init", "--vi-sigma-init", VALUE),
    ("marginal.algorithms", "--marginal", LIST),
    ("marginal.steps", "--steps", VALUE),
    ("marginal.burnin", "--marginal-burnin", VALUE),
    ("marginal.ladder", "--ladder", VALUE),
    ("marginal.ladder_alpha", "--ladder-alpha", VALUE),
    ("marginal.ladder_beta", "--ladder-beta", VALUE),
    ("marginal.samples", "--is-samples", VALUE),
    ("bootstrap.replicates", "--bootstrap", VALUE),
    ("bootstrap.reset", "--bootstrap-continue", OFF),
    ("bootstrap.include_observed", "--bootstrap-observed", ON),
    ("bootstrap.compact", "--no-bootstrap-compact", OFF),
    ("bootstrap.threads", "--bootstrap-threads", VALUE),
)

#: Spec fields no flag sets, because no builder reads them: the data type
#: follows from the substitution model.
UNEXPRESSED = ("data.datatype", "data.genetic_code", "variational.mu_init")


def to_argv(spec: AnalysisSpec) -> List[str]:
    """The ``phyforge`` arguments that parse back into ``spec``."""
    defaults = parse([spec.method, "-i", "-", "-t", "-"])
    argv = [spec.method, "-i", str(spec.data.alignment), "-t", str(spec.data.tree)]
    for path, flag, how in OPTIONS:
        value, default = _get(spec, path), _get(defaults, path)
        if value == default:
            continue
        if how == VALUE:
            argv += [flag, _text(value)]
        elif how == ON and value:
            argv.append(flag)
        elif how == OFF and not value:
            argv.append(flag)
        elif how == TRISTATE and value is not None:
            argv.append(flag if value else flag.replace("--", "--no-", 1))
        elif how == LIST and value:
            argv += [flag, *value]
    argv += _parameter_flags(spec)
    return argv


def to_command(spec: AnalysisSpec, output: str = None) -> str:
    """``to_argv`` as one shell-quoted line, with ``-o output`` if given."""
    argv = ["phyforge", *to_argv(spec)]
    if output:
        argv += ["-o", output]
    return shlex.join(argv)


def _parameter_flags(spec: AnalysisSpec) -> List[str]:
    argv: List[str] = []
    for id_, override in spec.parameters.items():
        if not override.estimate:
            fix = (
                id_
                if override.init is None
                else f"{id_}={format_values(override.init)}"
            )
            argv += ["--fix", fix]
        elif override.init is not None:
            argv += ["--init", f"{id_}={format_values(override.init)}"]
        if override.prior is not None:
            argv += ["--prior", f"{id_}={override.prior}"]
        if override.mu is not None:
            argv += ["--vi-mu", f"{id_}={format_values(override.mu)}"]
        if override.sigma is not None:
            argv += ["--vi-sigma", f"{id_}={format_values(override.sigma)}"]
        if override.weight is not None:
            argv += ["--weight", f"{id_}={format_values(override.weight)}"]
    return argv


def _get(spec: Any, path: str) -> Any:
    for name in path.split("."):
        spec = getattr(spec, name)
    return spec


def _text(value: Any) -> str:
    if isinstance(value, float):
        return format_values(value)
    return str(value)
