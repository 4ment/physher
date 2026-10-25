"""Writing a spec back as a command line, and reading it in again."""

import dataclasses
import shlex

import pytest

from phyforge import AnalysisSpec
from phyforge.cli.args import parse
from phyforge.cli.argv import OPTIONS, UNEXPRESSED, to_argv, to_command
from phyforge.spec import ParameterOverride, PriorChoice

#: A value other than the CLI default for every option that is not a number
#: or a switch, which the test cannot make up for itself.
SAMPLES = {
    "stem": "run",
    "data.init_parsimony": False,
    "data.init_parsimony_min": 0.001,
    "substitution.kind": "GTR",
    "substitution.frequencies": "equal",
    "site.distribution": "weibull",
    "site.quadrature": "gausslaguerre",
    "site.rate_parameterization": "rate_ratios",
    "site.cat.assignment": "argmax",
    "site.cat.prior": "npmle",
    "site.cat.probe": 50,
    "clock.kind": "ncln",
    "clock.rate": 0.002,
    "clock.dates.source": "csv",
    "clock.dates.path": "dates file.csv",
    "clock.dates.regex": r"\|(\d+)$",
    "clock.dates.date_format": "yyyy-MM-dd",
    "coalescent.kind": "skygrid",
    "coalescent.grid": 12,
    "coalescent.cutoff": 30.5,
    "priors.branch_lengths": "none",
    "priors.clock": "exponential",
    "optimizer.algorithm": "bfgs",
    "mcmc.cpo_burnin": 7,
    "mcmc.lppd_burnin": 8,
    "variational.checkpoint": "checkpoint.json",
    "variational.update": "sgd",
    "marginal.algorithms": ("ss", "ps"),
    "marginal.burnin": 9,
    "marginal.ladder": "uniform",
}


def get(spec, path):
    for name in path.split("."):
        spec = getattr(spec, name)
    return spec


def with_value(spec, path, value):
    """A copy of ``spec`` with the field at ``path`` set to ``value``."""
    head, _, rest = path.partition(".")
    if not rest:
        return dataclasses.replace(spec, **{head: value})
    return dataclasses.replace(
        spec, **{head: with_value(getattr(spec, head), rest, value)}
    )


def different(path, default):
    if path in SAMPLES:
        return SAMPLES[path]
    if isinstance(default, bool):
        return not default
    if isinstance(default, int):
        return default + 3
    if isinstance(default, float):
        return default * 2 + 0.25
    raise AssertionError(f"add a sample value for {path}")


def leaf_paths(spec, prefix=""):
    for field in dataclasses.fields(spec):
        value = getattr(spec, field.name)
        path = f"{prefix}{field.name}"
        if dataclasses.is_dataclass(value):
            yield from leaf_paths(value, f"{path}.")
        else:
            yield path


def defaults():
    return parse(["ml", "-i", "aln.fa", "-t", "tree.nwk"])


def test_the_defaults_need_no_options():
    assert to_argv(defaults()) == ["ml", "-i", "aln.fa", "-t", "tree.nwk"]


def test_every_field_has_a_flag():
    covered = {path for path, _, _ in OPTIONS}
    covered |= set(UNEXPRESSED) | {
        "method",
        "data.alignment",
        "data.tree",
        "parameters",
    }
    assert set(leaf_paths(AnalysisSpec())) <= covered


@pytest.mark.parametrize("path", [path for path, _, _ in OPTIONS])
def test_each_option_round_trips(path):
    spec = with_value(defaults(), path, different(path, get(defaults(), path)))
    if path == "marginal.steps":
        # --steps on its own still means stepping stone; a ladder with no
        # estimator to read it is refused by validate anyway.
        spec = with_value(spec, "marginal.algorithms", ("ss",))
    assert parse(to_argv(spec)) == spec


@pytest.mark.parametrize("method", ["ml", "map", "mcmc", "advi"])
def test_the_method_round_trips(method):
    spec = dataclasses.replace(defaults(), method=method)
    assert parse(to_argv(spec)) == spec


def test_overrides_round_trip():
    spec = dataclasses.replace(
        defaults(),
        method="mcmc",
        parameters={
            "kappa": ParameterOverride(
                init=2.5,
                prior=PriorChoice("lognormal", {"mu": 1.0, "sigma": 1.25}),
                weight=10.0,
            ),
            "shape": ParameterOverride(estimate=False, init=0.5),
            "frequencies": ParameterOverride(estimate=False, init=[0.1, 0.2, 0.3, 0.4]),
            "rates": ParameterOverride(
                prior=PriorChoice(
                    "dirichlet", {"alpha": [1.0, 2.0, 1.0, 1.0, 2.0, 1.0]}
                )
            ),
            "theta": ParameterOverride(estimate=False),
            "tree.distances": ParameterOverride(prior=PriorChoice("none")),
        },
    )
    assert parse(to_argv(spec)) == spec


def test_variational_overrides_round_trip():
    spec = dataclasses.replace(
        defaults(),
        method="advi",
        parameters={"kappa": ParameterOverride(mu=0.1, sigma=[0.2])},
    )
    reparsed = parse(to_argv(spec))
    assert reparsed.parameters["kappa"].mu == 0.1
    assert reparsed.parameters["kappa"].sigma == 0.2


def test_floats_survive_exactly():
    spec = with_value(defaults(), "clock.rate_init", 0.1 + 0.2)
    assert parse(to_argv(spec)).clock.rate_init == 0.1 + 0.2


def test_the_command_is_shell_quoted():
    spec = dataclasses.replace(
        defaults(),
        method="map",
        parameters={
            "kappa": ParameterOverride(
                prior=PriorChoice("gamma", {"shape": 2, "rate": 1})
            )
        },
    )
    command = to_command(spec, output="out.json")
    assert shlex.split(command)[1:] == to_argv(spec) + ["-o", "out.json"]
    assert "'kappa=gamma(shape=2,rate=1)'" in command
