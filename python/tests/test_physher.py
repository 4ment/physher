"""End-to-end: does physher actually accept what we generate?

Skipped unless a physher binary is found, so the suite still runs without a
build.  Set PHYSHER_BIN to point at one explicitly.
"""

import json
import os
import shutil
import subprocess

import pytest

from phyforge import AnalysisSpec, build
from phyforge.spec import (
    BootstrapOptions,
    CATSpec,
    ClockSpec,
    CoalescentSpec,
    DataSpec,
    MarginalOptions,
    MCMCOptions,
    SiteModelSpec,
    SubstitutionSpec,
    VariationalOptions,
)

REPOSITORY = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def find_physher():
    override = os.environ.get("PHYSHER_BIN")
    if override:
        return override if os.path.exists(override) else None
    local = os.path.join(REPOSITORY, "build", "physher")
    return local if os.path.exists(local) else shutil.which("physher")


PHYSHER = find_physher()

pytestmark = pytest.mark.skipif(PHYSHER is None, reason="no physher binary")


def run(document, tmp_path, alignment, tree_file):
    """Run ``document`` in a throwaway directory next to copies of the data."""
    for source in (alignment, tree_file):
        shutil.copy(source, tmp_path)
    document = json.loads(
        json.dumps(document)
        .replace(alignment, os.path.basename(alignment))
        .replace(tree_file, os.path.basename(tree_file))
    )
    config = tmp_path / "config.json"
    config.write_text(json.dumps(document))
    result = subprocess.run(
        [PHYSHER, "config.json"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, result.stderr[-2000:]
    return result.stdout


def likelihood(output):
    for line in output.splitlines():
        if line.startswith("treelikelihood:"):
            return float(line.split(":")[1])
    raise AssertionError(f"no likelihood in output:\n{output}")


def spec_for(alignment, tree_file, **changes):
    return AnalysisSpec(data=DataSpec(alignment=alignment, tree=tree_file), **changes)


@pytest.mark.parametrize(
    "model,expected",
    [
        ("JC69", -4336.6),
        ("K80", -4207.7),
        ("F81", -4304.8),
        ("HKY", -4174.8),
        ("GTR", -4158.8),
    ],
)
def test_nested_models_reach_their_optimum(
    alignment, tree_file, tmp_path, model, expected
):
    """Each model should reach its known optimum on the fluA example."""
    document = build(
        spec_for(alignment, tree_file, substitution=SubstitutionSpec(kind=model))
    )
    assert likelihood(run(document, tmp_path, alignment, tree_file)) == pytest.approx(
        expected, abs=0.1
    )


def test_parsimony_start_reaches_the_same_optimum(alignment, tree_file, tmp_path):
    """A parsimony start is a warm start, not a different model."""
    document = build(
        AnalysisSpec(
            data=DataSpec(alignment=alignment, tree=tree_file, init_parsimony=True)
        )
    )
    assert likelihood(run(document, tmp_path, alignment, tree_file)) == pytest.approx(
        -4336.6, abs=0.1
    )


def test_gamma_rates_improve_the_fit(alignment, tree_file, tmp_path):
    # The repository's own GTR-G4-ML baseline is -4081.86.
    document = build(
        spec_for(
            alignment,
            tree_file,
            substitution=SubstitutionSpec(kind="GTR"),
            site=SiteModelSpec(categories=4),
        )
    )
    assert likelihood(run(document, tmp_path, alignment, tree_file)) == pytest.approx(
        -4081.86, abs=0.1
    )


@pytest.mark.parametrize("distribution", ["gamma", "weibull", "lognormal"])
def test_every_rate_distribution_runs(alignment, tree_file, tmp_path, distribution):
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(categories=4, distribution=distribution),
        )
    )
    assert likelihood(run(document, tmp_path, alignment, tree_file)) < -4000


@pytest.mark.parametrize(
    "parameterization",
    ["rate_increments", "rate_ratios", "rate_shape", "mean_contribution"],
)
def test_free_rates_reach_the_same_optimum(
    alignment, tree_file, tmp_path, parameterization
):
    """The four charts describe the same model, so they agree at the optimum.

    -4250.53 is what the JC69-D4*-ML.json examples reach.
    """
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=4,
                distribution="discrete",
                rate_parameterization=parameterization,
            ),
        )
    )
    assert likelihood(run(document, tmp_path, alignment, tree_file)) == pytest.approx(
        -4250.53, abs=0.1
    )


def test_em_reaches_the_free_rate_optimum(alignment, tree_file, tmp_path):
    """EM is a different route to the same mixture, so it lands in the same place.

    -4250.53 is what the coordinate sweep above reaches, and what
    JC69-D4-EM-ML.json records.
    """
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=4, distribution="discrete", free_rate_em=True
            ),
        )
    )
    assert likelihood(run(document, tmp_path, alignment, tree_file)) == pytest.approx(
        -4250.53, abs=0.1
    )


def test_cat_reports_a_mixture_likelihood(alignment, tree_file, tmp_path):
    """The comparable number, which JC69-CAT20-ML.json puts at -4257.31.

    Parsimony is turned off to stay comparable with that baseline, which
    starts flat.  CAT reads its per-pattern rate assignment off the starting
    branch lengths, so unlike the other models the start picks the optimum:
    from parsimony this run settles at -4265.2 instead.
    """
    document = build(
        AnalysisSpec(
            data=DataSpec(alignment=alignment, tree=tree_file, init_parsimony=False),
            site=SiteModelSpec(
                categories=20,
                distribution="cat",
                cat=CATSpec(assignment="posterior_mean", prior="npmle"),
            ),
        )
    )
    output = run(document, tmp_path, alignment, tree_file)
    line = next(line for line in output.splitlines() if line.startswith("CAT mixture"))
    assert float(line.split()[3]) == pytest.approx(-4257.31, abs=0.5)


@pytest.mark.parametrize("clock", ["strict", "ncln"])
def test_time_trees_run(alignment, tree_file, tmp_path, clock):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="map",
            substitution=SubstitutionSpec(kind="HKY"),
            clock=ClockSpec(kind=clock),
            coalescent=CoalescentSpec(kind="constant"),
        )
    )
    run(document, tmp_path, alignment, tree_file)


@pytest.mark.parametrize(
    "coalescent,extra",
    [
        ("constant", {}),
        ("exponential", {}),
        ("skyride", {}),
        ("skygrid", {"grid": 10, "cutoff": 20.0}),
        ("skyglide", {"grid": 10, "cutoff": 20.0}),
    ],
)
def test_every_coalescent_runs(alignment, tree_file, tmp_path, coalescent, extra):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="advi",
            clock=ClockSpec(kind="strict"),
            coalescent=CoalescentSpec(kind=coalescent, **extra),
            variational=VariationalOptions(iterations=200, samples=20),
        )
    )
    run(document, tmp_path, alignment, tree_file)


def test_mcmc_runs_and_samples(alignment, tree_file, tmp_path):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            stem="chain",
            substitution=SubstitutionSpec(kind="HKY"),
            clock=ClockSpec(kind="strict"),
            coalescent=CoalescentSpec(kind="constant"),
            mcmc=MCMCOptions(length=5000, every=1000),
        )
    )
    run(document, tmp_path, alignment, tree_file)
    lines = (tmp_path / "chain.log").read_text().splitlines()
    assert len(lines) > 2
    assert lines[0].split("\t")[:3] == ["iter", "joint", "treelikelihood"]


# --------------------------------------------------------------------------
# marginal likelihood
# --------------------------------------------------------------------------


def marginal(output, label):
    """The number physher printed for one estimator, by the line it names."""
    for line in output.splitlines():
        if line.startswith(label):
            return float(line.split(":")[1])
    raise AssertionError(f"no {label!r} in output:\n{output}")


def test_stepping_stone_runs_the_whole_ladder(alignment, tree_file, tmp_path):
    """Every rung is run and read back, and the two estimators agree roughly."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            stem="ss",
            substitution=SubstitutionSpec(kind="HKY"),
            mcmc=MCMCOptions(length=5000, every=50),
            marginal=MarginalOptions(algorithms=("ss", "ps"), steps=4),
        )
    )
    output = run(document, tmp_path, alignment, tree_file)
    # One trace per rung, named "0ss.log" .. "3ss.log" by physher.
    assert all((tmp_path / f"{i}ss.log").exists() for i in range(4))
    stepping = marginal(output, "Stepping stone marginal likelihood")
    path = marginal(output, "Path sampling marginal likelihood")
    assert stepping == pytest.approx(path, rel=0.1)


def test_the_harmonic_mean_reads_a_plain_chain(alignment, tree_file, tmp_path):
    """No ladder: one chain, and an estimator that passes over its log."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            stem="hm",
            substitution=SubstitutionSpec(kind="HKY"),
            mcmc=MCMCOptions(length=5000, every=50),
            marginal=MarginalOptions(algorithms=("shm",)),
        )
    )
    output = run(document, tmp_path, alignment, tree_file)
    assert (tmp_path / "hm.log").exists()
    # The harmonic mean is famously biased upwards, but it is bounded below by
    # nothing useful; all this asks is that it read a real trace.
    assert marginal(output, "Smoothed harmonic mean") < 0


def test_importance_sampling_bounds_the_elbo(alignment, tree_file, tmp_path):
    """log p(y) >= ELBO, so a converged approximation should show it."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="advi",
            substitution=SubstitutionSpec(kind="HKY"),
            variational=VariationalOptions(iterations=2000, samples=0),
            marginal=MarginalOptions(algorithms=("is",), samples=1000),
        )
    )
    output = run(document, tmp_path, alignment, tree_file)
    # The optimizer's last "<iteration>\t<elbo>  <rate>" line.
    last = [line for line in output.splitlines() if line.startswith("2000\t")][-1]
    elbo = float(last.split("\t")[1].split()[0])
    assert marginal(output, "Marginal likelihood using IS") > elbo


# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------


def test_the_observed_replicate_reproduces_the_estimate(alignment, tree_file, tmp_path):
    """include_observed refits on the real weights, so it lands on the MLE."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            stem="boot",
            substitution=SubstitutionSpec(kind="HKY"),
            bootstrap=BootstrapOptions(replicates=3, include_observed=True),
        )
    )
    output = run(document, tmp_path, alignment, tree_file)
    rows = (tmp_path / "boot.bootstrap.csv").read_text().splitlines()
    header, body = rows[0].split("\t"), [row.split("\t") for row in rows[1:]]
    assert len(body) == 4  # 3 replicates, plus the observed control
    column = header.index("treelikelihood")
    assert float(body[0][column]) == pytest.approx(likelihood(output), abs=1e-4)
    # A resampled alignment is a different dataset, so it fits differently.
    assert all(float(row[column]) != float(body[0][column]) for row in body[1:])


def test_bootstrap_replicates_vary(alignment, tree_file, tmp_path):
    document = build(
        spec_for(
            alignment,
            tree_file,
            stem="boot",
            bootstrap=BootstrapOptions(replicates=5),
        )
    )
    run(document, tmp_path, alignment, tree_file)
    rows = (tmp_path / "boot.bootstrap.csv").read_text().splitlines()
    assert len(rows) == 6  # a header, then one row per replicate
    trees = (tmp_path / "boot.bootstrap.trees").read_text().splitlines()
    assert len(trees) == 5
