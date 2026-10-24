"""The spec and its validation."""

import pytest

from phyforge import AnalysisSpec, check
from phyforge.spec import (
    BootstrapOptions,
    CATSpec,
    ClockSpec,
    CoalescentSpec,
    DataSpec,
    MarginalOptions,
    OptimizerOptions,
    SiteModelSpec,
    SubstitutionSpec,
)


def base(alignment, tree_file, **changes):
    return AnalysisSpec(data=DataSpec(alignment=alignment, tree=tree_file), **changes)


def test_round_trips_through_a_dict(alignment, tree_file):
    spec = base(alignment, tree_file, method="advi")
    assert AnalysisSpec.from_dict(spec.to_dict()) == spec


def test_valid_spec_has_no_problems(alignment, tree_file):
    assert check(base(alignment, tree_file)) == []


def test_missing_alignment_is_reported(tree_file):
    problems = check(AnalysisSpec(data=DataSpec(tree=tree_file)))
    assert [problem.field for problem in problems] == ["alignment"]


def test_parsimony_is_accepted_on_an_unrooted_tree(alignment, tree_file):
    spec = AnalysisSpec(
        data=DataSpec(alignment=alignment, tree=tree_file, init_parsimony=True)
    )
    assert check(spec) == []


def test_parsimony_is_refused_on_a_time_tree(alignment, tree_file):
    spec = AnalysisSpec(
        data=DataSpec(alignment=alignment, tree=tree_file, init_parsimony=True),
        clock=ClockSpec(kind="strict"),
    )
    assert [problem.field for problem in check(spec)] == ["init_parsimony"]


def test_parsimony_and_keep_branch_lengths_conflict(alignment, tree_file):
    spec = AnalysisSpec(
        data=DataSpec(
            alignment=alignment,
            tree=tree_file,
            init_parsimony=True,
            keep_branch_lengths=True,
        )
    )
    assert [problem.field for problem in check(spec)] == ["init_parsimony"]


def test_parsimony_is_on_by_default_on_an_unrooted_tree(alignment, tree_file):
    spec = AnalysisSpec(data=DataSpec(alignment=alignment, tree=tree_file))
    assert spec.data.init_parsimony is None
    assert spec.uses_parsimony is True


@pytest.mark.parametrize(
    "spec",
    [
        AnalysisSpec(data=DataSpec(keep_branch_lengths=True)),
        AnalysisSpec(clock=ClockSpec(kind="strict")),
        AnalysisSpec(data=DataSpec(init_parsimony=False)),
    ],
)
def test_the_default_stands_aside_where_parsimony_does_not_apply(spec):
    assert spec.uses_parsimony is False


def test_the_default_is_not_a_conflict(alignment, tree_file):
    # left unset, parsimony never argues with the flags it is refused next to
    conflicting = (
        AnalysisSpec(
            data=DataSpec(alignment=alignment, tree=tree_file, keep_branch_lengths=True)
        ),
        AnalysisSpec(
            data=DataSpec(alignment=alignment, tree=tree_file),
            method="map",
            clock=ClockSpec(kind="strict"),
        ),
    )
    for spec in conflicting:
        assert [
            problem.field
            for problem in check(spec)
            if problem.field.startswith("init_parsimony")
        ] == []


def test_parsimony_floor_needs_parsimony(alignment, tree_file):
    spec = AnalysisSpec(
        data=DataSpec(
            alignment=alignment,
            tree=tree_file,
            init_parsimony=False,
            init_parsimony_min=1e-8,
        )
    )
    assert [problem.field for problem in check(spec)] == ["init_parsimony_min"]


def test_parsimony_floor_must_not_be_negative(alignment, tree_file):
    spec = AnalysisSpec(
        data=DataSpec(
            alignment=alignment,
            tree=tree_file,
            init_parsimony=True,
            init_parsimony_min=-1.0,
        )
    )
    assert [problem.field for problem in check(spec)] == ["init_parsimony_min"]


def test_unknown_model_is_reported(alignment, tree_file):
    problems = check(
        base(alignment, tree_file, substitution=SubstitutionSpec(kind="TN93"))
    )
    assert problems[0].field == "model"


def test_coalescent_needs_a_clock(alignment, tree_file):
    problems = check(
        base(
            alignment,
            tree_file,
            method="mcmc",
            coalescent=CoalescentSpec(kind="constant"),
        )
    )
    assert "needs a clock" in str(problems[0])


def test_skygrid_needs_a_grid_and_cutoff(alignment, tree_file):
    problems = check(
        base(
            alignment,
            tree_file,
            method="mcmc",
            clock=ClockSpec(kind="strict"),
            coalescent=CoalescentSpec(kind="skygrid"),
        )
    )
    assert {problem.field for problem in problems} == {"grid", "cutoff"}


def test_codon_models_are_refused(alignment, tree_file):
    problems = check(
        base(alignment, tree_file, substitution=SubstitutionSpec(kind="GY94"))
    )
    assert problems[0].field == "model"


def test_free_rates_refuse_the_gradient_methods(alignment, tree_file):
    """No analytic gradient, and physher does not guard against one."""
    site = SiteModelSpec(categories=4, distribution="discrete")
    problems = check(base(alignment, tree_file, method="advi", site=site))
    assert problems[0].field == "method"
    problems = check(
        base(
            alignment,
            tree_file,
            site=site,
            optimizer=OptimizerOptions(algorithm="bfgs"),
        )
    )
    assert problems[0].field == "algorithm"


def test_em_refuses_what_physher_declines_to_run(alignment, tree_file):
    """Each of these makes physher skip the EM step rather than fail."""
    problems = check(
        base(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=4,
                distribution="discrete",
                rate_parameterization="rate_shape",
                free_rate_em=True,
            ),
        )
    )
    assert [problem.field for problem in problems] == ["em"]

    problems = check(
        base(
            alignment,
            tree_file,
            method="mcmc",
            site=SiteModelSpec(
                categories=4, distribution="discrete", free_rate_em=True
            ),
        )
    )
    assert "em" in {problem.field for problem in problems}

    problems = check(
        base(
            alignment,
            tree_file,
            clock=ClockSpec(kind="strict"),
            site=SiteModelSpec(
                categories=4, distribution="discrete", free_rate_em=True
            ),
        )
    )
    assert "em" in {problem.field for problem in problems}


def test_em_needs_a_mixture_to_split(alignment, tree_file):
    problems = check(
        base(
            alignment,
            tree_file,
            site=SiteModelSpec(categories=4, distribution="gamma", free_rate_em=True),
        )
    )
    assert [problem.field for problem in problems] == ["em"]


def test_em_accepts_the_free_rate_default(alignment, tree_file):
    assert (
        check(
            base(
                alignment,
                tree_file,
                site=SiteModelSpec(
                    categories=4, distribution="discrete", free_rate_em=True
                ),
            )
        )
        == []
    )


def test_ordered_free_rates_have_no_invariant_form(alignment, tree_file):
    problems = check(
        base(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=4,
                distribution="discrete",
                invariant=True,
                rate_parameterization="rate_ratios",
            ),
        )
    )
    assert problems[0].field == "invariant"


def test_cat_is_maximum_likelihood_only(alignment, tree_file):
    problems = check(
        base(
            alignment,
            tree_file,
            method="mcmc",
            site=SiteModelSpec(categories=20, distribution="cat"),
        )
    )
    assert [problem.field for problem in problems] == ["method"]


def test_cat_needs_nucleotides(alignment, tree_file):
    problems = check(
        base(
            alignment,
            tree_file,
            substitution=SubstitutionSpec(kind="LG"),
            site=SiteModelSpec(categories=20, distribution="cat"),
        )
    )
    assert "nucleotides only" in str(problems[0])


def test_cat_prior_is_npmle_or_a_shape(alignment, tree_file):
    def problems_for(prior):
        return check(
            base(
                alignment,
                tree_file,
                site=SiteModelSpec(
                    categories=20, distribution="cat", cat=CATSpec(prior=prior)
                ),
            )
        )

    assert problems_for("npmle") == []
    assert problems_for("3.0") == []
    assert problems_for("gamma")[0].field == "cat_prior"


def test_a_probe_grid_only_means_something_to_the_posterior_mean(alignment, tree_file):
    problems = check(
        base(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=20,
                distribution="cat",
                cat=CATSpec(assignment="argmax", probe=100),
            ),
        )
    )
    assert problems[0].field == "cat_probe"


def test_cpo_needs_a_posterior_sample(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.mcmc.cpo = True
    assert [problem.field for problem in check(spec)] == ["cpo"]


def test_cpo_and_stepping_stone_are_refused(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.mcmc.cpo = True
    spec.marginal = MarginalOptions(algorithms=("ss",), steps=10)
    assert [problem.field for problem in check(spec)] == ["cpo"]


def test_cpo_burnin_cannot_swallow_the_chain(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.mcmc.cpo = True
    spec.mcmc.cpo_burnin = spec.mcmc.pattern_samples()
    assert [problem.field for problem in check(spec)] == ["cpo_burnin"]


def test_cpo_off_is_never_checked(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.mcmc.cpo_burnin = -1
    assert check(spec) == []


def test_lppd_needs_a_posterior_sample(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.mcmc.lppd = True
    assert [problem.field for problem in check(spec)] == ["lppd"]


def test_lppd_burnin_cannot_swallow_the_chain(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.mcmc.lppd = True
    spec.mcmc.lppd_burnin = spec.mcmc.pattern_samples()
    assert [problem.field for problem in check(spec)] == ["lppd_burnin"]


def test_lppd_and_cpo_are_checked_independently(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.mcmc.cpo = True
    spec.mcmc.lppd = True
    spec.marginal = MarginalOptions(algorithms=("ss",), steps=10)
    assert [problem.field for problem in check(spec)] == ["cpo", "lppd"]


def test_lppd_off_is_never_checked(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.mcmc.lppd_burnin = -1
    assert check(spec) == []


# --------------------------------------------------------------------------
# marginal likelihood
# --------------------------------------------------------------------------


def test_no_marginal_estimator_is_not_a_problem(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    assert check(spec) == []


def test_unknown_estimator_is_refused(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.marginal = MarginalOptions(algorithms=("gss",))
    assert [problem.field for problem in check(spec)] == ["marginal"]


def test_stepping_stone_needs_a_ladder(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.marginal = MarginalOptions(algorithms=("ss",))
    assert [problem.field for problem in check(spec)] == ["steps"]


def test_a_ladder_nothing_reads_is_refused(alignment, tree_file):
    """Every chain but the last would be thrown away."""
    spec = base(alignment, tree_file, method="mcmc")
    spec.marginal = MarginalOptions(algorithms=("hm",), steps=8)
    assert [problem.field for problem in check(spec)] == ["steps"]


def test_the_harmonic_mean_needs_no_ladder(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.marginal = MarginalOptions(algorithms=("shm",))
    assert check(spec) == []


def test_chain_estimators_need_a_chain(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.marginal = MarginalOptions(algorithms=("ss",), steps=8)
    assert [problem.field for problem in check(spec)] == ["marginal"]


def test_importance_sampling_needs_an_approximation(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.marginal = MarginalOptions(algorithms=("is",))
    assert [problem.field for problem in check(spec)] == ["marginal"]


def test_importance_sampling_runs_under_advi(alignment, tree_file):
    spec = base(alignment, tree_file, method="advi")
    spec.marginal = MarginalOptions(algorithms=("is",))
    assert check(spec) == []


def test_marginal_burnin_cannot_swallow_the_chain(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.marginal = MarginalOptions(
        algorithms=("hm",), burnin=spec.mcmc.logged_samples()
    )
    assert [problem.field for problem in check(spec)] == ["marginal_burnin"]


def test_the_ladder_shape_is_checked(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.marginal = MarginalOptions(algorithms=("ss",), steps=4, ladder_alpha=0)
    assert [problem.field for problem in check(spec)] == ["ladder_alpha"]


# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------


def test_no_bootstrap_is_not_a_problem(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    assert check(spec) == []


def test_the_bootstrap_needs_a_point_estimate(alignment, tree_file):
    spec = base(alignment, tree_file, method="mcmc")
    spec.bootstrap = BootstrapOptions(replicates=100)
    assert [problem.field for problem in check(spec)] == ["bootstrap"]


@pytest.mark.parametrize("method", ["ml", "map"])
def test_the_bootstrap_runs_after_either_point_estimate(alignment, tree_file, method):
    spec = base(alignment, tree_file, method=method)
    spec.bootstrap = BootstrapOptions(replicates=100)
    assert check(spec) == []


def test_negative_replicates_are_refused(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.bootstrap = BootstrapOptions(replicates=-1)
    assert [problem.field for problem in check(spec)] == ["bootstrap"]


@pytest.mark.parametrize("threads", [0, -2])
def test_threads_below_one_are_refused(alignment, tree_file, threads):
    spec = base(alignment, tree_file, method="ml")
    spec.bootstrap = BootstrapOptions(replicates=10, threads=threads)
    assert [problem.field for problem in check(spec)] == ["bootstrap"]


def test_threads_need_independent_replicates(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.bootstrap = BootstrapOptions(replicates=10, threads=4, reset=False)
    assert [problem.field for problem in check(spec)] == ["bootstrap"]


def test_threads_need_compacted_pattern_sets(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.bootstrap = BootstrapOptions(replicates=10, threads=4, compact=False)
    assert [problem.field for problem in check(spec)] == ["bootstrap"]


def test_threads_are_valid_with_the_defaults(alignment, tree_file):
    spec = base(alignment, tree_file, method="ml")
    spec.bootstrap = BootstrapOptions(replicates=10, threads=4)
    assert check(spec) == []
