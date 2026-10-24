"""Building documents, checked against the schema physher enforces."""

import pytest

from phyforge import AnalysisSpec, build, make_context
from phyforge.build.analysis import build as build_document
from phyforge.spec import (
    BootstrapOptions,
    CATSpec,
    ClockSpec,
    CoalescentSpec,
    DataSpec,
    MarginalOptions,
    MCMCOptions,
    OptimizerOptions,
    SiteModelSpec,
    SubstitutionSpec,
)


def spec_for(alignment, tree_file, **changes):
    return AnalysisSpec(data=DataSpec(alignment=alignment, tree=tree_file), **changes)


def find(node, predicate):
    """Depth-first search for the first node satisfying ``predicate``."""
    if isinstance(node, dict):
        if predicate(node):
            return node
        for value in node.values():
            found = find(value, predicate)
            if found is not None:
                return found
    elif isinstance(node, list):
        for value in node:
            found = find(value, predicate)
            if found is not None:
                return found
    return None


def of_type(type_name):
    return lambda node: node.get("type") == type_name


def walk(node):
    if isinstance(node, dict):
        yield node
        for value in node.values():
            yield from walk(value)
    elif isinstance(node, list):
        for value in node:
            yield from walk(value)


# --------------------------------------------------------------------------
# schema rules physher enforces in json_validate
# --------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["ml", "map", "mcmc", "advi"])
def test_parameters_never_use_the_forbidden_value_keys(alignment, tree_file, method):
    # A coalescent is a prior, so it only belongs on the Bayesian methods.
    coalescent = CoalescentSpec(kind="constant" if method != "ml" else None)
    document = build(
        spec_for(
            alignment,
            tree_file,
            method=method,
            substitution=SubstitutionSpec(kind="GTR"),
            site=SiteModelSpec(categories=4),
            clock=ClockSpec(kind="strict"),
            coalescent=coalescent,
        )
    )
    for node in walk(document):
        if node.get("type") == "parameter":
            assert "value" not in node and "values" not in node


def test_models_use_kind_not_model(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            substitution=SubstitutionSpec(kind="HKY"),
            clock=ClockSpec(kind="strict"),
            coalescent=CoalescentSpec(kind="constant"),
        )
    )
    for type_name in ("substitutionmodel", "branchmodel", "coalescent"):
        node = find(document, of_type(type_name))
        assert "kind" in node and "model" not in node


TREE_LENGTH_COLUMN = {"ref": "@tree", "quantity": "treeLength"}


def trace_columns(document):
    return find(document, lambda node: "columns" in node)["columns"]


def test_ml_logs_the_tree_length_of_an_unrooted_tree(alignment, tree_file):
    for method in ("ml", "map"):
        document = build(spec_for(alignment, tree_file, method=method))
        assert TREE_LENGTH_COLUMN in trace_columns(document)


def test_a_time_tree_does_not_log_a_tree_length(alignment, tree_file):
    document = build(spec_for(alignment, tree_file, clock=ClockSpec(kind="strict")))
    assert TREE_LENGTH_COLUMN not in trace_columns(document)


def test_only_the_optimizers_log_a_tree_length(alignment, tree_file):
    """MCMC and ADVI traces are left as they were."""
    for method in ("mcmc", "advi"):
        document = build(spec_for(alignment, tree_file, method=method))
        for node in walk(document):
            assert TREE_LENGTH_COLUMN not in node.get("columns", [])


def test_loggers_use_columns_not_models(alignment, tree_file):
    document = build(spec_for(alignment, tree_file, method="mcmc"))
    for node in walk(document):
        if node.get("type") == "logger":
            assert "models" not in node and "x" not in node
            # a logger is either a tree logger or a column logger, never both
            assert ("tree" in node) != ("columns" in node)


# --------------------------------------------------------------------------
# model structure
# --------------------------------------------------------------------------


def test_unrooted_tree_has_branch_lengths(alignment, tree_file):
    document = build(spec_for(alignment, tree_file))
    tree = find(document, of_type("tree"))
    assert tree["branch_lengths"]["dimension"] == 2 * 69 - 3
    assert "time" not in tree


def test_parsimony_initializer_references_the_site_pattern(alignment, tree_file):
    # unset rather than True: the initializer is what an unrooted run gets by
    # default
    document = build(spec_for(alignment, tree_file))
    tree = find(document, of_type("tree"))
    assert tree["initializer"] == {
        "algorithm": "parsimony",
        "sitepattern": "@patterns",
    }
    # the referenced site pattern has to be declared before the tree reads it
    likelihood = find(document, of_type("treelikelihood"))
    keys = list(likelihood)
    assert keys.index("sitepattern") < keys.index("tree")
    assert likelihood["sitepattern"]["id"] == "patterns"


def test_parsimony_initializer_carries_a_floor(alignment, tree_file):
    document = build(
        AnalysisSpec(
            data=DataSpec(
                alignment=alignment,
                tree=tree_file,
                init_parsimony=True,
                init_parsimony_min=1e-8,
            )
        )
    )
    assert find(document, of_type("tree"))["initializer"]["min"] == 1e-8


def test_no_initializer_when_parsimony_is_turned_off(alignment, tree_file):
    document = build(
        AnalysisSpec(
            data=DataSpec(alignment=alignment, tree=tree_file, init_parsimony=False)
        )
    )
    assert "initializer" not in find(document, of_type("tree"))


@pytest.mark.parametrize(
    "changes",
    [
        {"data": {"keep_branch_lengths": True}},
        {"clock": ClockSpec(kind="strict")},
    ],
)
def test_no_initializer_where_parsimony_does_not_apply(alignment, tree_file, changes):
    data = DataSpec(alignment=alignment, tree=tree_file, **changes.pop("data", {}))
    document = build(AnalysisSpec(data=data, method="mcmc", **changes))
    assert "initializer" not in find(document, of_type("tree"))


def test_time_tree_is_reparameterized(alignment, tree_file):
    document = build(
        spec_for(alignment, tree_file, method="mcmc", clock=ClockSpec(kind="strict"))
    )
    tree = find(document, of_type("tree"))
    assert tree["time"] is True
    assert tree["proportions"]["x"]["dimension"] == 69 - 2
    # both height parameters are transformed, whatever the method
    assert tree["root_height"]["x"]["type"] == "parameter"


def test_time_tree_carries_dates(alignment, tree_file):
    document = build(
        spec_for(alignment, tree_file, method="mcmc", clock=ClockSpec(kind="strict"))
    )
    taxa = find(document, of_type("Taxa"))["taxa"]
    assert len(taxa) == 69
    assert taxa[0]["attributes"]["date"] == 1981.0


def test_jc69_frequencies_are_pinned(alignment, tree_file):
    document = build(spec_for(alignment, tree_file))
    frequencies = find(document, of_type("Simplex"))
    assert frequencies["x"] == [0.25] * 4


def test_gtr_estimates_rates_and_frequencies(alignment, tree_file):
    context = make_context(
        spec_for(alignment, tree_file, substitution=SubstitutionSpec(kind="GTR"))
    )
    build_document(context.spec, context)
    roles = {free.role for free in context.free_parameters}
    assert {"rates", "frequencies", "branch_lengths"} <= roles


def test_k80_is_emitted_as_hky_with_equal_frequencies(alignment, tree_file):
    # physher's own "k80" kind does not optimize; K80 is HKY with pinned
    # frequencies, so that is what gets written.
    document = build(
        spec_for(alignment, tree_file, substitution=SubstitutionSpec(kind="K80"))
    )
    model = find(document, of_type("substitutionmodel"))
    assert model["kind"] == "hky"
    assert model["frequencies"]["x"] == [0.25] * 4


def test_lognormal_site_model_uses_scale(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(categories=4, distribution="lognormal"),
        )
    )
    site_model = find(document, of_type("sitemodel"))
    assert "scale" in site_model and "shape" not in site_model


@pytest.mark.parametrize(
    "parameterization,key,dimension,simplex",
    [
        ("rate_increments", "rate_increments", 4, False),
        ("rate_ratios", "rate_ratios", 3, False),
        ("rate_shape", "rate_shape", 4, True),
        ("mean_contribution", "mean_contribution", 4, True),
    ],
)
def test_free_rates_match_their_parameterization(
    alignment, tree_file, parameterization, key, dimension, simplex
):
    """physher checks the parameter against the key that named it."""
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
    site_model = find(document, of_type("sitemodel"))
    assert site_model["distribution"] == "discrete"
    rates = site_model[key]
    assert (rates["type"] == "Simplex") is simplex
    if simplex:
        assert rates["x"]["dimension"] == dimension - 1
    else:
        assert rates["dimension"] == dimension
    # every chart needs the weights, and only the ratios need a top rate
    assert site_model["proportions"]["x"]["dimension"] == 3
    assert ("top_rate" in site_model) is (parameterization == "rate_ratios")


def test_invariant_free_rates_lengthen_only_the_weights(alignment, tree_file):
    # The invariant class sits at rate 0 and contributes nothing to the
    # weighted mean, so it is one more element of the simplex but not of the
    # increments.
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(categories=4, distribution="discrete", invariant=True),
        )
    )
    site_model = find(document, of_type("sitemodel"))
    assert site_model["invariant"] is True
    assert site_model["rate_increments"]["dimension"] == 4
    assert site_model["proportions"]["x"]["dimension"] == 4
    assert "proportion_invariant" not in site_model


def test_cat_reads_its_categories_from_the_site_pattern(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(categories=20, distribution="cat"),
        )
    )
    site_model = find(document, of_type("sitemodel"))
    assert site_model["sitepattern"] == "&patterns"
    assert site_model["categories"] == 20
    assert site_model["parameters"]["dimension"] == 20
    assert "distribution" not in site_model


def test_cat_rates_are_not_free_parameters(alignment, tree_file):
    """The cat estimator writes them, so nothing else may move them."""
    spec = spec_for(
        alignment, tree_file, site=SiteModelSpec(categories=20, distribution="cat")
    )
    context = make_context(spec)
    document = build_document(spec, context)
    assert [free.role for free in context.free_parameters] == ["branch_lengths"]
    # they are logged all the same, by asking the site model for them
    columns = document["physher"][1]["columns"]
    assert {"ref": "@sitemodel", "quantity": "rates"} in columns


def test_em_replaces_the_free_rate_sweep(alignment, tree_file):
    """EM moves the rates and the weights, so nothing else may sweep them."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=4, distribution="discrete", free_rate_em=True
            ),
        )
    )
    optimizer = document["physher"][0]
    assert [entry["algorithm"] for entry in optimizer["list"]] == ["serial", "em"]
    em = optimizer["list"][-1]
    assert em["model"] == "@treelikelihood"
    # the parameters are still built and still logged; only who moves them
    # changes
    assert "parameters" not in em
    site_model = find(document, of_type("sitemodel"))
    assert "rate_increments" in site_model and "proportions" in site_model


def test_em_leaves_the_other_sweeps_alone(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            substitution=SubstitutionSpec(kind="HKY"),
            site=SiteModelSpec(
                categories=4, distribution="discrete", free_rate_em=True
            ),
        )
    )
    optimizer = document["physher"][0]
    assert [entry["algorithm"] for entry in optimizer["list"]] == [
        "serial",
        "serial",
        "serial",
        "em",
    ]


def test_cat_schedules_a_reassignment_after_the_branch_lengths(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=20,
                distribution="cat",
                cat=CATSpec(assignment="posterior_mean", prior="npmle"),
            ),
        )
    )
    optimizer = document["physher"][0]
    assert [entry["algorithm"] for entry in optimizer["list"]] == ["serial", "cat"]
    assert optimizer["min"] == 2
    cat = optimizer["list"][1]
    assert cat["model"] == "@treelikelihood"
    assert cat["assignment"] == "posterior_mean" and cat["prior"] == "npmle"


def test_cat_reports_a_comparable_likelihood_last(alignment, tree_file):
    document = build(
        spec_for(
            alignment, tree_file, site=SiteModelSpec(categories=20, distribution="cat")
        )
    )
    report = document["physher"][-1]
    assert report["type"] == "cat" and report["assign"] is False

    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=20, distribution="cat", cat=CATSpec(mixture=False)
            ),
        )
    )
    assert not any(node.get("type") == "cat" for node in document["physher"])


def test_cat_prior_keeps_a_gamma_shape_a_number(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            site=SiteModelSpec(
                categories=20,
                distribution="cat",
                cat=CATSpec(assignment="argmax", prior="3"),
            ),
        )
    )
    assert document["physher"][0]["list"][1]["prior"] == 3.0


def test_free_rates_are_given_a_proper_prior(alignment, tree_file):
    # Both ordered charts are flat along theta -> c theta, so a flat prior
    # would leave an improper posterior.
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="map",
            site=SiteModelSpec(categories=4, distribution="discrete"),
        )
    )
    priors = {node["id"]: node for node in document["prior"]["distributions"]}
    assert priors["prior.site_rates"]["distribution"] == "exponential"
    assert priors["prior.site_proportions"]["distribution"] == "dirichlet"


def test_ncln_rates_span_every_branch(alignment, tree_file):
    document = build(
        spec_for(alignment, tree_file, method="map", clock=ClockSpec(kind="ncln"))
    )
    branch_model = find(document, of_type("branchmodel"))
    assert branch_model["rate"]["dimension"] == 2 * 69 - 2
    assert "location" in branch_model and "scale" in branch_model


# --------------------------------------------------------------------------
# method-specific structure
# --------------------------------------------------------------------------


def test_ml_has_no_prior(alignment, tree_file):
    document = build(spec_for(alignment, tree_file))
    assert "prior" not in document and "joint" not in document


def test_map_targets_the_joint(alignment, tree_file):
    document = build(spec_for(alignment, tree_file, method="map"))
    assert document["joint"]["distributions"] == ["@treelikelihood", "@prior"]
    assert document["physher"][0]["target"] == "@joint"


def test_map_has_no_jacobian(alignment, tree_file):
    # A MAP estimate is reported in the constrained space, so adding the
    # Jacobian would move the reported optimum.
    document = build(
        spec_for(alignment, tree_file, method="map", clock=ClockSpec(kind="strict"))
    )
    assert "jacobian" not in document


def test_mcmc_covers_every_free_parameter(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        method="mcmc",
        substitution=SubstitutionSpec(kind="GTR"),
        clock=ClockSpec(kind="strict"),
        coalescent=CoalescentSpec(kind="constant"),
    )
    context = make_context(spec)
    document = build_document(spec, context)
    operators = document["physher"][0]["operators"]
    moved = {operator["x"] for operator in operators}
    assert moved == {f"&{free.id}" for free in context.free_parameters}


def test_operator_weights_are_the_parameter_dimension(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        method="mcmc",
        substitution=SubstitutionSpec(kind="GTR"),
        clock=ClockSpec(kind="strict"),
    )
    context = make_context(spec)
    document = build_document(spec, context)
    dimensions = {free.id: free.dimension for free in context.free_parameters}
    weights = {
        operator["x"].lstrip("&"): operator["weight"]
        for operator in document["physher"][0]["operators"]
    }
    assert weights == dimensions
    # not all one: the branch rates and the clock rate cannot share a weight
    assert len(set(weights.values())) > 1


def test_cpo_logs_the_patterns_it_reads_back(alignment, tree_file):
    spec = spec_for(alignment, tree_file, method="mcmc")
    spec.mcmc.cpo = True
    document = build(spec)
    cpo = find(document, of_type("cpo"))
    logger = find(
        document,
        lambda node: any(
            column.get("quantity") == "patternLogLikelihood"
            for column in node.get("columns", [])
            if isinstance(column, dict)
        ),
    )
    assert logger["file"] == cpo["filename"]
    assert logger["columns"][0]["ref"] == cpo["model"]
    # a tenth of the 1000 rows a million-iteration chain logs every 1000
    assert cpo["burnin"] == 100


def test_lppd_logs_the_patterns_it_reads_back(alignment, tree_file):
    spec = spec_for(alignment, tree_file, method="mcmc")
    spec.mcmc.lppd = True
    document = build(spec)
    lppd = find(document, of_type("predictive"))
    logger = find(
        document,
        lambda node: any(
            column.get("quantity") == "patternLogLikelihood"
            for column in node.get("columns", [])
            if isinstance(column, dict)
        ),
    )
    assert logger["file"] == lppd["filename"]
    assert logger["columns"][0]["ref"] == lppd["model"]
    assert lppd["burnin"] == 100


def test_cpo_and_lppd_share_one_pattern_logger(alignment, tree_file):
    spec = spec_for(alignment, tree_file, method="mcmc")
    spec.mcmc.cpo = True
    spec.mcmc.lppd = True
    document = build(spec)
    loggers = [
        node
        for node in walk(document)
        if any(
            isinstance(column, dict)
            and column.get("quantity") == "patternLogLikelihood"
            for column in node.get("columns", [])
        )
    ]
    assert len(loggers) == 1
    cpo = find(document, of_type("cpo"))
    lppd = find(document, of_type("predictive"))
    assert cpo["filename"] == lppd["filename"] == loggers[0]["file"]


def test_without_cpo_nothing_logs_patterns(alignment, tree_file):
    document = build(spec_for(alignment, tree_file, method="mcmc"))
    assert find(document, of_type("cpo")) is None
    assert find(document, of_type("predictive")) is None
    assert not [
        node
        for node in walk(document)
        if any(
            isinstance(column, dict)
            and column.get("quantity") == "patternLogLikelihood"
            for column in node.get("columns", [])
        )
    ]


def test_mcmc_jacobian_lists_only_transformed_parameters(alignment, tree_file):
    spec = spec_for(alignment, tree_file, method="mcmc", clock=ClockSpec(kind="strict"))
    context = make_context(spec)
    document = build_document(spec, context)
    expected = [
        f"&{free.id}" for free in context.free_parameters if free.is_transformed
    ]
    assert document["jacobian"]["parameters"] == expected


def test_advi_has_one_factor_per_free_parameter(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        method="advi",
        substitution=SubstitutionSpec(kind="GTR"),
        clock=ClockSpec(kind="strict"),
    )
    context = make_context(spec)
    document = build_document(spec, context)
    factors = document["variational"]["distributions"]
    assert len(factors) == len(context.free_parameters)
    # every factor sits on an unconstrained value
    approximated = {factor["x"] for factor in factors}
    assert approximated == {f"&{free.free_id}" for free in context.free_parameters}


def test_advi_simplex_factor_is_one_dimension_smaller(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="advi",
            substitution=SubstitutionSpec(kind="GTR"),
        )
    )
    factor = find(
        document["variational"],
        lambda node: node.get("x") == "&frequencies.unconstrained",
    )
    assert factor["parameters"]["mu"]["dimension"] == 3


def test_gradient_optimizer_moves_unconstrained_parameters(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        substitution=SubstitutionSpec(kind="GTR"),
        optimizer=OptimizerOptions(algorithm="bfgs"),
    )
    context = make_context(spec)
    document = build_document(spec, context)
    moved = document["physher"][0]["parameters"]
    assert "&rates.unconstrained" in moved


def test_declaration_order_precedes_reference(alignment, tree_file):
    """physher builds top-level nodes in order, so ids must exist when used."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="advi",
            clock=ClockSpec(kind="strict"),
            coalescent=CoalescentSpec(kind="constant"),
        )
    )
    keys = [key for key in document if key != "physher"]
    assert keys.index("treelikelihood") < keys.index("prior")
    assert keys.index("prior") < keys.index("joint")
    assert keys.index("joint") < keys.index("variational")
    assert keys.index("variational") < keys.index("elbo")


# --------------------------------------------------------------------------
# marginal likelihood
# --------------------------------------------------------------------------


def test_no_marginal_estimator_by_default(alignment, tree_file):
    document = build(spec_for(alignment, tree_file, method="mcmc"))
    types = [node["type"] for node in document["physher"]]
    assert "mmcmc" not in types
    assert "marginallikelihood" not in types


def test_stepping_stone_wraps_the_chain_in_a_ladder(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            marginal=MarginalOptions(algorithms=("ss",), steps=8),
        )
    )
    ladder, estimator = document["physher"]
    assert ladder["type"] == "mmcmc"
    assert ladder["mcmc"]["type"] == "mcmc"
    assert estimator["type"] == "marginallikelihood"
    assert estimator["algorithm"] == "ss"


def test_the_two_ladders_agree(alignment, tree_file):
    """The chains and the estimator must be built at the same temperatures."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            marginal=MarginalOptions(
                algorithms=("ss",), steps=8, ladder="beta", ladder_alpha=0.25
            ),
        )
    )
    ladder, estimator = document["physher"]
    keys = ("steps", "distribution", "alpha", "beta")
    assert {key: ladder[key] for key in keys} == {key: estimator[key] for key in keys}
    assert ladder["alpha"] == 0.25


def test_a_uniform_ladder_carries_no_beta_shapes(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            marginal=MarginalOptions(algorithms=("ps",), steps=4, ladder="uniform"),
        )
    )
    ladder = document["physher"][0]
    assert ladder["distribution"] == "uniform"
    assert "alpha" not in ladder


def test_each_estimator_gets_its_own_node(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            marginal=MarginalOptions(algorithms=("ss", "ps", "ps2"), steps=4),
        )
    )
    estimators = [
        node for node in document["physher"] if node["type"] == "marginallikelihood"
    ]
    assert [node["algorithm"] for node in estimators] == ["ss", "ps", "ps2"]
    assert len({node["id"] for node in estimators}) == 3


def test_the_harmonic_mean_needs_no_ladder(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            marginal=MarginalOptions(algorithms=("hm",)),
        )
    )
    chain, estimator = document["physher"]
    assert chain["type"] == "mcmc"
    assert estimator["algorithm"] == "hm"
    assert "steps" not in estimator


def test_the_estimator_reads_the_chain_log(alignment, tree_file):
    """The file it names and the column it looks for are the chain's own."""
    spec = spec_for(
        alignment,
        tree_file,
        method="mcmc",
        marginal=MarginalOptions(algorithms=("hm",)),
    )
    document = build(spec)
    chain, estimator = document["physher"]
    log = find(chain["log"], lambda node: node.get("file", "").endswith(".log"))
    assert estimator["file"] == log["file"]
    assert estimator["treelikelihood"] == "treelikelihood"


def test_marginal_burnin_counts_logged_rows(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="mcmc",
            mcmc=MCMCOptions(length=10000, every=100),
            marginal=MarginalOptions(algorithms=("hm",)),
        )
    )
    assert document["physher"][1]["burnin"] == 10


def test_importance_sampling_reweights_the_approximation(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            method="advi",
            marginal=MarginalOptions(algorithms=("is",), samples=500),
        )
    )
    sampler = document["physher"][-1]
    assert sampler["type"] == "vbis"
    assert sampler["model"] == "@joint"
    assert sampler["distribution"] == "@variational"
    assert sampler["samples"] == 500
    # The joint carries a log Jacobian term the self-normalized form drops.
    assert sampler["normalize"] is False


# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------


def test_no_bootstrap_by_default(alignment, tree_file):
    document = build(spec_for(alignment, tree_file))
    assert not [node for node in document["physher"] if node["type"] == "bootstrap"]


def test_the_bootstrap_follows_the_loggers(alignment, tree_file):
    """What the run reports is the fit on the observed data, not a replicate."""
    document = build(
        spec_for(alignment, tree_file, bootstrap=BootstrapOptions(replicates=50))
    )
    types = [node["type"] for node in document["physher"]]
    assert types[-1] == "bootstrap"
    assert types[0] == "optimizer"


def test_replicates_are_refitted_the_same_way(alignment, tree_file):
    document = build(
        spec_for(alignment, tree_file, bootstrap=BootstrapOptions(replicates=10))
    )
    optimizer, bootstrap = document["physher"][0], document["physher"][-1]
    replicate = bootstrap["optimizer"]
    assert replicate["algorithm"] == optimizer["algorithm"]
    assert [step["algorithm"] for step in replicate["list"]] == [
        step["algorithm"] for step in optimizer["list"]
    ]
    # Distinct ids: physher registers every optimizer it reads by id.
    assert replicate["id"] != optimizer["id"]
    assert not {step["id"] for step in replicate["list"]} & {
        step["id"] for step in optimizer["list"]
    }


def test_the_bootstrap_restarts_from_the_estimate(alignment, tree_file):
    """Every free parameter is snapshotted; a meta schedule declares none."""
    spec = spec_for(alignment, tree_file, bootstrap=BootstrapOptions(replicates=10))
    context = make_context(spec)
    document = build_document(spec, context)
    bootstrap = document["physher"][-1]
    assert bootstrap["parameters"] == [
        f"&{free.free_id}" for free in context.free_parameters
    ]
    assert bootstrap["reset"] is True


def test_bootstrap_options_reach_the_node(alignment, tree_file):
    document = build(
        spec_for(
            alignment,
            tree_file,
            bootstrap=BootstrapOptions(
                replicates=7,
                reset=False,
                include_observed=True,
                compact=False,
            ),
        )
    )
    bootstrap = document["physher"][-1]
    assert bootstrap["replicates"] == 7
    assert bootstrap["reset"] is False
    assert bootstrap["include_observed"] is True
    assert bootstrap["compact"] is False


def test_one_thread_is_left_to_the_default(alignment, tree_file):
    document = build(
        spec_for(alignment, tree_file, bootstrap=BootstrapOptions(replicates=10))
    )
    assert "threads" not in document["physher"][-1]


def test_bootstrap_threads_reach_the_node(alignment, tree_file):
    document = build(
        spec_for(
            alignment, tree_file, bootstrap=BootstrapOptions(replicates=10, threads=4)
        )
    )
    assert document["physher"][-1]["threads"] == 4


def test_the_bootstrap_writes_files_of_its_own(alignment, tree_file):
    """Its loggers must not overwrite the ones the point estimate wrote."""
    document = build(
        spec_for(
            alignment,
            tree_file,
            stem="run",
            bootstrap=BootstrapOptions(replicates=10),
        )
    )
    files = {
        node["file"] for node in document["physher"][-1]["loggers"] if "file" in node
    }
    assert files == {"run.bootstrap.csv", "run.bootstrap.trees"}
    assert "run.tree" not in files
