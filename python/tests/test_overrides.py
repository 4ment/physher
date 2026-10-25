"""Per-parameter overrides: starting values, fixing, priors, advi and mcmc."""

import math

import pytest

from phyforge import AnalysisSpec, SpecError, build, check, describe
from phyforge.build.context import simplex_to_unconstrained
from phyforge.spec import (
    ClockSpec,
    CoalescentSpec,
    DataSpec,
    ParameterOverride,
    PriorChoice,
    SiteModelSpec,
    SubstitutionSpec,
)


def spec_for(alignment, tree_file, overrides=None, **changes):
    changes.setdefault("substitution", SubstitutionSpec(kind="HKY"))
    return AnalysisSpec(
        data=DataSpec(alignment=alignment, tree=tree_file),
        parameters=overrides or {},
        **changes,
    )


def node_with_id(node, id_):
    if isinstance(node, dict):
        if node.get("id") == id_:
            return node
        children = node.values()
    elif isinstance(node, list):
        children = node
    else:
        return None
    for child in children:
        found = node_with_id(child, id_)
        if found is not None:
            return found
    return None


def priors_of(document):
    return {d["x"]: d for d in document["prior"]["distributions"] if "x" in d}


def problem_fields(spec):
    _, problems = describe(spec)
    return [problem.field for problem in problems]


# --------------------------------------------------------------------------
# the text form of a prior
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text",
    ["gamma(shape=2,rate=10)", "dirichlet(alpha=[1,1,2,2])", "oneonx", "none"],
)
def test_a_prior_round_trips_through_its_text(text):
    assert str(PriorChoice.parse(text)) == text


def test_a_prior_parses_spaces_and_case():
    assert PriorChoice.parse(" Gamma( shape = 2 , rate=0.5 )") == PriorChoice(
        "gamma", {"shape": 2.0, "rate": 0.5}
    )


@pytest.mark.parametrize("text", ["gamma(shape)", "gamma(shape=x)", "(a=1)"])
def test_a_malformed_prior_is_refused(text):
    with pytest.raises(ValueError):
        PriorChoice.parse(text)


def test_overrides_round_trip_through_a_dict(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {
            "kappa": ParameterOverride(
                init=2.0, prior=PriorChoice("lognormal", {"mu": 1, "sigma": 1.25})
            ),
            "shape": ParameterOverride(estimate=False, init=0.5),
        },
        method="map",
    )
    assert AnalysisSpec.from_dict(spec.to_dict()) == spec


def test_an_override_reads_a_prior_written_as_text():
    override = ParameterOverride.from_dict({"prior": "exponential(mean=2)"})
    assert override.prior == PriorChoice("exponential", {"mean": 2.0})


# --------------------------------------------------------------------------
# describe
# --------------------------------------------------------------------------


def test_describe_lists_every_parameter_with_its_defaults(alignment, tree_file):
    spec = spec_for(alignment, tree_file, method="map", site=SiteModelSpec(4))
    rows, problems = describe(spec)
    assert problems == []
    by_id = {row.id: row for row in rows}
    assert set(by_id) == {"tree.distances", "frequencies", "kappa", "shape"}
    assert by_id["shape"].prior == PriorChoice("exponential", {"lambda": 1.0})
    assert by_id["kappa"].prior is None
    assert by_id["frequencies"].default_init == [0.25] * 4
    assert all(row.prior_editable for row in rows)


def test_describe_marks_structural_priors(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        method="mcmc",
        clock=ClockSpec(kind="strict"),
        coalescent=CoalescentSpec(kind="skyride"),
    )
    rows, _ = describe(spec)
    by_id = {row.id: row for row in rows}
    assert not by_id["theta"].prior_editable
    assert not by_id["proportions"].prior_editable
    assert by_id["precision"].prior_editable


def test_describe_still_lists_a_fixed_parameter(alignment, tree_file):
    spec = spec_for(
        alignment, tree_file, {"kappa": ParameterOverride(estimate=False, init=3.0)}
    )
    rows, _ = describe(spec)
    kappa = next(row for row in rows if row.id == "kappa")
    assert not kappa.estimated
    assert kappa.init == 3.0


def test_describe_lists_nothing_for_an_invalid_spec(tree_file):
    rows, problems = describe(AnalysisSpec(data=DataSpec(tree=tree_file)))
    assert rows == []
    assert [problem.field for problem in problems] == ["data.alignment"]


def test_describe_can_ignore_overrides_of_absent_parameters(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"shape": ParameterOverride(init=0.5)},
        substitution=SubstitutionSpec(kind="HKY"),
    )
    assert problem_fields(spec) == ["parameters.shape"]
    assert describe(spec, strict=False)[1] == []


# --------------------------------------------------------------------------
# starting values
# --------------------------------------------------------------------------


def test_an_init_override_wins_over_the_shorthand(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"kappa": ParameterOverride(init=4.0)},
        substitution=SubstitutionSpec(kind="HKY", kappa_init=2.0),
    )
    assert node_with_id(build(spec), "kappa")["x"] == 4.0


def test_an_init_override_moves_the_unconstrained_value(alignment, tree_file):
    spec = spec_for(
        alignment, tree_file, {"kappa": ParameterOverride(init=4.0)}, method="advi"
    )
    document = build(spec)
    assert node_with_id(document, "kappa.unconstrained")["x"] == pytest.approx(
        math.log(4.0)
    )
    # and the variational location follows it
    assert node_with_id(document, "mu.kappa")["x"] == pytest.approx(math.log(4.0))


def test_a_vector_init_gives_every_element(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"site.rates": ParameterOverride(init=[1.0, 2.0, 3.0])},
        site=SiteModelSpec(categories=3, distribution="discrete"),
    )
    assert node_with_id(build(spec), "site.rates")["x"] == [1.0, 2.0, 3.0]


def test_a_simplex_init_is_written_as_proportions(alignment, tree_file):
    frequencies = [0.1, 0.2, 0.3, 0.4]
    spec = spec_for(
        alignment,
        tree_file,
        {"frequencies": ParameterOverride(init=frequencies)},
        method="mcmc",
    )
    assert node_with_id(build(spec), "frequencies")["x"] == pytest.approx(frequencies)


def test_a_transformed_simplex_init_goes_through_stick_breaking(alignment, tree_file):
    frequencies = [0.1, 0.2, 0.3, 0.4]
    spec = spec_for(
        alignment, tree_file, {"frequencies": ParameterOverride(init=frequencies)}
    )
    inner = node_with_id(build(spec), "frequencies.unconstrained")
    assert inner["x"] == pytest.approx(simplex_to_unconstrained(frequencies))


def test_stick_breaking_inverts_physher_transform():
    # _inverse_transform_simplex_stan, written out
    values = [0.1, 0.2, 0.3, 0.4]
    y = simplex_to_unconstrained(values)
    stick, x = 1.0, []
    for k, yk in enumerate(y):
        x.append(stick / (1 + math.exp(-(yk - math.log(len(y) - k)))))
        stick -= x[-1]
    x.append(stick)
    assert x == pytest.approx(values)
    assert simplex_to_unconstrained([0.25] * 4) == pytest.approx([0.0] * 3)


@pytest.mark.parametrize(
    "init",
    [0.0, -1.0, [1.0, 2.0]],
    ids=["boundary", "negative", "wrong-length"],
)
def test_a_bad_init_is_a_problem(alignment, tree_file, init):
    spec = spec_for(alignment, tree_file, {"kappa": ParameterOverride(init=init)})
    assert problem_fields(spec) == ["parameters.kappa.init"]
    with pytest.raises(SpecError):
        build(spec)


@pytest.mark.parametrize("init", [[0.5, 0.5, 0.5, 0.5], [0.5, 0.5]])
def test_a_simplex_init_must_be_proportions(alignment, tree_file, init):
    spec = spec_for(alignment, tree_file, {"frequencies": ParameterOverride(init=init)})
    assert problem_fields(spec) == ["parameters.frequencies.init"]


def test_branch_length_init_needs_parsimony_off(alignment, tree_file):
    overrides = {"tree.distances": ParameterOverride(init=0.05)}
    spec = spec_for(alignment, tree_file, overrides)
    assert problem_fields(spec) == ["parameters.tree.distances.init"]
    spec.data.init_parsimony = False
    assert problem_fields(spec) == []


def test_every_problem_is_reported_at_once(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {
            "kappa": ParameterOverride(init=-1.0),
            "frequencies": ParameterOverride(init=[1.0]),
            "nothing": ParameterOverride(init=1.0),
        },
    )
    assert problem_fields(spec) == [
        "parameters.frequencies.init",
        "parameters.kappa.init",
        "parameters.nothing",
    ]


# --------------------------------------------------------------------------
# fixing
# --------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["ml", "map", "mcmc", "advi"])
def test_a_fixed_parameter_is_left_alone(alignment, tree_file, method):
    spec = spec_for(
        alignment,
        tree_file,
        {"shape": ParameterOverride(estimate=False, init=0.5)},
        method=method,
        site=SiteModelSpec(4),
    )
    document = build(spec)
    references = [
        node
        for node in _strings(document["physher"])
        + _strings(document.get("variational", {}))
        + _strings(document.get("prior", {}))
        if node.strip("&@").startswith("shape")
    ]
    assert references == []
    shape = node_with_id(document, "shape")
    assert shape is not None
    assert shape.get("x") == 0.5 or node_with_id(shape, "shape.unconstrained")[
        "x"
    ] == pytest.approx(math.log(0.5))


def test_a_fixed_simplex_keeps_its_proportions(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"frequencies": ParameterOverride(estimate=False, init=[0.1, 0.2, 0.3, 0.4])},
        method="mcmc",
    )
    document = build(spec)
    assert node_with_id(document, "frequencies")["x"] == pytest.approx(
        [0.1, 0.2, 0.3, 0.4]
    )
    assert "&frequencies" not in priors_of(document)


def test_fixing_theta_drops_the_gmrf(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"theta": ParameterOverride(estimate=False)},
        method="map",
        clock=ClockSpec(kind="strict"),
        coalescent=CoalescentSpec(kind="skyride"),
    )
    document = build(spec)
    distributions = [d.get("distribution") for d in document["prior"]["distributions"]]
    assert "gmrf" not in distributions


def test_a_fixed_parameter_takes_nothing_else(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"kappa": ParameterOverride(estimate=False, prior=PriorChoice("oneonx"))},
        method="map",
    )
    assert [problem.field for problem in check(spec)] == ["parameters.kappa.prior"]


def _strings(node):
    if isinstance(node, str):
        return [node]
    if isinstance(node, dict):
        return [s for value in node.values() for s in _strings(value)]
    if isinstance(node, list):
        return [s for value in node for s in _strings(value)]
    return []


# --------------------------------------------------------------------------
# priors
# --------------------------------------------------------------------------


def test_a_prior_override_replaces_the_default(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"shape": ParameterOverride(prior=PriorChoice.parse("gamma(shape=2,rate=4)"))},
        method="map",
        site=SiteModelSpec(4),
    )
    prior = priors_of(build(spec))["&shape"]
    assert prior["id"] == "prior.shape"
    assert prior["distribution"] == "gamma"
    assert prior["parameters"] == {"shape": 2.0, "rate": 4.0}


def test_a_parameter_without_a_default_can_be_given_a_prior(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {
            "kappa": ParameterOverride(
                prior=PriorChoice.parse("lognormal(mu=1,sigma=1.25)")
            )
        },
        method="map",
    )
    prior = priors_of(build(spec))["&kappa"]
    assert prior["id"] == "prior.kappa"
    assert prior["parameters"] == {"mu": 1.0, "sigma": 1.25}


def test_none_removes_a_prior(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"tree.distances": ParameterOverride(prior=PriorChoice("none"))},
        method="map",
    )
    assert "&tree.distances" not in priors_of(build(spec))


def test_ctmcscale_is_given_the_tree(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"kappa": ParameterOverride(prior=PriorChoice("ctmcscale"))},
        method="map",
    )
    assert priors_of(build(spec))["&kappa"]["tree"] == "@tree"


def test_a_vector_hyperparameter_must_fit_the_parameter(alignment, tree_file):
    prior = PriorChoice("dirichlet", {"alpha": [1.0, 2.0]})
    spec = spec_for(
        alignment,
        tree_file,
        {"frequencies": ParameterOverride(prior=prior)},
        method="map",
    )
    assert problem_fields(spec) == ["parameters.frequencies.prior"]
    prior.parameters["alpha"] = [1.0, 2.0, 3.0, 4.0]
    assert problem_fields(spec) == []


@pytest.mark.parametrize(
    "prior, field",
    [
        ("beta(alpha=1,beta=1)", "parameters.kappa.prior"),  # wrong kind
        ("gamma(shape=1)", "parameters.kappa.prior"),  # missing rate
        ("gamma(shape=1,rate=-1)", "parameters.kappa.prior"),  # not positive
        ("student(df=3)", "parameters.kappa.prior"),  # unknown
    ],
)
def test_an_unsuitable_prior_is_a_problem(alignment, tree_file, prior, field):
    spec = spec_for(
        alignment,
        tree_file,
        {"kappa": ParameterOverride(prior=PriorChoice.parse(prior))},
        method="map",
    )
    assert problem_fields(spec) == [field]


def test_a_structural_prior_cannot_be_replaced(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {
            "proportions": ParameterOverride(
                prior=PriorChoice("beta", {"alpha": 1, "beta": 1})
            )
        },
        method="map",
        clock=ClockSpec(kind="strict"),
        coalescent=CoalescentSpec(kind="constant"),
    )
    assert problem_fields(spec) == ["parameters.proportions.prior"]


def test_ml_takes_no_prior(alignment, tree_file):
    spec = spec_for(
        alignment, tree_file, {"kappa": ParameterOverride(prior=PriorChoice("oneonx"))}
    )
    assert [problem.field for problem in check(spec)] == ["parameters.kappa.prior"]


# --------------------------------------------------------------------------
# advi and mcmc
# --------------------------------------------------------------------------


def test_variational_location_and_scale_can_be_set(alignment, tree_file):
    spec = spec_for(
        alignment,
        tree_file,
        {"kappa": ParameterOverride(mu=1.5, sigma=0.05)},
        method="advi",
    )
    document = build(spec)
    assert node_with_id(document, "mu.kappa")["x"] == 1.5
    assert node_with_id(document, "sigma.kappa")["x"] == 0.05


def test_a_simplex_scale_counts_unconstrained_values(alignment, tree_file):
    overrides = {"frequencies": ParameterOverride(sigma=[0.1, 0.2, 0.3, 0.4])}
    spec = spec_for(alignment, tree_file, overrides, method="advi")
    assert problem_fields(spec) == ["parameters.frequencies.sigma"]
    overrides["frequencies"].sigma = [0.1, 0.2, 0.3]
    assert node_with_id(build(spec), "sigma.frequencies")["x"] == [0.1, 0.2, 0.3]


def test_an_operator_weight_can_be_set(alignment, tree_file):
    spec = spec_for(
        alignment, tree_file, {"kappa": ParameterOverride(weight=10)}, method="mcmc"
    )
    operators = build(spec)["physher"][0]["operators"]
    weights = {op["id"]: op["weight"] for op in operators if op["x"] == "&kappa"}
    assert weights == {"op.scaler.kappa": 10, "op.slider.kappa": 10}


@pytest.mark.parametrize(
    "override, method, field",
    [
        (ParameterOverride(mu=1.0), "map", "parameters.kappa.mu"),
        (ParameterOverride(sigma=0.0), "advi", "parameters.kappa.sigma"),
        (ParameterOverride(weight=2), "advi", "parameters.kappa.weight"),
        (ParameterOverride(weight=0), "mcmc", "parameters.kappa.weight"),
    ],
)
def test_method_options_need_their_method(
    alignment, tree_file, override, method, field
):
    spec = spec_for(alignment, tree_file, {"kappa": override}, method=method)
    assert [problem.field for problem in check(spec)] == [field]
