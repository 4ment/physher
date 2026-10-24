"""The CLI, which should do nothing but build a spec and print JSON."""

import json

from phyforge.cli import main, parse


def test_arguments_become_a_spec(alignment, tree_file):
    spec = parse(
        [
            "mcmc",
            "-i",
            alignment,
            "-t",
            tree_file,
            "-m",
            "GTR",
            "-C",
            "4",
            "-c",
            "strict",
        ]
    )
    assert spec.method == "mcmc"
    assert spec.substitution.kind == "GTR"
    assert spec.site.categories == 4
    assert spec.clock.kind == "strict"


def test_parsimony_flags_reach_the_spec(alignment, tree_file):
    spec = parse(
        [
            "ml",
            "-i",
            alignment,
            "-t",
            tree_file,
            "--init-parsimony",
            "--init-parsimony-min",
            "1e-8",
        ]
    )
    assert spec.data.init_parsimony is True
    assert spec.data.init_parsimony_min == 1e-8


def test_parsimony_defaults_to_unset(alignment, tree_file):
    spec = parse(["ml", "-i", alignment, "-t", tree_file])
    assert spec.data.init_parsimony is None
    assert spec.data.init_parsimony_min is None
    # unset means on wherever it applies, which an unrooted tree does
    assert spec.uses_parsimony is True


def test_parsimony_can_be_turned_off(alignment, tree_file):
    spec = parse(["ml", "-i", alignment, "-t", tree_file, "--no-init-parsimony"])
    assert spec.data.init_parsimony is False
    assert spec.uses_parsimony is False


def test_writes_a_document(alignment, tree_file, tmp_path, capsys):
    output = tmp_path / "config.json"
    code = main(["ml", "-i", alignment, "-t", tree_file, "-o", str(output)])
    assert code == 0
    document = json.loads(output.read_text())
    assert "treelikelihood" in document and "physher" in document


def test_prints_to_stdout(alignment, tree_file, capsys):
    assert main(["ml", "-i", alignment, "-t", tree_file]) == 0
    assert "treelikelihood" in json.loads(capsys.readouterr().out)


def test_reports_invalid_spec_without_writing(tree_file, capsys):
    code = main(["ml", "-i", "missing.fa", "-t", tree_file])
    assert code == 2
    assert "does not exist" in capsys.readouterr().err


def test_site_model_flag_covers_cat_and_free_rates(alignment, tree_file):
    spec = parse(
        [
            "ml",
            "-i",
            alignment,
            "-t",
            tree_file,
            "--site-model",
            "cat",
            "-C",
            "20",
            "--cat-prior",
            "npmle",
        ]
    )
    assert spec.site.is_cat and spec.site.cat.prior == "npmle"

    spec = parse(
        [
            "ml",
            "-i",
            alignment,
            "-t",
            tree_file,
            "-D",
            "discrete",
            "-C",
            "4",
            "--rate-parameterization",
            "mean_contribution",
        ]
    )
    assert spec.site.is_free_rates
    assert spec.site.rate_parameterization == "mean_contribution"
    assert not spec.site.uses_em

    spec = parse(
        [
            "ml",
            "-i",
            alignment,
            "-t",
            tree_file,
            "-D",
            "discrete",
            "-C",
            "4",
            "--cat-em",
        ]
    )
    assert spec.site.uses_em


def test_cpo_reaches_the_document(alignment, tree_file, tmp_path):
    output = tmp_path / "config.json"
    code = main(
        [
            "mcmc",
            "-i",
            alignment,
            "-t",
            tree_file,
            "--cpo",
            "--cpo-burnin",
            "7",
            "-o",
            str(output),
        ]
    )
    assert code == 0
    document = json.loads(output.read_text())
    cpo = [node for node in document["physher"] if node["type"] == "cpo"]
    assert len(cpo) == 1
    assert cpo[0]["burnin"] == 7


def test_lppd_reaches_the_document(alignment, tree_file, tmp_path):
    output = tmp_path / "config.json"
    code = main(
        [
            "mcmc",
            "-i",
            alignment,
            "-t",
            tree_file,
            "--lppd",
            "--lppd-burnin",
            "7",
            "-o",
            str(output),
        ]
    )
    assert code == 0
    document = json.loads(output.read_text())
    lppd = [node for node in document["physher"] if node["type"] == "predictive"]
    assert len(lppd) == 1
    assert lppd[0]["burnin"] == 7
    assert not [node for node in document["physher"] if node["type"] == "cpo"]


def test_steps_alone_still_means_stepping_stone(alignment, tree_file):
    """--steps was the only way to ask for a marginal likelihood; keep it so."""
    spec = parse(["mcmc", "-i", alignment, "-t", tree_file, "--steps", "16"])
    assert spec.marginal.algorithms == ("ss",)
    assert spec.marginal.steps == 16


def test_marginal_takes_several_estimators(alignment, tree_file):
    spec = parse(
        [
            "mcmc",
            "-i",
            alignment,
            "-t",
            tree_file,
            "--marginal",
            "ss",
            "ps",
            "--steps",
            "8",
            "--ladder",
            "uniform",
            "--marginal-burnin",
            "5",
        ]
    )
    assert spec.marginal.algorithms == ("ss", "ps")
    assert spec.marginal.ladder == "uniform"
    assert spec.marginal.burnin == 5


def test_marginal_is_off_without_the_flag(alignment, tree_file):
    spec = parse(["mcmc", "-i", alignment, "-t", tree_file])
    assert spec.marginal.algorithms == ()
    assert not spec.marginal.enabled


def test_importance_sampling_reaches_the_document(alignment, tree_file, tmp_path):
    output = tmp_path / "config.json"
    code = main(
        [
            "advi",
            "-i",
            alignment,
            "-t",
            tree_file,
            "--marginal",
            "is",
            "--is-samples",
            "250",
            "-o",
            str(output),
        ]
    )
    assert code == 0
    document = json.loads(output.read_text())
    sampler = [node for node in document["physher"] if node["type"] == "vbis"]
    assert len(sampler) == 1
    assert sampler[0]["samples"] == 250


def test_bootstrap_flags_become_a_spec(alignment, tree_file):
    spec = parse(
        [
            "ml",
            "-i",
            alignment,
            "-t",
            tree_file,
            "-b",
            "200",
            "--bootstrap-continue",
            "--bootstrap-observed",
        ]
    )
    assert spec.bootstrap.replicates == 200
    assert spec.bootstrap.reset is False
    assert spec.bootstrap.include_observed is True
    assert spec.bootstrap.compact is True
    assert spec.bootstrap.threads == 1


def test_bootstrap_threads_flag(alignment, tree_file):
    spec = parse(
        ["ml", "-i", alignment, "-t", tree_file, "-b", "20", "--bootstrap-threads", "6"]
    )
    assert spec.bootstrap.threads == 6


def test_bootstrap_reaches_the_document(alignment, tree_file, tmp_path):
    output = tmp_path / "config.json"
    code = main(["ml", "-i", alignment, "-t", tree_file, "-b", "20", "-o", str(output)])
    assert code == 0
    document = json.loads(output.read_text())
    bootstrap = [node for node in document["physher"] if node["type"] == "bootstrap"]
    assert len(bootstrap) == 1
    assert bootstrap[0]["replicates"] == 20


def test_an_unusable_marginal_request_is_reported(alignment, tree_file, capsys):
    code = main(["ml", "-i", alignment, "-t", tree_file, "--marginal", "hm"])
    assert code != 0
    assert "marginal" in capsys.readouterr().err
