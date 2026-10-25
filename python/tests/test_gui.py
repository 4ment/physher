"""The Streamlit front end, driven headless through ``AppTest``.

Skipped without the ``gui`` extra.  Files are given as paths rather than
uploads, which ``AppTest`` cannot simulate; the upload path is the same code
once the bytes are on disk.
"""

import os
import shlex

import pytest

pytest.importorskip("streamlit")

from streamlit.testing.v1 import AppTest  # noqa: E402

from phyforge.cli.args import parse  # noqa: E402

APP = os.path.join(os.path.dirname(__file__), "..", "src", "phyforge", "gui", "app.py")


@pytest.fixture
def app(alignment, tree_file):
    at = AppTest.from_file(APP, default_timeout=60).run()
    at.text_input(key="local:alignment").input(alignment).run()
    at.text_input(key="local:tree").input(tree_file).run()
    return at


def command(at):
    """The phyforge command the Output tab shows, as argv."""
    lines = [block.value for block in at.code if block.value.startswith("phyforge ")]
    assert lines, [error.value for error in at.error]
    argv = shlex.split(lines[0])[1:]
    return argv[: argv.index("-o")]


def no_exceptions(at):
    assert [exception.value for exception in at.exception] == []


def test_the_app_starts_empty():
    at = AppTest.from_file(APP, default_timeout=60).run()
    no_exceptions(at)
    assert len(at.tabs) == 8


def test_loaded_data_gives_the_default_command(app, alignment, tree_file):
    no_exceptions(app)
    assert command(app) == ["ml", "-i", alignment, "-t", tree_file]


def test_the_command_builds_what_the_gui_built(app, alignment, tree_file):
    app.radio(key="w:method").set_value("mcmc").run()
    app.selectbox(key="w:substitution.kind").set_value("HKY").run()
    app.selectbox(key="w:@rate_model").set_value("gamma").run()
    app.text_input(key="p:kappa:init").input("2.5").run()
    app.checkbox(key="p:shape:estimate").uncheck().run()
    app.selectbox(key="p:kappa:prior.dist").set_value("lognormal").run()
    app.text_input(key="p:kappa:weight").input("10").run()
    no_exceptions(app)

    spec = parse(command(app))
    assert spec.method == "mcmc"
    assert spec.site.categories == 4
    assert spec.parameters["kappa"].init == 2.5
    assert spec.parameters["kappa"].weight == 10.0
    assert str(spec.parameters["kappa"].prior) == "lognormal(mu=0,sigma=1)"
    assert spec.parameters["shape"].estimate is False


def test_a_bad_value_is_shown_and_blocks_the_build(app):
    app.radio(key="w:method").set_value("map").run()
    app.selectbox(key="w:substitution.kind").set_value("HKY").run()
    app.text_input(key="p:kappa:init").input("-1").run()
    errors = [error.value for error in app.error]
    assert any("kappa must start inside" in error for error in errors)
    assert not any(block.value.startswith("phyforge ") for block in app.code)


def test_switching_models_keeps_an_override(app):
    app.selectbox(key="w:substitution.kind").set_value("HKY").run()
    app.text_input(key="p:kappa:init").input("3").run()
    app.selectbox(key="w:substitution.kind").set_value("GTR").run()
    assert "--init" not in command(app)
    app.selectbox(key="w:substitution.kind").set_value("HKY").run()
    assert "kappa=3" in command(app)


def test_a_time_tree_with_a_skygrid(app):
    app.checkbox(key="w:@phyforge.time_tree").check().run()
    app.radio(key="w:method").set_value("advi").run()
    app.selectbox(key="w:coalescent.kind").set_value("skygrid").run()
    assert any("grid points" in error.value for error in app.error)
    app.number_input(key="w:coalescent.grid").set_value(10).run()
    app.number_input(key="w:coalescent.cutoff").set_value(20.0).run()
    app.text_input(key="p:rate:sigma").input("0.05").run()
    no_exceptions(app)
    spec = parse(command(app))
    assert spec.clock.kind == "strict"
    assert spec.coalescent.kind == "skygrid"
    assert spec.coalescent.grid == 10
    assert spec.parameters["rate"].sigma == 0.05


def test_root_to_tip_runs_on_dated_tips(app):
    app.checkbox(key="w:@phyforge.time_tree").check().run()
    app.button(key="r2t_compute").click().run()
    no_exceptions(app)
    metrics = {metric.label: metric.value for metric in app.metric}
    assert metrics["n tips"] == "69"


SESSION_SCRIPT = """
import streamlit as st
from phyforge.gui import session, state, uploads
from phyforge.spec import ParameterOverride

state.init_state()
if "dumped" not in st.session_state:
    state.values()["method"] = "map"
    state.values()["substitution.kind"] = "HKY"
    state.overrides()["kappa"] = ParameterOverride(init=2.0)
    st.session_state[state.TIME_TREE] = True
    uploads.store("alignment", "aln.fa", b">a\\nACGT\\n>b\\nACGA\\n")
    uploads.use_local("tree", TREE)
    st.session_state["dumped"] = session.dump()
    # forget everything, then read it back
    for key in list(st.session_state):
        if key != "dumped":
            del st.session_state[key]
    state.init_state()
    session.load(st.session_state["dumped"])
"""


def test_a_session_round_trips(tree_file):
    at = AppTest.from_string(
        f"TREE = {tree_file!r}\n" + SESSION_SCRIPT, default_timeout=60
    ).run()
    no_exceptions(at)
    values = at.session_state["phyforge.values"]
    assert values["method"] == "map"
    assert values["substitution.kind"] == "HKY"
    assert at.session_state["phyforge.overrides"]["kappa"].init == 2.0
    assert at.session_state["phyforge.time_tree"] is True
    files = at.session_state["phyforge.files"]
    # the upload's content travelled in the session file; the path did not
    assert files["alignment"]["name"] == "aln.fa"
    with open(files["alignment"]["path"]) as fp:
        assert fp.read().startswith(">a")
    assert files["tree"] == {"name": tree_file, "path": tree_file, "upload": False}
    assert values["data.tree"] == tree_file
