import os

import pytest

EXAMPLES = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "examples", "fluA")
)


@pytest.fixture
def alignment():
    return os.path.join(EXAMPLES, "fluA.fa")


@pytest.fixture
def tree_file():
    return os.path.join(EXAMPLES, "fluA-rooted.nxs")
