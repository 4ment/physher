"""Root-to-tip regression and the newick parser it reads trees with."""

import random

import pytest

np = pytest.importorskip("numpy")

from phyforge.gui.root_to_tip import root_to_tip  # noqa: E402
from phyforge.io.newick import parse_newick  # noqa: E402


def test_the_parser_keeps_labels_lengths_and_shape():
    tree = parse_newick("((a:1,'b c':2.5)90:0.5,[comment]d:3e-1);")
    assert [leaf.name for leaf in tree.leaves()] == ["a", "b c", "d"]
    assert [leaf.length for leaf in tree.leaves()] == [1.0, 2.5, 0.3]
    clade = tree.children[0]
    assert clade.name == "90" and clade.length == 0.5
    assert clade.children[0].parent is clade


def test_a_perfect_clock_is_recovered():
    # Tips sampled at 2000, 2010 and 2020 on a clock of 0.5 per year from a
    # root in 1990, written unrooted so the root has to be found.
    newick = "(a:5,b:10,c:15);"
    dates = {"a": 2000.0, "b": 2010.0, "c": 2020.0}
    result = root_to_tip(newick, dates)
    assert result["r2"] == pytest.approx(1.0)
    assert result["slope"] == pytest.approx(0.5)
    assert result["tmrca"] == pytest.approx(1990.0)


def _random_tree(n, rng):
    nodes = [f"t{i}:{rng.uniform(0.01, 0.5):.4f}" for i in range(n)]
    while len(nodes) > 1:
        a, b = sorted(rng.sample(range(len(nodes)), 2), reverse=True)
        x, y = nodes.pop(a), nodes.pop(b)
        nodes.append(f"({x},{y}):{rng.uniform(0.01, 0.5):.4f}")
    return nodes[0].rsplit(":", 1)[0] + ";"


def _brute_force(newick, dates, criterion, points=401):
    """Every branch at ``points`` positions, distances recomputed from scratch."""
    tree = parse_newick(newick)
    leaves = [leaf for leaf in tree.leaves() if leaf.name in dates]
    x = np.array([dates[leaf.name] for leaf in leaves])

    def ancestors(node):
        found, depth = {}, 0.0
        while node is not None:
            found[id(node)] = depth
            depth += node.length or 0.0
            node = node.parent
        return found

    def distance(u, v):
        above, depth = ancestors(u), 0.0
        while id(v) not in above:
            depth += v.length or 0.0
            v = v.parent
        return depth + above[id(v)]

    best = np.inf
    for node in tree.preorder():
        if node.parent is None or not node.length:
            continue
        below = {id(leaf) for leaf in node.leaves()}
        base = np.array([distance(node, leaf) for leaf in leaves])
        sign = np.array([1.0 if id(leaf) in below else -1.0 for leaf in leaves])
        for t in np.linspace(0, node.length, points):
            y = base + t * sign
            slope, intercept = np.polyfit(x, y, 1)
            ssr = np.sum((y - slope * x - intercept) ** 2)
            score = ssr if criterion == "ssr" else ssr / np.sum((y - y.mean()) ** 2)
            best = min(best, score)
    return best


@pytest.mark.parametrize("criterion", ["r2", "ssr"])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_the_closed_form_is_at_least_as_good_as_a_grid(criterion, seed):
    rng = random.Random(seed)
    n = rng.randint(5, 10)
    newick = _random_tree(n, rng)
    dates = {f"t{i}": rng.uniform(2000, 2020) for i in range(n)}
    result = root_to_tip(newick, dates, criterion)
    y, x = np.array(result["root_to_tip"]), np.array(result["dates"])
    ssr = np.sum((y - result["slope"] * x - result["intercept"]) ** 2)
    score = ssr if criterion == "ssr" else 1 - result["r2"]
    assert score <= _brute_force(newick, dates, criterion) + 1e-9


def test_too_few_dated_tips_is_refused():
    with pytest.raises(ValueError, match="3 dated tips"):
        root_to_tip("(a:1,b:1,c:1);", {"a": 2000.0, "b": 2001.0})
