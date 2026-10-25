"""Root-to-tip regression, searching every branch for the best root.

The same analysis TempEst does: regress each tip's distance from a root on
its date, and put the root wherever that fits a clock best.  For a root a
fraction of the way up the branch above ``c``, the distance to tip ``i`` is

    y_i(t) = d(c, i) + t * s_i,    s_i = +1 below c, -1 elsewhere,

which is linear in ``t``.  So the residual sum of squares of the regression is
a quadratic in ``t``, the coefficient of determination a ratio of two, and the
best point on each branch has a closed form: no line search, and one pass of
numpy per branch.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, TypedDict

import numpy as np

from ..io.newick import Node, parse_newick

__all__ = ["CRITERIA", "RootToTipResult", "root_to_tip"]

#: What the root search optimizes.
CRITERIA = {
    "r2": "r² (maximise the coefficient of determination)",
    "ssr": "SSR (minimise the sum of squared residuals)",
}


class RootToTipResult(TypedDict):
    leaf_names: List[str]
    dates: List[float]
    root_to_tip: List[float]
    slope: float
    intercept: float
    r2: float
    tmrca: Optional[float]
    n: int
    criterion: str


def root_to_tip(
    newick: str, dates: Dict[str, float], criterion: str = "r2"
) -> RootToTipResult:
    """Find the best root for a clock and return the regression there.

    Tips without a date are left out of the fit; there have to be at least
    three that have one, with at least two distinct dates.
    """
    if criterion not in CRITERIA:
        raise ValueError(f"unknown criterion {criterion!r}")
    tree = parse_newick(newick)
    leaves = [leaf for leaf in tree.leaves() if leaf.name in dates]
    if len(leaves) < 3:
        raise ValueError(f"need at least 3 dated tips, not {len(leaves)}")
    x = np.array([float(dates[leaf.name]) for leaf in leaves])
    if np.ptp(x) == 0:
        raise ValueError("every tip has the same date, so there is no slope to fit")
    if any(node.length is None for node in tree.preorder() if node is not tree):
        raise ValueError("the tree has no branch lengths to measure distances with")

    distances, below = _distances(tree, leaves)
    residual = _residual_maker(x)

    best: Optional[Tuple[float, np.ndarray]] = None
    for node in tree.preorder():
        if node.parent is None or not node.length or node.length <= 0:
            continue
        a = distances[id(node)]
        b = np.where(below[id(node)], 1.0, -1.0)
        t, score = _best_point(a, b, node.length, residual, criterion)
        if best is None or score < best[0]:
            best = (score, a + t * b)
    if best is None:
        raise ValueError("no branch of positive length to place a root on")

    y = best[1]
    slope, intercept = np.polyfit(x, y, 1)
    ssr = float(np.sum((y - (slope * x + intercept)) ** 2))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ssr / sst if sst > 0 else 0.0
    return RootToTipResult(
        leaf_names=[leaf.name for leaf in leaves],
        dates=x.tolist(),
        root_to_tip=y.tolist(),
        slope=float(slope),
        intercept=float(intercept),
        r2=min(1.0, max(0.0, r2)),
        tmrca=float(-intercept / slope) if slope != 0 else None,
        n=len(leaves),
        criterion=criterion,
    )


def _distances(
    tree: Node, leaves: List[Node]
) -> Tuple[Dict[int, np.ndarray], Dict[int, np.ndarray]]:
    """Distance from every node to every dated tip, and which tips it subtends."""
    index = {id(leaf): i for i, leaf in enumerate(leaves)}
    below: Dict[int, np.ndarray] = {}
    for node in reversed(list(tree.preorder())):
        mask = np.zeros(len(leaves), dtype=bool)
        if id(node) in index:
            mask[index[id(node)]] = True
        for child in node.children:
            mask |= below[id(child)]
        below[id(node)] = mask

    # From the root down: moving onto a child brings the tips below it closer
    # by the branch length and takes every other tip further away.
    depth = {id(tree): 0.0}
    for node in tree.preorder():
        if node.parent is not None:
            depth[id(node)] = depth[id(node.parent)] + (node.length or 0.0)
    distances = {id(tree): np.array([depth[id(leaf)] for leaf in leaves])}
    for node in tree.preorder():
        if node.parent is None:
            continue
        length = node.length or 0.0
        sign = np.where(below[id(node)], -1.0, 1.0)
        distances[id(node)] = distances[id(node.parent)] + sign * length
    return distances, below


def _residual_maker(x: np.ndarray):
    """The map from a vector to its residuals once regressed on ``x``."""
    design = np.column_stack([np.ones_like(x), x])
    projection = design @ np.linalg.pinv(design)
    return lambda v: v - projection @ v


def _best_point(a, b, length, residual, criterion) -> Tuple[float, float]:
    """Where on a branch of ``length`` the criterion is best, and its value.

    The residuals at ``t`` are ``Ma + t Mb``; centring gives the deviations the
    total sum of squares is taken over, ``Ca + t Cb``.
    """
    ma, mb = residual(a), residual(b)
    ssr = (mb @ mb, 2 * ma @ mb, ma @ ma)
    if criterion == "ssr":
        candidates = [0.0, length]
        if ssr[0] > 0:
            candidates.append(-ssr[1] / (2 * ssr[0]))
        objective = _quadratic(*ssr)
    else:
        ca, cb = a - a.mean(), b - b.mean()
        sst = (cb @ cb, 2 * ca @ cb, ca @ ca)
        candidates = [0.0, length] + _ratio_stationary_points(ssr, sst)
        numerator, denominator = _quadratic(*ssr), _quadratic(*sst)

        def objective(t):
            total = denominator(t)
            # 1 - r2: the ratio itself, so that smaller is better here too
            return numerator(t) / total if total > 0 else np.inf

    points = [t for t in candidates if 0.0 <= t <= length]
    values = [objective(t) for t in points]
    best = int(np.argmin(values))
    return points[best], float(values[best])


def _quadratic(p, q, r):
    return lambda t: p * t * t + q * t + r


def _ratio_stationary_points(numerator, denominator) -> List[float]:
    """Roots of the derivative of (A t^2 + B t + C) / (D t^2 + E t + F)."""
    a, b, c = numerator
    d, e, f = denominator
    roots = np.roots([a * e - b * d, 2 * (a * f - c * d), b * f - c * e])
    return [float(root.real) for root in roots if abs(root.imag) < 1e-12]
