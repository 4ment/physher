"""The ids phyforge gives the objects it creates.

Collected in one place because every layer needs to agree on them: the
builders declare them, the loggers and operators reference them, and a GUI
displays them.  Instantiate ``Ids`` to rename everything at once.
"""

from __future__ import annotations

from dataclasses import dataclass

UNCONSTRAINED_SUFFIX = ".unconstrained"


@dataclass(frozen=True)
class Ids:
    taxa: str = "taxa"
    alignment: str = "seqs"
    site_pattern: str = "patterns"
    site_model: str = "sitemodel"
    substitution_model: str = "substmodel"
    tree: str = "tree"
    branch_model: str = "branchmodel"
    tree_likelihood: str = "treelikelihood"
    prior: str = "prior"
    joint: str = "joint"
    jacobian: str = "jacobian"
    variational: str = "variational"
    elbo: str = "elbo"
    grad_elbo: str = "gradelbo"

    # parameters
    branch_lengths: str = "tree.distances"
    root_height: str = "root.height"
    proportions: str = "proportions"
    rate: str = "rate"
    branch_rates: str = "branch_rates"
    clock_location: str = "clock.location"
    clock_scale: str = "clock.scale"
    frequencies: str = "frequencies"
    rates: str = "rates"
    kappa: str = "kappa"
    shape: str = "shape"
    scale: str = "scale"
    proportion_invariant: str = "proportion_invariant"
    #: free category rates, and the weights that go with them.  Distinct from
    #: `proportions`, which is the tree's node-height reparameterization.
    site_rates: str = "site.rates"
    site_proportions: str = "site.proportions"
    top_rate: str = "site.top_rate"
    #: the empirical CAT rates, which the cat estimator writes rather than an
    #: optimizer, and so are never registered as a free parameter
    cat_rates: str = "catrates"
    theta: str = "theta"
    growth: str = "growth"
    precision: str = "precision"

    def unconstrained(self, id_: str) -> str:
        return f"{id_}{UNCONSTRAINED_SUFFIX}"
