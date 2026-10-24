"""The branch model -- how a molecular clock maps time to substitutions."""

from __future__ import annotations

import math

from ..nodes import model_ref, parameter
from .context import POSITIVE, BuildContext
from .parameters import free_parameter

ROLE_RATE = "rate"
ROLE_BRANCH_RATES = "branch_rates"
ROLE_CLOCK_SCALE = "clock_scale"

#: The fixed reference lognormal the ``ncln`` clock draws its raw per-branch
#: rates from (``_arbitrary_calculate_rates`` in ``branchmodel.c``).  Writing
#: it down here lets the prior match the transform exactly.
NCLN_REFERENCE_MU = -0.5 * math.log(2.0)
NCLN_REFERENCE_SIGMA = math.sqrt(math.log(2.0))


def build(context: BuildContext) -> dict:
    """Build the ``branchmodel`` node.  Only called for a time tree."""
    spec = context.spec.clock
    ids = context.ids
    node = {
        "id": ids.branch_model,
        "type": "branchmodel",
        "kind": spec.kind,
        "tree": model_ref(ids.tree),
    }

    if spec.kind == "ncln":
        return _uncorrelated(context, node)

    if spec.rate is None:
        node["rate"] = free_parameter(
            context, ids.rate, spec.rate_init, role=ROLE_RATE, kind=POSITIVE
        )
    else:
        # Pinned: not registered, so no operator or variational factor touches it.
        node["rate"] = parameter(ids.rate, spec.rate, lower=0)
    return node


def _uncorrelated(context: BuildContext, node: dict) -> dict:
    """The non-centred uncorrelated lognormal clock.

    ``rate`` holds one raw draw per branch on a *fixed* reference lognormal;
    ``location`` scales them to the mean rate and ``scale`` sets how far they
    spread.  Keeping the per-branch variation in the fixed reference scale is
    what makes the parameterization non-centred, so ``location`` and ``scale``
    stay identifiable when the data say little about individual branches.
    """
    spec = context.spec.clock
    ids = context.ids

    node["rate"] = free_parameter(
        context,
        ids.branch_rates,
        1.0,
        role=ROLE_BRANCH_RATES,
        kind=POSITIVE,
        dimension=context.rooted_branch_count,
        log=False,
    )
    if spec.rate is None:
        node["location"] = free_parameter(
            context, ids.clock_location, spec.rate_init, role=ROLE_RATE, kind=POSITIVE
        )
    else:
        node["location"] = parameter(ids.clock_location, spec.rate, lower=0)
    node["scale"] = free_parameter(
        context,
        ids.clock_scale,
        spec.ncln_scale_init,
        role=ROLE_CLOCK_SCALE,
        kind=POSITIVE,
    )
    return node
