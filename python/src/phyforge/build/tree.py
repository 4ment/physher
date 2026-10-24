"""The tree node, and the taxa block that carries tip dates."""

from __future__ import annotations

from ..io.trees import read_newick
from ..nodes import model_ref
from .context import POSITIVE, UNIT, BuildContext
from .parameters import free_parameter

ROLE_BRANCH_LENGTHS = "branch_lengths"
ROLE_ROOT_HEIGHT = "root_height"
ROLE_PROPORTIONS = "proportions"


def build(context: BuildContext) -> dict:
    """Build the ``tree`` node for the spec.

    A time tree is reparameterized as a root height plus one proportion per
    remaining internal node, always through unconstrained parameters: that
    reparameterization is what makes the height ordering constraint implicit
    rather than something a sampler has to respect.
    """
    spec = context.spec
    ids = context.ids
    node = {"id": ids.tree, "type": "tree"}

    if spec.clock.is_time_tree:
        node["time"] = True
        node["root_height"] = free_parameter(
            context,
            ids.root_height,
            spec.clock.root_height_init,
            role=ROLE_ROOT_HEIGHT,
            kind=POSITIVE,
            transformed=True,
        )
        node["proportions"] = free_parameter(
            context,
            ids.proportions,
            0.5,
            role=ROLE_PROPORTIONS,
            kind=UNIT,
            dimension=context.internal_count,
            transformed=True,
            log=False,
        )
        node["taxa"] = build_taxa(context)
    else:
        node["branch_lengths"] = free_parameter(
            context,
            ids.branch_lengths,
            0.1,
            role=ROLE_BRANCH_LENGTHS,
            kind=POSITIVE,
            dimension=context.branch_count,
            log=False,
        )

    if spec.data.embed_tree:
        node["newick"] = read_newick(spec.data.tree)
    else:
        node["file"] = spec.data.tree

    if spec.data.keep_branch_lengths:
        node["keep_branch_lengths"] = True
    if spec.uses_parsimony:
        # physher overwrites the branch lengths of the topology loaded above
        # with their parsimony estimates.  The site pattern is referenced
        # rather than repeated: the tree node is built after the likelihood's
        # own sitepattern, so the id is already registered by the time
        # physher reads this.
        initializer = {
            "algorithm": "parsimony",
            "sitepattern": model_ref(ids.site_pattern),
        }
        if spec.data.init_parsimony_min is not None:
            initializer["min"] = spec.data.init_parsimony_min
        node["initializer"] = initializer
    return node


def build_taxa(context: BuildContext) -> dict:
    """The taxa block, with a ``date`` attribute per taxon for a time tree."""
    ids = context.ids
    taxa = []
    for name in context.taxa:
        taxon = {"id": name, "type": "Taxon"}
        if context.dates is not None:
            taxon["attributes"] = {"date": context.dates[name]}
        taxa.append(taxon)
    return {"id": ids.taxa, "type": "Taxa", "taxa": taxa}
