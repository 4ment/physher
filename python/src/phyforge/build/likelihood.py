"""The tree likelihood, and the alignment it is computed from."""

from __future__ import annotations

from ..io.fasta import read_fasta
from ..spec import SUBSTITUTION_MODELS
from . import branchmodel, sitemodel, substitution, tree
from .context import BuildContext


def build(context: BuildContext) -> dict:
    """Build the ``treelikelihood`` node and everything under it."""
    spec = context.spec
    ids = context.ids
    node = {
        "id": ids.tree_likelihood,
        "type": "treelikelihood",
        "sitepattern": build_site_pattern(context),
        "sitemodel": sitemodel.build(context),
        "substitutionmodel": substitution.build(context),
        "tree": tree.build(context),
    }
    if spec.data.tip_states:
        node["tipstates"] = True
    if not spec.data.sse:
        node["sse"] = False
    if spec.clock.is_time_tree:
        node["branchmodel"] = branchmodel.build(context)
        # The height reparameterization's Jacobian is added once, as its own
        # term, so it is not double counted here.
        node["include_jacobian"] = False
    return node


def build_site_pattern(context: BuildContext) -> dict:
    """Build the ``sitepattern`` node wrapping the alignment."""
    spec = context.spec
    ids = context.ids
    datatype = SUBSTITUTION_MODELS[spec.substitution.kind.upper()]

    alignment = {"id": ids.alignment, "type": "alignment"}
    if spec.data.embed_alignment:
        alignment["sequences"] = {
            sequence.taxon: sequence.sequence
            for sequence in read_fasta(spec.data.alignment)
        }
    else:
        alignment["file"] = spec.data.alignment

    return {
        "id": ids.site_pattern,
        "type": "sitepattern",
        "datatype": datatype,
        "alignment": alignment,
    }
