"""Assembling a whole physher document from a spec.

physher reads the top-level object in order, building every entry except
``physher`` (the run list) and ``init`` as a model, so a node may only
reference ids declared before it.  The ordering below reflects that.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from ..io.dates import dates_from_csv, dates_from_names
from ..io.fasta import read_fasta
from ..io.trees import taxa_from_tree
from ..nodes import compound, model_ref
from ..spec import AnalysisSpec
from ..validate import validate
from . import actions, likelihood, priors, variational
from .context import BuildContext

__all__ = ["build", "make_context"]


def build(spec: AnalysisSpec, context: Optional[BuildContext] = None) -> Dict:
    """Build the physher JSON document for ``spec``.

    Raises :class:`phyforge.validate.SpecError` if the spec is inconsistent,
    so a caller never gets a config that cannot run.
    """
    validate(spec)
    if context is None:
        context = make_context(spec)

    document: Dict[str, object] = {}
    ids = context.ids

    # The likelihood first: it declares the tree, the substitution and site
    # models, and with them most of the free parameters everything else maps
    # over.
    tree_likelihood = likelihood.build(context)
    document[ids.tree_likelihood] = tree_likelihood

    variational_parameters: List[str] = []

    if spec.uses_priors:
        document[ids.prior] = priors.build(context)
        members: List[object] = [
            model_ref(ids.tree_likelihood),
            model_ref(ids.prior),
        ]
        if _needs_jacobian(context):
            document[ids.jacobian] = variational.build_jacobian(context)
            members.append(model_ref(ids.jacobian))
        if spec.clock.is_time_tree:
            # The tree joins the joint so its stored state takes part in the
            # propose/accept cycle.
            members.append(model_ref(ids.tree))
        document[ids.joint] = compound(ids.joint, members)

    if spec.method == "advi":
        approximation, variational_parameters = variational.build(context)
        document[ids.variational] = approximation
        for bound in variational.build_bounds(context):
            document[bound["id"]] = bound

    document["physher"] = actions.build(context, variational_parameters)
    return document


def make_context(spec: AnalysisSpec) -> BuildContext:
    """Read the input files a spec refers to and set up the build state."""
    taxa = _read_taxa(spec)
    dates = _read_dates(spec, taxa) if spec.clock.is_time_tree else None
    return BuildContext(spec=spec, taxa=taxa, dates=dates)


def _read_taxa(spec: AnalysisSpec) -> List[str]:
    if spec.data.alignment is not None:
        return [sequence.taxon for sequence in read_fasta(spec.data.alignment)]
    return taxa_from_tree(spec.data.tree)


def _read_dates(spec: AnalysisSpec, taxa: List[str]) -> Dict[str, float]:
    dates_spec = spec.clock.dates
    if dates_spec.source == "zero":
        return {taxon: 0.0 for taxon in taxa}
    if dates_spec.source == "csv":
        dates = dates_from_csv(dates_spec.path, dates_spec.date_format)
        missing = [taxon for taxon in taxa if taxon not in dates]
        if missing:
            raise KeyError(
                f"{len(missing)} taxa are missing from {dates_spec.path}; "
                f"first is {missing[0]!r}"
            )
        return {taxon: dates[taxon] for taxon in taxa}
    return dates_from_names(taxa, dates_spec.regex, dates_spec.date_format)


def _needs_jacobian(context: BuildContext) -> bool:
    """Whether the joint needs an explicit log Jacobian determinant term.

    Only for the methods that treat the unconstrained values as the random
    variables.  A MAP estimate is reported in the constrained space, so adding
    the Jacobian there would move the reported optimum.
    """
    if context.spec.method not in ("mcmc", "advi"):
        return False
    return any(free.is_transformed for free in context.free_parameters)
