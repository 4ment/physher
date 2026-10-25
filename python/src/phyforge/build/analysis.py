"""Assembling a whole physher document from a spec.

physher reads the top-level object in order, building every entry except
``physher`` (the run list) and ``init`` as a model, so a node may only
reference ids declared before it.  The ordering below reflects that.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from ..io.dates import dates_from_csv, dates_from_names
from ..io.fasta import read_fasta
from ..io.trees import taxa_from_tree
from ..nodes import compound, model_ref
from ..spec import AnalysisSpec, PriorChoice, Values
from ..validate import SpecError, check, validate
from . import actions, likelihood, priors, variational
from .context import BuildContext
from .overrides import check_unused

__all__ = ["build", "describe", "make_context", "ParameterInfo"]


def build(spec: AnalysisSpec, context: Optional[BuildContext] = None) -> Dict:
    """Build the physher JSON document for ``spec``.

    Raises :class:`phyforge.validate.SpecError` if the spec is inconsistent,
    so a caller never gets a config that cannot run.
    """
    validate(spec)
    if context is None:
        context = make_context(spec)
    document = _assemble(spec, context)
    check_unused(context)
    if context.problems:
        raise context.problems[0]
    return document


@dataclass
class ParameterInfo:
    """One parameter of an analysis, as a front end shows it for editing.

    Each ``default_*`` is what the builders choose on their own; the field
    without the prefix is what the document actually carries once the spec's
    overrides are applied.  ``mu``/``sigma`` are only set under advi and
    ``weight`` only under mcmc.
    """

    id: str
    role: str
    kind: str
    dimension: int
    #: number of unconstrained values: one fewer than ``dimension`` for a
    #: transformed simplex
    free_dimension: int
    estimated: bool
    init: Values
    default_init: Values
    prior: Optional[PriorChoice]
    default_prior: Optional[PriorChoice]
    #: why the prior cannot be overridden, when it is part of the model
    prior_note: Optional[str]
    mu: Optional[Values]
    sigma: Optional[Values]
    weight: Optional[float]

    @property
    def prior_editable(self) -> bool:
        return self.prior_note is None


def describe(
    spec: AnalysisSpec, *, strict: bool = True
) -> Tuple[List[ParameterInfo], List[SpecError]]:
    """Every parameter of ``spec``'s analysis, and every problem with the spec.

    This is a dry run of :func:`build` that reports rather than raises: a
    front end calls it to list what can be edited, and to show every problem
    at once.  With ``strict=False`` an override naming a parameter this
    analysis does not have is not a problem -- a GUI keeps those around for
    when the model that has the parameter is selected again.  Nothing is
    listed while the spec itself is invalid, since then there is no model.
    """
    problems = check(spec)
    if problems:
        return [], problems
    context = make_context(spec)
    _assemble(spec, context)
    if strict:
        check_unused(context)
    rows = [
        ParameterInfo(
            id=free.id,
            role=free.role,
            kind=free.kind,
            dimension=free.dimension,
            free_dimension=free.free_dimension,
            estimated=free.estimated,
            init=free.init,
            default_init=free.default_init,
            prior=free.prior,
            default_prior=free.default_prior,
            prior_note=free.prior_note,
            mu=free.mu,
            sigma=free.sigma,
            weight=free.weight,
        )
        for free in context.all_parameters
    ]
    return rows, list(context.problems)


def _assemble(spec: AnalysisSpec, context: BuildContext) -> Dict:
    """Every node of the document, recording override problems on ``context``."""
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
