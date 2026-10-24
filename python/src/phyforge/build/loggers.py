"""Trace and tree loggers."""

from __future__ import annotations

from typing import List, Optional

from ..nodes import model_ref, param_ref
from . import sitemodel
from .context import BuildContext


def tree_length_column(context: BuildContext) -> Optional[dict]:
    """The tree-length column, or None when the tree carries no branch lengths.

    Only worth logging for an unrooted tree, whose branch lengths are the free
    parameters and are otherwise left out of the trace. A time tree's branches
    are a rate times a height, both of which already have columns of their own.
    """
    if context.spec.clock.is_time_tree:
        return None
    return {"ref": model_ref(context.ids.tree), "quantity": "treeLength"}


def columns(
    context: BuildContext,
    targets: Optional[List[str]] = None,
    *,
    tree_length: bool = False,
) -> List:
    """The default set of logged columns: the target densities, then parameters.

    Vector parameters that are only a nuisance -- branch lengths, node height
    proportions -- are left out; the tree logger records those.
    """
    result: List = list(targets or [])
    if tree_length:
        column = tree_length_column(context)
        if column is not None:
            result.append(column)
    for free in context.free_parameters:
        if free.log:
            result.append(param_ref(free.id))
    # A model may hold quantities no free parameter stands for, such as the
    # empirical rates CAT reads off the alignment.
    result.extend(sitemodel.log_columns(context))
    return result


def trace(
    id_: str,
    context: BuildContext,
    *,
    file: Optional[str] = None,
    every: Optional[int] = None,
    targets: Optional[List[str]] = None,
    tree_length: bool = False,
) -> dict:
    node = {
        "id": id_,
        "type": "logger",
        "columns": columns(context, targets, tree_length=tree_length),
    }
    if file is not None:
        node["file"] = file
    if every is not None:
        node["every"] = every
    return node


def tree(
    id_: str,
    context: BuildContext,
    *,
    file: Optional[str] = None,
    every: Optional[int] = None,
) -> dict:
    node = {"id": id_, "type": "logger", "tree": model_ref(context.ids.tree)}
    if file is not None:
        node["file"] = file
    if every is not None:
        node["every"] = every
    return node


def patterns(
    id_: str,
    context: BuildContext,
    *,
    file: str,
    every: Optional[int] = None,
) -> dict:
    """A logger writing one column per site pattern log likelihood.

    The cpo and lppd estimators read its columns back by name, so this logger
    carries nothing else: the names they look for are built from the tree
    likelihood's id. A wide format keeps the exponentiated per-pattern means
    from losing their small terms to rounding.
    """
    node: dict = {
        "id": id_,
        "type": "logger",
        "file": file,
        "columns": [
            {
                "ref": model_ref(context.ids.tree_likelihood),
                "quantity": "patternLogLikelihood",
                "format": "%.10e",
            }
        ],
    }
    if every is not None:
        node["every"] = every
    return node


def target_refs(context: BuildContext) -> List[str]:
    """References to the densities worth logging for this spec."""
    ids = context.ids
    if context.spec.uses_priors:
        return [
            model_ref(ids.joint),
            model_ref(ids.tree_likelihood),
            model_ref(ids.prior),
        ]
    return [model_ref(ids.tree_likelihood)]
