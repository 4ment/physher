"""Applying per-parameter overrides where they need the model to judge them.

``validate`` checks what an override says on its own; whether a prior is a
density on the parameter it names, or a vector hyperparameter has the right
length, is only known once the builders have declared that parameter.  Those
checks live here, and record problems on the context rather than raising, so a
caller hears about every override at once.
"""

from __future__ import annotations

from typing import List, Optional

from ..nodes import distribution, model_ref, param_ref
from ..spec import NO_PRIOR, PRIOR_DISTRIBUTIONS, PriorChoice, Values
from ..validate import SpecError
from .context import BuildContext, FreeParameter

__all__ = ["PriorWriter", "check_unused", "vector_length_problem"]


class PriorWriter:
    """Writes the prior on each free parameter, the override's or the default.

    Every free parameter goes through :meth:`write` or :meth:`structural`
    exactly once; :meth:`remaining` then gives the ones no model builder has a
    default for -- kappa, say -- the chance to carry an overridden prior too.
    """

    def __init__(self, context: BuildContext):
        self.context = context
        self.seen: set = set()

    def write(
        self,
        free: Optional[FreeParameter],
        id_: str,
        default: Optional[PriorChoice],
    ) -> List[dict]:
        """The prior node on ``free``, or nothing if it ends up without one."""
        if free is None:
            return []
        self.seen.add(free.id)
        free.default_prior = default
        choice = self.context.override(free.id).prior or default
        if choice is not None and choice.distribution == NO_PRIOR:
            choice = None
        if choice is not None and not self._suits(free, choice):
            choice = default
        free.prior = choice
        if choice is None:
            return []
        extra = {}
        if choice.distribution == "ctmcscale":
            # The CTMC reference prior scales itself by the tree's length.
            extra["tree"] = model_ref(self.context.ids.tree)
        return [
            distribution(
                id_,
                choice.distribution,
                param_ref(free.id),
                parameters=dict(choice.parameters) or None,
                **extra,
            )
        ]

    def structural(self, free: Optional[FreeParameter], note: str) -> None:
        """Mark the prior on ``free`` as part of the model, not a choice."""
        if free is None:
            return
        self.seen.add(free.id)
        free.prior_note = note
        if self.context.override(free.id).prior is not None:
            self.context.problem(
                f"the prior on {free.id} is {note}, which cannot be replaced",
                free.id,
                "prior",
            )

    def remaining(self) -> List[dict]:
        """Priors on the free parameters no builder wrote a default for."""
        priors = []
        for free in list(self.context.free_parameters):
            if free.id not in self.seen:
                priors.extend(self.write(free, f"prior.{free.id}", None))
        return priors

    def _suits(self, free: FreeParameter, choice: PriorChoice) -> bool:
        family = PRIOR_DISTRIBUTIONS[choice.distribution]
        if free.kind not in family.kinds:
            self.context.problem(
                f"{choice.distribution} is a density on a {' or '.join(family.kinds)} "
                f"parameter, and {free.id} is {free.kind}",
                free.id,
                "prior",
            )
            return False
        for key, value in choice.parameters.items():
            problem = vector_length_problem(value, free.dimension)
            if problem:
                self.context.problem(
                    f"{choice.distribution} {key} for {free.id}: {problem}",
                    free.id,
                    "prior",
                )
                return False
        return True


def vector_length_problem(value: Values, dimension: int) -> Optional[str]:
    """Why ``value`` cannot be given to a parameter of ``dimension``, if it cannot."""
    if isinstance(value, (list, tuple)) and len(value) not in (1, dimension):
        return f"expected 1 or {dimension} value(s), not {len(value)}"
    return None


def check_unused(context: BuildContext) -> None:
    """Record an override that names no parameter of this analysis."""
    known = [free.id for free in context.all_parameters]
    for id_ in context.spec.parameters:
        if id_ not in known:
            context.problems.append(
                SpecError(
                    f"this analysis has no parameter {id_!r} (it has "
                    f"{', '.join(known)})",
                    f"parameters.{id_}",
                )
            )
