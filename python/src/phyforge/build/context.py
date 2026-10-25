"""Build-time state shared by the model builders.

The important piece is the free-parameter registry.  Rather than walking the
finished JSON to guess which parameters are estimated -- fragile, and blind to
transforms -- each builder *declares* the free parameters it creates.  The
downstream builders (MCMC operators, variational distributions, the Jacobian
term, loggers) are then straightforward maps over that list.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from ..ids import Ids
from ..spec import (  # noqa: F401 -- re-exported for the builders
    POSITIVE,
    REAL,
    SIMPLEX,
    UNIT,
    AnalysisSpec,
    ParameterOverride,
    PriorChoice,
    Values,
)
from ..validate import SpecError

# REAL: the whole real line, untransformed.  POSITIVE: log/exp transform.
# UNIT: (0, 1), logit transform.  SIMPLEX: stick-breaking, one fewer free
# element.


@dataclass
class FreeParameter:
    """A quantity the analysis estimates."""

    id: str
    role: str
    kind: str = REAL
    dimension: int = 1
    #: id of the unconstrained parameter backing this one, when transformed
    unconstrained_id: Optional[str] = None
    #: starting value in the unconstrained space, used to centre a
    #: variational approximation
    unconstrained_init: Values = 0.0
    #: whether this parameter should appear in trace logs by default
    log: bool = True
    #: starting value in the constrained space, as written into the document
    init: Optional[Values] = None
    #: the starting value the builder would have used without an override
    default_init: Optional[Values] = None
    #: False for a parameter an override holds fixed
    estimated: bool = True
    #: the prior the builder wrote, and the one it would have written without
    #: an override; None means no prior
    prior: Optional[PriorChoice] = None
    default_prior: Optional[PriorChoice] = None
    #: set when the prior is part of the model's structure (a coalescent, a
    #: GMRF, a reference distribution) and so cannot be overridden
    prior_note: Optional[str] = None
    #: the variational normal's location and scale, once advi has placed it
    mu: Optional[Values] = None
    sigma: Optional[Values] = None
    #: the weight of each of its MCMC operators, once mcmc has built them
    weight: Optional[float] = None

    @property
    def free_id(self) -> str:
        """The id an optimizer or variational distribution attaches to.

        Gradients accumulate in the leaves of the autodiff graph, so a
        transformed parameter is driven through its unconstrained backing.
        """
        return self.unconstrained_id or self.id

    @property
    def is_transformed(self) -> bool:
        return self.unconstrained_id is not None

    @property
    def free_dimension(self) -> int:
        """Number of free values -- a K-simplex has K-1 of them."""
        if self.kind == SIMPLEX and self.is_transformed:
            return self.dimension - 1
        return self.dimension


@dataclass
class BuildContext:
    """Everything the builders need beyond the spec itself."""

    spec: AnalysisSpec
    taxa: List[str]
    dates: Optional[Dict[str, float]] = None
    ids: Ids = field(default_factory=Ids)
    free_parameters: List[FreeParameter] = field(default_factory=list)
    #: parameters an override holds at their starting value; kept apart so
    #: nothing that maps over the free ones has to skip them
    fixed_parameters: List[FreeParameter] = field(default_factory=list)
    #: problems with the overrides, found while the builders apply them
    problems: List[SpecError] = field(default_factory=list)

    @property
    def taxa_count(self) -> int:
        return len(self.taxa)

    @property
    def branch_count(self) -> int:
        """Number of branch lengths in an unrooted tree."""
        return 2 * self.taxa_count - 3

    @property
    def rooted_branch_count(self) -> int:
        """Number of branches in a rooted tree: every node except the root."""
        return 2 * self.taxa_count - 2

    @property
    def internal_count(self) -> int:
        """Number of non-root internal nodes, i.e. the proportions dimension."""
        return self.taxa_count - 2

    def register(self, parameter: FreeParameter) -> FreeParameter:
        if parameter.estimated:
            self.free_parameters.append(parameter)
        else:
            self.fixed_parameters.append(parameter)
        return parameter

    def override(self, id_: str) -> ParameterOverride:
        return self.spec.override(id_)

    def problem(self, message: str, id_: str, attribute: str) -> None:
        self.problems.append(SpecError(message, f"parameters.{id_}.{attribute}"))

    @property
    def all_parameters(self) -> List[FreeParameter]:
        """Every parameter an override could name, estimated or fixed."""
        return self.free_parameters + self.fixed_parameters

    def by_id(self, id_: str) -> Optional[FreeParameter]:
        for parameter in self.all_parameters:
            if parameter.id == id_:
                return parameter
        return None

    def by_role(self, role: str) -> Optional[FreeParameter]:
        for parameter in self.free_parameters:
            if parameter.role == role:
                return parameter
        return None


def to_unconstrained(value: Values, kind: str) -> Values:
    """Map a constrained starting value into the unconstrained space.

    A vector maps element by element, except a simplex, which goes through
    the stick-breaking transform as a whole.
    """
    if kind == SIMPLEX:
        return simplex_to_unconstrained(value)
    if isinstance(value, (list, tuple)):
        return [_scalar_to_unconstrained(v, kind) for v in value]
    return _scalar_to_unconstrained(value, kind)


def simplex_to_unconstrained(values: Sequence[float]) -> List[float]:
    """Inverse of physher's stick-breaking map (``_transform_simplex_stan``).

    ``K`` proportions become ``K - 1`` unconstrained values; equal proportions
    map to zeros, which is where a simplex starts when none is given.
    """
    dimension = len(values) - 1
    result, used = [], 0.0
    for i in range(dimension):
        z = values[i] / (1.0 - used)
        result.append(math.log(z / (1.0 - z)) + math.log(dimension - i))
        used += values[i]
    return result


def _scalar_to_unconstrained(value: float, kind: str) -> float:
    if kind == POSITIVE:
        if value <= 0.0:
            raise ValueError(f"a positive parameter cannot start at {value}")
        return math.log(value)
    if kind == UNIT:
        if not 0.0 < value < 1.0:
            raise ValueError(f"a unit parameter cannot start at {value}")
        return math.log(value / (1.0 - value))
    return value
