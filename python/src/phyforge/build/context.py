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
from typing import Dict, List, Optional

from ..ids import Ids
from ..spec import AnalysisSpec

#: A parameter constrained to the whole real line.
REAL = "real"
#: Constrained to be positive; log/exp transform.
POSITIVE = "positive"
#: Constrained to (0, 1); logit transform.
UNIT = "unit"
#: A simplex; stick-breaking transform with one fewer free element.
SIMPLEX = "simplex"


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
    unconstrained_init: float = 0.0
    #: whether this parameter should appear in trace logs by default
    log: bool = True

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
        self.free_parameters.append(parameter)
        return parameter

    def by_role(self, role: str) -> Optional[FreeParameter]:
        for parameter in self.free_parameters:
            if parameter.role == role:
                return parameter
        return None


def to_unconstrained(value: float, kind: str) -> float:
    """Map a constrained starting value into the unconstrained space."""
    if kind == POSITIVE:
        if value <= 0.0:
            raise ValueError(f"a positive parameter cannot start at {value}")
        return math.log(value)
    if kind == UNIT:
        if not 0.0 < value < 1.0:
            raise ValueError(f"a unit parameter cannot start at {value}")
        return math.log(value / (1.0 - value))
    return value
