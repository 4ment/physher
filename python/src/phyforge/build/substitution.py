"""The substitution model node."""

from __future__ import annotations

from ..nodes import simplex
from ..spec import SUBSTITUTION_MODELS
from .context import POSITIVE, BuildContext
from .parameters import free_parameter, free_simplex

#: Number of states per data type, which sizes the frequencies simplex.
STATE_COUNT = {"nucleotide": 4, "aa": 20}

#: Models whose rate matrix is symmetric with equal frequencies by definition.
EQUAL_FREQUENCY_MODELS = ("JC69", "K80", "SYM")

#: physher's own "k80" kind does not optimize: its kappa runs away and the
#: likelihood diverges (-56457 against -4207 for the equivalent model on the
#: fluA example).  K80 is HKY with the frequencies pinned at 1/4, so it is
#: emitted that way instead -- which is what examples/fluA/K80-ELBO.json does.
KIND_SUBSTITUTIONS = {"k80": "hky"}

ROLE_FREQUENCIES = "frequencies"
ROLE_RATES = "rates"
ROLE_KAPPA = "kappa"


def state_count(datatype: str) -> int:
    try:
        return STATE_COUNT[datatype]
    except KeyError:
        raise ValueError(f"unsupported data type {datatype!r}") from None


def build(context: BuildContext) -> dict:
    """Build the ``substitutionmodel`` node for the spec's model."""
    spec = context.spec.substitution
    datatype = SUBSTITUTION_MODELS[spec.kind.upper()]
    states = state_count(datatype)
    ids = context.ids

    kind = spec.normalized_kind
    node = {
        "id": ids.substitution_model,
        "type": "substitutionmodel",
        "kind": KIND_SUBSTITUTIONS.get(kind, kind),
        "datatype": datatype,
    }

    node["frequencies"] = _frequencies(context, spec, states)

    if spec.has_free_rates:
        # The six exchangeabilities of a reversible nucleotide model, on the
        # simplex because only their ratios matter.
        node["rates"] = free_simplex(context, ids.rates, 6, role=ROLE_RATES)
    elif spec.has_kappa:
        node["rates"] = {
            "kappa": free_parameter(
                context,
                ids.kappa,
                spec.kappa_init,
                role=ROLE_KAPPA,
                kind=POSITIVE,
            )
        }
    return node


def _frequencies(context: BuildContext, spec, states: int) -> dict:
    ids = context.ids
    if spec.kind.upper() in EQUAL_FREQUENCY_MODELS or spec.frequencies == "equal":
        # Pinned: spelled out as values so nothing registers it as free.
        return simplex(ids.frequencies, [1.0 / states] * states)
    return free_simplex(context, ids.frequencies, states, role=ROLE_FREQUENCIES)
