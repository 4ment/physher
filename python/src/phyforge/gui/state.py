"""Session state: the spec being edited, and what it currently amounts to.

The spec lives in ``st.session_state`` as a flat mapping from dotted path
(``site.categories``) to value, seeded from the *CLI's* defaults so that the
command shown in the Output tab only lists what the user changed.  Widgets are
bound to a path (see :mod:`phyforge.gui.widgets`) and write back to it from an
``on_change`` callback, which runs before the rerun -- so by the time a page is
drawn, :func:`analyse` has already seen every edit.

Keeping values outside the widgets matters because Streamlit forgets a
widget's state on any run in which it is not drawn: the clock options would
reset every time the tree was made unrooted and rooted again.

Per-parameter overrides are kept for every id the user ever edited, even ones
the current model does not have; only those it has are passed to the
builders, so switching from HKY to GTR and back restores ``kappa``'s.
"""

from __future__ import annotations

import dataclasses
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import streamlit as st

from ..build import ParameterInfo, describe
from ..cli.args import parse
from ..spec import AnalysisSpec, ParameterOverride
from ..validate import SpecError

VALUES = "phyforge.values"
OVERRIDES = "phyforge.overrides"
FILES = "phyforge.files"
#: GUI-only switches that are not spec fields
TIME_TREE = "phyforge.time_tree"
DATATYPE = "phyforge.datatype"
#: id of the spec revision every widget key was last seeded from; bumped to
#: throw every widget's state away when a session is loaded
REVISION = "phyforge.revision"


def cli_defaults() -> AnalysisSpec:
    """The spec ``phyforge ml -i ALN -t TREE`` builds, with no data yet."""
    spec = parse(["ml", "-i", "-", "-t", "-"])
    spec.data.alignment = None
    spec.data.tree = None
    # The GUI toggles the clock with its own switch and keeps the choice of
    # clock for when it is turned on again.
    spec.clock.kind = "strict"
    return spec


def flatten(data: Dict[str, Any], prefix: str = "") -> Dict[str, Any]:
    flat: Dict[str, Any] = {}
    for key, value in data.items():
        if isinstance(value, dict) and key != "parameters":
            flat.update(flatten(value, f"{prefix}{key}."))
        else:
            flat[f"{prefix}{key}"] = value
    return flat


def unflatten(flat: Dict[str, Any]) -> Dict[str, Any]:
    nested: Dict[str, Any] = {}
    for path, value in flat.items():
        target = nested
        *heads, last = path.split(".")
        for head in heads:
            target = target.setdefault(head, {})
        target[last] = value
    return nested


def init_state() -> None:
    state = st.session_state
    if VALUES not in state:
        defaults = cli_defaults().to_dict()
        defaults.pop("parameters")
        state[VALUES] = flatten(defaults)
        state[OVERRIDES] = {}
        state[FILES] = {}
        state[TIME_TREE] = False
        state[DATATYPE] = "nucleotide"
        state[REVISION] = 0


def values() -> Dict[str, Any]:
    return st.session_state[VALUES]


def overrides() -> Dict[str, ParameterOverride]:
    return st.session_state[OVERRIDES]


def files() -> Dict[str, Dict[str, Any]]:
    """Input files by role (alignment, tree, dates): original name and path."""
    return st.session_state[FILES]


def is_time_tree() -> bool:
    return bool(st.session_state[TIME_TREE])


def spec_from_values(
    flat: Dict[str, Any],
    parameters: Dict[str, ParameterOverride],
    time_tree: bool,
) -> AnalysisSpec:
    nested = unflatten(flat)
    if not time_tree:
        nested["clock"]["kind"] = None
        nested["coalescent"]["kind"] = None
    nested["parameters"] = {
        id_: dataclasses.asdict(override) for id_, override in parameters.items()
    }
    return AnalysisSpec.from_dict(nested)


@dataclass
class Analysis:
    """What the current settings amount to, worked out once per run."""

    #: the spec with only the overrides the model has a parameter for
    spec: AnalysisSpec
    rows: List[ParameterInfo]
    problems: List[SpecError]
    #: ids of overrides kept for parameters this model does not have
    inactive: List[str]
    by_field: Dict[str, List[str]] = field(default_factory=dict)
    #: fields a widget showed the problem for, so the rest are listed once
    shown: set = field(default_factory=set)

    def problems_for(self, path: str) -> List[str]:
        self.shown.add(path)
        return self.by_field.get(path, [])

    def row(self, id_: str) -> Optional[ParameterInfo]:
        return next((row for row in self.rows if row.id == id_), None)

    def rows_for(self, roles) -> List[ParameterInfo]:
        return [row for row in self.rows if row.role in roles]


def analyse() -> Analysis:
    everything = spec_from_values(values(), overrides(), is_time_tree())
    try:
        rows, problems = describe(everything, strict=False)
    except (OSError, ValueError, KeyError) as error:
        # describe validates first, so this is a file that exists but cannot
        # be read as what it claims to be.
        rows, problems = [], [SpecError(str(error))]
    known = {row.id for row in rows}
    active = {id_: o for id_, o in overrides().items() if id_ in known}
    inactive = [id_ for id_ in overrides() if id_ not in known] if rows else []
    spec = dataclasses.replace(everything, parameters=active)
    by_field: Dict[str, List[str]] = defaultdict(list)
    for problem in problems:
        by_field[problem.field or ""].append(ValueError.__str__(problem))
    return Analysis(spec, rows, problems, inactive, dict(by_field))


def set_override(id_: str, **changes: Any) -> None:
    """Change fields of one override, dropping it once it is back to defaults."""
    current = overrides().get(id_, ParameterOverride())
    updated = dataclasses.replace(current, **changes)
    if updated.is_empty:
        overrides().pop(id_, None)
    else:
        overrides()[id_] = updated


def reset_widgets() -> None:
    """Forget every widget's state, so each is seeded again from the values."""
    for key in list(st.session_state):
        if isinstance(key, str) and key.startswith(("w:", "p:", "local:")):
            del st.session_state[key]
    st.session_state[REVISION] += 1
