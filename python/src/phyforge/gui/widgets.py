"""Widgets bound to spec fields and to per-parameter overrides.

Every bound widget has the key ``w:<path>`` or ``p:<id>:<attribute>``, is
seeded from the stored value when Streamlit has no state for it, and writes
back from ``on_change``.  The problems :func:`phyforge.describe` reports for a
field are shown right under the widget that sets it.
"""

from __future__ import annotations

from typing import Any, Callable, Iterable, List, Optional, Sequence

import streamlit as st

from ..build import ParameterInfo
from ..spec import (
    NO_PRIOR,
    PRIOR_DISTRIBUTIONS,
    PriorChoice,
    format_values,
    parse_values,
)
from . import state
from .state import Analysis

#: Starting hyperparameters when a distribution is picked, so that choosing
#: one gives a valid prior straight away.
DEFAULT_HYPERPARAMETERS = {
    "lambda": 1.0,
    "mean": 1.0,
    "shape": 1.0,
    "rate": 1.0,
    "scale": 1.0,
    "mu": 0.0,
    "sigma": 1.0,
    "stdev": 1.0,
    "tau": 1.0,
    "alpha": 1.0,
    "beta": 1.0,
    "a": 1.0,
    "b": 1.0,
    "location": 0.0,
}

DEFAULT = "default"


# --------------------------------------------------------------------------
# spec fields
# --------------------------------------------------------------------------


def bind(
    widget: Callable[..., Any],
    path: str,
    label: str,
    analysis: Analysis,
    *,
    encode: Callable[[Any], Any] = lambda value: value,
    decode: Callable[[Any], Any] = lambda value: value,
    container: Any = None,
    **kwargs: Any,
) -> Any:
    """Draw ``widget`` for the spec field at ``path`` and return its value."""
    key = f"w:{path}"
    if key not in st.session_state:
        st.session_state[key] = encode(state.values()[path])

    def changed() -> None:
        state.values()[path] = decode(st.session_state[key])

    target = container or st
    raw = getattr(target, widget.__name__)(label, key=key, on_change=changed, **kwargs)
    show_problems(analysis, path, container=target)
    return decode(raw)


def bind_session(
    widget: Callable[..., Any],
    name: str,
    label: str,
    *,
    container: Any = None,
    **kwargs: Any,
) -> Any:
    """Like :func:`bind`, for a GUI-only setting kept at ``session_state[name]``."""
    key = f"w:@{name}"
    if key not in st.session_state:
        st.session_state[key] = st.session_state[name]

    def changed() -> None:
        st.session_state[name] = st.session_state[key]

    target = container or st
    return getattr(target, widget.__name__)(label, key=key, on_change=changed, **kwargs)


def select(
    path: str,
    label: str,
    options: Sequence[Any],
    analysis: Analysis,
    **kwargs: Any,
) -> Any:
    """A selectbox for ``path`` whose stored value falls back to the first option.

    The options of one field can depend on another (the substitution models
    on the data type), so a stored value may no longer be offered.
    """
    options = list(options)
    key = f"w:{path}"
    if state.values()[path] not in options:
        state.values()[path] = options[0]
        st.session_state.pop(key, None)
    return bind(st.selectbox, path, label, analysis, options=options, **kwargs)


def optional_text(path: str, label: str, analysis: Analysis, **kwargs: Any) -> Any:
    """A text input for an optional string field: blank means None."""
    return bind(
        st.text_input,
        path,
        label,
        analysis,
        encode=lambda value: "" if value is None else str(value),
        decode=lambda value: value.strip() or None,
        **kwargs,
    )


def show_problems(analysis: Analysis, path: str, container: Any = None) -> None:
    for message in analysis.problems_for(path):
        (container or st).error(message, icon="⚠️")


# --------------------------------------------------------------------------
# per-parameter overrides
# --------------------------------------------------------------------------


def short(values: Any) -> str:
    """A compact rendering of a value for placeholders and captions."""
    if values is None:
        return ""
    if isinstance(values, list):
        if len(set(values)) == 1:
            return f"{values[0]:.6g} × {len(values)}"
        shown = ", ".join(f"{value:.4g}" for value in values[:4])
        return shown + (", …" if len(values) > 4 else "")
    return f"{values:.6g}"


def _values_input(
    container: Any,
    row: ParameterInfo,
    attribute: str,
    label: str,
    placeholder: str,
    analysis: Analysis,
    *,
    scalar: bool = False,
) -> None:
    """A text input for a vector-or-scalar override field; blank clears it."""
    key = f"p:{row.id}:{attribute}"
    error_key = f"{key}:error"
    if key not in st.session_state:
        current = getattr(state.overrides().get(row.id), attribute, None)
        st.session_state[key] = "" if current is None else format_values(current)

    def changed() -> None:
        text = st.session_state[key].strip()
        st.session_state.pop(error_key, None)
        if not text:
            state.set_override(row.id, **{attribute: None})
            return
        try:
            value = parse_values(text)
            if scalar and isinstance(value, list):
                raise ValueError("expected a single number")
        except ValueError as error:
            st.session_state[error_key] = str(error)
            return
        state.set_override(row.id, **{attribute: value})

    container.text_input(
        label,
        key=key,
        on_change=changed,
        placeholder=placeholder,
        label_visibility="collapsed",
    )
    if error_key in st.session_state:
        container.error(st.session_state[error_key], icon="⚠️")
    show_problems(analysis, f"parameters.{row.id}.{attribute}", container=container)


def _estimate_checkbox(container: Any, row: ParameterInfo, analysis: Analysis) -> None:
    key = f"p:{row.id}:estimate"
    if key not in st.session_state:
        override = state.overrides().get(row.id)
        st.session_state[key] = True if override is None else override.estimate

    def changed() -> None:
        state.set_override(row.id, estimate=st.session_state[key])

    container.checkbox(
        "estimate",
        key=key,
        on_change=changed,
        help="Untick to hold this parameter at its starting value",
    )


def _label(container: Any, row: ParameterInfo) -> None:
    container.markdown(f"`{row.id}`")
    container.caption(f"{row.kind} · {row.dimension}")


def init_table(
    analysis: Analysis, roles: Iterable[str], title: str = "Parameters"
) -> None:
    """Starting value and estimate/fix, one row per parameter of ``roles``."""
    rows = analysis.rows_for(set(roles))
    if not rows:
        return
    st.subheader(title)
    st.caption(
        "Leave a starting value blank for the default shown. A vector takes one "
        "number for every element, or a comma-separated list."
    )
    header = st.columns([2, 4, 1])
    header[0].markdown("**Parameter**")
    header[1].markdown("**Starting value**")
    header[2].markdown("**Estimate**")
    for row in rows:
        columns = st.columns([2, 4, 1])
        _label(columns[0], row)
        _values_input(
            columns[1],
            row,
            "init",
            "start",
            f"default: {short(row.default_init)}",
            analysis,
        )
        _estimate_checkbox(columns[2], row, analysis)
        show_problems(analysis, f"parameters.{row.id}.estimate")


def prior_table(analysis: Analysis) -> None:
    """The prior of every estimated parameter, editable where it is a choice."""
    rows = [row for row in analysis.rows if row.estimated]
    if not rows:
        st.info("Load an alignment and a tree to list the parameters.")
        return
    for row in rows:
        with st.container(border=True):
            columns = st.columns([2, 3, 5])
            _label(columns[0], row)
            if not row.prior_editable:
                columns[1].markdown(f"_{row.prior_note}_")
                columns[2].caption("Part of the model; not a separate choice.")
                show_problems(analysis, f"parameters.{row.id}.prior")
                continue
            _prior_editor(columns[1], columns[2], row, analysis)
            show_problems(analysis, f"parameters.{row.id}.prior")


def _prior_editor(
    left: Any, right: Any, row: ParameterInfo, analysis: Analysis
) -> None:
    base = f"p:{row.id}:prior"
    override = state.overrides().get(row.id)
    current: Optional[PriorChoice] = override.prior if override else None

    options = [DEFAULT, NO_PRIOR] + [
        name for name, family in PRIOR_DISTRIBUTIONS.items() if row.kind in family.kinds
    ]
    default_text = str(row.default_prior) if row.default_prior else "no prior"
    if f"{base}.dist" not in st.session_state:
        st.session_state[f"{base}.dist"] = (
            current.distribution if current is not None else DEFAULT
        )

    def commit() -> None:
        _commit_prior(row.id, base)

    distribution = left.selectbox(
        "distribution",
        options,
        key=f"{base}.dist",
        on_change=commit,
        format_func=lambda name: (
            f"default: {default_text}" if name == DEFAULT else name
        ),
        label_visibility="collapsed",
    )
    if distribution in (DEFAULT, NO_PRIOR):
        if distribution == NO_PRIOR:
            right.caption("No prior: the parameter is left improper.")
        return

    family = PRIOR_DISTRIBUTIONS[distribution]
    sets = family.hyperparameters
    if len(sets) > 1:
        if f"{base}.set" not in st.session_state:
            chosen = sets[0]
            if current is not None and current.distribution == distribution:
                chosen = next(
                    (keys for keys in sets if set(keys) == set(current.parameters)),
                    sets[0],
                )
            st.session_state[f"{base}.set"] = ",".join(chosen)
        left.selectbox(
            "parameterization",
            [",".join(keys) for keys in sets],
            key=f"{base}.set",
            on_change=commit,
            label_visibility="collapsed",
        )
    keys = _chosen_keys(base, sets)
    if not keys:
        right.caption("Takes no hyperparameters.")
        return
    inputs = right.columns(len(keys))
    for column, name in zip(inputs, keys):
        key = f"{base}.{name}"
        if key not in st.session_state:
            value = None
            if current is not None and current.distribution == distribution:
                value = current.parameters.get(name)
            if value is None:
                value = DEFAULT_HYPERPARAMETERS.get(name, 1.0)
            st.session_state[key] = format_values(value)
        column.text_input(name, key=key, on_change=commit)
    if f"{base}:error" in st.session_state:
        right.error(st.session_state[f"{base}:error"], icon="⚠️")


def _chosen_keys(base: str, sets) -> List[str]:
    if len(sets) == 1:
        return list(sets[0])
    chosen = st.session_state.get(f"{base}.set") or ",".join(sets[0])
    return [key for key in chosen.split(",") if key]


def _commit_prior(id_: str, base: str) -> None:
    """Turn the prior widgets of one parameter into its override."""
    st.session_state.pop(f"{base}:error", None)
    distribution = st.session_state[f"{base}.dist"]
    if distribution == DEFAULT:
        state.set_override(id_, prior=None)
        return
    if distribution == NO_PRIOR:
        state.set_override(id_, prior=PriorChoice(NO_PRIOR))
        return
    sets = PRIOR_DISTRIBUTIONS[distribution].hyperparameters
    names = _chosen_keys(base, sets)
    if set(names) not in [set(keys) for keys in sets]:
        names = list(sets[0])
        st.session_state[f"{base}.set"] = ",".join(names)
    parameters = {}
    for name in names:
        text = st.session_state.get(f"{base}.{name}")
        if text is None:
            text = format_values(DEFAULT_HYPERPARAMETERS.get(name, 1.0))
        try:
            parameters[name] = parse_values(text)
        except ValueError:
            st.session_state[f"{base}:error"] = f"{name} must be a number"
            return
    state.set_override(id_, prior=PriorChoice(distribution, parameters))


def variational_table(analysis: Analysis) -> None:
    """advi: where each factor of the mean-field approximation starts."""
    rows = [row for row in analysis.rows if row.estimated]
    if not rows:
        return
    st.caption(
        "One normal per parameter, on the unconstrained value; a simplex of K "
        "proportions has K − 1 of them. Blank keeps the default shown."
    )
    header = st.columns([2, 3, 3])
    header[0].markdown("**Parameter**")
    header[1].markdown("**μ**")
    header[2].markdown("**σ**")
    for row in rows:
        override = state.overrides().get(row.id)
        columns = st.columns([2, 3, 3])
        _label(columns[0], row)
        mu_default = row.mu if override is None or override.mu is None else None
        sigma_default = (
            row.sigma if override is None or override.sigma is None else None
        )
        _values_input(
            columns[1], row, "mu", "mu", f"default: {short(mu_default)}", analysis
        )
        _values_input(
            columns[2],
            row,
            "sigma",
            "sigma",
            f"default: {short(sigma_default)}",
            analysis,
        )


def weight_table(analysis: Analysis) -> None:
    """mcmc: how often each parameter's operators are picked."""
    rows = [row for row in analysis.rows if row.estimated]
    if not rows:
        return
    st.caption(
        "Operators are picked in proportion to their weight; the default is "
        "the parameter's dimension, which evens out the moves per element."
    )
    for row in rows:
        columns = st.columns([2, 3, 3])
        _label(columns[0], row)
        _values_input(
            columns[1],
            row,
            "weight",
            "weight",
            f"default: {row.dimension}",
            analysis,
            scalar=True,
        )
