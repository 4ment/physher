"""The ``phyforge`` command."""

from __future__ import annotations

import argparse
import json
import sys
from typing import List, Optional, Sequence

from ..build import ParameterInfo, build, describe
from ..validate import SpecError, check
from .args import make_parser, to_spec


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["gui"]:
        from ..gui.launch import launch

        return launch(argv[1:])

    parser = make_parser()
    arguments = parser.parse_args(argv)
    try:
        spec = to_spec(arguments)
    except (argparse.ArgumentTypeError, ValueError, TypeError, OSError) as error:
        parser.error(str(error))

    problems = check(spec)
    if problems:
        for problem in problems:
            print(f"phyforge: {problem}", file=sys.stderr)
        return 2

    if arguments.list_parameters:
        rows, problems = describe(spec)
        print(format_parameters(rows, spec.method))
        for problem in problems:
            print(f"phyforge: {problem}", file=sys.stderr)
        return 2 if problems else 0

    try:
        document = build(spec)
    except (SpecError, ValueError, KeyError, OSError) as error:
        print(f"phyforge: {error}", file=sys.stderr)
        return 1

    text = json.dumps(document, indent=2)
    if arguments.output:
        with open(arguments.output, "w") as fp:
            fp.write(text + "\n")
    else:
        print(text)
    return 0


def format_parameters(rows: List[ParameterInfo], method: str) -> str:
    """The ``--list-parameters`` table: one row per parameter, as overridden."""
    headers = ["id", "kind", "dim", "estimate", "init", "prior"]
    if method == "advi":
        headers += ["mu", "sigma"]
    if method == "mcmc":
        headers.append("weight")
    table = [headers]
    for row in rows:
        if row.prior_note:
            prior = f"({row.prior_note})"
        else:
            prior = "-" if row.prior is None else str(row.prior)
        line = [
            row.id,
            row.kind,
            str(row.dimension),
            "yes" if row.estimated else "fixed",
            _short(row.init),
            prior if method != "ml" else "-",
        ]
        if method == "advi":
            line += [_short(row.mu), _short(row.sigma)]
        if method == "mcmc":
            line.append(_short(row.weight))
        table.append(line)
    widths = [max(len(line[i]) for line in table) for i in range(len(headers))]
    return "\n".join(
        "  ".join(cell.ljust(width) for cell, width in zip(line, widths)).rstrip()
        for line in table
    )


def _short(values) -> str:
    """A value for the table: vectors that repeat one number are shown once."""
    if values is None:
        return "-"
    if not isinstance(values, list):
        return f"{values:.6g}"
    if len(set(values)) == 1:
        return f"{values[0]:.6g} (x{len(values)})"
    shown = ",".join(f"{value:.6g}" for value in values[:5])
    return shown + (",..." if len(values) > 5 else "")


if __name__ == "__main__":
    sys.exit(main())
