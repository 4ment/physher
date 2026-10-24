"""The ``phyforge`` command."""

from __future__ import annotations

import json
import sys
from typing import Optional, Sequence

from ..build import build
from ..validate import SpecError, check
from .args import make_parser, to_spec


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = make_parser()
    arguments = parser.parse_args(argv)
    spec = to_spec(arguments)

    problems = check(spec)
    if problems:
        for problem in problems:
            print(f"phyforge: {problem}", file=sys.stderr)
        return 2

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


if __name__ == "__main__":
    sys.exit(main())
