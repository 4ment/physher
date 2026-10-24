"""Reading taxon names out of a tree file.

The extension is not trusted -- ``examples/fluA/fluA-rooted.nxs`` is plain
newick -- so the format is sniffed from the first non-blank line.
"""

from __future__ import annotations

import re
from enum import Enum
from typing import List


class TreeFormat(str, Enum):
    NEWICK = "newick"
    NEXUS = "nexus"


def sniff_format(path: str) -> TreeFormat:
    with open(path, "r") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            return (
                TreeFormat.NEXUS
                if line.upper().startswith("#NEXUS")
                else TreeFormat.NEWICK
            )
    raise ValueError(f"{path} is empty")


def read_newick(path: str) -> str:
    """Return the first tree in ``path`` as a newick string.

    For a nexus file this is the tree definition with any TRANSLATE table
    already substituted back in, so the result always carries real taxon names.
    """
    if sniff_format(path) == TreeFormat.NEWICK:
        with open(path, "r") as fp:
            text = fp.read().strip()
        return _first_tree(text)

    with open(path, "r") as fp:
        text = fp.read()

    translation = _read_translate_block(text)
    match = re.search(
        r"^\s*tree\s+[^=]*=\s*(?:\[[^\]]*\]\s*)*(\(.*?;)",
        text,
        re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )
    if match is None:
        raise ValueError(f"no tree found in {path}")
    newick = match.group(1).strip()
    if translation:
        newick = _apply_translation(newick, translation)
    return newick


def taxa_from_tree(path: str) -> List[str]:
    """The tip labels of the first tree in ``path``, in newick order."""
    return _tip_labels(read_newick(path))


def _first_tree(text: str) -> str:
    end = text.find(";")
    return text if end == -1 else text[: end + 1]


def _read_translate_block(text: str) -> dict:
    match = re.search(r"\btranslate\b(.*?);", text, re.IGNORECASE | re.DOTALL)
    if match is None:
        return {}
    translation = {}
    for entry in match.group(1).split(","):
        parts = entry.split(None, 1)
        if len(parts) == 2:
            translation[parts[0].strip()] = _unquote(parts[1].strip())
    return translation


def _apply_translation(newick: str, translation: dict) -> str:
    def substitute(match: "re.Match") -> str:
        label = match.group(1)
        return match.group(0).replace(label, translation.get(label, label), 1)

    return re.sub(r"(?<=[(,])\s*([^(),:;\[\]]+)", substitute, newick)


def _unquote(label: str) -> str:
    label = label.strip()
    if len(label) > 1 and label[0] == label[-1] and label[0] in "'\"":
        return label[1:-1].replace("''", "'")
    return label


def _tip_labels(newick: str) -> List[str]:
    """Tip labels of a newick string.

    A label is a tip when it is preceded by ``(`` or ``,`` -- an internal node
    label follows the closing ``)`` of its clade instead.
    """
    labels: List[str] = []
    index = 0
    length = len(newick)
    previous = ""
    while index < length:
        char = newick[index]
        if char == "'":
            end = index + 1
            while end < length:
                if newick[end] == "'":
                    if end + 1 < length and newick[end + 1] == "'":
                        end += 2
                        continue
                    break
                end += 1
            if previous in "(,":
                labels.append(newick[index + 1 : end].replace("''", "'"))
            index = end + 1
            previous = "label"
            continue
        if char == "[":
            index = newick.find("]", index)
            index = length if index == -1 else index + 1
            continue
        if char in "(),:;":
            previous = char
            index += 1
            continue
        if char.isspace():
            index += 1
            continue
        end = index
        while end < length and newick[end] not in "(),:;[":
            end += 1
        if previous in "(,":
            labels.append(newick[index:end].strip())
        index = end
        previous = "label"
    if not labels:
        raise ValueError("could not find any tip label in the tree")
    return labels
