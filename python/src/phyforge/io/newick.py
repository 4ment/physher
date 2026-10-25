"""A newick parser that keeps the tree, for the analyses that need its shape.

:mod:`phyforge.io.trees` only ever needs the tip labels, so it scans for them;
root-to-tip regression needs the topology and the branch lengths as well.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, List, Optional

__all__ = ["Node", "parse_newick"]


@dataclass(eq=False)
class Node:
    name: str = ""
    #: length of the branch above this node; None when the newick gives none
    length: Optional[float] = None
    children: List["Node"] = field(default_factory=list)
    parent: Optional["Node"] = None

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def preorder(self) -> Iterator["Node"]:
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def leaves(self) -> List["Node"]:
        return [node for node in self.preorder() if node.is_leaf]

    def neighbours(self) -> List["Node"]:
        """Adjacent nodes, the tree read as unrooted."""
        return self.children + ([self.parent] if self.parent is not None else [])


def parse_newick(text: str) -> Node:
    """Parse the first tree in ``text``; comments in ``[...]`` are skipped."""
    parser = _Parser(text)
    root = parser.subtree()
    parser.skip()
    if parser.peek() not in (";", ""):
        raise ValueError(f"unexpected {parser.peek()!r} at position {parser.index}")
    return root


class _Parser:
    def __init__(self, text: str):
        self.text = text
        self.index = 0

    def peek(self) -> str:
        return self.text[self.index] if self.index < len(self.text) else ""

    def skip(self) -> None:
        """Skip whitespace and bracketed comments."""
        while self.index < len(self.text):
            char = self.text[self.index]
            if char.isspace():
                self.index += 1
            elif char == "[":
                end = self.text.find("]", self.index)
                self.index = len(self.text) if end == -1 else end + 1
            else:
                break

    def subtree(self) -> Node:
        self.skip()
        node = Node()
        if self.peek() == "(":
            self.index += 1
            while True:
                child = self.subtree()
                child.parent = node
                node.children.append(child)
                self.skip()
                char = self.peek()
                self.index += 1
                if char == ")":
                    break
                if char != ",":
                    raise ValueError(
                        f"expected ',' or ')' at position {self.index - 1}"
                    )
        node.name = self.label()
        self.skip()
        if self.peek() == ":":
            self.index += 1
            self.skip()
            node.length = float(self.token())
        return node

    def label(self) -> str:
        self.skip()
        if self.peek() == "'":
            end = self.index + 1
            parts = []
            while end < len(self.text):
                if self.text[end] == "'":
                    if self.text[end + 1 : end + 2] == "'":
                        parts.append("'")
                        end += 2
                        continue
                    break
                parts.append(self.text[end])
                end += 1
            self.index = end + 1
            return "".join(parts)
        return self.token().strip()

    def token(self) -> str:
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in "(),:;[":
            self.index += 1
        return self.text[start : self.index]
