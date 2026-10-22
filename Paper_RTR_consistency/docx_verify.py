#!/usr/bin/env python3
"""Verify a tracked-change docx by simulating Accept All and Reject All.

Reject All must reproduce the original manuscript exactly; Accept All must give
the intended revised text. Anything else means the revisions are malformed.
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
w = lambda t: f"{{{W}}}{t}"  # noqa: E731


def paragraphs(path: Path, mode: str) -> list[str]:
    """mode: 'accept', 'reject', or 'markup'."""
    root = ET.fromstring(zipfile.ZipFile(path).read("word/document.xml"))
    out = []
    for p in root.iter(w("p")):
        if mode == "accept":
            ppr = p.find(w("pPr"))
            if ppr is not None:
                rpr = ppr.find(w("rPr"))
                if rpr is not None and rpr.find(w("del")) is not None:
                    continue  # paragraph mark deleted -> paragraph goes away
        if mode == "reject":
            ppr = p.find(w("pPr"))
            if ppr is not None:
                rpr = ppr.find(w("rPr"))
                if rpr is not None and rpr.find(w("ins")) is not None:
                    continue  # paragraph mark inserted -> never existed
        buf = []
        for node in p.iter():
            if node.tag == w("t"):
                parent_ins = any(a.tag == w("ins") for a in ancestors(p, node))
                if mode == "reject" and parent_ins:
                    continue
                buf.append(node.text or "")
            elif node.tag == w("delText"):
                if mode == "reject":
                    buf.append(node.text or "")
        out.append("".join(buf))
    return out


def ancestors(root: ET.Element, target: ET.Element) -> list[ET.Element]:
    parent = {c: p for p in root.iter() for c in p}
    chain, node = [], target
    while node in parent:
        node = parent[node]
        chain.append(node)
    return chain


def main() -> None:
    orig, edited = Path(sys.argv[1]), Path(sys.argv[2])
    base = paragraphs(orig, "markup")
    rej = paragraphs(edited, "reject")
    acc = paragraphs(edited, "accept")

    print(f"original paragraphs: {len(base)}")
    print(f"after Reject All:    {len(rej)}")
    print(f"after Accept All:    {len(acc)}")

    bad = 0
    if len(base) != len(rej):
        print("\n!! Reject All changes the paragraph count")
        bad += 1
    else:
        for i, (a, b) in enumerate(zip(base, rej)):
            if a != b:
                bad += 1
                print(f"\n!! p{i} differs after Reject All")
                print(f"   original: {a[:150]!r}")
                print(f"   rejected: {b[:150]!r}")
                if bad > 4:
                    break
    print("\nReject All reproduces the original exactly."
          if bad == 0 else f"\n{bad} mismatch(es) after Reject All.")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
