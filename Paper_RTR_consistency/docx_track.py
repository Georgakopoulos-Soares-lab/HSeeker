#!/usr/bin/env python3
"""Minimal Word tracked-change toolkit for the HSeeker revision.

Operates directly on word/document.xml. Every edit is emitted as a real
<w:ins>/<w:del> revision attributed to one author, so Word offers Accept and
Reject and leaves clean text behind. No character colour is applied: the colour
of a tracked change belongs to the reader's Word settings, not to the file.

Paragraph indices are positions in the full <w:p> sequence of the document,
counting empty paragraphs, so they match the numbering used in REVISION_PLAN.md.
"""
from __future__ import annotations

import copy
import re
import shutil
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
w = lambda t: f"{{{W}}}{t}"  # noqa: E731

AUTHOR = "HSeeker revision"
DATE = "2026-10-15T00:00:00Z"
# Inserted text carries this colour explicitly, so it stays blue after Accept All
# and not only while the revision marks are displayed. To return the accepted
# manuscript to ordinary black, select all and set the font colour to Automatic.
INS_COLOUR = "4A90D9"

PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
for prefix, uri in (("w", W), ("r", R), ("a", A), ("pic", PIC), ("wp", WP), ("", PKG)):
    ET.register_namespace(prefix, uri)


class Doc:
    def __init__(self, path: Path):
        self.src = Path(path)
        self.zin = zipfile.ZipFile(self.src)
        self.tree = ET.fromstring(self.zin.read("word/document.xml"))
        self.body = self.tree.find(w("body"))
        self._id = 9000
        self.log: list[str] = []
        # Snapshot the paragraph elements once. All public methods address
        # paragraphs by their index in the SUBMITTED document, so inserting or
        # deleting never renumbers a later edit.
        self._orig = list(self.tree.iter(w("p")))
        self.rels = ET.fromstring(self.zin.read("word/_rels/document.xml.rels"))
        self.new_media: dict[str, bytes] = {}

    # ---- helpers ---------------------------------------------------------
    def paras(self) -> list[ET.Element]:
        return list(self.tree.iter(w("p")))

    def para(self, i: int) -> ET.Element:
        return self._orig[i]

    def text(self, i: int) -> str:
        return "".join(t.text or "" for t in self.para(i).iter(w("t")))

    def nid(self) -> str:
        self._id += 1
        return str(self._id)

    def _rev_attrs(self) -> dict:
        return {w("id"): self.nid(), w("author"): AUTHOR, w("date"): DATE}

    @staticmethod
    def _parent_map(root: ET.Element) -> dict:
        return {c: p for p in root.iter() for c in p}

    # ---- run surgery -----------------------------------------------------
    @staticmethod
    def _runs_with_text(p: ET.Element) -> list[ET.Element]:
        """Runs carrying a w:t, in document order, descending into hyperlinks.

        Citation superscripts in this manuscript sit inside w:hyperlink, so a
        direct-children-only walk would compute wrong text offsets.
        """
        out = []
        for child in p:
            if child.tag == w("r") and child.find(w("t")) is not None:
                out.append(child)
            elif child.tag == w("hyperlink"):
                out.extend(r for r in child
                           if r.tag == w("r") and r.find(w("t")) is not None)
        return out

    def _holder(self, p: ET.Element, run: ET.Element) -> ET.Element:
        """The element a run actually hangs from (the paragraph, or a hyperlink)."""
        for child in p:
            if child is run:
                return p
            if child.tag == w("hyperlink") and any(r is run for r in child):
                return child
        raise KeyError("run not found in paragraph")

    def _anchor(self, p: ET.Element, run: ET.Element) -> int:
        """Index within p of the top-level child that contains this run."""
        holder = self._holder(p, run)
        return list(p).index(run if holder is p else holder)

    def _make_run(self, template: ET.Element | None, text: str) -> ET.Element:
        r = ET.Element(w("r"))
        if template is not None:
            rpr = template.find(w("rPr"))
            if rpr is not None:
                r.append(copy.deepcopy(rpr))
        t = ET.SubElement(r, w("t"))
        t.text = text
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        return r

    def _to_deleted(self, run: ET.Element) -> ET.Element:
        """Copy a run with w:t rewritten as w:delText."""
        d = copy.deepcopy(run)
        for t in d.findall(w("t")):
            t.tag = w("delText")
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        return d

    @staticmethod
    def _colourise(run: ET.Element, hexval: str) -> None:
        """Force a character colour on a run, keeping rPr first as the schema requires."""
        rpr = run.find(w("rPr"))
        if rpr is None:
            rpr = ET.Element(w("rPr"))
            run.insert(0, rpr)
        for old in rpr.findall(w("color")):
            rpr.remove(old)
        col = ET.Element(w("color"))
        col.set(w("val"), hexval)
        rpr.insert(0, col)

    def _wrap_ins(self, runs: list[ET.Element]) -> ET.Element:
        ins = ET.Element(w("ins"), self._rev_attrs())
        for r in runs:
            self._colourise(r, INS_COLOUR)
            ins.append(r)
        return ins

    def _wrap_del(self, runs: list[ET.Element]) -> ET.Element:
        dl = ET.Element(w("del"), self._rev_attrs())
        for r in runs:
            dl.append(self._to_deleted(r))
        return dl

    # ---- public edits ----------------------------------------------------
    def replace_text(self, i: int, old: str, new: str, label: str = "") -> None:
        """Replace one occurrence of `old` inside paragraph i, as a tracked change.

        Handles the substring spanning several runs, and preserves the character
        formatting of the run where the match starts.
        """
        p = self.para(i)
        runs = self._runs_with_text(p)
        spans, pos = [], 0
        for r in runs:
            txt = "".join(t.text or "" for t in r.findall(w("t")))
            spans.append((pos, pos + len(txt), r, txt))
            pos += len(txt)
        full = "".join(s[3] for s in spans)
        at = full.find(old)
        if at < 0:
            raise KeyError(f"p{i}: {old[:60]!r} not found")
        end = at + len(old)

        touched = [s for s in spans if s[0] < end and s[1] > at]
        first = touched[0]
        holder = self._holder(p, first[2])
        if any(self._holder(p, s[2]) is not holder for s in touched):
            raise ValueError(f"p{i}: match spans a hyperlink boundary; "
                             f"split the edit: {old[:40]!r}")
        anchor_idx = list(holder).index(first[2])

        head = first[3][: at - first[0]]
        tail = touched[-1][3][end - touched[-1][0]:]

        rebuilt: list[ET.Element] = []
        if head:
            rebuilt.append(self._make_run(first[2], head))
        mid = [self._make_run(s[2], s[3][max(at, s[0]) - s[0]: min(end, s[1]) - s[0]])
               for s in touched]
        rebuilt.append(self._wrap_del(mid))
        if new:
            rebuilt.append(self._wrap_ins([self._make_run(first[2], new)]))
        if tail:
            rebuilt.append(self._make_run(touched[-1][2], tail))

        for s in touched:
            holder.remove(s[2])
        for off, el in enumerate(rebuilt):
            holder.insert(anchor_idx + off, el)
        self.log.append(f"p{i} replace{' · ' + label if label else ''}: "
                        f"{old[:48]!r} -> {new[:48]!r}")

    def replace_paragraph(self, i: int, new: str, label: str = "") -> None:
        """Delete every run in paragraph i and insert `new` in its place."""
        p = self.para(i)
        runs = self._runs_with_text(p)
        if not runs:
            raise KeyError(f"p{i}: no text runs")
        tmpl = runs[0]
        anchor = self._anchor(p, runs[0])
        for r in runs:
            self._holder(p, r).remove(r)
        for h in [c for c in p if c.tag == w("hyperlink") and len(c) == 0]:
            p.remove(h)
        p.insert(anchor, self._wrap_del(runs))
        p.insert(anchor + 1, self._wrap_ins([self._make_run(tmpl, new)]))
        self.log.append(f"p{i} replace whole{' · ' + label if label else ''}")

    def delete_paragraph(self, i: int, label: str = "") -> None:
        """Mark every run deleted and mark the paragraph mark itself deleted."""
        p = self.para(i)
        runs = self._runs_with_text(p)
        if runs:
            anchor = self._anchor(p, runs[0])
            for r in runs:
                self._holder(p, r).remove(r)
            for h in [c for c in p if c.tag == w("hyperlink") and len(c) == 0]:
                p.remove(h)
            p.insert(anchor, self._wrap_del(runs))
        ppr = p.find(w("pPr"))
        if ppr is None:
            ppr = ET.Element(w("pPr"))
            p.insert(0, ppr)
        rpr = ppr.find(w("rPr"))
        if rpr is None:
            rpr = ET.SubElement(ppr, w("rPr"))
        rpr.insert(0, ET.Element(w("del"), self._rev_attrs()))
        self.log.append(f"p{i} delete{' · ' + label if label else ''}")

    def insert_paragraph_after(self, i: int, text: str, style_from: int | None = None,
                               label: str = "") -> None:
        """Insert a new paragraph after i, whole-paragraph tracked insertion."""
        p = self.para(i)
        parent = self._parent_map(self.tree)[p]
        at = list(parent).index(p)

        np = ET.Element(w("p"))
        src_ppr = self.para(style_from if style_from is not None else i).find(w("pPr"))
        ppr = copy.deepcopy(src_ppr) if src_ppr is not None else ET.Element(w("pPr"))
        for stale in ppr.findall(w("rPr")):
            ppr.remove(stale)
        rpr = ET.SubElement(ppr, w("rPr"))
        rpr.append(ET.Element(w("ins"), self._rev_attrs()))  # paragraph mark inserted
        np.append(ppr)

        tmpl_runs = self._runs_with_text(
            self.para(style_from if style_from is not None else i))
        tmpl = tmpl_runs[0] if tmpl_runs else None
        np.append(self._wrap_ins([self._make_run(tmpl, text)]))
        parent.insert(at + 1, np)
        self.log.append(f"after p{i} insert{' · ' + label if label else ''}: {text[:52]!r}")

    # ---- images ----------------------------------------------------------
    def swap_image(self, i: int, png: Path, label: str = "", n: int = 1) -> None:
        """Point the nth image in paragraph i at new bytes, keeping its size."""
        p = self.para(i)
        blips = list(p.iter(f"{{{A}}}blip"))
        if len(blips) < n:
            raise KeyError(f"p{i}: no image {n} (has {len(blips)})")
        blip = blips[n - 1]
        rid = blip.get(f"{{{R}}}embed")
        target = None
        for rel in self.rels:
            if rel.get("Id") == rid:
                target = rel.get("Target")
        if target is None:
            raise KeyError(f"p{i}: relationship {rid} missing")
        self.new_media["word/" + target] = Path(png).read_bytes()
        self.log.append(f"p{i} image <- {Path(png).name}{' · ' + label if label else ''}")

    def insert_image_after(self, i: int, png: Path, template: int,
                           label: str = "") -> None:
        """Insert a new image paragraph after i, copying the frame of `template`.

        Keeps the template's display width and rescales the height to the new
        image's aspect ratio, so the figure does not come in distorted.
        """
        import struct

        src = self.para(template)
        drawing = next((c for c in src.iter(w("drawing"))), None)
        if drawing is None:
            raise KeyError(f"p{template}: no drawing to copy")

        data = Path(png).read_bytes()
        if data[:8] != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"{png} is not a PNG")
        iw, ih = struct.unpack(">II", data[16:24])

        # new media part + relationship
        n = 1
        while f"word/media/hseeker_fig{n}.png" in self.new_media:
            n += 1
        target = f"media/hseeker_fig{n}.png"
        self.new_media["word/" + target] = data
        rid = f"rIdHS{n}"
        rel = ET.SubElement(self.rels, "{http://schemas.openxmlformats.org/package/2006/relationships}Relationship")
        rel.set("Id", rid)
        rel.set("Type", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image")
        rel.set("Target", target)

        np = ET.Element(w("p"))
        src_ppr = src.find(w("pPr"))
        ppr = copy.deepcopy(src_ppr) if src_ppr is not None else ET.Element(w("pPr"))
        for stale in ppr.findall(w("rPr")):
            ppr.remove(stale)
        rpr = ET.SubElement(ppr, w("rPr"))
        rpr.append(ET.Element(w("ins"), self._rev_attrs()))
        np.append(ppr)

        run = ET.Element(w("r"))
        dr = copy.deepcopy(drawing)
        for blip in dr.iter(f"{{{A}}}blip"):
            blip.set(f"{{{R}}}embed", rid)
        for ext in dr.iter(f"{{{WP}}}extent"):
            cx = int(ext.get("cx"))
            ext.set("cy", str(int(round(cx * ih / iw))))
        for ext in dr.iter(f"{{{A}}}ext"):
            if ext.get("cx") and ext.get("cy"):
                cx = int(ext.get("cx"))
                ext.set("cy", str(int(round(cx * ih / iw))))
        for dp in dr.iter(f"{{{WP}}}docPr"):
            dp.set("id", str(9500 + n))
            dp.set("name", f"HSeeker figure {n}")
        run.append(dr)
        np.append(self._wrap_ins([run]))

        parent = self._parent_map(self.tree)[self.para(i)]
        parent.insert(list(parent).index(self.para(i)) + 1, np)
        self.log.append(f"after p{i} insert image {Path(png).name}"
                        f"{' · ' + label if label else ''}")

    def insert_table_after(self, i: int, header: list[str], rows: list[list[str]],
                           label: str = "") -> None:
        """Insert a table after paragraph i as a tracked insertion."""
        def cell(text: str, bold: bool) -> ET.Element:
            tc = ET.Element(w("tc"))
            tcpr = ET.SubElement(tc, w("tcPr"))
            ET.SubElement(tcpr, w("tcW")).set(w("type"), "auto")
            p = ET.SubElement(tc, w("p"))
            ppr = ET.SubElement(p, w("pPr"))
            rpr0 = ET.SubElement(ppr, w("rPr"))
            rpr0.append(ET.Element(w("ins"), self._rev_attrs()))
            r = ET.Element(w("r"))
            if bold:
                rpr = ET.SubElement(r, w("rPr"))
                ET.SubElement(rpr, w("b"))
            t = ET.SubElement(r, w("t"))
            t.text = text
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
            p.append(self._wrap_ins([r]))
            return tc

        tbl = ET.Element(w("tbl"))
        pr = ET.SubElement(tbl, w("tblPr"))
        ET.SubElement(pr, w("tblStyle")).set(w("val"), "TableGrid")
        ET.SubElement(pr, w("tblW")).set(w("type"), "auto")
        borders = ET.SubElement(pr, w("tblBorders"))
        for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
            e = ET.SubElement(borders, w(edge))
            e.set(w("val"), "single"); e.set(w("sz"), "4")
            e.set(w("space"), "0"); e.set(w("color"), "999999")
        for ri, cells in enumerate([header] + rows):
            tr = ET.SubElement(tbl, w("tr"))
            trpr = ET.SubElement(tr, w("trPr"))
            trpr.append(ET.Element(w("ins"), self._rev_attrs()))  # row inserted
            for text in cells:
                tr.append(cell(text, bold=(ri == 0)))

        anchor = self.para(i)
        parent = self._parent_map(self.tree)[anchor]
        at = list(parent).index(anchor) + 1
        parent.insert(at, tbl)
        # Word requires a paragraph after a table
        tail = ET.Element(w("p"))
        ppr = ET.SubElement(tail, w("pPr"))
        rpr = ET.SubElement(ppr, w("rPr"))
        rpr.append(ET.Element(w("ins"), self._rev_attrs()))
        parent.insert(at + 1, tail)
        self.log.append(f"after p{i} insert table {len(rows)}x{len(header)}"
                        f"{' · ' + label if label else ''}")

    # ---- save ------------------------------------------------------------
    def save(self, out: Path) -> None:
        out = Path(out)
        xml = ET.tostring(self.tree, encoding="UTF-8", xml_declaration=True)
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zo:
            for item in self.zin.infolist():
                if item.filename == "word/document.xml":
                    zo.writestr(item, xml)
                elif item.filename == "word/_rels/document.xml.rels":
                    zo.writestr(item, ET.tostring(self.rels, encoding="UTF-8",
                                                  xml_declaration=True))
                elif item.filename in self.new_media:
                    zo.writestr(item, self.new_media[item.filename])
                else:
                    zo.writestr(item, self.zin.read(item.filename))
            for name, blob in self.new_media.items():
                if name not in self.zin.namelist():
                    zo.writestr(name, blob)
        print(f"wrote {out}  ({len(self.log)} tracked edits)")
