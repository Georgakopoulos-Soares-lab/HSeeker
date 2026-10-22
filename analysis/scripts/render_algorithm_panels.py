#!/usr/bin/env python3
"""Render the algorithm panels of main.tex as images.

No LaTeX is available on this machine, so this parses the subset of the
algorithmic environment that main.tex actually uses and typesets it directly:
caption, rules, numbered lines, keyword emphasis and right-aligned comments.

Output: analysis/results/plots/algorithm{1..4}.png
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
TEX = ROOT / "main.tex"
OUT = ROOT / "analysis" / "results" / "plots"

SYMBOLS = {
    r"\gets": "\u2190", r"\leftarrow": "\u2190", r"\in": "\u2208",
    r"\ge": "\u2265", r"\le": "\u2264", r"\ell": "l",
    r"\sigma": "\u03c3", r"\tau": "\u03c4", r"\delta": "\u03b4",
    r"\cdot": "\u00b7", r"\bot": "\u22a5", r"\sum": "\u03a3",
    r"\lfloor": "\u230a", r"\rfloor": "\u230b",
    r"\triangleright": "\u25b7", r"\times": "\u00d7", r"\neq": "\u2260",
    r"\infty": "\u221e", r"\approx": "\u2248", r"\pm": "\u00b1",
    r"\max": "max", r"\min": "min", r"\log": "log",
}
STRIP_WRAPPERS = ("mathrm", "mathit", "mathcal", "text", "textit",
                  "textbf", "textsc", "textstyle", "mathbf", "emph")
KEYWORDS = ("if", "then", "else if", "else", "end if", "for", "end for",
            "while", "end while", "do", "return", "procedure", "end procedure",
            "Require:", "Ensure:")


def strip_wrappers(s: str) -> str:
    changed = True
    while changed:
        changed = False
        for cmd in STRIP_WRAPPERS:
            i = s.find("\\" + cmd + "{")
            if i < 0:
                continue
            j, depth = i + len(cmd) + 2, 1
            while j < len(s) and depth:
                depth += (s[j] == "{") - (s[j] == "}")
                j += 1
            inner = s[i + len(cmd) + 2: j - 1]
            if cmd == "textsc":
                inner = inner.capitalize()
            s = s[:i] + inner + s[j:]
            changed = True
    return s


LABELS: dict[str, int] = {}


def _cases(m: re.Match) -> str:
    """Render a LaTeX cases block inline: value if condition; ...; value otherwise."""
    parts = []
    for row in re.split(r"\\\\", m.group(1)):
        row = row.strip()
        if not row:
            continue
        val, _, cond = row.partition("&")
        parts.append(f"{val.strip()} {cond.strip()}" if cond.strip() else val.strip())
    return "; ".join(parts)


def detex(s: str) -> str:
    s = re.sub(r"\\ref\{([^}]*)\}",
               lambda m: str(LABELS.get(m.group(1), "")), s)
    s = re.sub(r"\\label\{[^}]*\}", "", s)
    s = s.replace(r"\{", "\x01").replace(r"\}", "\x02")
    s = re.sub(r"\\begin\{cases\}(.*?)\\end\{cases\}", _cases, s, flags=re.S)
    s = strip_wrappers(s)
    for k, v in sorted(SYMBOLS.items(), key=lambda kv: -len(kv[0])):
        s = s.replace(k, v)
    s = s.replace(r"\bigl", "").replace(r"\bigr", "").replace(r"\,", " ")
    s = s.replace(r"\;", " ").replace(r"\!", "").replace(r"\\", " ")
    s = s.replace("$", "").replace("~", " ")
    s = re.sub(r"_\{([^}]*)\}", r"_\1", s)
    s = re.sub(r"\^\{([^}]*)\}", r"^\1", s)
    s = s.replace("{", "").replace("}", "")
    s = s.replace("``", '"').replace("''", '"')
    s = s.replace(r"\_", "_")
    s = re.sub(r"\\[A-Za-z]+", "", s)
    s = s.replace("\x01", "{").replace("\x02", "}")
    return re.sub(r"\s+", " ", s).strip()


def parse(block: str) -> tuple[str, list[tuple[int, str, str, bool]]]:
    """Return (caption, [(indent, text, comment, numbered)])."""
    cap = re.search(r"\\caption\{(.*?)\}\s*\n", block, re.S)
    caption = detex(cap.group(1)) if cap else ""
    body = block[block.index("\\begin{algorithmic}"):]
    body = body[body.index("\n") + 1: body.rindex("\\end{algorithmic}")]

    lines, indent = [], 0
    cmds = ("State", "If", "ElsIf", "Else", "EndIf", "For", "EndFor", "While",
            "EndWhile", "Return", "Procedure", "EndProcedure", "Statex",
            "Require", "Ensure")
    for raw in re.split(r"(?=\\(?:" + "|".join(cmds) + r")\b)", body):
        raw = raw.strip()
        if not raw:
            continue
        comment = ""
        m = re.search(r"\\Comment\{(.*)\}\s*$", raw, re.S)
        if m:
            comment = detex(m.group(1))
            raw = raw[:m.start()].rstrip()

        def emit(txt, num=True, d=0):
            if not txt.strip():
                return
            budget = 92 - (indent + d) * 3
            if comment:
                budget -= min(len(comment) + 4, 46)
            if len(txt) <= budget:
                lines.append((indent + d, txt, comment, num))
                return
            words, cur, out = txt.split(" "), "", []
            for wd in words:
                if cur and len(cur) + 1 + len(wd) > budget:
                    out.append(cur); cur = wd
                else:
                    cur = f"{cur} {wd}".strip()
            if cur:
                out.append(cur)
            for k, piece in enumerate(out):
                lines.append((indent + d + (0 if k == 0 else 1), piece,
                              comment if k == 0 else "", num if k == 0 else False))

        if raw.startswith("\\Require"):
            emit("Require: " + detex(raw[8:]), False)
        elif raw.startswith("\\Ensure"):
            emit("Ensure: " + detex(raw[7:]), False)
        elif raw.startswith("\\Statex"):
            emit(detex(raw[7:]), False)
        elif raw.startswith("\\State"):
            rest = raw[6:].strip()
            if rest.startswith("\\Return"):
                emit("return " + detex(rest[7:]))
            else:
                emit(detex(rest))
        elif raw.startswith("\\If"):
            emit("if " + detex(raw[3:]) + " then"); indent += 1
        elif raw.startswith("\\ElsIf"):
            indent -= 1; emit("else if " + detex(raw[6:]) + " then"); indent += 1
        elif raw.startswith("\\Else"):
            indent -= 1; emit("else"); indent += 1
        elif raw.startswith("\\EndIf"):
            indent -= 1; emit("end if")
        elif raw.startswith("\\For"):
            emit("for " + detex(raw[4:]) + " do"); indent += 1
        elif raw.startswith("\\EndFor"):
            indent -= 1; emit("end for")
        elif raw.startswith("\\While"):
            emit("while " + detex(raw[6:]) + " do"); indent += 1
        elif raw.startswith("\\EndWhile"):
            indent -= 1; emit("end while")
        elif raw.startswith("\\Procedure"):
            a = re.match(r"\\Procedure\{([^}]*)\}\{(.*)\}", raw, re.S)
            emit(f"procedure {a.group(1)}({detex(a.group(2))})" if a
                 else "procedure"); indent += 1
        elif raw.startswith("\\EndProcedure"):
            indent -= 1; emit("end procedure")
        elif raw.startswith("\\Return"):
            emit("return " + detex(raw[7:]))
        else:
            emit(detex(raw))
    return caption, lines


def render(caption: str, lines, path: Path, number: int) -> None:
    fs, lh = 8.6, 0.235                    # fontsize, line height (inches)
    cw = fs * 0.6022 / 72                  # DejaVu Sans Mono advance, inches
    width = 7.4
    n = len(lines)
    height = 0.70 + n * lh

    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, width); ax.set_ylim(0, height)

    x0, top = 0.12, height - 0.12
    ax.plot([x0, width - 0.12], [top, top], color="black", lw=1.4)
    ax.text(x0, top - 0.17, f"Algorithm {number}: {caption}", fontsize=fs + 0.9,
            fontweight="bold", va="center", ha="left")
    ax.plot([x0, width - 0.12], [top - 0.32, top - 0.32], color="black", lw=0.9)

    y = top - 0.32 - lh
    num = 0
    for ind, text, comment, numbered in lines:
        if numbered:
            num += 1
            ax.text(x0 + 0.26, y, f"{num:>2}:", fontsize=fs - 0.9,
                    family="monospace", color="#666666", va="center", ha="right")
        tx = x0 + 0.38 + ind * 0.22
        low = text.lower()
        kw = next((k for k in KEYWORDS if low.startswith(k)), None) if numbered \
            else next((k for k in ("Require:", "Ensure:") if text.startswith(k)), None)
        if kw:
            ax.text(tx, y, text[:len(kw)], fontsize=fs, family="monospace",
                    fontweight="bold", va="center", ha="left")
            rest = text[len(kw):]
            if rest:
                ax.text(tx + len(kw) * cw, y, rest, fontsize=fs,
                        family="monospace", va="center", ha="left")
        else:
            style = "italic" if not numbered and text.startswith("\u25b7") else "normal"
            ax.text(tx, y, text, fontsize=fs, family="monospace",
                    style=style, va="center", ha="left",
                    color="#333333" if style == "italic" else "black")
        if comment:
            ax.text(width - 0.16, y, "\u25b7 " + comment, fontsize=fs - 1.0,
                    family="monospace", color="#555555", va="center", ha="right")
        y -= lh

    ax.plot([x0, width - 0.12], [y + lh / 2, y + lh / 2], color="black", lw=1.4)
    fig.savefig(path, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print(f"  Algorithm {number}: {len(lines)} lines -> {path.name}")


def main() -> None:
    src = TEX.read_text()
    blocks = re.findall(r"\\begin\{algorithm\}.*?\\end\{algorithm\}", src, re.S)
    for i, b in enumerate(blocks, 1):
        for lab in re.findall(r"\\label\{([^}]*)\}", b):
            LABELS[lab] = i
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(blocks)} algorithm panels in {TEX.name}")
    for i, b in enumerate(blocks, 1):
        caption, lines = parse(b)
        caption = caption.replace("Thermodynamic Stability Scoring",
                                  "Heuristic Stability Scoring")
        render(caption, lines, OUT / f"algorithm{i}.png", i)


if __name__ == "__main__":
    main()
