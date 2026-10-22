#!/usr/bin/env python3
"""Assemble the revised Figure 1.

Panel A is the redrawn candidate-identification flowchart. Panel B is the
submitted biological schematic with three edits inside its middle row, made on
the raster because no vector source for that artwork exists in the repository:

  1. the "Sensitivity Analysis / AUC" inset is removed, that analysis having
     been withdrawn;
  2. the detection tree it sat above is moved up into the space, so the panel
     loses a crowded corner rather than gaining a hole;
  3. "H-DNA Stability" is relabelled, so the figure does not reintroduce the
     free-energy reading of the score.

If the original BioRender file turns up, prefer editing that and re-exporting.
"""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
PLOTS = ROOT / "analysis" / "results" / "plots"
SUBMITTED = Path("/tmp/claude-903286/-scratch-10899-kimopro-HSeeker/"
                 "2ace8cff-952f-44f6-83eb-34d39660be40/scratchpad/Figure1.png")

# Regions in the submitted 1591 x 2048 render of Figure 1.
PANEL_B_TOP = 622          # first row of the biological schematic
# The inset, the checklist icon beside it and the arrow that pointed at it all
# belong to the withdrawn analysis and go together.
CLEAR_BOXES = [(82, 1128, 522, 1428),    # AUC chart, axis labels, checklist icon
               (478, 1200, 578, 1276)]   # arrow that pointed at the inset
AUC_BOX = (82, 1128, 522, 1428)
TREE_BOX = (45, 1432, 500, 1576)     # detection tree, moved up
LABEL_BOX = (608, 1134, 900, 1200)   # "H-DNA Stability" dashed label

FONTS = ("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
         "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")


def font(size: int) -> ImageFont.FreeTypeFont:
    for f in FONTS:
        if Path(f).exists():
            return ImageFont.truetype(f, size)
    return ImageFont.load_default()


def dashed_rect(d: ImageDraw.ImageDraw, box, colour, width=3, dash=11, gap=8):
    x0, y0, x1, y1 = box
    for x in range(x0, x1, dash + gap):
        d.line([(x, y0), (min(x + dash, x1), y0)], fill=colour, width=width)
        d.line([(x, y1), (min(x + dash, x1), y1)], fill=colour, width=width)
    for y in range(y0, y1, dash + gap):
        d.line([(x0, y), (x0, min(y + dash, y1))], fill=colour, width=width)
        d.line([(x1, y), (x1, min(y + dash, y1))], fill=colour, width=width)


def build_panel_b() -> Image.Image:
    src = Image.open(SUBMITTED).convert("RGB")
    b = src.crop((0, PANEL_B_TOP, src.width, src.height))
    off = PANEL_B_TOP

    tree = b.crop((TREE_BOX[0], TREE_BOX[1] - off, TREE_BOX[2], TREE_BOX[3] - off))

    d = ImageDraw.Draw(b)
    for bx in CLEAR_BOXES:
        d.rectangle((bx[0], bx[1] - off, bx[2], bx[3] - off), fill="white")
    d.rectangle((TREE_BOX[0] - 8, TREE_BOX[1] - off - 4, TREE_BOX[2] + 40,
                 TREE_BOX[3] - off + 26), fill="white")
    # put the tree back, vertically centred in the space the inset vacated
    new_top = (AUC_BOX[1] - off) + ((AUC_BOX[3] - AUC_BOX[1]) +
                                    (TREE_BOX[3] - TREE_BOX[1])) // 2 - tree.height // 2
    b.paste(tree, (TREE_BOX[0], new_top))

    # relabel the stability box
    d.rectangle((LABEL_BOX[0], LABEL_BOX[1] - off, LABEL_BOX[2] + 60,
                 LABEL_BOX[3] - off), fill="white")
    lx0, ly0 = LABEL_BOX[0] - 6, LABEL_BOX[1] - off - 2
    lx1, ly1 = LABEL_BOX[2] + 58, LABEL_BOX[3] - off - 2
    dashed_rect(d, (lx0, ly0, lx1, ly1), (110, 118, 180))
    f = font(26)
    text = "Stability score\n(heuristic ranking)"
    tb = d.multiline_textbbox((0, 0), text, font=f, align="center", spacing=6)
    d.multiline_text(((lx0 + lx1) // 2 - (tb[2] - tb[0]) // 2,
                      (ly0 + ly1) // 2 - (tb[3] - tb[1]) // 2 - 4),
                     text, font=f, fill="black", align="center", spacing=6)
    return b


def main() -> None:
    panel_a = Image.open(PLOTS / "figure1_panelA_workflow.png").convert("RGB")
    panel_b = build_panel_b()

    width = max(panel_a.width, panel_b.width)
    if panel_a.width != width:
        panel_a = panel_a.resize(
            (width, round(panel_a.height * width / panel_a.width)), Image.LANCZOS)
    if panel_b.width != width:
        panel_b = panel_b.resize(
            (width, round(panel_b.height * width / panel_b.width)), Image.LANCZOS)

    gap = 26
    out = Image.new("RGB", (width, panel_a.height + gap + panel_b.height), "white")
    out.paste(panel_a, (0, 0))
    out.paste(panel_b, (0, panel_a.height + gap))

    d = ImageDraw.Draw(out)
    d.text((14, panel_a.height + gap + 6), "B", font=font(46), fill="black")

    out.save(PLOTS / "figure1_workflow.png")
    out.save(PLOTS / "figure1_workflow.pdf", "PDF", resolution=300)
    print(f"wrote {PLOTS / 'figure1_workflow.png'}  {out.size[0]}x{out.size[1]}")


if __name__ == "__main__":
    main()
