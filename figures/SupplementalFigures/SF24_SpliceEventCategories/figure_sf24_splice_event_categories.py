"""SF24 — Isopair splice event category schematic (12 events).

Twelve sub-panels arranged 3 rows × 4 cols. Each sub-panel is a two-track
exon schematic: Reference isoform (top) vs Comparator isoform (bottom),
using the same coordinate examples as the canonical Isopair glossary
(`results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_mashr.Rmd`,
§ Event Type Glossary).

The 12 events, grouped as in the Isopair Rmd:
  - Transcript-end variants:   Alt_TSS, Alt_TES
  - Alternative splice sites:  A5SS, A3SS
  - Exon skip / miss:          SE, Missing_Internal
  - Intron retention (partial + full + boundary variants):
                               Partial_IR_5, Partial_IR_3,
                               IR, IR_diff_5, IR_diff_3, IR_diff_5_3

The intent is diagrammatic — not to render actual isoforms — so coordinates
match the glossary's illustrative examples exactly, ensuring the paper's
schematic is consistent with the source-of-truth definitions.

Style: matplotlib rendered with ggplot-mimic theme for structural consistency
with SF1-SF23 (grey panel bg, white gridlines), but with axes hidden — this
is a category diagram, not a data plot.
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

HERE = Path(__file__).resolve().parent
LIB  = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    assert_text_within_canvas,
    panel_letter,
    NMD_COLOR,
    CONTROL_COLOR,
    REF_COLOR,
    STRIP_BG,
    TITLE_C,
    AXIS_C,
    PANEL_BG,
    BODY_FS,
    HEADER_FS,
)

apply_ggplot_rcparams()

# Reference (top track) uses REF_COLOR (grey); comparator (bottom track) uses
# a distinct accent so the reader can immediately tell which is which. The
# accent is NMD_COLOR (orange) — matching the Isopair pair-analysis
# convention of NMD/Susceptible = orange used throughout SF25-SF27, §4
# figures.
COMP_COLOR = NMD_COLOR
EXON_H = 0.28   # exon block height in axes-y units
Y_REF  = 0.68   # reference row center
Y_COMP = 0.28   # comparator row center

# Coordinate examples from the Rmd glossary (§ Event Type Glossary chunks).
# Each entry: (ref_start_list, ref_end_list, comp_start_list, comp_end_list,
#              subtitle_note). Coordinates are illustrative genomic positions
#              on the + strand (5' → 3' left → right).
EVENTS = [
    # ── Transcript-end variants ──
    dict(
        code="Alt_TSS", title="Alt TSS",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[250, 400, 700], comp_e=[350, 550, 850],
        note="Comparator starts downstream",
    ),
    dict(
        code="Alt_TES", title="Alt TES",
        ref_s=[100, 400, 700], ref_e=[200, 550, 900],
        comp_s=[100, 400, 700], comp_e=[200, 550, 800],
        note="Comparator ends earlier",
    ),
    # ── Alternative splice sites ──
    dict(
        code="A5SS", title="A5SS",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 400, 700], comp_e=[200, 600, 850],
        note="Same acceptor, different donor",
    ),
    dict(
        code="A3SS", title="A3SS",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 350, 700], comp_e=[200, 550, 850],
        note="Same donor, different acceptor",
    ),
    # ── Exon skip / miss ──
    dict(
        code="SE", title="Skipped Exon (SE)",
        ref_s=[100, 400, 600, 800], ref_e=[200, 500, 700, 900],
        comp_s=[100, 400, 800],     comp_e=[200, 500, 900],
        note="Comparator skips one internal exon",
    ),
    dict(
        code="Missing_Internal", title="Missing Internal",
        ref_s=[100, 350, 550, 750, 950], ref_e=[200, 450, 650, 850, 1050],
        comp_s=[100, 950],               comp_e=[200, 1050],
        note="Comparator skips multiple internal exons",
    ),
    # ── Partial intron retention ──
    dict(
        code="Partial_IR_5", title="Partial IR (5′)",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 400, 700], comp_e=[200, 620, 850],
        note="Exon extends into intron from 5′ side",
    ),
    dict(
        code="Partial_IR_3", title="Partial IR (3′)",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 400, 640], comp_e=[200, 550, 850],
        note="Exon extends into intron from 3′ side",
    ),
    # ── Full intron retention (aggregate row) ──
    dict(
        code="IR", title="Full IR",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 400],       comp_e=[200, 850],
        note="Comparator retains whole intron",
    ),
    dict(
        code="IR_diff_5", title="IR diff 5′",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 430],       comp_e=[200, 850],
        note="Retaining exon starts after ref exon",
    ),
    dict(
        code="IR_diff_3", title="IR diff 3′",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 400],       comp_e=[200, 820],
        note="Retaining exon ends before ref exon",
    ),
    dict(
        code="IR_diff_5_3", title="IR diff 5′ + 3′",
        ref_s=[100, 400, 700], ref_e=[200, 550, 850],
        comp_s=[100, 430],       comp_e=[200, 820],
        note="Retaining exon falls short at both ends",
    ),
]


def draw_track(ax, y_center, starts, ends, color, x_min, x_max, label_text):
    """Draw a single isoform exon track: connecting intron line + exon blocks."""
    # Intron line from first exon start to last exon end
    ax.plot([starts[0], ends[-1]], [y_center, y_center],
            color=AXIS_C, linewidth=0.9, zorder=2)
    # Exon rectangles
    for s, e in zip(starts, ends):
        ax.add_patch(mpatches.Rectangle(
            (s, y_center - EXON_H / 2), e - s, EXON_H,
            facecolor=color, edgecolor=TITLE_C, linewidth=0.7,
            zorder=3,
        ))
    # Row label to the left of the track
    ax.text(x_min - (x_max - x_min) * 0.02, y_center, label_text,
            ha="right", va="center", fontsize=BODY_FS - 1,
            color=TITLE_C, fontweight="normal")


def draw_event_panel(ax, event):
    """Draw one 2-track schematic in `ax`. No numeric axes."""
    ax.set_facecolor(PANEL_BG)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])

    x_all = event["ref_s"] + event["ref_e"] + event["comp_s"] + event["comp_e"]
    x_min = min(x_all)
    x_max = max(x_all)
    pad_l = (x_max - x_min) * 0.22  # room for "Reference" / "Comparator" labels
    pad_r = (x_max - x_min) * 0.05

    ax.set_xlim(x_min - pad_l, x_max + pad_r)
    ax.set_ylim(0.0, 1.05)

    draw_track(ax, Y_REF,  event["ref_s"],  event["ref_e"],
               REF_COLOR,  x_min, x_max, "Reference")
    draw_track(ax, Y_COMP, event["comp_s"], event["comp_e"],
               COMP_COLOR, x_min, x_max, "Comparator")

    # Facet header (event title, grey strip) — placed inside the panel to
    # avoid needing extra headroom above every one of 12 sub-panels
    ax.text(
        0.5, 0.98, event["title"],
        transform=ax.transAxes,
        ha="center", va="top",
        fontsize=HEADER_FS - 2, fontweight="bold",
        color=TITLE_C,
        bbox=dict(facecolor=STRIP_BG, edgecolor="none", boxstyle="square,pad=0.4"),
        zorder=5,
    )
    # Subtitle note below the tracks
    ax.text(
        0.5, 0.02, event["note"],
        transform=ax.transAxes,
        ha="center", va="bottom",
        fontsize=BODY_FS - 2, color=TITLE_C, style="italic",
        zorder=5,
    )


def build_figure():
    fig, axes = plt.subplots(3, 4, figsize=(15.0, 8.4))
    fig.subplots_adjust(left=0.03, right=0.99, top=0.93, bottom=0.05,
                        wspace=0.18, hspace=0.35)
    for ax, event in zip(axes.flat, EVENTS):
        draw_event_panel(ax, event)
    # Section header
    fig.text(
        0.5, 0.975,
        "Twelve splice event categories detected by Isopair",
        ha="center", va="center",
        fontsize=HEADER_FS + 2, fontweight="bold",
        color=TITLE_C,
    )
    # Legend strip below header explaining color mapping
    fig.text(
        0.5, 0.945,
        "Reference (top, grey) = dominant non-NMD isoform · "
        "Comparator (bottom, orange) = paired NMD/Control isoform · "
        "boxes = exons · lines = introns",
        ha="center", va="center",
        fontsize=BODY_FS, color=TITLE_C, style="italic",
    )
    return fig


def main():
    fig = build_figure()
    assert_text_within_canvas(fig)
    out_png = HERE / "figure_sf24_splice_event_categories.png"
    out_pdf = HERE / "figure_sf24_splice_event_categories.pdf"
    fig.savefig(out_png, dpi=200, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    print(f"wrote {out_png.name} and {out_pdf.name}")
    print(f"  n events = {len(EVENTS)}")


if __name__ == "__main__":
    main()
