"""SF41 — Raw GC content across the AUG and stop-codon windows.

Two-panel supplemental figure showing mean rolling-GC content of transcript
sequences per position, comparing NMD-susceptible vs Control classes, across
the AUG and stop-codon windows the deep-learning model consumes.

The paper's claim (§5, at SF41): GC content differentiates NMD from Control
transcripts in the stop-codon window but not in the AUG window, consistent
with PTCs turning exon sequence (GC-rich) into 3'UTR sequence (GC-poor)
downstream of a premature stop. Post-stop, Control transcripts show the
expected sharp drop in GC content (real 3'UTR); NMD transcripts stay
GC-rich because their "post-stop" sequence is still coding.

Data source: results_4ct/gc_content_across_{atg,stop}_window_atg500_stop500.tsv
(from 09_export_gc_content.py --branch {atg,stop}). 50-bp sliding window
with 10-bp step. Shading is ±1 standard error of the mean.

Style: ggplot-mimic theme (grey panel + white gridlines) to match SF1-SF23.
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LIB  = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    facet_header,
    assert_text_within_canvas,
    NMD_COLOR,
    CONTROL_COLOR,
    TITLE_C,
    BODY_FS,
    HEADER_FS,
    panel_letter,
)

apply_ggplot_rcparams()

DATA = HERE / "data"


def load(branch):
    # ref-AUG-only cohort (see 09d_export_gc_content_refaug_only.py in the
    # model repo): NMD isoforms restricted to is_ref_cds=1 on ORF0; Controls
    # unchanged. Drops the ~23% of NMD test transcripts that fall through to
    # the TD2 fallback, whose "stop" would be a downstream non-PTC codon.
    df = pd.read_csv(
        DATA / f"gc_content_across_{branch}_window_refaug_only_atg500_stop500.tsv",
        sep="\t",
    )
    return df.sort_values(["class", "rel_mid"]).reset_index(drop=True)


def draw_panel(ax, df, *, codon_label):
    style_axes_ggplot(ax)
    ax.axvline(0, color="#555555", linewidth=0.8, linestyle=":", zorder=2)
    for cls, color, label in [
        ("NMD",     NMD_COLOR,     "NMD susceptible"),
        ("Control", CONTROL_COLOR, "Control"),
    ]:
        sub = df[df["class"] == cls]
        x   = sub["rel_mid"].values
        y   = sub["mean_gc"].values
        se  = sub["se_gc"].values
        ax.fill_between(x, y - se, y + se, color=color, alpha=0.15, zorder=3)
        ax.plot(x, y, color=color, linewidth=1.8, label=label, zorder=4)
    ax.set_xlabel(f"Position relative to {codon_label} (nt)", fontsize=BODY_FS)
    ax.set_ylabel("Mean GC content (50-nt window)", fontsize=BODY_FS)
    ax.set_xlim(df["rel_mid"].min(), df["rel_mid"].max())
    ax.set_ylim(0.35, 0.60)
    return ax


def main():
    atg  = load("atg")
    stop = load("stop")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.82, bottom=0.12, wspace=0.06)

    draw_panel(axes[0], atg,  codon_label="start codon")
    draw_panel(axes[1], stop, codon_label="stop codon")
    axes[1].set_ylabel("")

    facet_header(axes[0], "Start (AUG) window", height=0.08)
    facet_header(axes[1], "Stop codon window",   height=0.08)

    panel_letter(axes[0], "A", x=-0.07, y=1.14)
    panel_letter(axes[1], "B", x=-0.07, y=1.14)

    axes[1].legend(loc="upper right", frameon=True, framealpha=0.9,
                    edgecolor="none", facecolor="white", fontsize=BODY_FS)

    assert_text_within_canvas(fig)

    out_png = HERE / "figure_sf41_gc_content_stop_window.png"
    out_pdf = HERE / "figure_sf41_gc_content_stop_window.pdf"
    fig.savefig(out_png, dpi=200, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    print(f"wrote {out_png.name} and {out_pdf.name}")
    print(f"  ATG rows: {len(atg)}   STOP rows: {len(stop)}")


if __name__ == "__main__":
    main()
