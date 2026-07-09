"""SF36 — Absolute Shapley values across the AUG and stop-codon windows.

Two-panel supplemental figure showing per-position |SHAP| aggregated across
all input channels, comparing NMD-susceptible vs Control transcript classes,
across the two 500-nt sequence windows the DL model consumes.

Data source: results_4ct/shap_profile_{atg,stop}_joint_atg500_stop500.tsv
(from 09b_export_subgroup_profiles.py, 5-run pooled DeepSHAP).

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


def load_window(branch):
    df = pd.read_csv(DATA / f"shap_profile_{branch}_joint_atg500_stop500.tsv", sep="\t")
    # Aggregate across channels: total |SHAP| per position, per class.
    agg = df.groupby("relative_position", as_index=False)[["nmd_abs_shap", "ctrl_abs_shap"]].sum()
    agg = agg.sort_values("relative_position").reset_index(drop=True)
    return agg


def smooth(series, window=11):
    """Rolling-mean smoothing — the per-position raw signal has trinucleotide-
    frame periodicity that reads as noise at figure scale; a small rolling
    average clarifies the class-level pattern without hiding real features."""
    return series.rolling(window=window, center=True, min_periods=1).mean()


def draw_panel(ax, agg, *, codon_label):
    style_axes_ggplot(ax)
    ax.plot(agg["relative_position"], smooth(agg["nmd_abs_shap"]),
            color=NMD_COLOR, linewidth=1.6, label="NMD susceptible", zorder=3)
    ax.plot(agg["relative_position"], smooth(agg["ctrl_abs_shap"]),
            color=CONTROL_COLOR, linewidth=1.6, label="Control", zorder=3)
    ax.axvline(0, color="#555555", linewidth=0.8, linestyle=":", zorder=2)
    ax.set_xlabel(f"Position relative to {codon_label} (nt)", fontsize=BODY_FS)
    ax.set_ylabel("Total |SHAP|", fontsize=BODY_FS)
    ax.set_xlim(agg["relative_position"].min(), agg["relative_position"].max())
    ax.set_ylim(bottom=0)
    return ax


def main():
    atg = load_window("atg")
    stop = load_window("stop")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.82, bottom=0.12, wspace=0.20)

    draw_panel(axes[0], atg, codon_label="start codon")
    draw_panel(axes[1], stop, codon_label="stop codon")

    facet_header(axes[0], "Start (AUG) window", height=0.08)
    facet_header(axes[1], "Stop codon window",   height=0.08)

    panel_letter(axes[0], "A", x=-0.07, y=1.14)
    panel_letter(axes[1], "B", x=-0.07, y=1.14)

    # Legend in the top-right of the second panel
    axes[1].legend(loc="upper right", frameon=True, framealpha=0.9,
                    edgecolor="none", facecolor="white", fontsize=BODY_FS)

    assert_text_within_canvas(fig)

    out_png = HERE / "figure_sf36_shap_across_windows.png"
    out_pdf = HERE / "figure_sf36_shap_across_windows.pdf"
    fig.savefig(out_png, dpi=200, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    print(f"wrote {out_png.name} and {out_pdf.name}")
    print(f"  n positions ATG: {len(atg)}   STOP: {len(stop)}")


if __name__ == "__main__":
    main()
