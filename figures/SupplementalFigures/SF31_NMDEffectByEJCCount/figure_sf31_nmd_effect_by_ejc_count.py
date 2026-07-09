"""SF31 — NMD effect size by downstream EJC count (PTC+ comparators).

Violin + box distributions of mashr posterior-mean logFC (SMG1i vs DMSO,
averaged across AT2/LAE/FB/MV) as a function of the comparator's number of
downstream EJCs. Restricted to PTC-positive NMD comparators (has_ptc == TRUE)
from SF30's population. EJC counts ≥ 7 collapsed into a "7+" bin.

Data source: `sf31_ejc_count_logfc.tsv` (produced by data_export.R alongside
SF30's TSV). Matches the population and derivation of the Rmd chunk
`goal1-fig3-ejc-boxplot` in `05_final_report_mashr.Rmd`
§ "NMD Strength by Downstream EJC Count".
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
    assert_text_within_canvas,
    CONTROL_COLOR,
    AXIS_C,
    TITLE_C,
    BODY_FS,
    HEADER_FS,
)

apply_ggplot_rcparams()

DATA = HERE / "data"


def main():
    df = pd.read_csv(DATA / "sf31_ejc_count_logfc.tsv", sep="\t")
    df = df.dropna(subset=["mean_logFC", "n_downstream_ejcs"])

    # 7+ bin
    df["ejc_bin"] = np.where(df["n_downstream_ejcs"] >= 7, "7+",
                              df["n_downstream_ejcs"].astype(int).astype(str))
    order = [str(i) for i in range(1, 7)] + ["7+"]
    df = df[df["ejc_bin"].isin(order)].copy()

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    fig.subplots_adjust(left=0.09, right=0.97, top=0.86, bottom=0.14)
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    data = [df.loc[df["ejc_bin"] == b, "mean_logFC"].to_numpy() for b in order]
    positions = np.arange(len(order))

    parts = ax.violinplot(data, positions=positions, widths=0.75,
                           showmeans=False, showmedians=False, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(CONTROL_COLOR)
        body.set_edgecolor(AXIS_C)
        body.set_alpha(0.55)
        body.set_zorder(3)

    # Overlay narrow boxplots
    ax.boxplot(data, positions=positions, widths=0.14,
                showfliers=False, patch_artist=True,
                boxprops=dict(facecolor="white", edgecolor=TITLE_C, linewidth=1.0),
                medianprops=dict(color=TITLE_C, linewidth=1.4),
                whiskerprops=dict(color=TITLE_C, linewidth=1.0),
                capprops=dict(color=TITLE_C, linewidth=1.0),
                zorder=5)

    # n = labels above each violin
    ymax = df["mean_logFC"].quantile(0.98)
    y_label = ymax + 0.35
    for i, b in enumerate(order):
        n = len(df.loc[df["ejc_bin"] == b])
        ax.text(positions[i], y_label, f"n = {n:,}",
                ha="center", va="bottom", fontsize=BODY_FS - 1, color=TITLE_C)

    # Median line
    medians = [np.median(d) for d in data]
    ax.plot(positions, medians, color="#c0392b", linewidth=1.2,
            marker="o", markersize=4, zorder=6, label="Median trend")
    ax.legend(loc="lower right", frameon=True, facecolor="white",
              edgecolor="none", fontsize=BODY_FS - 1)

    ax.set_xticks(positions)
    ax.set_xticklabels(order, fontsize=BODY_FS, color=TITLE_C)
    ax.set_xlabel("Number of downstream EJCs (comparator)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Mean logFC (mashr posterior mean)",
                   fontsize=BODY_FS, color=TITLE_C)
    y_hi = max(df["mean_logFC"].quantile(0.995), y_label + 0.6)
    y_lo = df["mean_logFC"].quantile(0.005) - 0.2
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlim(-0.6, len(order) - 0.4)

    fig.text(
        0.5, 0.94,
        f"NMD effect size by downstream EJC count (PTC+ comparators, n = {len(df):,})",
        ha="center", va="center",
        fontsize=HEADER_FS, fontweight="bold", color=TITLE_C,
    )

    assert_text_within_canvas(fig)

    out_png = HERE / "figure_sf31_nmd_effect_by_ejc_count.png"
    out_pdf = HERE / "figure_sf31_nmd_effect_by_ejc_count.pdf"
    fig.savefig(out_png, dpi=200, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    print(f"wrote {out_png.name} and {out_pdf.name}")
    print(f"  medians: {[round(m,2) for m in medians]}")
    print(f"  n_by_bin: {[len(d) for d in data]}")


if __name__ == "__main__":
    main()
