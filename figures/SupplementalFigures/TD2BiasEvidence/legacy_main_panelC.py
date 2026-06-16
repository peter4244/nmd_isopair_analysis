"""
Figure 4 Panel B — TD2-selected CDS length vs reference-AUG ORF length.

For the 900 "Recovered PTC" pairs (effectively_ptc), comparing the
TD2-selected CDS length on the comparator with the ref-AUG-traced ORF
length on the same comparator. TD2 picks a longer ORF in **100% of
pairs** (literal "always longer"), with TD2 median = 3,168 nt vs
ref-AUG median = 372 nt — an 8.5× ratio of medians.

Data:
  ./data/panelC_td2_vs_refaug_length.tsv
    Columns: comparator_isoform_id, ORF_source, length_nt
    Source: ref_atg_analysis.rds + cds.rds. Exported by ./data_export.R.

Style: publications overlay — 6×4 cell at 1.5:1, no in-panel title,
no bbox_inches="tight" on save, validator wired into main().
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14

C_TD2_FILL  = "#67a9cf"; C_TD2_LINE  = "#2b8cbe"
C_REF_FILL  = "#ef8a62"; C_REF_LINE  = "#d6604d"
C_TITLE = "#222222"
C_AXIS = "#555555"


def load_data():
    return pd.read_csv(HERE.parent / "data" / "panelC_td2_vs_refaug_length.tsv", sep="\t")


def build_figure(df):
    fig, ax = plt.subplots(figsize=(6, 4))

    # Log-10 lengths
    x_grid = np.linspace(np.log10(50), np.log10(15000), 400)
    peaks = {}

    for label, fill, line in [
        ("Reference ATG ORF", C_REF_FILL, C_REF_LINE),
        ("TD2-selected CDS",  C_TD2_FILL, C_TD2_LINE),
    ]:
        vals = np.log10(df.loc[df["ORF_source"] == label, "length_nt"].to_numpy())
        kde = gaussian_kde(vals, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(10 ** x_grid, y, alpha=0.45, color=fill, linewidth=0)
        ax.plot(10 ** x_grid, y, color=line, linewidth=1.6)
        peak_i = int(np.argmax(y))
        peaks[label] = (10 ** x_grid[peak_i], y[peak_i], line, len(vals))

    ax.set_xscale("log")
    ax.set_xlim(50, 15000)
    y_max = max(p[1] for p in peaks.values())
    ax.set_ylim(0, y_max * 1.22)
    ax.set_xlabel("ORF length (nt, log scale)", fontsize=BODY_FS, color=C_TITLE)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=BODY_FS, colors=C_TITLE)

    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(C_AXIS)
    ax.spines["bottom"].set_linewidth(0.8)

    # Inline curve labels above each KDE peak (avoids legend-over-curve overlap)
    for label, (xp, yp, line_c, n) in peaks.items():
        ax.text(xp, yp + y_max * 0.06,
                f"{label}\n(n={n:,})",
                ha="center", va="bottom",
                fontsize=BODY_FS - 1, color=line_c, fontweight="bold")

    fig.subplots_adjust(left=0.04, right=0.97, bottom=0.18, top=0.96)
    return fig


def main():
    df = load_data()
    fig = build_figure(df)

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure4_panelC_td2_vs_refaug_length.pdf", facecolor="white")
    fig.savefig(out_dir / "figure4_panelC_td2_vs_refaug_length.png", dpi=300, facecolor="white")
    print("Saved: figure4_panelC_td2_vs_refaug_length.pdf and .png")


if __name__ == "__main__":
    main()
