"""
Figure 3 Panel D — Stop-codon distance from last EJC (NMD vs Control).

Density of distance from the comparator's stop codon to its last exon-exon
junction. NMD comparators concentrate to the RIGHT of the 50-nt PTC threshold
(stops downstream of the last EJC); Control comparators concentrate LEFT
(normal stops in the last exon).

Stop source: ref-AUG-traced stop where ref AUG was found in the comparator
(categories effectively_ptc, no_downstream_ejc, truncated_no_ejc; ~86% of NMD
and ~68% of Control). TD2 stop position elsewhere.

Data:
  ./data/panelD_stop_codon_distance.tsv
    Columns: comparator_isoform_id, comparison ("NMD" | "Control"),
             distance (int, nt), stop_source ("ref_aug_traced" | "td2"),
             category (ref-AUG class or NA)
    Source: cod_c2 / cod_c4 (canonical 2,496-pair gene-matched coding-coding
            set) + structures.rds (last EJC tx position) + ref_atg_analysis
            (ref-AUG stop) + ptc.rds (TD2 fallback distance).
    Exported by ./data_export.R.

Sample size: 2,496 / 2,496 (same denominator as Panels B, C, E).

X-axis clipped to [-1000, 1500] nt; clipped counts printed to stdout for
methodology file disclosure.

Style: HEADER_FS=18 bold, BODY_FS=14 regular. NMD coral, Control light blue.
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

C_NMD_FILL = "#ef8a62"
C_NMD_LINE = "#d6604d"
C_CTRL_FILL = "#67a9cf"
C_CTRL_LINE = "#2b8cbe"
C_TITLE = "#222222"
C_AXIS = "#555555"
C_THRESHOLD = "#cc0000"

X_MIN, X_MAX = -1000, 1500


def load_data():
    return pd.read_csv(HERE.parent / "data" / "panelD_stop_codon_distance.tsv", sep="\t")


def build_figure(df):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    x_grid = np.linspace(X_MIN, X_MAX, 800)
    n_clipped = {}
    for label, fill, line in [
        ("NMD", C_NMD_FILL, C_NMD_LINE),
        ("Control", C_CTRL_FILL, C_CTRL_LINE),
    ]:
        vals = df.loc[df["comparison"] == label, "distance"].dropna().to_numpy()
        kde = gaussian_kde(vals, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(x_grid, y, alpha=0.45, color=fill, linewidth=0,
                        label=f"{label} (n={len(vals):,})")
        ax.plot(x_grid, y, color=line, linewidth=1.6)
        n_clipped[label] = int(((vals < X_MIN) | (vals > X_MAX)).sum())

    # PTC threshold
    ax.axvline(50, color=C_THRESHOLD, linestyle="--", linewidth=1.2, alpha=0.7)
    ymax = ax.get_ylim()[1]
    ax.text(50 + 40, ymax * 0.96, "50-nt PTC\nthreshold",
            ha="left", va="top", fontsize=BODY_FS, color=C_THRESHOLD)

    ax.set_title("Stop codon to last EJC", fontsize=HEADER_FS,
                 fontweight="bold", color=C_TITLE, loc="left", pad=12)
    ax.set_xlabel("Distance from last EJC (nt)", fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylabel("Density", fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="both", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(bottom=0)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    legend = ax.legend(fontsize=BODY_FS, loc="upper left", frameon=False)
    for text in legend.get_texts():
        text.set_color(C_TITLE)

    print(f"  Clipped at x∈[{X_MIN}, {X_MAX}]: NMD={n_clipped['NMD']}, Control={n_clipped['Control']}")

    fig.tight_layout()
    return fig


def main():
    df = load_data()
    fig = build_figure(df)
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure3_panelD_stop_codon_distance.pdf",
                bbox_inches="tight", facecolor="white")
    fig.savefig(out_dir / "figure3_panelD_stop_codon_distance.png",
                dpi=300, bbox_inches="tight", facecolor="white")
    print("Saved: figure3_panelD_stop_codon_distance.pdf and .png")


if __name__ == "__main__":
    main()
