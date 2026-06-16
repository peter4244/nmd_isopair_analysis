"""
Figure 4 Panel C — Kozak PWM score: reference AUG vs TD2-predicted AUG.

For 1,032 NMD+/PTC− pairs where the reference AUG is exonic in the
comparator AND TD2 picked a different AUG, the Kozak PWM (log-odds)
score is computed at both positions via `Isopair::scoreKozakPWM`.
Reference AUGs are stronger initiation contexts: ref median = 0.88
vs TD2 median = −0.31, paired Wilcoxon p = 2.5×10⁻⁵⁹.

Data:
  ./data/panelD_kozak.tsv
    Columns: comparator_isoform_id, ATG_source, score
    Source: ref_atg_analysis.rds + transcript FASTA + Isopair PWM.
    Exported by ./data_export.R.

Style: publications overlay — 6×4 cell at 1.5:1, no in-panel title,
no bbox_inches="tight" on save, validator wired into main().
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14

C_REF = "#4393c3"   # Reference ATG: reference-blue (matches Fig 3 schematic ref-track)
C_TD2 = "#ef8a62"   # TD2-predicted ATG: NMD coral (the "bad" choice in this story)
C_TITLE = "#222222"
C_AXIS = "#555555"


def load_data():
    return pd.read_csv(HERE.parent / "data" / "panelD_kozak.tsv", sep="\t")


def build_figure(df):
    fig, ax = plt.subplots(figsize=(6, 4))

    sources = ["Reference ATG", "TD2-predicted ATG"]
    palette = {"Reference ATG": C_REF, "TD2-predicted ATG": C_TD2}
    data = [df.loc[df["ATG_source"] == s, "score"].to_numpy() for s in sources]
    n_per_group = {s: len(d) for s, d in zip(sources, data)}

    sns.violinplot(data=df, x="ATG_source", y="score", order=sources,
                   hue="ATG_source", palette=palette, legend=False,
                   cut=0, inner="quart", density_norm="width",
                   linewidth=0.9, ax=ax)
    for patch in ax.collections:
        patch.set_edgecolor(C_TITLE); patch.set_alpha(0.78)
    for line in ax.lines:
        line.set_color(C_TITLE); line.set_linewidth(1.0)

    # Significance bracket above
    ymax = max(np.max(d) for d in data)
    ymin = min(np.min(d) for d in data)
    y_bracket = ymax + 0.6
    ax.plot([0, 0, 1, 1], [y_bracket - 0.15, y_bracket, y_bracket, y_bracket - 0.15],
            color=C_TITLE, linewidth=0.8)
    ax.text(0.5, y_bracket + 0.1,
            "paired Wilcoxon\np = 2.5×10$^{-59}$",
            ha="center", va="bottom",
            fontsize=BODY_FS - 1, color=C_TITLE)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{s}\n(n={n_per_group[s]:,})" for s in sources],
                       fontsize=BODY_FS - 1, color=C_TITLE)
    ax.set_xlabel("")
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_ylabel("Kozak PWM score (log-odds)", fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylim(bottom=ymin - 0.3, top=y_bracket + 1.4)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.20, top=0.96)
    return fig


def main():
    df = load_data()
    fig = build_figure(df)

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure4_panelD_kozak.pdf", facecolor="white")
    fig.savefig(out_dir / "figure4_panelD_kozak.png", dpi=300, facecolor="white")
    print("Saved: figure4_panelD_kozak.pdf and .png")


if __name__ == "__main__":
    main()
