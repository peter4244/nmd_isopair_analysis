"""
Figure 5 Panel E — SHAP nucleotide logo across ±20 nt around the stop codon.

Mirrors old Fig 6F ("SHAP values at stop codon"): two stacked logos —
NMD on top, Control on bottom — with per-position A/C/G/T letters scaled
by |SHAP|. Model has learned UGA as the leakiest stop and the importance
of the +4 position for readthrough.

Data:
  ./data/panelE_stop_profile.tsv (same schema as Panel D).
  Source: shap_profile_stop_joint_atg500_stop500.tsv (4ct).

Style: publications overlay — 6×4 cell at 1.5:1.
"""

import sys
from pathlib import Path

import logomaker as lm
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

BODY_FS = 14
LABEL_FS = 11
ZOOM = 20
CODON_LABEL = "stop"
NUC_COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}

C_TITLE = "#222222"
C_AXIS  = "#555555"
C_CODON = "#444444"


def load_matrix(value_col):
    df = pd.read_csv(HERE.parent / "data" / "panelE_stop_profile.tsv", sep="\t")
    df = df[df["channel"].isin(["A", "C", "G", "T"])]
    df = df[df["relative_position"].between(-ZOOM, ZOOM)]
    mat = (
        df.pivot(index="relative_position", columns="channel", values=value_col)
          .reindex(columns=["A", "C", "G", "T"])
          .sort_index()
    )
    return mat


def draw_logo(ax, mat, title):
    lm.Logo(mat, ax=ax, color_scheme=NUC_COLORS, baseline_width=0)
    ax.axvline(0, color=C_CODON, linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_ylabel("|SHAP|", fontsize=LABEL_FS, color=C_TITLE)
    ax.set_xlim(-ZOOM - 0.5, ZOOM + 0.5)
    ax.text(0.01, 0.94, title, transform=ax.transAxes,
            ha="left", va="top", fontsize=LABEL_FS, fontweight="bold",
            color=C_TITLE)
    ax.tick_params(axis="both", labelsize=LABEL_FS - 1, colors=C_TITLE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.7)


def build_figure():
    fig, axes = plt.subplots(2, 1, figsize=(6, 4), sharex=True)

    nmd_mat = load_matrix("nmd_abs_shap")
    ctrl_mat = load_matrix("ctrl_abs_shap")
    ymax = max(nmd_mat.values.sum(axis=1).max(),
               ctrl_mat.values.sum(axis=1).max()) * 1.05

    draw_logo(axes[0], nmd_mat, "NMD")
    draw_logo(axes[1], ctrl_mat, "Control")

    for a in axes:
        a.set_ylim(0, ymax)

    axes[1].set_xlabel(f"Position relative to {CODON_LABEL}", fontsize=BODY_FS, color=C_TITLE)
    axes[0].tick_params(axis="x", labelbottom=False)

    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.16, top=0.96, hspace=0.10)
    return fig


def main():
    fig = build_figure()

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelE_stop_logo.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelE_stop_logo.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelE_stop_logo.pdf and .png")


if __name__ == "__main__":
    main()
