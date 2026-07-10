"""
SF38 — Attention distribution across the model's five candidate ORFs,
NMD susceptible vs non-NMD.

Renders natively from the current 4CT model's per-isoform attention weights
(`uorf_attention_predictions.tsv`, columns attn_0…attn_4). Replaces the
earlier legacy-PNG placeholder that carried stale 10-ORF ranks, "SQANTI CDS"
subtitle, and a pre-4CT-switch n = 15,584 cohort.

Population: held-out test split (chr 1/3/5/7), excluding the "test_paralog"
partition per the training-time paralog dedup convention.

Panels:
  (A) Attention weight per candidate ORF rank, by class (boxplot).
  (B) Shannon entropy of the per-isoform attention vector across the five
      candidates, KDE per class.
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    render_and_validate,
    docx_body_fs,
    panel_letter,
    TITLE_C,
    AXIS_C,
)
apply_ggplot_rcparams()

NATIVE_W = 11.0
BODY_FS = docx_body_fs(NATIVE_W)
LABEL_FS = BODY_FS - 2

DATA_ATTN = Path("/Users/petecastaldi/claude_projects/NMD_orf_model_v5_4ct/"
                 "results_4ct/uorf_attention_predictions.tsv")
DATA_LABELS = Path("/Users/petecastaldi/claude_projects/NMD_orf_model_v5_4ct/"
                   "results_4ct/predictions_atg500_stop500.tsv")

# Main-paper Fig 4 / Fig 5G palette
COLOR_NMD  = "#ef8a62"   # coral
COLOR_CTRL = "#92c5de"   # light blue

N_ORFS = 5

# Test split = chr 1 / 3 / 5 / 7 (train-time held-out set). We take labels
# from predictions_atg500_stop500.tsv (canonical training-time source, also
# used by SF36 / SF37 / SF40 / SF41) rather than the attention TSV's own
# label column — 4 isoforms disagree between the two, and standardising on
# one source lets all SF36-SF42 cite the same 2,268 / 7,863 test-set n's.
TEST_CHRS = {"chr1", "chr3", "chr5", "chr7"}


def load():
    attn = pd.read_csv(DATA_ATTN, sep="\t")
    labels = pd.read_csv(DATA_LABELS, sep="\t")
    labels = labels[labels["chr"].isin(TEST_CHRS)][["isoform_id", "label"]]
    # Drop the attention TSV's own label; use labels from the canonical file.
    attn = attn.drop(columns=[c for c in ["label"] if c in attn.columns])
    df = attn.merge(labels, on="isoform_id", how="inner")
    return df


def shannon_entropy(row):
    p = np.array([row[f"attn_{k}"] for k in range(N_ORFS)], dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-(p * np.log2(p)).sum())


def render_panel_A(ax, df):
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    positions_ctrl = np.arange(N_ORFS) * 3 - 0.55
    positions_nmd  = np.arange(N_ORFS) * 3 + 0.55

    data_ctrl = [df.loc[df["label"] == 0, f"attn_{k}"].to_numpy() for k in range(N_ORFS)]
    data_nmd  = [df.loc[df["label"] == 1, f"attn_{k}"].to_numpy() for k in range(N_ORFS)]

    def draw_boxes(data, positions, color):
        ax.boxplot(
            data, positions=positions, widths=0.85,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor=color, edgecolor=TITLE_C, linewidth=0.9),
            medianprops=dict(color=TITLE_C, linewidth=1.4),
            whiskerprops=dict(color=TITLE_C, linewidth=0.9),
            capprops=dict(color=TITLE_C, linewidth=0.9),
            zorder=3,
        )

    draw_boxes(data_ctrl, positions_ctrl, COLOR_CTRL)
    draw_boxes(data_nmd,  positions_nmd,  COLOR_NMD)

    ax.set_xticks(np.arange(N_ORFS) * 3)
    ax.set_xticklabels([str(k) for k in range(N_ORFS)],
                        fontsize=BODY_FS, color=TITLE_C)
    ax.set_xlim(-2, (N_ORFS - 1) * 3 + 2)
    ax.set_ylim(-0.03, 1.03)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Candidate ORF rank (0 = main ORF)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Attention weight",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=TITLE_C)


def render_panel_B(ax, df):
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    df["entropy"] = df.apply(shannon_entropy, axis=1)
    e_nmd  = df.loc[df["label"] == 1, "entropy"].to_numpy()
    e_ctrl = df.loc[df["label"] == 0, "entropy"].to_numpy()

    max_ent = float(np.log2(N_ORFS))
    x_grid = np.linspace(0, max_ent, 300)

    for e, color in [(e_nmd, COLOR_NMD), (e_ctrl, COLOR_CTRL)]:
        kde = gaussian_kde(e, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(x_grid, y, alpha=0.55, color=color,
                        edgecolor=TITLE_C, linewidth=0.9, zorder=3)

    ax.set_xlim(0, max_ent * 1.02)
    ax.set_xlabel("Per-isoform attention entropy (bits)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Density", fontsize=BODY_FS, color=TITLE_C)
    ax.tick_params(axis="both", labelsize=BODY_FS, colors=TITLE_C)


def build_figure():
    df = load()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(NATIVE_W, 4.8))
    fig.subplots_adjust(left=0.09, right=0.95, top=0.80, bottom=0.22,
                        wspace=0.32)

    render_panel_A(axA, df)
    render_panel_B(axB, df)

    panel_letter(axA, "A", x=-0.14, y=1.22)
    panel_letter(axB, "B", x=-0.14, y=1.22)

    # Shared legend in the bottom margin — both panels share the same two
    # series; per-panel legends overlapped the boxes/density. Sample sizes
    # live in the caption, so the legend carries only the class colors.
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_NMD,  edgecolor=TITLE_C,
                      linewidth=0.9, label="NMD susceptible"),
        plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_CTRL, edgecolor=TITLE_C,
                      linewidth=0.9, label="Control"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=BODY_FS, bbox_to_anchor=(0.5, 0.01))

    return fig, df


def main():
    fig, df = build_figure()
    render_and_validate(fig, HERE.parent / "figure_s_attention_distribution",
                        native_width_in=NATIVE_W)
    print(f"Saved figure_s_attention_distribution.{{pdf,png}}")
    print(f"  NMD:     n = {int((df['label'] == 1).sum()):,}")
    print(f"  Control: n = {int((df['label'] == 0).sum()):,}")


if __name__ == "__main__":
    main()
