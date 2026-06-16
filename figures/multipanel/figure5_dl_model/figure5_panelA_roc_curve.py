"""
Figure 5 Panel A — ROC curve of the best NMD-prediction model.

Mirrors the orf_model_report fig1d_roc_curve content: best-model
predictions on the chr-1/3/5/7 paralog-free holdout, AUC + AUPRC
annotated.

Data:
  ./data/panelA_predictions.tsv  — per-test-transcript (label, prob)
  ./data/panelA_metrics.tsv      — best-model AUC, AUPRC, n_test, n_nmd

Style: publications overlay — 6×4 cell at 1.5:1, no in-panel title,
no bbox_inches="tight" on save, validator wired into main().
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc as roc_auc_score

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14

# Curves: NMD-coral, with a subtle diagonal reference line.
C_ROC = "#d6604d"
C_DIAG = "#888888"
C_TITLE = "#222222"
C_AXIS = "#555555"


def load_data():
    preds = pd.read_csv(HERE.parent / "data" / "panelA_predictions.tsv", sep="\t")
    metrics = pd.read_csv(HERE.parent / "data" / "panelA_metrics.tsv", sep="\t").iloc[0]
    return preds, metrics


def build_figure(preds, metrics):
    fig, ax = plt.subplots(figsize=(6, 4))

    fpr, tpr, _ = roc_curve(preds["label"], preds["prob"])
    auc_val = roc_auc_score(fpr, tpr)

    ax.plot(fpr, tpr, color=C_ROC, linewidth=2.0, solid_capstyle="round")
    ax.plot([0, 1], [0, 1], color=C_DIAG, linestyle="--", linewidth=1.0, alpha=0.7)

    ax.fill_between(fpr, tpr, color=C_ROC, alpha=0.08)

    # Annotation block at lower-right with AUC + AUPRC + n
    n_test = int(metrics["n_test"])
    n_nmd = int(metrics["n_nmd"])
    n_ctrl = n_test - n_nmd
    auprc = float(metrics["auprc"])
    txt = (
        f"AUC = {auc_val:.3f}\n"
        f"AUPRC = {auprc:.3f}\n"
        f"n = {n_test:,} test transcripts\n"
        f"({n_nmd:,} NMD | {n_ctrl:,} Control)"
    )
    ax.text(
        0.97, 0.06, txt,
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=BODY_FS, color=C_TITLE,
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("False positive rate", fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylabel("True positive rate", fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="both", labelsize=BODY_FS, colors=C_TITLE)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.16, top=0.96)
    return fig


def main():
    preds, metrics = load_data()
    fig = build_figure(preds, metrics)

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelA_roc_curve.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelA_roc_curve.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelA_roc_curve.pdf and .png")


if __name__ == "__main__":
    main()
