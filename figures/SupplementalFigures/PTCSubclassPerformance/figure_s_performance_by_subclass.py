"""
Supplemental Figure — Deep-learning NMD-call discrimination, stratified by
PTC subclass within the n = 1,166 ref-AUG-traceable subset.

Companion to the §5 manuscript sentence:

    "However, overall predictive performance was substantially lower in
     NMD+/PTC- isoforms (SFx)."

Two panels at the held-out chr-1/3/5/7 paralog-free test split:

  (A) ROC curves for two contrasts:
        — NMD+/PTC+ retained  vs  Control
        — NMD+/PTC- retained  vs  Control
      AUC and AUPRC annotated.

  (B) Per-isoform predicted NMD probability stratified by subclass
      (strip + boxplot overlay), with the 0.5 decision threshold marked.

KEY NUMBERS (from data/predprob_by_subclass_n1166.tsv):
  NMD+/PTC+ vs Control:  AUC = 0.955, AUPRC = 0.940   (n=255 vs n=276)
  NMD+/PTC- vs Control:  AUC = 0.737, AUPRC = 0.258   (n=30  vs n=276)
  Mean predicted prob:   Control 0.19, NMD+/PTC+ 0.85, NMD+/PTC- 0.39.

Scope note: test-set predictions only (chr-1/3/5/7, paralog-free). NMD+/
PTC- has n = 30 at the test split.

Data: data/predprob_by_subclass_n1166.tsv (built by sibling
PTCSubclassBranchSHAP/data_export.R).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

try:
    from validate_figure_layout import validate_multipanel_layout
except Exception:
    validate_multipanel_layout = None

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

DATA = HERE.parent / "data" / "predprob_by_subclass_n1166.tsv"

GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC- retained", "Control"]
GROUP_LABEL = {
    "NMD+/PTC+":          "NMD+/PTC+",
    "NMD+/PTC- retained": "NMD+/PTC−",
    "Control":            "Control",
}
# Canonical NMD subgroup palette — verbatim from orf_model_report_v5.Rmd:48-50
# and project_nmd_figure_palette memory. Do NOT swap for matplotlib defaults;
# the same group must read in the same color across the whole paper.
GROUP_COLOR = {
    "NMD+/PTC+":          "#ef8a62",   # peach (= binary NMD)
    "NMD+/PTC- retained": "#d95f02",   # orange (Rmd "NMD PTC-, ref ATG retained")
    "Control":            "#67a9cf",   # blue (= binary Control)
}


def load_split():
    df = pd.read_csv(DATA, sep="\t")
    ctrl = df[df["group"] == "Control"]
    ptcp = df[df["group"] == "NMD+/PTC+"]
    ptcn = df[df["group"] == "NMD+/PTC- retained"]
    return df, ctrl, ptcp, ptcn


def contrast_metrics(pos: pd.DataFrame, neg: pd.DataFrame):
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    p = np.concatenate([pos["prob"].values, neg["prob"].values])
    fpr, tpr, _ = roc_curve(y, p)
    auc = roc_auc_score(y, p)
    auprc = average_precision_score(y, p)
    return fpr, tpr, auc, auprc


def build_figure():
    df, ctrl, ptcp, ptcn = load_split()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.2, 3.8))
    fig.subplots_adjust(left=0.085, right=0.99, bottom=0.16, top=0.88, wspace=0.28)

    # ── Panel A — ROC curves
    fpr1, tpr1, auc1, aupr1 = contrast_metrics(ptcp, ctrl)
    fpr2, tpr2, auc2, aupr2 = contrast_metrics(ptcn, ctrl)
    axA.plot(fpr1, tpr1, color=GROUP_COLOR["NMD+/PTC+"], lw=2.0,
             label=f"NMD+/PTC+ vs Ctrl\nAUC = {auc1:.2f}, AUPRC = {aupr1:.2f}\n(n = {len(ptcp)} vs {len(ctrl)})")
    axA.plot(fpr2, tpr2, color=GROUP_COLOR["NMD+/PTC- retained"], lw=2.0,
             label=f"NMD+/PTC− vs Ctrl\nAUC = {auc2:.2f}, AUPRC = {aupr2:.2f}\n(n = {len(ptcn)} vs {len(ctrl)})")
    axA.plot([0, 1], [0, 1], "--", color="#666666", lw=0.8)
    axA.set_xlim(-0.02, 1.02)
    axA.set_ylim(-0.02, 1.02)
    axA.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axA.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axA.set_xlabel("False positive rate", fontsize=9)
    axA.set_ylabel("True positive rate", fontsize=9)
    axA.set_title("Discrimination vs Control", fontsize=10, pad=4)
    axA.tick_params(axis="both", labelsize=8)
    axA.legend(loc="lower right", fontsize=7.5, frameon=False, handlelength=2.0,
               handletextpad=0.5, borderaxespad=0.4, labelspacing=0.9)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.set_aspect("equal")

    # ── Panel B — predicted probability per subgroup
    rng = np.random.default_rng(seed=1)
    xs = []
    ys = []
    cs = []
    positions = list(range(len(GROUP_ORDER)))
    for i, g in enumerate(GROUP_ORDER):
        vals = df.loc[df["group"] == g, "prob"].values
        jitter = rng.uniform(-0.18, 0.18, size=len(vals))
        xs.extend(i + jitter)
        ys.extend(vals)
        cs.extend([GROUP_COLOR[g]] * len(vals))
    axB.scatter(xs, ys, c=cs, s=10, alpha=0.55, edgecolor="none",
                zorder=2, rasterized=True)

    # Boxplot overlay
    bdata = [df.loc[df["group"] == g, "prob"].values for g in GROUP_ORDER]
    bp = axB.boxplot(
        bdata, positions=positions, widths=0.45,
        patch_artist=True, showfliers=False,
        boxprops=dict(facecolor="none", edgecolor="black", linewidth=1.0),
        medianprops=dict(color="black", linewidth=1.4),
        whiskerprops=dict(color="black", linewidth=1.0),
        capprops=dict(color="black", linewidth=1.0),
        zorder=3,
    )

    axB.axhline(0.5, color="#aa6600", lw=0.9, ls="--", zorder=1)
    axB.text(2.45, 0.51, "decision\nthreshold (0.5)", color="#aa6600",
             fontsize=7, ha="right", va="bottom", style="italic")
    axB.set_xticks(positions)
    axB.set_xticklabels([
        f"{GROUP_LABEL[g]}\n(n = {int((df['group'] == g).sum()):,})"
        for g in GROUP_ORDER
    ], fontsize=8)
    axB.set_xlim(-0.6, len(GROUP_ORDER) - 0.4)
    axB.set_ylim(-0.03, 1.03)
    axB.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axB.set_ylabel("Model NMD probability  P(NMD | x)", fontsize=9)
    axB.set_title("Predicted probability per subclass", fontsize=10, pad=4)
    axB.tick_params(axis="y", labelsize=8)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    # Panel labels
    axA.text(-0.13, 1.02, "A", transform=axA.transAxes,
             fontsize=14, fontweight="bold", va="bottom", ha="left")
    axB.text(-0.13, 1.02, "B", transform=axB.transAxes,
             fontsize=14, fontweight="bold", va="bottom", ha="left")

    # Footnote
    fig.text(
        0.5, 0.025,
        f"AUC drops from {auc1:.2f} (NMD+/PTC+ vs Ctrl) to {auc2:.2f} "
        "(NMD+/PTC− vs Ctrl) — held-out test set, chr-1/3/5/7 paralog-free.",
        ha="center", va="bottom", fontsize=8, style="italic",
    )
    return fig, (axA, axB)


def main():
    fig, axes = build_figure()
    if validate_multipanel_layout is not None:
        try:
            validate_multipanel_layout(fig)
        except Exception as e:
            print(f"[validate_multipanel_layout] non-fatal: {e}")
    out = HERE.parent / "figure_s_performance_by_subclass"
    fig.savefig(f"{out}.pdf", facecolor="white")
    fig.savefig(f"{out}.png", dpi=300, facecolor="white")
    print(f"Saved: {out}.pdf and .png")


if __name__ == "__main__":
    main()
