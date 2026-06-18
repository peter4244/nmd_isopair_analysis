"""
Supplemental Figure — Branch-level KernelSHAP stratified by PTC subclass
across the FULL n = 1,166 ref-AUG-traceable subset (not test-set-restricted).

Style: matches Figure 5 Panel C exactly — branch order Structural / Stop /
ATG, same colour palette, percent label above each bar, mean |SHAP|
inside each bar, no in-panel title. Three side-by-side mini-panels (one
per subgroup) share a common y-axis so heights compare directly.

Companion to the §5 manuscript sentence:

    "We also examined the Shapley scores for the START and STOP windows
     stratified by isoform PTC subclass, and we observed that the START
     window information is roughly three times more important in the
     NMD+/PTC- than in the NMD+/PTC+ isoforms (SFx)."

KEY NUMBER (full cohort, no test-set restriction):
  Mean within-isoform ATG-branch share:
    NMD+/PTC+         (n=1,016) →  9.2%
    NMD+/PTC- retained (n=  95) → 20.4%
    Control           (n=1,107) → 11.3%
  Ratio NMD+/PTC- vs NMD+/PTC+ ≈ 2.21×.

The Panel-C-style bar labels show share-of-subgroup-total (sum of
subgroup mean |SHAP|), giving:
    NMD+/PTC+   Structural 61.9% / Stop 29.1% / ATG  9.0%
    NMD+/PTC-   Structural 52.6% / Stop 29.2% / ATG 18.1%
    Control     Structural 62.3% / Stop 27.1% / ATG 10.6%

Manuscript text reads "roughly three times more important" — the actual
share-of-total ratio is 2.21× (full cohort) or 2.20× (test-set only).
Find/replace pair to revise the manuscript prose is in README.md.

Data: data/branch_shap_by_subclass_n1166.tsv (built by data_export.R).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

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

# ── Panel C styling constants (mirrored verbatim from figure5_panelC_branch_importance.py) ──
HEADER_FS = 18
BODY_FS = 14

BRANCH_COLORS = {
    "Structural": "#d95f02",
    "Stop":       "#FF6B6B",
    "ATG":        "#4ECDC4",
}
C_TITLE = "#222222"
C_AXIS  = "#555555"

BRANCH_ORDER = ["Structural", "Stop", "ATG"]

DATA = HERE.parent / "data" / "branch_shap_by_subclass_n1166.tsv"

GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC- retained", "Control"]
GROUP_LABEL = {
    "NMD+/PTC+":          "NMD+/PTC+",
    "NMD+/PTC- retained": "NMD+/PTC−",
    "Control":            "Control",
}


def per_group_summary(df: pd.DataFrame, group: str):
    """Return a DataFrame matching Figure 5 Panel C's per-branch table:
       columns branch, mean_abs_shap, pct (within-subgroup share-of-total),
       with branch order = Structural / Stop / ATG."""
    sub = df[df["group"] == group]
    means = (
        sub.groupby("branch")["abs_shap"].mean().reindex(BRANCH_ORDER).reset_index()
    )
    total = means["abs_shap"].sum()
    means["pct"] = 100 * means["abs_shap"] / total
    means.columns = ["branch", "mean_abs_shap", "pct"]
    means["n"] = int((df["group"] == group).shape[0]) if False else \
        int(df.loc[df["group"] == group, "comparator_isoform_id"].nunique())
    return means


def draw_panel(ax, summary: pd.DataFrame, group_label: str, n: int, y_top: float):
    """Render one Panel-C-style mini-panel on `ax` for one subgroup.

    Y-axis is within-subgroup share of total mean |SHAP|, in %. The mean
    |SHAP| absolute value is retained inside each bar (white bold) so the
    Panel-C dual-label pattern is preserved.
    """
    branches = summary["branch"].tolist()
    pcts = summary["pct"].to_numpy()                # within-subgroup share, %
    vals = summary["mean_abs_shap"].to_numpy()      # absolute mean |SHAP|
    colors = [BRANCH_COLORS[b] for b in branches]

    x = range(len(branches))
    ax.bar(x, pcts, color=colors, edgecolor="none", width=0.62)

    # Percent ABOVE each bar (bold black, 1 decimal).
    for i, (v, p) in enumerate(zip(vals, pcts)):
        ax.text(i, p + y_top * 0.025, f"{p:.1f}%",
                ha="center", va="bottom",
                fontsize=BODY_FS, fontweight="bold", color=C_TITLE)
        # Mean |SHAP| INSIDE each bar (white bold), or above if bar is too
        # short to host a centered white label legibly.
        if p > y_top * 0.10:
            ax.text(i, p / 2, f"{v:.2f}",
                    ha="center", va="center",
                    fontsize=BODY_FS - 2, fontweight="bold", color="white")
        else:
            ax.text(i, p + y_top * 0.09, f"{v:.2f}",
                    ha="center", va="bottom",
                    fontsize=BODY_FS - 2, fontweight="bold", color=C_TITLE)

    ax.set_xticks(list(x))
    ax.set_xticklabels(branches, fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)
    ax.set_ylim(0, y_top)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    # Subgroup header (above the axes — replaces Panel C's nonexistent title
    # while still tagging which subgroup this mini-panel describes).
    ax.text(0.5, 1.04,
            f"{group_label}",
            transform=ax.transAxes,
            ha="center", va="bottom",
            fontsize=HEADER_FS, fontweight="bold", color=C_TITLE)
    ax.text(0.5, 0.97,
            f"n = {n:,}",
            transform=ax.transAxes,
            ha="center", va="bottom",
            fontsize=BODY_FS - 2, color=C_AXIS)


def build_figure():
    df = pd.read_csv(DATA, sep="\t")
    summaries = {g: per_group_summary(df, g) for g in GROUP_ORDER}
    n_by_group = {
        g: int(df.loc[df["group"] == g, "comparator_isoform_id"].nunique())
        for g in GROUP_ORDER
    }
    # Y-axis is now percentage (within-subgroup share). The largest bar
    # across all subgroups is the Control Structural at ~62%; cap at 75%
    # to leave headroom for the percent-above-bar labels without colliding
    # with the subgroup header.
    y_top = 75.0

    # Three Panel-C-style mini-panels side by side. Width per panel
    # matches Panel C's 6×4 sizing; total width 18 in.
    fig, axes = plt.subplots(1, 3, figsize=(18, 4.5), sharey=True)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.12, top=0.83,
                        wspace=0.20)

    for ax, group in zip(axes, GROUP_ORDER):
        draw_panel(ax, summaries[group], GROUP_LABEL[group],
                   n_by_group[group], y_top)

    # Y-axis label only on leftmost panel (matches Panel C).
    axes[0].set_ylabel("Within-subgroup share of |SHAP|  (%)",
                       fontsize=BODY_FS, color=C_TITLE)
    # Tick marks every 20% — clean for a 0–75% range.
    for ax in axes:
        ax.set_yticks([0, 20, 40, 60])

    # Panel labels (top-left of each axes)
    for ax, letter in zip(axes, ["A", "B", "C"]):
        ax.text(-0.04, 1.13, letter, transform=ax.transAxes,
                fontsize=HEADER_FS + 2, fontweight="bold",
                va="bottom", ha="left", color=C_TITLE)

    # Footnote
    pct_atg_pos = summaries["NMD+/PTC+"].loc[
        summaries["NMD+/PTC+"]["branch"] == "ATG", "pct"].iloc[0]
    pct_atg_neg = summaries["NMD+/PTC- retained"].loc[
        summaries["NMD+/PTC- retained"]["branch"] == "ATG", "pct"].iloc[0]
    fig.text(
        0.5, 0.018,
        f"ATG-branch share of subgroup total: NMD+/PTC− = {pct_atg_neg:.1f}% "
        f"vs NMD+/PTC+ = {pct_atg_pos:.1f}%  →  {pct_atg_neg/pct_atg_pos:.1f}× higher in PTC−"
        f"  (n = 1,016 / 95 / 1,107 isoforms; full n = 1,166 subset).",
        ha="center", va="bottom", fontsize=11, style="italic", color=C_TITLE,
    )
    return fig, axes


def main():
    fig, axes = build_figure()
    if validate_multipanel_layout is not None:
        try:
            validate_multipanel_layout(fig)
        except Exception as e:
            print(f"[validate_multipanel_layout] non-fatal: {e}")
    out = HERE.parent / "figure_s_branch_shap_by_subclass"
    fig.savefig(f"{out}.pdf", facecolor="white")
    fig.savefig(f"{out}.png", dpi=300, facecolor="white")
    print(f"Saved: {out}.pdf and .png")


if __name__ == "__main__":
    main()
