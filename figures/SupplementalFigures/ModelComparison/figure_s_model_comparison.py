"""
Supplemental Figure — Head-to-head comparison of our NMD predictor against
two variant-level NMD prediction models (NMDetective-B, Lindeboom 2019;
NMDEP-rule baseline, Saadat 2025) on the same cohort of long-read-observed
isoforms, using the mashr posterior-mean log fold-change under SMG1 inhibition
as the continuous gold standard.

Two panels:

  (A) Predicted score vs gold-standard log fold-change, one subplot per
      model, points coloured by PTC subclass (NMD+/PTC+, NMD+/PTC−, Control).
      Spearman correlation annotated per panel.

  (B) Within-subgroup Spearman correlation per model (grouped bar chart;
      bars = subclass × model). Shows where each model's predictions hold
      up vs break down.

Data: code/nmd_predictor_comparison/per_isoform_scores_2026.6.20.tsv
      code/nmd_predictor_comparison/metrics_summary_2026.6.20.tsv

Style: canonical NMD palette (peach/yellow/blue for PTC+/PTC−/Control);
       Arial; validator-clean.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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

ANALYSIS = HERE.parents[3] / "code" / "nmd_predictor_comparison"
SCORES = ANALYSIS / "per_isoform_scores_2026.6.20.tsv"
METRICS = ANALYSIS / "metrics_summary_2026.6.20.tsv"

# Canonical palette (project_nmd_figure_palette)
GROUP_LABEL = {
    "NMD+/PTC+": "NMD+/PTC+",
    "NMD+/PTC-": "NMD+/PTC−",
    "Control":   "Control",
}
GROUP_COLOR = {
    "NMD+/PTC+": "#ef8a62",
    "NMD+/PTC-": "#d95f02",
    "Control":   "#67a9cf",
}
SUBCLASS_ORDER = ["NMD+/PTC+", "NMD+/PTC-", "Control"]

MODEL_LABEL = {
    "nmdetective_b": "NMDetective-B\n(Lindeboom 2019)",
    "nmdep_rule":    "NMDEP rule baseline\n(Saadat 2025)",
    "our_model":     "Our model\n(this work)",
}
MODEL_ORDER = ["nmdetective_b", "nmdep_rule", "our_model"]
MODEL_SCORE_COL = {
    "nmdetective_b": "nmdetective_b_score",
    "nmdep_rule":    "nmdep_rule_score",
    "our_model":     "our_model_prob",
}


def load():
    s = pd.read_csv(SCORES, sep="\t")
    m = pd.read_csv(METRICS, sep="\t")
    # Head-to-head intersection: isoforms with all three model scores
    mask = (
        s["nmdetective_b_score"].notna()
        & s["nmdep_rule_score"].notna()
        & s["our_model_prob"].notna()
    )
    return s[mask].copy(), m


def panel_a_scatter(ax, df, model_key):
    """One subplot in Panel A: predicted score vs gold-standard log-FC."""
    col = MODEL_SCORE_COL[model_key]
    for sub in SUBCLASS_ORDER:
        d = df[df["subclass"] == sub]
        ax.scatter(
            d[col], d["mashr_posterior_mean_logfc"],
            color=GROUP_COLOR[sub], alpha=0.6, s=18,
            edgecolor="none", label=GROUP_LABEL[sub],
        )
    # Pooled Spearman
    sp = df[[col, "mashr_posterior_mean_logfc"]].corr(method="spearman").iloc[0, 1]
    ax.text(
        0.04, 0.95, f"Spearman = {sp:.2f}",
        transform=ax.transAxes, fontsize=10, color="#222222",
        va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#999", lw=0.5),
    )
    ax.set_title(MODEL_LABEL[model_key], fontsize=11, color="#222222", pad=4)
    ax.set_xlabel("Predicted NMD score", fontsize=10, color="#222222")
    if model_key == "our_model":
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    else:
        ax.set_xticks([0.0, 0.2, 0.4, 0.6])
    ax.tick_params(axis="both", labelsize=9, colors="#222222")
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)
    for sd in ("bottom", "left"):
        ax.spines[sd].set_color("#555555")
        ax.spines[sd].set_linewidth(0.8)


def panel_b_bars(ax, metrics, head_to_head_n):
    """Grouped bar chart of per-subclass Spearman within head-to-head set."""
    rows = []
    for model in MODEL_ORDER:
        for sub in SUBCLASS_ORDER:
            key = f"head-to-head:{sub}"
            r = metrics[(metrics["model"] == model) & (metrics["stratum"] == key)]
            if r.empty:
                continue
            rows.append({
                "model":   MODEL_LABEL[model].replace("\n", " "),
                "subclass": sub,
                "spearman": r.iloc[0]["spearman"],
                "n":         int(r.iloc[0]["n"]),
            })
    d = pd.DataFrame(rows)

    x = np.arange(len(MODEL_ORDER))
    w = 0.27
    for i, sub in enumerate(SUBCLASS_ORDER):
        vals_by_model = []
        for m in MODEL_ORDER:
            row = d[(d["model"] == MODEL_LABEL[m].replace("\n", " "))
                    & (d["subclass"] == sub)]
            vals_by_model.append(row.iloc[0]["spearman"] if not row.empty else np.nan)
        bars = ax.bar(
            x + (i - 1) * w, vals_by_model, width=w,
            color=GROUP_COLOR[sub], edgecolor="black", linewidth=0.5,
            label=f"{GROUP_LABEL[sub]} (n = {int(d.loc[d['subclass']==sub, 'n'].iloc[0])})",
        )
        # Bar labels
        for rect, v in zip(bars, vals_by_model):
            if np.isnan(v):
                continue
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                v + (0.03 if v >= 0 else -0.06),
                f"{v:.2f}",
                ha="center", va="bottom" if v >= 0 else "top",
                fontsize=8, color="#222222",
            )

    ax.axhline(0, color="#555555", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["NMDetective-B\n(Lindeboom 2019)",
         "NMDEP rule baseline\n(Saadat 2025)",
         "Our model\n(this work)"],
        fontsize=9, color="#222222",
    )
    ax.set_ylabel("Spearman correlation with\nmashr log fold-change (SMG1i)",
                  fontsize=10, color="#222222")
    ax.set_ylim(-0.35, 0.75)
    ax.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6])
    ax.tick_params(axis="y", labelsize=9, colors="#222222")
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)
    for sd in ("bottom", "left"):
        ax.spines[sd].set_color("#555555")
        ax.spines[sd].set_linewidth(0.8)
    ax.legend(
        loc="upper left", fontsize=9, frameon=False,
        handlelength=1.5, handletextpad=0.5, borderaxespad=0.4,
        title=f"Within-subclass (n = {head_to_head_n} total)",
        title_fontsize=9,
    )


def build_figure():
    df, metrics = load()
    n_hh = len(df)

    fig = plt.figure(figsize=(16, 5.5))
    gs_outer = fig.add_gridspec(
        nrows=1, ncols=2,
        width_ratios=[3.0, 1.55],
        wspace=0.38, left=0.055, right=0.985, bottom=0.17, top=0.84,
    )
    gs_a = gs_outer[0, 0].subgridspec(nrows=1, ncols=3, wspace=0.45)
    gs = [gs_a[0, 0], gs_a[0, 1], gs_a[0, 2], gs_outer[0, 1]]

    axes_a = [fig.add_subplot(gs[i]) for i in range(3)]
    ax_b = fig.add_subplot(gs[3])

    # Panel A: scatter per model
    for ax, m in zip(axes_a, MODEL_ORDER):
        panel_a_scatter(ax, df, m)
    axes_a[0].set_ylabel(
        "Gold standard: mashr posterior mean\nlog fold-change under SMG1i",
        fontsize=10, color="#222222",
    )
    # Shared legend at the bottom of the first subplot to avoid overlap
    axes_a[0].legend(
        loc="lower right", fontsize=8, frameon=False,
        handlelength=1.0, handletextpad=0.4, borderaxespad=0.3,
        title="Subclass", title_fontsize=8,
    )

    # Panel B: stratified bars
    panel_b_bars(ax_b, metrics, n_hh)

    # Panel labels
    for ax, letter in zip([axes_a[0], ax_b], ["A", "B"]):
        ax.text(
            -0.13, 1.06, letter, transform=ax.transAxes,
            fontsize=15, fontweight="bold", va="bottom", ha="left",
            color="#222222",
        )

    fig.text(
        0.5, 0.025,
        f"Comparison on the {n_hh}-isoform intersection of the n = 1,166 ref-AUG-traceable "
        "cohort with the deep-learning model's H5 universe (full-cohort predictions: "
        "train + val + test).",
        ha="center", va="bottom", fontsize=9, color="#555555", style="italic",
    )
    return fig, axes_a + [ax_b]


def main():
    fig, axes = build_figure()
    if validate_multipanel_layout is not None:
        try:
            validate_multipanel_layout(fig)
        except Exception as e:
            print(f"[validate_multipanel_layout] non-fatal: {e}")
    out = HERE.parent / "figure_s_model_comparison"
    fig.savefig(f"{out}.pdf", facecolor="white")
    fig.savefig(f"{out}.png", dpi=300, facecolor="white")
    print(f"Saved: {out}.pdf and .png")


if __name__ == "__main__":
    main()
