"""
Supplemental Figure — TransDecoder2 (TD2) ORF-call bias on the §4 universe.

Three observations of TD2's PTC-avoidance behavior at two complementary scopes,
both on the same isoform universe as Figure 4 Panels C/D (2×3 layout,
18"×10" landscape):

Row 1 (broad scope, n = 1,166) — 1:1 gene-matched NMD comparator pairs in the
Figure 4 Section C universe: ENST-reference, ref-AUG-traceable categories,
TD2 CDS available, reference AUG exonic, re-intersected on
(gene_id, reference_isoform_id) with the Control side after the category
filter. 50% have TD2 == ref AUG (trivially match), 50% have TD2 ≠ ref AUG
(contribute the bias signal).

  A  TD2-selected CDS ORF length vs Reference-AUG-anchored ORF length (KDE)
  B  Paired Kozak PWM at the reference AUG vs at the TD2-predicted ATG (violin)
  C  TD2 ATG position relative to the reference AUG (histogram, spike at 0 = TD2==ref)

Row 2 (occult-PTC subset, n = 492) — the subset of Row 1 in which reference-AUG
tracing revealed a downstream-EJC stop indicative of a premature termination
codon that TD2's own CDS call missed (`effectively_ptc` ∩ `original_ptc ==
FALSE`). In this subset TD2 ATG ≠ ref AUG by construction; the bias is most
pronounced.

  D  TD2 vs ref-AUG ORF length on n = 492 occult-PTC subset (same as A)
  E  Paired Kozak PWM on the same n = 492 (same as B)
  F  TD2 ATG position relative to ref AUG on the same n = 492 (same as C)

Headline comparison: in the n=1,166 broad scope TD2 picks 4.69× longer ORFs
on average and 50% of pairs have TD2 agreeing with the reference; in the
occult-PTC subset TD2 picks 7.73× longer ORFs and 99.0% are downstream of the
reference. The bias is real across the full population but most extreme in
the population where TD2's PTC-avoidance is operationally engaged.

Data (in ./data/):
  panelA_td2_vs_refaug_length.tsv          # broad (n=1,166)
  panelB_kozak.tsv                          # broad (n=1,166)
  panelC_td2_position.tsv                   # broad (n=1,166)
  panelD_td2_vs_refaug_length_occult.tsv    # occult-PTC (n=492)
  panelE_kozak_occult.tsv                   # occult-PTC (n=492)
  panelF_td2_position_occult.tsv            # occult-PTC (n=492)
  summary.tsv
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    assert_text_within_canvas,
    panel_letter,
    TITLE_C,
    AXIS_C,
    PANEL_LETTER_FS,
)

apply_ggplot_rcparams()

BODY_FS = 12
LABEL_FS = 11
C_TITLE = TITLE_C
C_AXIS = AXIS_C
C_LBL = TITLE_C

C_REF_FILL  = "#ef8a62"; C_REF_LINE  = "#d6604d"
C_TD2_FILL  = "#67a9cf"; C_TD2_LINE  = "#2b8cbe"
C_REF_ATG = "#4393c3"
C_TD2_ATG = "#ef8a62"
C_HIST_FILL = "#67a9cf"
C_HIST_LINE = "#2b8cbe"
C_REF_LINE_VLINE = "#d6604d"


def render_length_kde(ax, tsv_name, n_label=None):
    style_axes_ggplot(ax, xgrid=True, ygrid=False)
    df = pd.read_csv(HERE.parent / "data" / tsv_name, sep="\t")
    x_grid = np.linspace(np.log10(50), np.log10(15000), 400)
    peaks = {}
    for label, fill, line in [
        ("Reference ATG ORF", C_REF_FILL, C_REF_LINE),
        ("TD2-selected CDS",  C_TD2_FILL, C_TD2_LINE),
    ]:
        vals = np.log10(df.loc[df["ORF_source"] == label, "length_nt"].to_numpy())
        kde = gaussian_kde(vals, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(10 ** x_grid, y, alpha=0.55, color=fill, linewidth=0, zorder=3)
        ax.plot(10 ** x_grid, y, color=line, linewidth=1.6, zorder=4)
        peak_i = int(np.argmax(y))
        peaks[label] = (10 ** x_grid[peak_i], y[peak_i], line, len(vals))
    ax.set_xscale("log")
    ax.set_xlim(50, 15000)
    y_max = max(p[1] for p in peaks.values())
    ax.set_ylim(0, y_max * 1.28)
    ax.set_xlabel("ORF length (nt, log scale)", fontsize=BODY_FS, color=C_TITLE)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=BODY_FS, colors=C_TITLE)
    for label, (xp, yp, line_c, n) in peaks.items():
        ax.text(xp, yp + y_max * 0.07,
                f"{label}\n(n={n:,})",
                ha="center", va="bottom",
                fontsize=BODY_FS - 1, color=line_c, fontweight="bold",
                zorder=5)


def render_kozak_violin(ax, tsv_name, pvalue_text):
    style_axes_ggplot(ax, xgrid=False, ygrid=True)
    df = pd.read_csv(HERE.parent / "data" / tsv_name, sep="\t")
    sources = ["Reference ATG", "TD2-predicted ATG"]
    palette = {"Reference ATG": C_REF_ATG, "TD2-predicted ATG": C_TD2_ATG}
    n_per = {s: int((df["ATG_source"] == s).sum()) for s in sources}

    # Directional pair count: how often is ref Kozak > TD2 Kozak?
    ref = df[df["ATG_source"] == "Reference ATG"][["comparator_isoform_id", "score"]]\
            .rename(columns={"score": "ref_score"})
    td2 = df[df["ATG_source"] == "TD2-predicted ATG"][["comparator_isoform_id", "score"]]\
            .rename(columns={"score": "td2_score"})
    paired = ref.merge(td2, on="comparator_isoform_id")
    n_pairs = len(paired)
    n_ref_higher = int((paired["ref_score"] > paired["td2_score"]).sum())
    pct_ref_higher = 100 * n_ref_higher / n_pairs

    sns.violinplot(data=df, x="ATG_source", y="score", order=sources,
                   hue="ATG_source", palette=palette, legend=False,
                   cut=0, inner="quart", density_norm="width",
                   linewidth=0.9, ax=ax)
    for patch in ax.collections:
        patch.set_edgecolor(C_TITLE); patch.set_alpha(0.85); patch.set_zorder(3)
    for line in ax.lines:
        line.set_color(C_TITLE); line.set_linewidth(1.0); line.set_zorder(4)

    data = [df.loc[df["ATG_source"] == s, "score"].to_numpy() for s in sources]
    ymax = max(np.max(d) for d in data); ymin = min(np.min(d) for d in data)
    y_br = ymax + 0.6
    ax.plot([0, 0, 1, 1], [y_br - 0.15, y_br, y_br, y_br - 0.15],
            color=C_TITLE, linewidth=0.8)
    ax.text(0.5, y_br + 0.1, f"paired Wilcoxon\n{pvalue_text}",
            ha="center", va="bottom", fontsize=BODY_FS - 2, color=C_TITLE)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{s}\n(n={n_per[s]:,})" for s in sources],
                       fontsize=LABEL_FS, color=C_TITLE)
    ax.set_xlabel("")
    ax.set_ylabel("Kozak PWM score (log-odds)", fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)
    ax.set_ylim(bottom=ymin - 0.3, top=y_br + 1.4)

    # Directional-pairs annotation in the lower-left of the violin panel
    ax.text(0.03, 0.04,
            f"ref AUG Kozak > TD2 Kozak\nin {n_ref_higher:,} / {n_pairs:,} pairs ({pct_ref_higher:.1f}%)",
            transform=ax.transAxes, fontsize=BODY_FS - 3, color=C_TITLE,
            ha="left", va="bottom",
            bbox=dict(facecolor="white", edgecolor=C_AXIS, alpha=0.9,
                      pad=3, linewidth=0.4))


def render_position_hist(ax, tsv_name):
    style_axes_ggplot(ax, xgrid=True, ygrid=True)
    df = pd.read_csv(HERE.parent / "data" / tsv_name, sep="\t")
    n = len(df)
    n_eq   = int((df["distance_nt"] == 0).sum())
    n_down = int((df["distance_nt"] > 0).sum())
    pct_eq   = 100 * n_eq / n
    pct_down = 100 * n_down / n
    pos = df.loc[df["distance_nt"] > 0, "distance_nt"]
    median_pos = float(np.median(pos)) if len(pos) > 0 else 0.0

    x_min, x_max = -500, 2500
    clipped = int(((df["distance_nt"] < x_min) | (df["distance_nt"] > x_max)).sum())
    vals = df["distance_nt"].clip(x_min, x_max)
    bins = np.linspace(x_min, x_max, 51)
    counts, edges, _ = ax.hist(vals, bins=bins, color=C_HIST_FILL,
                                edgecolor=C_HIST_LINE, linewidth=0.6, alpha=0.90,
                                zorder=3)

    ymax = counts.max() * 1.42
    ax.set_xlim(x_min, x_max); ax.set_ylim(0, ymax)
    ax.axvline(0, color=C_REF_LINE_VLINE, linewidth=1.4, linestyle="--", zorder=3)
    if median_pos > 0:
        ax.axvline(median_pos, color=C_TITLE, linewidth=1.0, linestyle=":")

    label_x = x_max - 80
    ax.text(label_x, ymax * 0.96, "ref AUG = 0",
            color=C_REF_LINE_VLINE, fontsize=BODY_FS - 2,
            ha="right", va="top", fontweight="bold")
    line_y = 0.86
    if n_eq > 0:
        ax.text(label_x, ymax * line_y,
                f"TD2 == ref AUG: {n_eq:,} / {n:,} ({pct_eq:.0f}%)",
                color=C_TITLE, fontsize=BODY_FS - 3, ha="right", va="top")
        line_y -= 0.08
    ax.text(label_x, ymax * line_y,
            f"TD2 downstream: {n_down:,} / {n:,} ({pct_down:.1f}%)",
            color=C_TITLE, fontsize=BODY_FS - 3, ha="right", va="top")
    line_y -= 0.08
    ax.text(label_x, ymax * line_y,
            f"median downstream offset: +{median_pos:.0f} nt",
            color=C_TITLE, fontsize=BODY_FS - 3, ha="right", va="top")

    if clipped > 0:
        ax.text(label_x, ymax * (line_y - 0.10),
                f"({clipped} outliers > +2,500 nt clipped)",
                fontsize=BODY_FS - 4, color=C_AXIS,
                ha="right", va="top", style="italic")

    ax.set_xlabel("TD2 ATG position relative to reference AUG (nt)",
                  fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylabel("Number of NMD comparator isoforms",
                  fontsize=BODY_FS - 1, color=C_TITLE)
    ax.tick_params(axis="both", labelsize=BODY_FS - 1, colors=C_TITLE)


def build_figure():
    fig = plt.figure(figsize=(18, 10))
    from matplotlib.gridspec import GridSpec
    # Leave bands above each row for a row-title; bands big enough so panel
    # labels (A/B/C and D/E/F) and the row title don't stack on top of each other.
    gs = GridSpec(2, 3, figure=fig,
                  left=0.06, right=0.99, top=0.86, bottom=0.07,
                  wspace=0.34, hspace=0.95)
    axes = np.array([[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)])

    render_length_kde(axes[0, 0], "panelA_td2_vs_refaug_length.tsv")
    render_kozak_violin(axes[0, 1], "panelB_kozak.tsv", "p = 1.2×10$^{-38}$")
    render_position_hist(axes[0, 2], "panelC_td2_position.tsv")

    render_length_kde(axes[1, 0], "panelD_td2_vs_refaug_length_occult.tsv")
    render_kozak_violin(axes[1, 1], "panelE_kozak_occult.tsv", "p = 2.6×10$^{-38}$")
    render_position_hist(axes[1, 2], "panelF_td2_position_occult.tsv")

    letters = [["A", "B", "C"], ["D", "E", "F"]]
    for ri in range(2):
        for ci in range(3):
            panel_letter(axes[ri, ci], letters[ri][ci], x=-0.20, y=1.04)

    # No row-group titles — caption identifies the two cohorts (Yul-style).
    return fig


def main():
    fig = build_figure()
    # Systematic canvas-bounds guard (see feedback_nmd_figure_canvas_bounds_guard).
    assert_text_within_canvas(fig)
    out = HERE.parent
    fig.savefig(out / "figure_s_td2_bias.pdf", facecolor="white")
    fig.savefig(out / "figure_s_td2_bias.png", dpi=300, facecolor="white")
    print("Saved: figure_s_td2_bias.{pdf,png}")


if __name__ == "__main__":
    main()
