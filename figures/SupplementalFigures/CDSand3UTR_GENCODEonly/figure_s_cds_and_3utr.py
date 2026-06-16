"""
Supplemental Figure — CDS length and two complementary 3'UTR-length measures
across NMD+/PTC+, NMD+/PTC-, and Control isoforms at TWO scopes.

2×3 layout (18"×8.5"):
  Row 1 — Section A (n = 72/118/190, all-3-ENST + coding-CDS;
                     own-GENCODE-stop classifier + 50-nt rule):
    A  CDS length (own GENCODE CDS; nt, log10)
    B  3'UTR length, translation-based (own stop → tx end; nt, log10)
    C  3'UTR length, non-PTC-stop based (walk past PTCs; nt, log10)

  Row 2 — Section C (n = 1050/116/1166, ENST-reference + ref-AUG-traceable,
                     re-intersected on (gene, reference) after category filter;
                     mechanism_class_4 → 3-group merged):
    D  CDS length (ref-AUG-projected ORF; nt, log10)
    E  3'UTR length, translation-based (ref-AUG-projected stop → tx end)
    F  3'UTR length, non-PTC-stop based (walk past PTCs)

The Section A row is fully GENCODE-anchored (every comparator has a curated
CDS). The Section C row uses ref-AUG-projected measurements so the comparator's
ORF/3'UTR are not subject to TD2's PTC-avoidance bias (Figure 4 / §4).

Data files (in ./data/):
  features_190_subset_long.tsv                  - row 1
  features_1166_subset_long.tsv                 - row 2
  panelA_cds_length_{descriptives,pairwise}.tsv
  panelC_utr3_translation_{descriptives,pairwise}.tsv
  panelD_utr3_non_ptc_stop_{descriptives,pairwise}.tsv
  panelD_cds_length_refaug_{descriptives,pairwise}.tsv
  panelE_utr3_translation_refaug_{descriptives,pairwise}.tsv
  panelF_utr3_non_ptc_stop_refaug_{descriptives,pairwise}.tsv
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

BODY_FS = 14
LABEL_FS = 12
PANEL_LABEL_FS = 22
ROW_TITLE_FS = 16
C_TITLE = "#222222"
C_AXIS = "#555555"
C_LBL = "#111111"

GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC-", "Control"]
COLORS = {"NMD+/PTC+": "#ef8a62", "NMD+/PTC-": "#fee08b", "Control": "#92c5de"}

# Row 1 — Section A (n=190), own-GENCODE-stop classifier
ROW_A = {
    "title": "Section A  ·  n = 72 / 118 / 190 (all-3-ENST + coding-CDS; own-GENCODE-stop)",
    "long_tsv": "features_190_subset_long.tsv",
    "panels": [
        {"letter": "A", "slug": "panelA_cds_length",
         "value_col": "cds_nt",
         "ylabel": "CDS length (nt, log10)"},
        {"letter": "B", "slug": "panelC_utr3_translation",
         "value_col": "utr3_to_tx_end_nt",
         "ylabel": "3'UTR length, translation-based\n(own stop → tx end; nt, log10)"},
        {"letter": "C", "slug": "panelD_utr3_non_ptc_stop",
         "value_col": "utr3_via_non_ptc_stop_nt",
         "ylabel": "3'UTR length, non-PTC-stop based\n(walk past PTCs; nt, log10)"},
    ],
}

# Row 2 — Section C (n=1,166), ref-AUG-projected classifier
ROW_C = {
    "title": "Section C  ·  n = 1,050 / 116 / 1,166 (ENST-ref + ref-AUG-traceable, re-intersected; ref-AUG-projected)",
    "long_tsv": "features_1166_subset_long.tsv",
    "panels": [
        {"letter": "D", "slug": "panelD_cds_length_refaug",
         "value_col": "cds_nt",
         "ylabel": "CDS length (ref-AUG ORF; nt, log10)"},
        {"letter": "E", "slug": "panelE_utr3_translation_refaug",
         "value_col": "utr3_to_tx_end_nt",
         "ylabel": "3'UTR length, translation-based\n(ref-AUG stop → tx end; nt, log10)"},
        {"letter": "F", "slug": "panelF_utr3_non_ptc_stop_refaug",
         "value_col": "utr3_via_non_ptc_stop_nt",
         "ylabel": "3'UTR length, non-PTC-stop based\n(walk past PTCs; nt, log10)"},
    ],
}


def sig_marker(p):
    if pd.isna(p): return "n.s."
    if p < 1e-4:   return "***"
    if p < 1e-3:   return "**"
    if p < 0.05:   return "*"
    return "n.s."


def render_panel(ax, long, desc, stats, value_col, ylabel):
    plot_df = long.copy()
    plot_df["y"] = np.log10(np.clip(plot_df[value_col], 1, None))
    plot_df = plot_df.dropna(subset=["y"])

    sns.violinplot(data=plot_df, x="group", y="y", order=GROUP_ORDER,
                   hue="group", palette=COLORS, legend=False,
                   cut=0, inner="quart", density_norm="width",
                   linewidth=0.9, ax=ax)
    for patch in ax.collections:
        patch.set_edgecolor(C_TITLE); patch.set_alpha(0.78)
    for line in ax.lines:
        line.set_color(C_TITLE); line.set_linewidth(1.0)

    data = [plot_df.loc[plot_df["group"] == g, "y"].to_numpy() for g in GROUP_ORDER]
    ymax = max(d.max() for d in data)
    y_top = ymax + 0.18
    spacing = 0.50
    pairs = [("NMD+/PTC-", "NMD+/PTC+", 0, 1, y_top),
             ("NMD+/PTC-", "Control",  1, 2, y_top + spacing),
             ("NMD+/PTC+", "Control",  0, 2, y_top + 2 * spacing)]
    for gx, gy, ix, iy, yy in pairs:
        row = stats[((stats["group_x"] == gx) & (stats["group_y"] == gy)) |
                    ((stats["group_x"] == gy) & (stats["group_y"] == gx))]
        if row.empty: continue
        p = row.iloc[0]["wilcox_p"]
        ax.plot([ix, ix, iy, iy], [yy - 0.04, yy, yy, yy - 0.04],
                color=C_TITLE, linewidth=0.6)
        ax.text((ix + iy) / 2, yy + 0.02, sig_marker(p),
                ha="center", va="bottom", fontsize=BODY_FS, color=C_TITLE)

    n_per   = {g: int(desc.loc[desc["group"] == g, "n"].iloc[0])      for g in GROUP_ORDER}
    med_per = {g: float(desc.loc[desc["group"] == g, "median"].iloc[0]) for g in GROUP_ORDER}

    def _fmt_med(m):
        # Integer when median ≥ 10 (the typical nt-length case); one decimal
        # for very small medians (e.g. longest-uORF medians of 0 nt or 15 nt).
        return f"{int(round(m)):,}" if m >= 10 else f"{m:.1f}"

    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(
        [f"{g}\nn = {n_per[g]:,}\nmed = {_fmt_med(med_per[g])} nt"
         for g in GROUP_ORDER],
        fontsize=LABEL_FS, color=C_TITLE)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)
    ax.set_ylim(bottom=0, top=y_top + 3 * spacing)

    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color(C_AXIS); ax.spines[s].set_linewidth(0.8)


def render_row(axes_row, row_spec):
    long = pd.read_csv(HERE.parent / "data" / row_spec["long_tsv"], sep="\t")
    for spec, ax in zip(row_spec["panels"], axes_row):
        desc  = pd.read_csv(HERE.parent / "data" / f"{spec['slug']}_descriptives.tsv", sep="\t")
        stats = pd.read_csv(HERE.parent / "data" / f"{spec['slug']}_pairwise.tsv", sep="\t")
        render_panel(ax, long, desc, stats, spec["value_col"], spec["ylabel"])
        # Panel letter — positioned ABOVE the panel, flush with the left spine
        # (axes-x = 0), NOT outside the panel to the left. Placing it OUTSIDE
        # (axes-x < 0) puts it in the same x-column as the ylabel text and
        # visually crowds the rotated ylabel — caught by validate_multipanel_layout's
        # cross-panel structural-text crowding check.
        ax.text(0.0, 1.04, spec["letter"], transform=ax.transAxes,
                fontsize=PANEL_LABEL_FS, fontweight="bold", color=C_LBL,
                ha="left", va="bottom")


def build_figure():
    # Wider canvas (taller too) leaves room for:
    #   - panel letters A/B/C above row 1 without clipping the figure top
    #   - panel letters D/E/F above row 2 without colliding with row 1's
    #     xtick labels (n=… lines)
    fig, axes = plt.subplots(2, 3, figsize=(18, 9.5))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.10,
                        wspace=0.30, hspace=0.70)

    render_row(axes[0], ROW_A)
    render_row(axes[1], ROW_C)

    return fig


def main():
    fig = build_figure()

    # Composite-level layout validator — catches clip-at-edge and
    # cross-panel text overlap that the per-panel single-axes validator
    # cannot see (it only sees one axes at a time). Run on every render.
    from validate_figure_layout import validate_multipanel_layout
    result = validate_multipanel_layout(fig, verbose=True)
    if result["summary"]["n_errors"] > 0:
        raise SystemExit(
            f"Layout validation failed with {result['summary']['n_errors']} "
            f"error(s). Fix the layout before saving — see report above.")

    out = HERE.parent
    fig.savefig(out / "figure_s_cds_and_3utr.pdf", facecolor="white")
    fig.savefig(out / "figure_s_cds_and_3utr.png", dpi=300, facecolor="white")
    print("Saved: figure_s_cds_and_3utr.{pdf,png}")


if __name__ == "__main__":
    main()
