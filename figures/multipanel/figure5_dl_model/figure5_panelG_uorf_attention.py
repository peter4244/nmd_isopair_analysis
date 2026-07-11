"""
Figure 5 Panel G — uORF attention proportion at the Subset 2 scope
(4-CT re-scope 2026-07-11: n=819 per arm; "n1166" filenames are the stable
subset-2 identifier).

The attention metric is computed under the **Path B-strict** rule
(TD2-bias-independent, structural-only identification of candidate uORFs):

    is_candidate_uorf =
      orf_length < 200 &                                # below GENCODE annotated
                                                        #   CDS floor of 270 nt
      frac_position < 0.5 &                             # start in 5' half
      (orf_start + orf_length) / tx_length < 0.6 &      # stop 5'-leaning
      kozak_score >= 1                                  # credible ATG

per-isoform uORF attention = sum(attention[is_candidate_uorf])

Subgroup structure at this scope (Subset 2 / mechanism_class_4 from Figure 4
Section C; ref-AUG-traceable, gene-matched 1:1):

    NMD+/PTC+          = category == "effectively_ptc"   (n=735 with model)
    NMD+/PTC- retained = category ∈ {no_downstream_ejc,
                                     truncated_no_ejc}    (n=54 with model)
    Control             = all 819 c4 comparators          (n=781 with model)

NMD+/PTC- *lost* (category ∈ ref_atg_lost / no_ref_isoform etc.) is
EXCLUDED from this scope by construction — those isoforms don't have a
projectable ref AUG and we therefore can't validate uORF calls against an
independent reference. See methodology.

Distribution is heavily zero-inflated for PTC+ and Control (they shouldn't
need uORFs to explain NMD trajectory). We render a violin overlaid with
strip plot so the zero mass is honest.

Data:
  ./data/panelG_uorf_attention_n1166.tsv  (per-isoform: group, v1, strict_attn)
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
C_TITLE = "#222222"
C_AXIS = "#555555"

GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC- retained", "Control"]
GROUP_DISPLAY = {
    "NMD+/PTC+": "NMD+/PTC+",
    "NMD+/PTC- retained": "NMD+/PTC−\nretained",
    "Control": "Control",
}
COLORS = {
    "NMD+/PTC+":          "#ef8a62",   # coral
    "NMD+/PTC- retained": "#fee08b",   # yellow
    "Control":            "#92c5de",   # light blue
}


from p_to_stars import p_to_stars as sig_marker  # noqa: E402 — unified 4-tier scheme


def lookup_p(stats, ga, gb):
    """Pull pre-computed Wilcoxon p from R-side TSV."""
    row = stats[((stats["group_x"] == ga) & (stats["group_y"] == gb)) |
                ((stats["group_x"] == gb) & (stats["group_y"] == ga))]
    if row.empty:
        return float("nan")
    return float(row.iloc[0]["wilcox_p"])


def build_figure():
    df = pd.read_csv(HERE.parent / "data" / "panelG_uorf_attention_n1166.tsv",
                     sep="\t")
    stats = pd.read_csv(HERE.parent / "data" / "panelG_uorf_attention_n1166_pairwise.tsv",
                        sep="\t")
    # Keep only isoforms with model coverage
    df = df[df["uorf_attention_frac"].notna()].copy()
    df = df[df["group"].isin(GROUP_ORDER)].copy()

    fig, ax = plt.subplots(figsize=(6, 4))

    sns.violinplot(
        data=df, x="group", y="strict_attn", order=GROUP_ORDER,
        hue="group", palette=COLORS, legend=False,
        cut=0, inner=None, density_norm="width",
        linewidth=0.9, ax=ax,
    )
    for patch in ax.collections:
        patch.set_edgecolor(C_TITLE)
        patch.set_alpha(0.55)

    # Strip overlay — point size + alpha scaled by group n for density balance
    n_per = df.groupby("group").size().to_dict()
    sns.stripplot(
        data=df, x="group", y="strict_attn", order=GROUP_ORDER,
        hue="group", palette=COLORS, legend=False,
        size=1.8, alpha=0.45, jitter=0.20,
        edgecolor=C_TITLE, linewidth=0.15, ax=ax,
    )

    # Median markers
    medians = df.groupby("group")["strict_attn"].median().to_dict()
    for i, g in enumerate(GROUP_ORDER):
        med = medians.get(g, 0.0)
        ax.hlines(med, i - 0.30, i + 0.30, colors=C_TITLE,
                  linewidth=1.8, zorder=5)

    # Significance brackets with Wilcoxon p
    ymax = float(df["strict_attn"].max())
    y_top = ymax * 1.10
    spacing = ymax * 0.16
    tick = ymax * 0.018
    text_off = ymax * 0.030

    pair_specs = [
        ("NMD+/PTC+", "NMD+/PTC- retained", 0, 1, y_top),
        ("NMD+/PTC- retained", "Control",   1, 2, y_top + spacing),
        ("NMD+/PTC+", "Control",            0, 2, y_top + 2 * spacing),
    ]
    for ga, gb, xa, xb, yy in pair_specs:
        p = lookup_p(stats, ga, gb)
        marker = sig_marker(p)
        ax.plot([xa, xa, xb, xb], [yy - tick, yy, yy, yy - tick],
                color=C_TITLE, linewidth=0.7)
        ax.text((xa + xb) / 2, yy + text_off, marker,
                ha="center", va="bottom", fontsize=BODY_FS, color=C_TITLE)

    # X tick labels with n
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(
        [f"{GROUP_DISPLAY[g]}\n(n={n_per.get(g,0):,})" for g in GROUP_ORDER],
        fontsize=LABEL_FS, color=C_TITLE,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Attention fraction on uORFs",
                  fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylim(-0.02, y_top + 2 * spacing + 4 * text_off)
    ax.tick_params(axis="y", labelsize=LABEL_FS, colors=C_TITLE)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    fig.subplots_adjust(left=0.15, right=0.97, bottom=0.20, top=0.96)
    return fig


def main():
    fig = build_figure()

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelG_uorf_attention.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelG_uorf_attention.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelG_uorf_attention.pdf and .png")


if __name__ == "__main__":
    main()
