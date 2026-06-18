"""
Supplemental Figure — Attention distribution: NMD vs Control.

Companion to the manuscript sentence in §5:

   "Attention analysis demonstrated that, while most of the attention was
   focused on ORF0 (the most likely CDS), attention was more broadly
   distributed for NMD than Control isoforms (SFx, panel A and B)."

Two panels:
  A — Attention weight per ORF rank (boxplots, NMD vs Control). Most
      mass is at rank 0 (the model's priority ORF — ref CDS > TD2 CDS >
      top-Kozak). NMD has a slightly lower median attention at rank 0
      and slightly more mass at ranks ≥1 vs Control.
  B — Shannon entropy of the per-isoform attention vector across the 5
      candidate ORFs. NMD's distribution is shifted toward higher entropy
      (more broadly distributed attention) compared to Control.

Both panels are extracted from the legacy deep_nmd_model report
(`fig5b_attention_by_rank.png`, `fig5a_entropy.png`). The underlying
computations don't depend on the stop-codon bug fix and are therefore
canonical as-is. When the cluster is back, the report's `fig5*` chunks
can be re-rendered natively against the latest model outputs.
"""

import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

LBL_FS = 18
LBL_COLOR = "#111111"

# Each panel is 2400×1500 native (aspect 1.6:1). Two side-by-side
# panels at cell 5.5 × 3.4 in (1.6:1).
CELL_W = 5.5
CELL_H = 3.4
FIG_W = CELL_W * 2     # 11.0
FIG_H = CELL_H         # 3.4

PANELS = {
    "A": "data/legacy_fig5b_attention_by_rank.png",
    "B": "data/legacy_fig5a_entropy.png",
}


def add_panel(fig, gs_slot, letter, paths):
    ax = fig.add_subplot(gs_slot)
    img = mpimg.imread(str(paths[letter]))
    ax.imshow(img, aspect="auto")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(
        0.0, 1.01, letter,
        transform=ax.transAxes,
        fontsize=LBL_FS, fontweight="bold", color=LBL_COLOR,
        ha="left", va="bottom",
    )


def build_figure():
    paths = {k: HERE.parent / v for k, v in PANELS.items()}

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    gs = GridSpec(
        1, 2, figure=fig,
        left=0.01, right=0.99, top=0.90, bottom=0.02,
        wspace=0.04,
    )
    add_panel(fig, gs[0, 0], "A", paths)
    add_panel(fig, gs[0, 1], "B", paths)
    return fig


def main():
    fig = build_figure()

    from validate_figure_layout import validate_multipanel_layout
    validate_multipanel_layout(fig, verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure_s_attention_distribution.pdf", facecolor="white")
    fig.savefig(out_dir / "figure_s_attention_distribution.png", dpi=300, facecolor="white")
    print(f"Saved: figure_s_attention_distribution.pdf and .png  ({FIG_W:.1f} × {FIG_H:.1f} in)")


if __name__ == "__main__":
    main()
