"""
SF37 — Stop codon usage in NMD susceptible vs non-NMD isoforms.

Native render (replaces the earlier legacy-PNG placeholder that showed DNA-form
labels TGA/TAA/TAG and pre-bug-fix frequencies). Within-class frequency of the
three canonical stop codons (UGA / UAA / UAG) at the main (rank-0) ORF of each
transcript in the deep-learning model's held-out test set.

Data: NMD_orf_model_v5_4ct/results_4ct/stop_codon_freq_by_class_sf37.tsv — the
canonical per-class frequencies from the `fig4a-stop-codon-freq` chunk of
orf_model_report_v5.Rmd (rank-0 ORF stop codon from selected_orfs.tsv joined to
test-set predictions; post-2026-04-30 stop-codon column bug fix). n = 2,268 NMD
susceptible + 7,863 Control = 10,131 test transcripts.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    docx_body_fs,
    render_and_validate,
    NMD_COLOR,
    CONTROL_COLOR,
    TITLE_C,
)
apply_ggplot_rcparams()

NATIVE_W = 6.5
BODY_FS = docx_body_fs(NATIVE_W)

DATA = Path("/Users/petecastaldi/claude_projects/NMD_orf_model_v5_4ct/"
            "results_4ct/stop_codon_freq_by_class_sf37.tsv")

CODON_ORDER = ["UGA", "UAA", "UAG"]


def main():
    df = pd.read_csv(DATA, sep="\t")

    def pcts(cls):
        s = df[df["class"] == cls].set_index("stop_codon")["pct"]
        return np.array([float(s[c]) for c in CODON_ORDER])

    nmd, ctrl = pcts("NMD"), pcts("Control")
    n_nmd = int(df.loc[df["class"] == "NMD", "class_total"].iloc[0])
    n_ctrl = int(df.loc[df["class"] == "Control", "class_total"].iloc[0])

    x = np.arange(len(CODON_ORDER))
    bw = 0.38

    fig, ax = plt.subplots(figsize=(NATIVE_W, 4.5))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.92, bottom=0.13)
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    b_nmd = ax.bar(x - bw / 2, nmd, width=bw, color=NMD_COLOR, edgecolor="white",
                   linewidth=0.6, zorder=3,
                   label=f"NMD susceptible (n = {n_nmd:,})")
    b_ctrl = ax.bar(x + bw / 2, ctrl, width=bw, color=CONTROL_COLOR,
                    edgecolor="white", linewidth=0.6, zorder=3,
                    label=f"Control (n = {n_ctrl:,})")

    for bars, vals in [(b_nmd, nmd), (b_ctrl, ctrl)]:
        for rect, v in zip(bars, vals):
            ax.text(rect.get_x() + rect.get_width() / 2, v + 1.0, f"{v:.1f}%",
                    ha="center", va="bottom", fontsize=BODY_FS - 2, color=TITLE_C)

    ax.set_xticks(x)
    ax.set_xticklabels(CODON_ORDER, fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylim(0, 66)
    ax.set_yticks([0, 20, 40, 60])
    ax.set_xlabel("Stop codon", fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Within-class frequency (%)", fontsize=BODY_FS, color=TITLE_C)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=TITLE_C)
    ax.legend(loc="upper right", frameon=True, facecolor="white",
              edgecolor="none", fontsize=BODY_FS - 2)

    render_and_validate(fig, HERE.parent / "figure_s_stop_codon_usage",
                        native_width_in=NATIVE_W)
    print(f"NMD n={n_nmd:,}, Control n={n_ctrl:,}")
    print(f"  UGA {nmd[0]}/{ctrl[0]}  UAA {nmd[1]}/{ctrl[1]}  UAG {nmd[2]}/{ctrl[2]}")


if __name__ == "__main__":
    main()
