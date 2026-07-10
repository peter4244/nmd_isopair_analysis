"""SF26 — Reference-isoform share of gene expression across pop_BC (n = 3,009 genes).

Standalone rebuild of the reference-share panel that was previously bundled
as Panel B of PairSetDescriptives. Split per Yul-era paper numbering so
the paper's SF26 reference resolves to one figure.

Value plotted per gene: reference-isoform DMSO 4-CT mean expression divided
by the parent gene's total non-NMD expression (both from the Isopair
pipeline; see Methods, pop_BC / Reference isoform).

Style: matplotlib rendered with ggplot-mimic theme (grey panel + white
gridlines) so the panel visually matches SF1-SF23. See
figures/lib/ggplot_style.py.
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LIB  = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    render_and_validate,
    docx_body_fs,
    BAR_COLOR,
    TITLE_C,
    HEADER_FS,
)

apply_ggplot_rcparams()

NATIVE_W = 6.5
BODY_FS = docx_body_fs(NATIVE_W)

DATA = HERE / "data"


def main():
    fra = pd.read_csv(DATA / "ref_expression_fraction.tsv", sep="\t")
    summary = pd.read_csv(DATA / "descriptives_summary.tsv", sep="\t").set_index("metric")["value"]
    n_genes = int(summary.loc["n_genes_in_pop_BC"])

    vals = (fra["ref_fraction_of_gene"].dropna() * 100).values  # → %
    med = float(np.median(vals))
    frac_ge_50 = float((vals >= 50).mean() * 100)

    fig, ax = plt.subplots(figsize=(NATIVE_W, 4.0))
    fig.subplots_adjust(left=0.11, right=0.95, top=0.88, bottom=0.14)

    style_axes_ggplot(ax)

    bins = np.arange(0, 105, 5)
    ax.hist(vals, bins=bins, color=BAR_COLOR, edgecolor="white", linewidth=0.6, zorder=3)

    ax.axvline(med, linestyle="--", color="#c0392b", linewidth=1.5, zorder=4)
    ax.axvline(50, linestyle=":",  color="#aa6600", linewidth=1.0, zorder=4)
    ymax = ax.get_ylim()[1]
    ax.text(
        med + 1.5, ymax * 0.94,
        f"median = {med:.1f}%",
        ha="left", va="top",
        fontsize=BODY_FS,
        color="#c0392b",
    )
    ax.text(
        50 - 1.5, ymax * 0.82,
        "50% threshold",
        ha="right", va="top",
        fontsize=BODY_FS,
        color="#aa6600",
    )

    ax.set_xlabel("Reference share of gene expression (%)", fontsize=BODY_FS)
    ax.set_ylabel("Number of genes", fontsize=BODY_FS)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])

    # No overall figure title — caption carries the title role (Yul-style).
    render_and_validate(fig, HERE / "figure_sf26_reference_share",
                        native_width_in=NATIVE_W)
    print(f"  n_genes = {n_genes:,}  median = {med:.1f}%  ≥50% = {frac_ge_50:.1f}%")


if __name__ == "__main__":
    main()
