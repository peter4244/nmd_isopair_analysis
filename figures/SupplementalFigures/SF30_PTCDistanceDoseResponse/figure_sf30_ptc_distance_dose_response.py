"""SF30 — NMD-response magnitude vs stop-codon distance to the last EJC.

Scatter of mashr posterior-mean logFC (SMG1i vs DMSO, averaged across the
four cell types AT2, LAE, FB, MV) against distance from the comparator's
stop codon to the last exon-junction complex. Loess overlay + 95% CI band.
Vertical dashed line at 50 nt marks the operational PTC threshold.

Data source: `sf30_ptc_distance_logfc.tsv` (produced by data_export.R from
the Isopair pipeline's canonical gene-matched C2 comparators + ptc.rds +
mashr_isoform_model_2026.3.10.rds). Matches the population and derivation
of the Rmd chunk `goal1-fig2-logfc-dist` in
`05_final_report_mashr.Rmd` § "NMD Response vs Stop Codon Distance".
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
LIB  = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    assert_text_within_canvas,
    NMD_COLOR,
    TITLE_C,
    BODY_FS,
    HEADER_FS,
)

apply_ggplot_rcparams()

DATA = HERE / "data"


def _tricube(u):
    u = np.abs(u)
    w = np.where(u < 1, (1 - u ** 3) ** 3, 0.0)
    return w


def loess_with_ci(x, y, frac=0.30, n_boot=200, seed=42, n_grid=200):
    """Simple local-linear (LOESS-like) fit with bootstrap 95% CI.

    Iterates over a fixed x-grid; at each grid point, fits a weighted linear
    regression to the nearest `k = ceil(frac * n)` observations, weighted by
    tricube of scaled distance. This reproduces the shape of ggplot's
    `geom_smooth(method="loess")` closely enough for a §4 supplemental
    figure — the point is the dose-response curve, not the exact smoother.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    k = max(int(np.ceil(frac * n)), 8)

    x_grid = np.linspace(np.quantile(x, 0.005), np.quantile(x, 0.995), n_grid)

    def _fit(xx, yy):
        out = np.empty_like(x_grid)
        for i, xg in enumerate(x_grid):
            d = np.abs(xx - xg)
            # local-window bandwidth = distance to k-th nearest
            if k >= len(xx):
                h = d.max() if d.max() > 0 else 1.0
            else:
                h = np.partition(d, k - 1)[k - 1]
                if h == 0:
                    h = np.max(d[d > 0]) if np.any(d > 0) else 1.0
            w = _tricube(d / h)
            if w.sum() == 0:
                out[i] = np.nan
                continue
            # Weighted linear regression via normal equations
            xw = xx - xg  # center on grid point
            sw   = w.sum()
            swx  = (w * xw).sum()
            swy  = (w * yy).sum()
            swxx = (w * xw * xw).sum()
            swxy = (w * xw * yy).sum()
            det = sw * swxx - swx * swx
            if det == 0:
                out[i] = swy / sw
            else:
                # intercept at xg
                out[i] = (swxx * swy - swx * swxy) / det
        return out

    fit_main = _fit(x, y)
    boots = np.empty((n_boot, n_grid))
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = _fit(x[idx], y[idx])
    lo = np.nanquantile(boots, 0.025, axis=0)
    hi = np.nanquantile(boots, 0.975, axis=0)
    return x_grid, fit_main, lo, hi


def main():
    df = pd.read_csv(DATA / "sf30_ptc_distance_logfc.tsv", sep="\t")
    df = df.dropna(subset=["ptc_distance", "mean_logFC"])

    # x-axis clip 1st/99th centile (matches Rmd)
    x_lo, x_hi = df["ptc_distance"].quantile([0.01, 0.99]).values
    plot_df = df[(df["ptc_distance"] >= x_lo) & (df["ptc_distance"] <= x_hi)].copy()
    n_plot = len(plot_df)

    rho, p = spearmanr(df["ptc_distance"], df["mean_logFC"])

    x = plot_df["ptc_distance"].to_numpy()
    y = plot_df["mean_logFC"].to_numpy()

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.86, bottom=0.14)
    style_axes_ggplot(ax)

    ax.scatter(x, y, s=6, alpha=0.18, color="#3a3a3a",
                edgecolors="none", zorder=3, rasterized=True)

    xg, yfit, lo, hi = loess_with_ci(x, y, frac=0.30)
    ax.fill_between(xg, lo, hi, color=NMD_COLOR, alpha=0.30, linewidth=0, zorder=4)
    ax.plot(xg, yfit, color=NMD_COLOR, linewidth=2.0, zorder=5)

    ax.axvline(50, linestyle="--", color="#c0392b", linewidth=1.0, zorder=2)
    ax.text(52, ax.get_ylim()[1] * 0.98, "PTC threshold (50 nt)",
            ha="left", va="top", fontsize=BODY_FS - 1, color="#c0392b")

    ax.set_xlabel("Distance from stop codon to last EJC (nt)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Mean logFC (mashr posterior mean)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_xlim(x_lo, x_hi)

    # Correlation annotation
    p_str = "< 10⁻¹⁶⁰" if p < 1e-160 else f"= {p:.2g}"
    ax.text(0.98, 0.04,
             f"Spearman ρ = {rho:.3f}\n"
             f"n = {n_plot:,} shown / {len(df):,} total\n"
             f"p {p_str}",
             transform=ax.transAxes, ha="right", va="bottom",
             fontsize=BODY_FS - 1, color=TITLE_C,
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.85,
                       boxstyle="square,pad=0.4"))

    fig.text(
        0.5, 0.94,
        "NMD response magnitude vs stop-codon distance to last EJC",
        ha="center", va="center",
        fontsize=HEADER_FS, fontweight="bold", color=TITLE_C,
    )

    assert_text_within_canvas(fig)

    out_png = HERE / "figure_sf30_ptc_distance_dose_response.png"
    out_pdf = HERE / "figure_sf30_ptc_distance_dose_response.pdf"
    fig.savefig(out_png, dpi=200, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    print(f"wrote {out_png.name} and {out_pdf.name}")
    print(f"  n comparators total: {len(df):,}")
    print(f"  n shown (1st-99th centile): {n_plot:,}")
    print(f"  Spearman rho = {rho:.4f}   p = {p:.2e}")


if __name__ == "__main__":
    main()
