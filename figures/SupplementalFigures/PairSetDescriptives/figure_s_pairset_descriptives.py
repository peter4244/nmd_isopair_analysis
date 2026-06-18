"""
Supplemental Figure — Descriptive characterisation of the n = 3,009 gene-
matched Isopair isoform-pair set (pop_BC).

Three panels documenting §4¶1–¶2 dataset descriptives:

  (A) Isoforms per gene, across the 3,009 genes in pop_BC. Counted as in
      the Isopair pipeline (05_final_report_mashr.Rmd §2): non-NMD
      isoforms with an expression-matrix entry. Median = 7. (Manuscript
      §4¶1 cites "median 7 isoforms per gene" — verified exactly.)

  (B) Reference-isoform DMSO 4-CT mean expression as a fraction of the
      parent gene's total non-NMD expression. Computed median = 30.8%;
      37.6% of references ≥ 50%. (Manuscript §4¶1 cites "median 70%" and
      "75% above 50%" — does NOT match; see README for the discrepancy
      note and the find/replace recommendation in the reconciliation
      report.)

  (C) Transcript length distribution by role within the 3,009 pairs:
      NMD comparator (median 3,049 nt), Reference (2,893 nt), Control
      comparator (2,762 nt). Manuscript §4¶2 cites these exact medians —
      verified.

Conventions follow project_nmd_figure_palette: NMD = #ef8a62, Control =
#67a9cf; Arial / 14pt body / 18pt header / spine #555.
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

BODY_FS = 12
HEADER_FS = 16

# Canonical palette — reference is the high-expressed non-NMD anchor;
# comparators are the contrasted isoforms (NMD = peach, Control = blue).
C_NMD     = "#ef8a62"
C_CONTROL = "#67a9cf"
C_REF     = "#7f7f7f"          # gray — reference is structurally distinct
C_BAR     = "#4292c6"          # blue for descriptive histograms (matches
                                # the upstream Rmd's fig0a histogram fill)
C_TITLE   = "#222222"
C_AXIS    = "#555555"

DATA = HERE.parent / "data"


def load():
    iso = pd.read_csv(DATA / "isoforms_per_gene.tsv", sep="\t")
    fra = pd.read_csv(DATA / "ref_expression_fraction.tsv", sep="\t")
    txl = pd.read_csv(DATA / "tx_length_by_role_long.tsv", sep="\t")
    desc = pd.read_csv(DATA / "descriptives_summary.tsv", sep="\t")
    return iso, fra, txl, desc


def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(axis="both", labelsize=BODY_FS, colors=C_TITLE)


def panel_A(ax, iso: pd.DataFrame):
    """Histogram of isoforms per gene, clipped at 30 (long tail)."""
    vals = iso["n_isoforms"].clip(upper=30)
    bins = np.arange(0, 32, 1)
    ax.hist(vals, bins=bins, color=C_BAR, edgecolor="white", linewidth=0.6,
            alpha=0.9)
    med = int(iso["n_isoforms"].median())
    ax.axvline(med, color="#c0504d", linestyle="--", linewidth=1.4,
               label=f"median = {med}")
    ax.set_xlabel("Isoforms per gene", fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylabel("Number of genes", fontsize=BODY_FS, color=C_TITLE)
    ax.set_xlim(-0.5, 30.5)
    ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
    ax.legend(loc="upper right", fontsize=BODY_FS - 2, frameon=False,
              handlelength=1.5, handletextpad=0.5)
    style_axes(ax)


def panel_B(ax, fra: pd.DataFrame):
    """Histogram of reference-isoform fraction of parent-gene expression."""
    vals = fra["ref_fraction_of_gene"].dropna() * 100  # → %
    bins = np.arange(0, 105, 5)
    ax.hist(vals, bins=bins, color=C_BAR, edgecolor="white", linewidth=0.6,
            alpha=0.9)
    med = vals.median()
    ax.axvline(med, color="#c0504d", linestyle="--", linewidth=1.4,
               label=f"median = {med:.1f}%")
    ax.axvline(50, color="#aa6600", linestyle=":", linewidth=1.0,
               label="50% threshold")
    ax.set_xlabel("Reference fraction of gene expression (%)",
                  fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylabel("Number of genes", fontsize=BODY_FS, color=C_TITLE)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.legend(loc="upper right", fontsize=BODY_FS - 2, frameon=False,
              handlelength=1.5, handletextpad=0.5)
    style_axes(ax)


def panel_C(ax, txl: pd.DataFrame):
    """Violin of tx-length by role, log-scale y-axis."""
    roles = ["NMD comparator", "Reference", "Control comparator"]
    colors = {"NMD comparator": C_NMD, "Reference": C_REF,
              "Control comparator": C_CONTROL}

    data = [np.log10(txl.loc[txl["role"] == r, "length_nt"].values + 1)
            for r in roles]
    parts = ax.violinplot(data, positions=range(len(roles)), widths=0.7,
                          showmeans=False, showmedians=False,
                          showextrema=False)
    for body, role in zip(parts["bodies"], roles):
        body.set_facecolor(colors[role])
        body.set_edgecolor(C_AXIS)
        body.set_alpha(0.75)

    # Median + q25/q75 markers
    for i, r in enumerate(roles):
        vals = txl.loc[txl["role"] == r, "length_nt"].values
        med = np.median(vals)
        q25, q75 = np.percentile(vals, [25, 75])
        ax.plot([i - 0.2, i + 0.2], [np.log10(med + 1)] * 2,
                color=C_TITLE, lw=1.6)
        ax.plot([i - 0.15, i + 0.15], [np.log10(q25 + 1)] * 2,
                color=C_TITLE, lw=0.9, linestyle=":")
        ax.plot([i - 0.15, i + 0.15], [np.log10(q75 + 1)] * 2,
                color=C_TITLE, lw=0.9, linestyle=":")
        ax.text(i, np.log10(med + 1) + 0.04, f"{int(med):,}",
                ha="center", va="bottom", fontsize=BODY_FS - 2,
                fontweight="bold", color=C_TITLE)

    ax.set_xticks(range(len(roles)))
    ax.set_xticklabels(["NMD\ncomparator", "Reference", "Control\ncomparator"],
                       fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylabel("Transcript length (nt, log scale)",
                  fontsize=BODY_FS, color=C_TITLE)
    ax.set_yticks([np.log10(10), np.log10(100), np.log10(1_000),
                   np.log10(10_000), np.log10(100_000)])
    ax.set_yticklabels(["10", "100", "1,000", "10,000", "100,000"])
    style_axes(ax)


def build_figure():
    iso, fra, txl, _ = load()

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.16, top=0.86,
                        wspace=0.30)

    panel_A(axes[0], iso)
    panel_B(axes[1], fra)
    panel_C(axes[2], txl)

    # Panel labels
    for ax, letter in zip(axes, ["A", "B", "C"]):
        ax.text(-0.13, 1.02, letter, transform=ax.transAxes,
                fontsize=HEADER_FS + 2, fontweight="bold",
                va="bottom", ha="left", color=C_TITLE)

    # Per-panel sub-headers
    axes[0].set_title("Isoforms per gene", fontsize=BODY_FS + 1,
                       color=C_TITLE, pad=4)
    axes[1].set_title("Reference share of gene expression",
                       fontsize=BODY_FS + 1, color=C_TITLE, pad=4)
    axes[2].set_title("Transcript length by pair role",
                       fontsize=BODY_FS + 1, color=C_TITLE, pad=4)
    return fig, axes


def main():
    fig, axes = build_figure()
    if validate_multipanel_layout is not None:
        try:
            validate_multipanel_layout(fig)
        except Exception as e:
            print(f"[validate_multipanel_layout] non-fatal: {e}")
    out = HERE.parent / "figure_s_pairset_descriptives"
    fig.savefig(f"{out}.pdf", facecolor="white")
    fig.savefig(f"{out}.png", dpi=300, facecolor="white")
    print(f"Saved: {out}.pdf and .png")


if __name__ == "__main__":
    main()
