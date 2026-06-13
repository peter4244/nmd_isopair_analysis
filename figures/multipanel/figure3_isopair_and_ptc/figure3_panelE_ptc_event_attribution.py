"""
Figure 3 Panel E — PTC-introducing event attribution.

Paired barplot: among events Isopair attributes as the cause of a PTC in
NMD comparators (NMD-causing), and among ALL splice events in Control
comparators (the background), what fraction is each event type?

Headline: skipped exon (SE) accounts for 52.8% of PTC-causing events vs
8.9% of Control events (~6× enrichment, p ~ 10⁻²²²). A3SS and A5SS are
also significantly enriched among PTC-causing events; terminal events
(Alt TSS, Alt TES) are depleted because they don't typically introduce
PTCs in the coding region.

Data:
  ./data/panelE_ptc_event_attribution.tsv
    Columns: event_type, n_ptc_events, pct_of_ptc, n_ctrl_events, pct_ctrl,
             enrichment, fisher_p, direction
    Source: isopair_wrapper/tables/table5b_ptc_event_enrichment.csv, copied
    by ./data_export.R.

Style: HEADER_FS=18 bold, BODY_FS=14 regular. PTC-causing (NMD-attributed)
coral, Control all-events light blue — matches Panels A/B/C/D.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14

C_NMD = "#ef8a62"
C_CTRL = "#67a9cf"
C_TITLE = "#222222"
C_AXIS = "#555555"

EVENT_LABELS = {
    "SE": "Skipped exon",
    "Alt_TSS": "Alt TSS",
    "Alt_TES": "Alt TES",
    "IR": "Intron retention",
    "A5SS": "A5SS",
    "A3SS": "A3SS",
    "Partial_IR_5": "Partial IR 5′",
    "Partial_IR_3": "Partial IR 3′",
    "IR_diff": "IR differs",
    "IR_diff_5": "IR differs 5′",
    "IR_diff_3": "IR differs 3′",
    "Missing_Internal": "Missing internal",
}


def sig_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def load_data():
    df = pd.read_csv(HERE.parent / "data" / "panelE_ptc_event_attribution.tsv", sep="\t")
    df = df.sort_values("pct_of_ptc", ascending=False).reset_index(drop=True)
    df["label"] = df["event_type"].map(EVENT_LABELS).fillna(df["event_type"])
    df["stars"] = df["fisher_p"].apply(sig_stars)
    return df


def build_figure(df):
    fig, ax = plt.subplots(figsize=(10, 6))

    n = len(df)
    x = np.arange(n)
    bw = 0.38

    ax.bar(x - bw / 2, df["pct_of_ptc"], bw,
           color=C_NMD, edgecolor="none", label="PTC-causing (NMD)")
    ax.bar(x + bw / 2, df["pct_ctrl"], bw,
           color=C_CTRL, edgecolor="none", label="All events (Control)")

    y_max = max(df["pct_of_ptc"].max(), df["pct_ctrl"].max())
    star_pad = y_max * 0.03
    for i, row in df.iterrows():
        top = max(row["pct_of_ptc"], row["pct_ctrl"]) + star_pad
        ax.text(i, top, row["stars"], ha="center", va="bottom",
                fontsize=BODY_FS, color=C_TITLE)

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=35, ha="right",
                       fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_ylabel("% of events", fontsize=BODY_FS, color=C_TITLE)

    # Pair-level PTC counts (constants from the ref-AUG effectively_ptc population
    # — see data_export.R pop_ptc_plus). Hard-coded here because the panelE TSV
    # is event-level, not pair-level.
    N_PTC_NMD = 1912
    N_PTC_CTRL = 288

    ax.text(0, 1.13, "PTC-causing event attribution",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=HEADER_FS, fontweight="bold", color=C_TITLE)
    ax.text(0, 1.02,
            f"PTCs identified: {N_PTC_NMD:,} NMD isoforms  |  {N_PTC_CTRL:,} Control isoforms",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=BODY_FS, color="#666666")

    ax.set_ylim(0, y_max * 1.15)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    legend = ax.legend(fontsize=BODY_FS, loc="upper right", frameon=False)
    for text in legend.get_texts():
        text.set_color(C_TITLE)

    fig.tight_layout()
    return fig


def main():
    df = load_data()
    fig = build_figure(df)
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure3_panelE_ptc_event_attribution.pdf",
                bbox_inches="tight", facecolor="white")
    fig.savefig(out_dir / "figure3_panelE_ptc_event_attribution.png",
                dpi=300, bbox_inches="tight", facecolor="white")
    print("Saved: figure3_panelE_ptc_event_attribution.pdf and .png")


if __name__ == "__main__":
    main()
