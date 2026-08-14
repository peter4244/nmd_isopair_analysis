"""MOCKUP — Figure 3D on the reference-AUG-traceable universe (n=819 per arm).

Not a deliverable. Exploratory look at what Panel D becomes when the scope moves from the
all-3-ENST GENCODE subset (130/130) to the reference-AUG-traceable cohort, with TD2's
anti-PTC bias removed by tracing the reference start codon into each comparator.

Same quantity, same convention, same styling as the shipped panel so the two are comparable:
  distance = last_ejc_tx_pos - stop_tx_pos ; positive = stop UPSTREAM of the last EJC.
The last-EJC arithmetic was validated at 260/260 exact agreement against the shipped panel data.

n=819 is a reconstruction of the published 833 cohort and is 14 pairs short of it; the PTC rate
agrees to 0.1 points (92.3% vs 92.2%). Labelled as such on the figure.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14
C_NMD_FILL, C_NMD_LINE = "#ef8a62", "#d6604d"
C_CTRL_FILL, C_CTRL_LINE = "#67a9cf", "#2b8cbe"
C_AXIS, C_THRESHOLD = "#555555", "#cc0000"
X_MIN, X_MAX = -1000, 1500

HERE = Path(__file__).resolve().parent
df = pd.read_csv(HERE / "panelD_819.tsv", sep="\t")

fig, ax = plt.subplots(figsize=(6, 4))
x_grid = np.linspace(X_MIN, X_MAX, 1000)

for label, fill, line in (("Control", C_CTRL_FILL, C_CTRL_LINE), ("NMD", C_NMD_FILL, C_NMD_LINE)):
    vals = df.loc[df["comparison"] == label, "distance"].to_numpy(float)
    kept = vals[(vals >= X_MIN) & (vals <= X_MAX)]
    kde = gaussian_kde(kept, bw_method="scott")
    y = kde(x_grid)
    ax.fill_between(x_grid, y, alpha=0.45, color=fill, linewidth=0)
    ax.plot(x_grid, y, color=line, linewidth=2,
            label=f"{label} (n={len(vals):,})")
    beyond = 100.0 * (vals >= 50).mean()
    print(f"  {label:8s} n={len(vals):4d}  clipped out={len(vals)-len(kept):3d}  "
          f">=+50nt {beyond:5.1f}%  median {np.median(vals):7.1f}")

ax.axvline(50, color=C_THRESHOLD, linestyle="--", linewidth=1.2, alpha=0.7)
ax.text(50, ax.get_ylim()[1] * 0.97, "  PTC ≥50 nt", color=C_THRESHOLD,
        fontsize=BODY_FS, va="top", ha="left")

ax.set_xlim(X_MIN, X_MAX)
ax.set_ylim(bottom=0)
ax.set_xlabel("Distance from last EJC (nt)", fontsize=HEADER_FS, color="#222222")
ax.set_yticks([])
ax.tick_params(axis="x", labelsize=BODY_FS, colors=C_AXIS)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.spines["bottom"].set_color(C_AXIS)

leg = ax.legend(fontsize=BODY_FS, loc="upper right", frameon=False, handlelength=1.4)
for t in leg.get_texts():
    t.set_color("#222222")

fig.tight_layout()
fig.savefig(HERE / "panelD_819_mockup.png", dpi=200, bbox_inches="tight")
print(f"\nwrote {HERE / 'panelD_819_mockup.png'}")
