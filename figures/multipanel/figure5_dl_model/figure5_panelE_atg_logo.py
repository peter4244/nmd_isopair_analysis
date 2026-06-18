"""
Figure 5 Panel E — DRAFT placeholder using the legacy Kozak-zoom NMD logo.

We don't have the per-isoform signed SHAP × input data locally (only
mean |SHAP| per channel per position). The legacy `motif_logo_atg_*.tsv`
that would let us re-render the signed logo from scratch was cleaned out
of results_4ct/ before this laptop's copy; the only signed-direction
artefact still on disk is the legacy rendered PNG:

  results/.../deep_nmd_model/figures/fig3a2_kozak_zoom.png   (2400 × 2400)

That image stacks NMD (top) and Control (bottom) sub-panels. We crop the
top half into ./data/legacy_kozak_zoom_nmd_top.png and embed it here as
a placeholder for the draft composite. The Kozak resonance reads cleanly
from this version — strong A/G at −3, C at −1, G at +4, all above the
zero line ('positive = toward NMD').

PLACEHOLDER. When the cluster comes back up:
  1. Re-run scripts/export_joint_motif_logos.py on Explorer to produce
     motif_logo_atg_joint_atg500_stop500.tsv with nmd_mean_contrib.
  2. scp that TSV to ./data/.
  3. Rewrite this script to draw the logo natively from that TSV (using
     logomaker with flip_below=True), restoring our Kozak band overlay
     + position annotations + Kozak-coordinate x-axis.

The legacy crop preserves the legacy axis labels (Position relative to
ATG / "T" not "U"); the final non-placeholder render will use AUG / U.

Companion Supplemental Figure shows the NMD-only and Control-only logos
side by side (the bottom half of the same legacy PNG is the Control
sub-panel).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

LEGACY_FULL = Path(
    "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/"
    "Version_6.0/isopair_wrapper/deep_nmd_model/figures/fig3a2_kozak_zoom.png"
)
PLACEHOLDER_CROP = HERE.parent / "data" / "legacy_kozak_zoom_nmd_top.png"


def ensure_crop():
    """Crop the top (NMD) sub-panel from the legacy PNG if not already cached."""
    if PLACEHOLDER_CROP.exists():
        return
    PLACEHOLDER_CROP.parent.mkdir(parents=True, exist_ok=True)
    im = Image.open(LEGACY_FULL)
    w, h = im.size
    top = im.crop((0, 0, w, h // 2))
    top.save(PLACEHOLDER_CROP)


def build_figure():
    ensure_crop()
    img = plt.imread(PLACEHOLDER_CROP)

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.imshow(img)
    ax.set_axis_off()
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
    return fig


def main():
    fig = build_figure()
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelE_atg_logo.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelE_atg_logo.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelE_atg_logo.pdf and .png  (DRAFT placeholder)")


if __name__ == "__main__":
    main()
