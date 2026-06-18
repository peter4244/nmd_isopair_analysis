"""
Figure 5 Panel F — DRAFT placeholder using the legacy stop-codon
signed-SHAP × input nucleotide logo (NMD samples).

Same construction as Panel E. The legacy "fig4a_stop_logo" is a
ggseqlogo render of mean signed SHAP × input around the stop codon
(NMD on top, Control on bottom). It was never written to the figures
folder as a standalone file, but it IS embedded as base64 in the
orf_model_report_v5.html. We extract it to:

    data/legacy_fig4a_stop_logo.png   (1920 × 1536, both NMD + Control)

and crop the TOP half (NMD only) into:

    data/legacy_stop_logo_nmd_top.png (1920 × 768)

which is then embedded here as the Panel F placeholder. Letters carry
proper signed direction; the visual style now matches Panel E.

PLACEHOLDER. When the cluster comes back up:
  1. Re-run scripts/export_joint_motif_logos.py on Explorer to produce
     motif_logo_stop_joint_atg500_stop500.tsv.
  2. scp that TSV to ./data/.
  3. Rewrite this script to draw the logo natively from that TSV (using
     logomaker with flip_below=True).

Legacy axis labels preserved as in the legacy render.
"""

import base64
import re
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

LEGACY_HTML = Path(
    "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/"
    "Version_6.0/isopair_wrapper/deep_nmd_model/orf_model_report_v5.html"
)
LEGACY_PNG_FULL = HERE.parent / "data" / "legacy_fig4a_stop_logo.png"
PLACEHOLDER_CROP = HERE.parent / "data" / "legacy_stop_logo_nmd_top.png"


def ensure_legacy_png():
    """Extract the fig4a_stop_logo base64 image from the legacy HTML if
    we haven't already."""
    if LEGACY_PNG_FULL.exists():
        return
    LEGACY_PNG_FULL.parent.mkdir(parents=True, exist_ok=True)
    html = LEGACY_HTML.read_text()
    idx = html.find("fig4a_stop_logo")
    sub = html[idx : idx + 6_000_000]
    m = re.search(r"data:image/png;base64,([A-Za-z0-9+/=]+)", sub)
    if not m:
        raise FileNotFoundError(
            "fig4a_stop_logo image not found in legacy HTML; cannot build placeholder."
        )
    LEGACY_PNG_FULL.write_bytes(base64.b64decode(m.group(1)))


def ensure_crop():
    if PLACEHOLDER_CROP.exists():
        return
    ensure_legacy_png()
    im = Image.open(LEGACY_PNG_FULL)
    w, h = im.size
    im.crop((0, 0, w, h // 2)).save(PLACEHOLDER_CROP)


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
    fig.savefig(out_dir / "figure5_panelF_stop_logo.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelF_stop_logo.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelF_stop_logo.pdf and .png  (DRAFT placeholder)")


if __name__ == "__main__":
    main()
