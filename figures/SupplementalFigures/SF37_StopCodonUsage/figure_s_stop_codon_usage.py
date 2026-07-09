"""
SF37 — Stop-codon usage in NMD susceptible vs non-NMD isoforms.

The figure shows within-class frequency of the three canonical stop
codons (TGA / TAA / TAG) at the priority ORF of the deep-learning
model's test cohort. Percentages: UGA (TGA) 55.9% vs 48.1% Control;
UAA 23.9% vs 29.9%; UAG 20.2% vs 22.0%.

Data source: the canonical render is `fig4a_stop_codon` from
`results/.../deep_nmd_model/orf_model_report_v5.html`. The legacy
rendered `fig3c_stop_identity.png` on disk is **stale** (pre–
2026-04-30 `scripts/patch_stop_codon.py` bug fix) and must not be
reused. We extract the post-fix v5 image from the HTML's base64 stream
and embed it here as a draft placeholder. The frequencies will be
recomputed natively when the cluster is back, but the directional
finding (UGA > NMD) is already canonical.

Placeholder pipeline (matches Figure 5 Panel F approach):
  1. Extract fig4a_stop_codon base64 from orf_model_report_v5.html
     into ./data/legacy_fig4a_stop_codon.png (1344 × 1152).
  2. Crop the top ~8.5% to drop the legacy subtitle (which is also
     dropped in the current Rmd source via the corresponding
     orf_model_report_v5.Rmd edit). Cropped PNG is stored at
     ./data/legacy_fig4a_stop_codon_nosubtitle.png.
  3. Embed the cropped version as the figure body via matplotlib imshow.
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

from ggplot_style import (
    apply_ggplot_rcparams,
    assert_text_within_canvas,
)
apply_ggplot_rcparams()

LEGACY_HTML = Path(
    "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/"
    "Version_6.0/isopair_wrapper/deep_nmd_model/orf_model_report_v5.html"
)
LEGACY_PNG_FULL = HERE.parent / "data" / "legacy_fig4a_stop_codon.png"
PLACEHOLDER_PNG = HERE.parent / "data" / "legacy_fig4a_stop_codon_nosubtitle.png"

# Fraction of the legacy PNG's vertical extent occupied by the subtitle
# band at the top (legend area underneath is part of the chart body, so
# we crop ONLY the subtitle line — empirically ~8.5%).
SUBTITLE_CROP_FRAC = 0.085


def ensure_extract():
    """Extract the post-bug-fix fig4a_stop_codon image from the v5 HTML
    if not already cached."""
    if LEGACY_PNG_FULL.exists():
        return
    LEGACY_PNG_FULL.parent.mkdir(parents=True, exist_ok=True)
    html = LEGACY_HTML.read_text()
    # First standalone fig4a_stop_codon (not fig4a_stop_codon_subgroup)
    m = next(re.finditer(r"fig4a_stop_codon(?!_)", html))
    sub = html[m.start() : m.start() + 6_000_000]
    img_m = re.search(r"data:image/png;base64,([A-Za-z0-9+/=]+)", sub)
    if img_m is None:
        raise FileNotFoundError(
            "fig4a_stop_codon image not found in v5 HTML; cannot build placeholder."
        )
    LEGACY_PNG_FULL.write_bytes(base64.b64decode(img_m.group(1)))


def ensure_crop():
    """Crop the top subtitle band to match the post-edit Rmd source
    (subtitle dropped in orf_model_report_v5.Rmd at lines ~1164–1168)."""
    if PLACEHOLDER_PNG.exists():
        return
    ensure_extract()
    im = Image.open(LEGACY_PNG_FULL)
    w, h = im.size
    im.crop((0, int(h * SUBTITLE_CROP_FRAC), w, h)).save(PLACEHOLDER_PNG)


def build_figure():
    ensure_crop()
    img = plt.imread(PLACEHOLDER_PNG)

    # Aspect after crop: 1344 / 1055 ≈ 1.274
    fig, ax = plt.subplots(figsize=(5.5, 4.3))
    ax.imshow(img)
    ax.set_axis_off()
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
    return fig


def main():
    fig = build_figure()
    assert_text_within_canvas(fig)
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure_s_stop_codon_usage.pdf", facecolor="white")
    fig.savefig(out_dir / "figure_s_stop_codon_usage.png", dpi=300, facecolor="white")
    print("Saved: figure_s_stop_codon_usage.pdf and .png  (DRAFT placeholder)")


if __name__ == "__main__":
    main()
