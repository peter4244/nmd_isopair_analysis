"""Matplotlib style helpers that mimic ggplot2's default theme.

Used by all §4 and §5 Supplemental Figures so they visually match the
R/ggplot-generated SF1–SF23 while keeping matplotlib as the tooling.

The three functions to call from a panel script:

  - `apply_ggplot_rcparams()` at import time to set font + fonttype.
  - `style_axes_ggplot(ax)` inside the panel render, right after building
    the axes. Applies the grey panel background and white gridlines.
  - `facet_header(ax, text)` to place a ggplot-style grey-strip header
    above an axes (matches the "AT2 / LAE / FB / MV" strip you see in
    SF1, SF3, SF19, SF22, etc.).

Paper-canonical colours are re-exported here so panels can `from
ggplot_style import CT_COLORS, NMD_COLOR, CONTROL_COLOR, ...` rather than
each rebuilding the palette.
"""

from __future__ import annotations
import matplotlib.pyplot as plt

# ── palette (matches project_nmd_figure_palette memory + Panel C py) ──
NMD_COLOR     = "#ef8a62"   # NMD comparator; also used for NMD/Susceptible class
CONTROL_COLOR = "#67a9cf"   # Control comparator; also used for non-NMD class
REF_COLOR     = "#7f7f7f"   # Reference isoform (structurally distinct)
BAR_COLOR     = "#4292c6"   # blue for descriptive histograms

# Cell-type palette used by the R/ggplot SF20 + SF21 pages.
CT_COLORS = {
    "AT2": "#F19797",   # peachy pink
    "LAE": "#A0C954",   # olive/green
    "FB":  "#3ABAB4",   # teal
    "MV":  "#B08CD3",   # lavender/purple
}
CT_ORDER = ["AT2", "LAE", "FB", "MV"]

# Structural colours for the theme itself.
PANEL_BG   = "#EBEBEB"   # ggplot theme_grey panel background
GRID_COLOR = "#FFFFFF"   # ggplot theme_grey major grid
STRIP_BG   = "#D9D9D9"   # facet strip background
TITLE_C    = "#222222"   # near-black for chart text
AXIS_C     = "#555555"   # spine + tick lines

# Font sizes: reader spec is two sizes only per figures_principles.md
HEADER_FS = 14
BODY_FS   = 11
PANEL_LETTER_FS = 18  # bold A/B/C labels top-left of each panel


def apply_ggplot_rcparams():
    """Set matplotlib rcParams to fonts/embedding compatible with SF1-SF23.

    Call this once at import time in a panel script (before creating any
    Figure). Sets Arial-family fonts and PDF/PS fonttype=42 so text is
    embedded as vector glyphs (searchable, non-rasterised).
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.titlesize": HEADER_FS,
        "axes.labelsize": BODY_FS,
        "xtick.labelsize": BODY_FS,
        "ytick.labelsize": BODY_FS,
        "legend.fontsize": BODY_FS,
    })


def style_axes_ggplot(ax, *, xgrid=True, ygrid=True, minor=False):
    """Apply the ggplot theme_grey look to `ax`.

    - Grey panel background (`PANEL_BG`).
    - White major gridlines behind data.
    - No visible spines (ggplot's default hides them; the panel bg
      provides the visual frame).
    - Tick marks removed, only labels kept.
    - Text drawn in `TITLE_C`.

    Call this immediately after creating the axes and before drawing.
    """
    ax.set_facecolor(PANEL_BG)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    # Grid behind data. axisbelow=True is critical so bars don't cover the grid.
    ax.set_axisbelow(True)
    ax.grid(
        axis="both" if (xgrid and ygrid) else ("x" if xgrid else "y"),
        which="major",
        color=GRID_COLOR,
        linewidth=0.8,
        zorder=0,
    )
    if minor:
        ax.minorticks_on()
        ax.grid(which="minor", color=GRID_COLOR, linewidth=0.4, zorder=0)
    # Suppress tick marks; keep labels in TITLE_C.
    ax.tick_params(
        axis="both",
        which="both",
        length=0,
        colors=TITLE_C,
        labelsize=BODY_FS,
    )
    # Axis labels — set colour only; leave the actual label text to caller.
    ax.xaxis.label.set_color(TITLE_C)
    ax.yaxis.label.set_color(TITLE_C)
    return ax


def facet_header(ax, text, *, height=0.06):
    """Place a ggplot-style facet strip header above an axes.

    Draws a grey rectangle sitting on top of the axes with `text` in
    bold black — mimics the strip labels you see above each facet in
    SF1 / SF3 / SF19 (e.g., "AT2", "LAE", "FB", "MV").

    `height` is the strip height in axes-relative coordinates
    (fraction of the plot area). Default 0.06 works for cells around
    2 inches tall; increase for bigger cells.
    """
    ax.text(
        0.5, 1.0 + height / 2,
        text,
        transform=ax.transAxes,
        ha="center", va="center",
        fontsize=HEADER_FS - 1,
        fontweight="bold",
        color=TITLE_C,
        bbox=dict(
            facecolor=STRIP_BG,
            edgecolor="none",
            boxstyle="square,pad=0.4",
        ),
    )


def panel_letter(ax, letter, *, x=-0.02, y=1.02):
    """Place a bold A/B/C letter in the top-left of a panel axes.

    Matches the sequential letter labels convention documented in
    ~/.claude/memory/figures_style_publications.md (composite letter
    label position).

    NB: this text sits ABOVE the axes (y > 1). If you also have a
    facet_header on the axes, increase y further (typically 1.10-1.14)
    AND make sure fig.subplots_adjust(top=...) leaves ≥ 0.15 headroom
    above the axes so the letter's ascender doesn't clip at the canvas
    top. Run assert_text_within_canvas(fig) before savefig to catch
    clipping automatically.
    """
    ax.text(
        x, y, letter,
        transform=ax.transAxes,
        fontsize=PANEL_LETTER_FS,
        fontweight="bold",
        color=TITLE_C,
        ha="left", va="bottom",
    )


def assert_text_within_canvas(fig, *, tolerance_px=1.0):
    """Raise AssertionError if any Text artist extends outside the figure canvas.

    This is the systematic guard for the "panel letter clipped at top"
    class of bug — panel letters and facet-header text sit at y > 1 in
    axes coords, and if fig.subplots_adjust(top=…) doesn't leave enough
    headroom their ascenders get chopped by the canvas boundary.

    Call this after all text is placed and BEFORE savefig. Raises with a
    detailed message pointing at the offending Text object so the user
    knows what to move.
    """
    # Force a draw so text bounding boxes are populated.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_bbox = fig.bbox  # in display pixels
    offenders = []
    for ax in fig.axes:
        for txt in list(ax.texts) + [ax.xaxis.label, ax.yaxis.label, ax.title]:
            s = txt.get_text() if hasattr(txt, "get_text") else ""
            if not s or not s.strip():
                continue
            try:
                bb = txt.get_window_extent(renderer=renderer)
            except Exception:
                continue
            if (bb.x0 < fig_bbox.x0 - tolerance_px
                    or bb.y0 < fig_bbox.y0 - tolerance_px
                    or bb.x1 > fig_bbox.x1 + tolerance_px
                    or bb.y1 > fig_bbox.y1 + tolerance_px):
                offenders.append(
                    f"  '{s}' extends outside canvas "
                    f"(text bbox x0={bb.x0:.0f} y0={bb.y0:.0f} x1={bb.x1:.0f} y1={bb.y1:.0f}; "
                    f"canvas x=[0..{fig_bbox.x1:.0f}] y=[0..{fig_bbox.y1:.0f}])"
                )
    if offenders:
        raise AssertionError(
            "Text elements would be clipped by the figure canvas — adjust "
            "fig.subplots_adjust(top=…) / panel_letter(y=…) / margins:\n"
            + "\n".join(offenders)
        )
