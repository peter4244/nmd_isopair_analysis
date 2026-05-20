#!/usr/bin/env python3
"""
validate_figure_layout.py

Validates spatial layout of matplotlib annotation-based figures (architecture
diagrams, flowcharts, schematics).  Checks for text-arrow collisions, text
overlap, crowding, centering, alignment, and padding.

Python port of validate_figure_layout.R — same 7 checks, adapted for
matplotlib's artist model.

Usage:
    from validate_figure_layout import validate_figure_layout
    fig, ax = plt.subplots(...)
    # ... build figure ...
    result = validate_figure_layout(fig, ax)

All thresholds are in data-coordinate units and can be tuned via arguments.
"""

import math
from dataclasses import dataclass, field
from typing import Optional

import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class TextElement:
    idx: int
    x: float
    y: float
    label: str
    fontsize: float
    ha: str
    va: str
    alpha: float
    bb_xmin: float = 0.0
    bb_xmax: float = 0.0
    bb_ymin: float = 0.0
    bb_ymax: float = 0.0


@dataclass
class SegmentElement:
    idx: int
    x: float
    y: float
    xend: float
    yend: float


@dataclass
class RectElement:
    idx: int
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    facecolor: Optional[str]
    alpha: float
    kind: str = 'rect'   # 'rect' / 'polygon' / 'circle'
    has_edge: bool = False   # outline-only shapes (no fill) are still "real"


@dataclass
class Issue:
    check: str
    severity: str  # ERROR, WARNING, INFO
    message: str
    suggestion: str = ""


# ---------------------------------------------------------------------------
# Extraction: pull texts, segments, and rects from a matplotlib Axes
# ---------------------------------------------------------------------------
def _extract_elements(fig, ax):
    """Extract texts, segments (arrows/lines), and rects from ax."""
    texts = []
    segments = []
    rects = []

    renderer = fig.canvas.get_renderer()
    inv_transform = ax.transData.inverted()

    # --- Texts: ax.texts + annotations (text part) ---
    text_artists = list(ax.texts)
    for child in ax.get_children():
        if isinstance(child, mtext.Annotation) and child not in text_artists:
            text_artists.append(child)

    for i, t in enumerate(text_artists):
        label = t.get_text()
        if not label or not label.strip():
            continue
        alpha = t.get_alpha()
        if alpha is None:
            alpha = 1.0

        pos = t.get_position()
        # For annotations the position is already in data coords if
        # xycoords='data' (the default for ax.annotate).  For ax.text
        # it's always data coords.
        x, y = pos

        fontsize = t.get_fontsize()
        ha = t.get_ha()
        va = t.get_va()

        te = TextElement(idx=i, x=x, y=y, label=label,
                         fontsize=fontsize, ha=ha, va=va, alpha=alpha)

        # Use renderer to get actual bounding box, then convert to data coords
        try:
            bbox_disp = t.get_window_extent(renderer)
            bbox_data = bbox_disp.transformed(inv_transform)
            te.bb_xmin = bbox_data.x0
            te.bb_xmax = bbox_data.x1
            te.bb_ymin = bbox_data.y0
            te.bb_ymax = bbox_data.y1
        except Exception:
            # Fallback: estimate from fontsize and character count
            te = _estimate_text_bbox(te, fig, ax)

        texts.append(te)

    # --- Segments: annotations with arrows, FancyArrowPatch, Line2D ---
    seg_idx = 0
    for child in ax.get_children():
        if isinstance(child, mtext.Annotation):
            arrow_props = child.arrowprops
            if arrow_props is not None:
                # Arrow goes from xytext (start) to xy (end)
                try:
                    xy = child.xy
                    xytext = child.xyann if hasattr(child, 'xyann') else child.get_position()
                    segments.append(SegmentElement(
                        idx=seg_idx, x=xytext[0], y=xytext[1],
                        xend=xy[0], yend=xy[1]))
                    seg_idx += 1
                except Exception:
                    pass

        elif isinstance(child, mpatches.FancyArrowPatch):
            # FancyArrowPatch stores positions differently
            try:
                # _posA_posB is set when using posA/posB constructor
                if hasattr(child, '_posA_posB') and child._posA_posB is not None:
                    posA, posB = child._posA_posB
                else:
                    path = child.get_path()
                    verts = path.vertices
                    posA = verts[0]
                    posB = verts[-1]
                segments.append(SegmentElement(
                    idx=seg_idx, x=posA[0], y=posA[1],
                    xend=posB[0], yend=posB[1]))
                seg_idx += 1
            except Exception:
                pass

        elif isinstance(child, mlines.Line2D):
            xdata = child.get_xdata()
            ydata = child.get_ydata()
            if len(xdata) >= 2:
                # Treat as segment from first to last point
                segments.append(SegmentElement(
                    idx=seg_idx, x=xdata[0], y=ydata[0],
                    xend=xdata[-1], yend=ydata[-1]))
                seg_idx += 1

    # --- Rects + Polygons + Circles: FancyBboxPatch, Rectangle, Polygon, Circle ---
    # Polygons (block arrows) and Circles (dots, markers) are also tracked so
    # the shape-shape overlap check can find arrow-vs-panel and dot-vs-rect
    # collisions. Existing text-rect-centering / text-rect-padding checks
    # filter to kind=='rect' so they don't get spurious hits on those.
    # Skip ax.patch (axes background) — lives in axes-coords.
    rect_idx = 0
    for child in ax.get_children():
        if child is ax.patch:
            continue
        kind = None
        if isinstance(child, (mpatches.FancyBboxPatch, mpatches.Rectangle)):
            kind = 'rect'
        elif isinstance(child, mpatches.Polygon):
            kind = 'polygon'
        elif isinstance(child, mpatches.Circle):
            kind = 'circle'

        if kind is not None:
            try:
                # Prefer logical bbox (no padding/stroke) when available;
                # fall back to the rendered window extent for Polygon/Circle.
                if kind == 'rect':
                    bb = child.get_bbox()
                    xmin, xmax = bb.x0, bb.x1
                    ymin, ymax = bb.y0, bb.y1
                else:
                    bbox_disp = child.get_window_extent(renderer)
                    bbox_data = bbox_disp.transformed(inv_transform)
                    xmin, xmax = bbox_data.x0, bbox_data.x1
                    ymin, ymax = bbox_data.y0, bbox_data.y1

                fc = child.get_facecolor()
                fc_alpha = fc[3] if len(fc) == 4 else 1.0
                facecolor = None if fc_alpha < 0.01 else str(fc)

                # An outline-only shape (no fill but visible stroke) is still
                # a "real" shape that other shapes can overlap.
                ec = child.get_edgecolor()
                ec_alpha = ec[3] if (ec is not None and len(ec) == 4) else 1.0
                lw = child.get_linewidth() or 0
                has_edge = ec_alpha > 0.05 and lw > 0.1

                patch_alpha = child.get_alpha()
                if patch_alpha is None:
                    patch_alpha = 1.0

                rects.append(RectElement(
                    idx=rect_idx,
                    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                    facecolor=facecolor, alpha=patch_alpha,
                    kind=kind, has_edge=has_edge))
                rect_idx += 1
            except Exception:
                pass

    return texts, segments, rects


def _estimate_text_bbox(te, fig, ax):
    """Fallback bbox estimation when renderer bbox fails."""
    fig_w, fig_h = fig.get_size_inches()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]

    mm_per_unit_x = (fig_w * 25.4) / x_range
    mm_per_unit_y = (fig_h * 25.4) / y_range

    sz_mm = te.fontsize * 0.353  # pt to mm
    lines = te.label.split('\n')
    n_lines = len(lines)
    max_nc = max(len(l) for l in lines)

    char_w_mm = sz_mm * 0.55
    char_h_mm = sz_mm * 1.2
    total_w_mm = max_nc * char_w_mm
    total_h_mm = n_lines * char_h_mm

    bbox_w = (total_w_mm / mm_per_unit_x) * 1.1
    bbox_h = (total_h_mm / mm_per_unit_y) * 1.1

    hj = {'center': 0.5, 'right': 1.0, 'left': 0.0}.get(te.ha, 0.5)
    vj = {'center': 0.5, 'top': 1.0, 'bottom': 0.0,
          'baseline': 0.0, 'center_baseline': 0.5}.get(te.va, 0.5)

    te.bb_xmin = te.x - hj * bbox_w
    te.bb_xmax = te.x + (1 - hj) * bbox_w
    te.bb_ymin = te.y - vj * bbox_h
    te.bb_ymax = te.y + (1 - vj) * bbox_h
    return te


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def _segment_intersects_rect(x1, y1, x2, y2, xmin, xmax, ymin, ymax):
    """Liang-Barsky line-rectangle intersection test."""
    dx = x2 - x1
    dy = y2 - y1
    p = [-dx, dx, -dy, dy]
    q = [x1 - xmin, xmax - x1, y1 - ymin, ymax - y1]

    t0, t1 = 0.0, 1.0
    for pk, qk in zip(p, q):
        if abs(pk) < 1e-12:
            if qk < 0:
                return False
        else:
            r = qk / pk
            if pk < 0:
                t0 = max(t0, r)
            else:
                t1 = min(t1, r)
            if t0 > t1:
                return False
    return True


def _point_to_segment_dist(px, py, x1, y1, x2, y2):
    """Minimum distance from point to line segment."""
    dx, dy = x2 - x1, y2 - y1
    len_sq = dx * dx + dy * dy
    if len_sq < 1e-12:
        return math.hypot(px - x1, py - y1)
    t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / len_sq))
    proj_x = x1 + t * dx
    proj_y = y1 + t * dy
    return math.hypot(px - proj_x, py - proj_y)


def _rect_to_segment_dist(rxmin, rxmax, rymin, rymax, sx1, sy1, sx2, sy2):
    """Min distance from rectangle boundary to line segment (sampled)."""
    mx, my = (rxmin + rxmax) / 2, (rymin + rymax) / 2
    pts = [
        (rxmin, rymin), (rxmax, rymin), (rxmin, rymax), (rxmax, rymax),
        (mx, rymin), (mx, rymax), (rxmin, my), (rxmax, my),
    ]
    dists = [_point_to_segment_dist(px, py, sx1, sy1, sx2, sy2)
             for px, py in pts]
    # Also check segment endpoints to rect boundary
    for ex, ey in [(sx1, sy1), (sx2, sy2)]:
        ddx = max(rxmin - ex, 0, ex - rxmax)
        ddy = max(rymin - ey, 0, ey - rymax)
        dists.append(math.hypot(ddx, ddy))
    return min(dists)


def _rects_overlap(xmin1, xmax1, ymin1, ymax1, xmin2, xmax2, ymin2, ymax2):
    return (xmin1 < xmax2 and xmax1 > xmin2 and
            ymin1 < ymax2 and ymax1 > ymin2)


# ---------------------------------------------------------------------------
# Report printer
# ---------------------------------------------------------------------------
def _print_report(result):
    s = result['summary']
    print(f"\n=== Figure Layout Validation ===")
    print(f"Extracted: {s['n_texts']} texts, {s['n_segments']} segments, "
          f"{s['n_rects']} rects\n")

    if s['n_errors'] > 0:
        print(f"ERRORS ({s['n_errors']}):")
        for i, iss in enumerate(result['errors'], 1):
            print(f"  [{i}] {iss.check}: {iss.message}")
            if iss.suggestion:
                print(f"      -> {iss.suggestion}")
        print()

    if s['n_warnings'] > 0:
        print(f"WARNINGS ({s['n_warnings']}):")
        for i, iss in enumerate(result['warnings'], 1):
            print(f"  [{i}] {iss.check}: {iss.message}")
            if iss.suggestion:
                print(f"      -> {iss.suggestion}")
        print()

    if s['n_info'] > 0:
        print(f"INFO ({s['n_info']}):")
        for i, iss in enumerate(result['info'], 1):
            print(f"  [{i}] {iss.check}: {iss.message}")
        print()

    if s['n_errors'] == 0 and s['n_warnings'] == 0:
        print("All checks passed.\n")
    elif s['n_errors'] == 0:
        print(f"No errors. {s['n_warnings']} warning(s) to review.\n")
    else:
        print(f"Result: {s['n_errors']} error(s) must be fixed.\n")


# ===========================================================================
# Main validation function
# ===========================================================================
def validate_figure_layout(
    fig,
    ax,
    *,
    collision_buffer=0.05,
    crowding_distance=0.4,
    centering_tolerance=0.3,
    alignment_y_tol=0.2,
    min_arrow_length=0.3,
    rect_padding_min=0.15,
    verbose=True,
):
    """
    Validate layout of a matplotlib annotation-based figure.

    Runs 8 checks:
      1. Text-segment collision (ERROR)
      2. Text-text overlap (ERROR)
      3. Text-segment proximity / crowding (WARNING)
      4. Text-rect centering (WARNING)
      5. Horizontal alignment (INFO)
      6. Arrow minimum length (WARNING)
      7. Text-rect padding / overflow (WARNING/ERROR)
      8. Shape-shape partial overlap (ERROR) — rects + polygons; skips
         intentional containment so background cards holding contents
         are still allowed.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes
    collision_buffer : float
        Extra buffer around text bbox for collision detection.
    crowding_distance : float
        Threshold for text-segment proximity warnings.
    centering_tolerance : float
        Max allowed horizontal offset of text from enclosing rect center.
    alignment_y_tol : float
        Tolerance for grouping texts into horizontal rows.
    min_arrow_length : float
        Minimum acceptable arrow/segment length.
    rect_padding_min : float
        Minimum padding between text bbox and enclosing rect edges.
    verbose : bool
        Print formatted report.

    Returns
    -------
    dict with keys: errors, warnings, info, summary
    """

    # Force a draw so renderers can compute extents
    fig.canvas.draw()

    errors = []
    warnings = []
    info = []

    # Phase A: Extract elements
    texts, segments, rects = _extract_elements(fig, ax)

    # --- Check 1: Text-segment collision (ERROR) ---
    if texts and segments:
        for te in texts:
            if te.alpha < 0.1:
                continue
            bb = (te.bb_xmin - collision_buffer, te.bb_xmax + collision_buffer,
                  te.bb_ymin - collision_buffer, te.bb_ymax + collision_buffer)
            for seg in segments:
                if _segment_intersects_rect(
                        seg.x, seg.y, seg.xend, seg.yend,
                        bb[0], bb[1], bb[2], bb[3]):
                    errors.append(Issue(
                        "text_segment_collision", "ERROR",
                        f"Text '{te.label[:40]}' at ({te.x:.1f},{te.y:.1f}) "
                        f"intersects segment ({seg.x:.1f},{seg.y:.1f})->"
                        f"({seg.xend:.1f},{seg.yend:.1f})",
                        "Move text away from the arrow path or reroute the segment"
                    ))

    # --- Check 2: Text-text overlap (ERROR) ---
    if len(texts) > 1:
        for i in range(len(texts)):
            if texts[i].alpha < 0.1:
                continue
            for j in range(i + 1, len(texts)):
                if texts[j].alpha < 0.1:
                    continue
                if _rects_overlap(
                        texts[i].bb_xmin, texts[i].bb_xmax,
                        texts[i].bb_ymin, texts[i].bb_ymax,
                        texts[j].bb_xmin, texts[j].bb_xmax,
                        texts[j].bb_ymin, texts[j].bb_ymax):
                    errors.append(Issue(
                        "text_text_overlap", "ERROR",
                        f"Text '{texts[i].label[:30]}' at ({texts[i].x:.1f},"
                        f"{texts[i].y:.1f}) overlaps '{texts[j].label[:30]}' "
                        f"at ({texts[j].x:.1f},{texts[j].y:.1f})",
                        "Increase vertical or horizontal spacing between these labels"
                    ))

    # --- Check 3: Text-segment proximity (WARNING) ---
    if texts and segments:
        for te in texts:
            if te.alpha < 0.1:
                continue
            for seg in segments:
                dist = _rect_to_segment_dist(
                    te.bb_xmin, te.bb_xmax, te.bb_ymin, te.bb_ymax,
                    seg.x, seg.y, seg.xend, seg.yend)
                if 0 < dist < crowding_distance:
                    warnings.append(Issue(
                        "text_segment_proximity", "WARNING",
                        f"Text '{te.label[:30]}' at ({te.x:.1f},{te.y:.1f}) "
                        f"is {dist:.2f} units from segment ({seg.x:.1f},"
                        f"{seg.y:.1f})->({seg.xend:.1f},{seg.yend:.1f})",
                        "Consider increasing spacing to avoid visual crowding"
                    ))

    # --- Check 4: Text-rect centering (WARNING) ---
    # Compute total axes area for filtering out large container rects
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    _axes_area = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0])
    _max_rect_area = _axes_area * 0.15  # skip rects > 15% of axes area

    if texts and rects:
        for te in texts:
            if te.alpha < 0.1:
                continue
            best_rect = None
            best_area = float('inf')
            for r in rects:
                if r.facecolor is None:
                    continue
                # Skip non-rect shapes (centering inside polygons isn't meaningful)
                if r.kind != 'rect':
                    continue
                # Skip low-alpha rects (background containers)
                if r.alpha < 0.5:
                    continue
                rect_area = (r.xmax - r.xmin) * (r.ymax - r.ymin)
                # Skip very large rects (axes frame, background panels)
                if rect_area > _max_rect_area:
                    continue
                inside_x = r.xmin <= te.x <= r.xmax
                inside_y = r.ymin <= te.y <= r.ymax
                if inside_x and inside_y:
                    area = (r.xmax - r.xmin) * (r.ymax - r.ymin)
                    if area < best_area:
                        best_area = area
                        best_rect = r
            if best_rect is not None:
                rx_center = (best_rect.xmin + best_rect.xmax) / 2
                offset = abs(te.x - rx_center)
                if offset > centering_tolerance:
                    warnings.append(Issue(
                        "text_rect_centering", "WARNING",
                        f"Text '{te.label[:30]}' at x={te.x:.1f} is {offset:.2f} "
                        f"units off-center from rect [{best_rect.xmin:.1f}, "
                        f"{best_rect.xmax:.1f}] (center={rx_center:.1f})",
                        f"Adjust text x to {rx_center:.1f}"
                    ))

    # --- Check 5: Horizontal alignment (INFO) ---
    if len(texts) > 2:
        ys = [(i, te.y) for i, te in enumerate(texts)]
        ys.sort(key=lambda t: t[1])
        groups = [[ys[0]]]
        for k in range(1, len(ys)):
            if abs(ys[k][1] - ys[k - 1][1]) <= alignment_y_tol:
                groups[-1].append(ys[k])
            else:
                groups.append([ys[k]])

        for g in groups:
            if len(g) < 3:
                continue
            mean_y = sum(t[1] for t in g) / len(g)
            for idx, y_val in g:
                if abs(y_val - mean_y) > alignment_y_tol:
                    te = texts[idx]
                    info.append(Issue(
                        "horizontal_alignment", "INFO",
                        f"Text '{te.label[:30]}' at y={te.y:.2f} deviates "
                        f"from group mean y={mean_y:.2f} (group of {len(g)} texts)",
                        f"Consider aligning y to {mean_y:.2f}"
                    ))

    # --- Check 6: Arrow minimum length (WARNING) ---
    if segments:
        for seg in segments:
            length = math.hypot(seg.xend - seg.x, seg.yend - seg.y)
            if length < min_arrow_length:
                warnings.append(Issue(
                    "arrow_min_length", "WARNING",
                    f"Segment ({seg.x:.1f},{seg.y:.1f})->({seg.xend:.1f},"
                    f"{seg.yend:.1f}) length={length:.2f} < minimum "
                    f"{min_arrow_length:.2f}",
                    "Increase spacing between connected elements or remove segment"
                ))

    # --- Check 7: Text-rect padding / overflow (WARNING/ERROR) ---
    if texts and rects:
        for te in texts:
            if te.alpha < 0.1:
                continue
            for r in rects:
                if r.facecolor is None:
                    continue
                # Padding is meaningful inside true rects, not block-arrow polygons
                if r.kind != 'rect':
                    continue

                inside = (te.bb_xmin >= r.xmin - 0.01 and
                          te.bb_xmax <= r.xmax + 0.01 and
                          te.bb_ymin >= r.ymin - 0.01 and
                          te.bb_ymax <= r.ymax + 0.01)

                center_inside = (r.xmin <= te.x <= r.xmax and
                                 r.ymin <= te.y <= r.ymax)

                if not inside and center_inside:
                    overflows = []
                    if te.bb_xmin < r.xmin:
                        overflows.append("left")
                    if te.bb_xmax > r.xmax:
                        overflows.append("right")
                    if te.bb_ymin < r.ymin:
                        overflows.append("bottom")
                    if te.bb_ymax > r.ymax:
                        overflows.append("top")
                    if overflows:
                        errors.append(Issue(
                            "text_rect_overflow", "ERROR",
                            f"Text '{te.label[:30]}' overflows rect "
                            f"[{r.xmin:.1f},{r.xmax:.1f}]x"
                            f"[{r.ymin:.1f},{r.ymax:.1f}] on: "
                            f"{', '.join(overflows)}",
                            "Reduce text size or widen the rect"
                        ))
                elif inside:
                    left_pad = te.bb_xmin - r.xmin
                    right_pad = r.xmax - te.bb_xmax
                    top_pad = r.ymax - te.bb_ymax
                    bottom_pad = te.bb_ymin - r.ymin
                    pads = {'left': left_pad, 'right': right_pad,
                            'top': top_pad, 'bottom': bottom_pad}
                    min_pad = min(pads.values())
                    if min_pad < rect_padding_min:
                        tight = min(pads, key=pads.get)
                        warnings.append(Issue(
                            "text_rect_padding", "WARNING",
                            f"Text '{te.label[:30]}' has only {min_pad:.2f} "
                            f"units padding on {tight} edge of rect "
                            f"[{r.xmin:.1f},{r.xmax:.1f}]x"
                            f"[{r.ymin:.1f},{r.ymax:.1f}]",
                            "Increase rect size or reduce text size"
                        ))

    # --- Check 8: Shape-shape partial overlap (ERROR) ---
    # Flags pairs of visible shapes (rects + polygons + circles) whose bboxes
    # overlap unless one fully contains the other (intentional containment).
    # A shape is considered "visible" if it has a fill OR a visible stroke, so
    # outline-only panels still get checked against arrows passing through.
    # Circle-circle pairs are skipped — overlapping circles are typically
    # decorative (smoke clouds, cell-cluster glyphs) and produce noise.
    def _shape_contained(inner, outer, tol=0.02):
        return (inner.xmin >= outer.xmin - tol and
                inner.xmax <= outer.xmax + tol and
                inner.ymin >= outer.ymin - tol and
                inner.ymax <= outer.ymax + tol)

    def _is_visible(r):
        return (r.facecolor is not None) or r.has_edge

    overlap_tol = 0.015   # tight: catch even sub-tenth-unit overlaps
    if len(rects) > 1:
        for i in range(len(rects)):
            r1 = rects[i]
            if not _is_visible(r1):
                continue
            for j in range(i + 1, len(rects)):
                r2 = rects[j]
                if not _is_visible(r2):
                    continue
                # Skip pairs where both are circles — usually decorative
                # (smoke clusters, cell nuclei) and not a layout problem.
                if r1.kind == 'circle' and r2.kind == 'circle':
                    continue
                ox = min(r1.xmax, r2.xmax) - max(r1.xmin, r2.xmin)
                oy = min(r1.ymax, r2.ymax) - max(r1.ymin, r2.ymin)
                if ox <= overlap_tol or oy <= overlap_tol:
                    continue
                if _shape_contained(r1, r2) or _shape_contained(r2, r1):
                    continue
                errors.append(Issue(
                    "shape_shape_overlap", "ERROR",
                    f"{r1.kind} [{r1.xmin:.2f},{r1.xmax:.2f}]x"
                    f"[{r1.ymin:.2f},{r1.ymax:.2f}] partially overlaps "
                    f"{r2.kind} [{r2.xmin:.2f},{r2.xmax:.2f}]x"
                    f"[{r2.ymin:.2f},{r2.ymax:.2f}] "
                    f"(overlap {ox:.2f} x {oy:.2f})",
                    "Reposition or resize so shapes either don't overlap, "
                    "or one fully contains the other"
                ))

    # --- Assemble and report ---
    result = {
        'errors': errors,
        'warnings': warnings,
        'info': info,
        'summary': {
            'n_texts': len(texts),
            'n_segments': len(segments),
            'n_rects': len(rects),
            'n_errors': len(errors),
            'n_warnings': len(warnings),
            'n_info': len(info),
        }
    }

    if verbose:
        _print_report(result)

    return result
