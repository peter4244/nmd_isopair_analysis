// transcript_viz.js — render an isoform's exon/CDS structure as an SVG track.
//
// IMPORTANT: This view is rendered in TRANSCRIPT coordinates, not genomic.
// Exons are drawn side-by-side with widths proportional to their lengths,
// separated by a small visual gap. This makes the structural composition
// (exon counts, CDS placement, UTR length) easy to compare across isoforms
// without long introns squashing the exons. The genomic span is shown in
// the meta grid above.
//
// Drawing conventions (Isopair-style):
//   - Each exon = filled rectangle; CDS portion taller than UTR portion.
//   - "Non-coding" badge replaces the ATG/STOP markers when no CDS annotated.
//   - ATG/STOP labels above the exons, with collision-avoidance.
//   - Transcript-coordinate scale at the bottom.
//   - 5' on left, 3' on right (always; transcript-coordinate orientation
//     does not depend on strand).

const COLORS = {
  cds:   "#4a5568",
  utr:   "#cbd5e0",
  axis:  "#9ca3af",
  text:  "#1f2937",
  muted: "#6b7280",
  nmd:   "#ef8a62",
};

const LAYOUT = {
  marginLeft: 30,
  marginRight: 24,
  marginTop: 32,        // headroom for ATG/STOP labels
  marginBottom: 36,     // room for tx-coordinate axis ticks
  exonGap: 4,           // px gap drawn between adjacent exons (transcript view)
  exonH: 18,
  utrH: 10,
  trackY: 40,
};

const NS = "http://www.w3.org/2000/svg";

/**
 * @param {HTMLElement} container — where to mount the SVG
 * @param {Object} iso — single isoform record from a per-gene shard
 */
export function renderTranscript(container, iso) {
  container.innerHTML = "";
  if (!iso || !iso.exon_starts || !iso.exon_ends) {
    container.innerHTML = '<p style="margin:1em;color:#6b7280;font-size:0.9em">No structural data for this isoform.</p>';
    return;
  }

  // Normalise (single-exon arrives as scalar after R round-trip)
  const exonStarts = [].concat(iso.exon_starts);
  const exonEnds   = [].concat(iso.exon_ends);
  const nExons     = exonStarts.length;
  if (nExons === 0) {
    container.innerHTML = '<p style="margin:1em;color:#6b7280;font-size:0.9em">No exon coordinates.</p>';
    return;
  }

  // Per-exon transcript-coordinate windows (1-based, 5'→3')
  // For minus strand we walk exons in reverse genomic order so that 5'→3'
  // along the transcript = left→right on the page.
  const order = (iso.strand === "-") ? d3range(nExons - 1, -1, -1) : d3range(0, nExons);
  const exons = order.map(i => {
    const gFrom = exonStarts[i];
    const gTo   = exonEnds[i];
    return { genomic_index: i, gFrom, gTo, length: gTo - gFrom + 1 };
  });
  // Assign tx coordinates (cumulative)
  let cursor = 1;
  for (const e of exons) {
    e.txFrom = cursor;
    e.txTo   = cursor + e.length - 1;
    cursor   = e.txTo + 1;
  }
  const txLen = cursor - 1;

  // Map a (genomic_index, genomic_position) to its tx coordinate
  const gToTx = (gIdx, gPos) => {
    const e = exons.find(x => x.genomic_index === gIdx);
    if (!e) return null;
    // Genomic offset within the exon — flipped for minus strand
    const off = (iso.strand === "-") ? (e.gTo - gPos) : (gPos - e.gFrom);
    return e.txFrom + off;
  };

  // CDS in tx coordinates (if annotated)
  let cdsTxFrom = null, cdsTxTo = null;
  const _cdsStartNum = (typeof iso.cds_start === "number") ? iso.cds_start : null;
  const _cdsStopNum  = (typeof iso.cds_stop  === "number") ? iso.cds_stop  : null;
  if (_cdsStartNum !== null && _cdsStopNum !== null) {
    // Find which exons contain cds_start and cds_stop
    let startTx = null, stopTx = null;
    for (let i = 0; i < nExons; i++) {
      if (_cdsStartNum >= exonStarts[i] && _cdsStartNum <= exonEnds[i]) {
        startTx = gToTx(i, _cdsStartNum);
      }
      if (_cdsStopNum >= exonStarts[i] && _cdsStopNum <= exonEnds[i]) {
        stopTx = gToTx(i, _cdsStopNum);
      }
    }
    if (startTx != null && stopTx != null) {
      cdsTxFrom = Math.min(startTx, stopTx);
      cdsTxTo   = Math.max(startTx, stopTx);
    }
  }

  // Layout
  const widthPx  = Math.max(container.clientWidth || 700, 400);
  const innerW   = widthPx - LAYOUT.marginLeft - LAYOUT.marginRight;
  const heightPx = LAYOUT.marginTop + LAYOUT.exonH + LAYOUT.marginBottom;

  // Reserve gap space, then scale exon widths to fill innerW
  const totalGapPx = (nExons - 1) * LAYOUT.exonGap;
  const drawableW  = Math.max(40, innerW - totalGapPx);
  const pxPerNt    = drawableW / txLen;

  // SVG
  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("width", widthPx);
  svg.setAttribute("height", heightPx);
  svg.setAttribute("viewBox", `0 0 ${widthPx} ${heightPx}`);

  // Compute screen-x for each exon
  let xCursor = LAYOUT.marginLeft;
  for (const e of exons) {
    e.px = xCursor;
    e.pxW = Math.max(2, e.length * pxPerNt);
    xCursor = e.px + e.pxW + LAYOUT.exonGap;
  }

  // 5'/3' edge labels — transcript-coord view always shows 5'→3' left-to-right,
  // regardless of genomic strand. Strand of the underlying genomic transcript
  // is noted in the top-right of the SVG.
  addText(svg, LAYOUT.marginLeft - 4, LAYOUT.trackY + LAYOUT.exonH / 2 + 4,
          "5'", { anchor: "end", size: 11, color: COLORS.muted });
  addText(svg, widthPx - LAYOUT.marginRight + 4, LAYOUT.trackY + LAYOUT.exonH / 2 + 4,
          "3'", { anchor: "start", size: 11, color: COLORS.muted });
  // Strand indicator + "transcript coordinates, direction always left→right"
  const strandTag = iso.strand === "-" ? "− strand" : "+ strand";
  addText(svg, widthPx - LAYOUT.marginRight, 12,
          `Transcript coords · 5' → 3' left→right · ${strandTag}`,
          { anchor: "end", size: 9, color: COLORS.muted });

  // Non-coding badge if no CDS available
  const hasCDS = (cdsTxFrom != null);
  if (!hasCDS && iso.coding_status === "non_coding") {
    addBadge(svg, widthPx / 2, LAYOUT.marginTop - 18, "non-coding");
  } else if (!hasCDS) {
    addBadge(svg, widthPx / 2, LAYOUT.marginTop - 18, "no annotated CDS");
  }

  // Draw exons (with CDS/UTR shading if applicable)
  for (const e of exons) {
    if (hasCDS) {
      // Three segments per exon: 5'UTR, CDS, 3'UTR
      const segs = [
        { txFrom: e.txFrom, txTo: Math.min(e.txTo, cdsTxFrom - 1), isCDS: false },
        { txFrom: Math.max(e.txFrom, cdsTxFrom), txTo: Math.min(e.txTo, cdsTxTo), isCDS: true },
        { txFrom: Math.max(e.txFrom, cdsTxTo + 1), txTo: e.txTo, isCDS: false },
      ].filter(s => s.txFrom <= s.txTo);
      for (const s of segs) {
        const x = e.px + (s.txFrom - e.txFrom) * pxPerNt;
        const w = Math.max(1.2, (s.txTo - s.txFrom + 1) * pxPerNt);
        const isCDS = s.isCDS;
        const r = document.createElementNS(NS, "rect");
        r.setAttribute("x", x);
        r.setAttribute("y", LAYOUT.trackY + (isCDS ? 0 : (LAYOUT.exonH - LAYOUT.utrH) / 2));
        r.setAttribute("width", w);
        r.setAttribute("height", isCDS ? LAYOUT.exonH : LAYOUT.utrH);
        r.setAttribute("fill", isCDS ? COLORS.cds : COLORS.utr);
        r.setAttribute("stroke", "#374151");
        r.setAttribute("stroke-width", "0.5");
        r.appendChild(makeTitle(
          `Exon ${e.genomic_index + 1}${isCDS ? " (CDS)" : " (UTR)"} · ${s.txTo - s.txFrom + 1} nt`
        ));
        svg.appendChild(r);
      }
    } else {
      const r = document.createElementNS(NS, "rect");
      r.setAttribute("x", e.px);
      r.setAttribute("y", LAYOUT.trackY + (LAYOUT.exonH - LAYOUT.utrH) / 2);
      r.setAttribute("width", e.pxW);
      r.setAttribute("height", LAYOUT.utrH);
      r.setAttribute("fill", COLORS.utr);
      r.setAttribute("stroke", "#374151");
      r.setAttribute("stroke-width", "0.5");
      r.appendChild(makeTitle(`Exon ${e.genomic_index + 1} · ${e.length} nt`));
      svg.appendChild(r);
    }

    // Exon number above
    const cx = e.px + e.pxW / 2;
    addText(svg, cx, LAYOUT.trackY - 4, String(e.genomic_index + 1),
            { anchor: "middle", size: 8, color: COLORS.muted });
  }

  // ATG/STOP markers (with overlap-aware label placement)
  if (hasCDS) {
    const atgX = txToPx(cdsTxFrom, exons, pxPerNt);
    const stopX = txToPx(cdsTxTo,   exons, pxPerNt);
    drawCdsMarker(svg, atgX, "ATG", LAYOUT.marginTop - 8, stopX - atgX < 50);
    drawCdsMarker(svg, stopX, "STOP", LAYOUT.marginTop - 8, stopX - atgX < 50);
  }

  // Transcript-coordinate axis (bottom)
  drawTxAxis(svg, exons, pxPerNt, widthPx, txLen);

  container.appendChild(svg);
}

// ── helpers ──
function d3range(a, b, step = 1) {
  const out = [];
  for (let i = a; (step > 0 ? i < b : i > b); i += step) out.push(i);
  return out;
}
function makeTitle(text) {
  const t = document.createElementNS(NS, "title");
  t.textContent = text;
  return t;
}
function addText(svg, x, y, str, { anchor = "start", size = 10, color = "#1f2937", bold = false } = {}) {
  const t = document.createElementNS(NS, "text");
  t.setAttribute("x", x);
  t.setAttribute("y", y);
  t.setAttribute("text-anchor", anchor);
  t.setAttribute("font-size", size);
  t.setAttribute("fill", color);
  t.setAttribute("font-family", "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif");
  if (bold) t.setAttribute("font-weight", "600");
  t.textContent = str;
  svg.appendChild(t);
  return t;
}
function addBadge(svg, cx, cy, label) {
  const padX = 6, padY = 2;
  const tmp = addText(svg, cx, cy + 3, label, { anchor: "middle", size: 10, color: "#fff" });
  const bb = tmp.getBBox ? tmp.getBBox() : { width: label.length * 6, height: 12 };
  svg.removeChild(tmp);
  const w = bb.width + padX * 2, h = bb.height + padY * 2;
  const rect = document.createElementNS(NS, "rect");
  rect.setAttribute("x", cx - w / 2);
  rect.setAttribute("y", cy - h / 2);
  rect.setAttribute("width", w);
  rect.setAttribute("height", h);
  rect.setAttribute("rx", 4);
  rect.setAttribute("fill", "#6b7280");
  svg.appendChild(rect);
  addText(svg, cx, cy + 3, label, { anchor: "middle", size: 10, color: "#fff", bold: true });
}
function txToPx(tx, exons, pxPerNt) {
  for (const e of exons) {
    if (tx >= e.txFrom && tx <= e.txTo) {
      return e.px + (tx - e.txFrom) * pxPerNt;
    }
  }
  return null;
}
function drawCdsMarker(svg, x, label, labelY, close) {
  // Vertical tick + label above. If markers are close together,
  // alternate label vertical offset to avoid overlap.
  const yShift = close && label === "STOP" ? 10 : 0;
  const tick = document.createElementNS(NS, "line");
  tick.setAttribute("x1", x); tick.setAttribute("y1", LAYOUT.trackY - 2);
  tick.setAttribute("x2", x); tick.setAttribute("y2", LAYOUT.trackY + LAYOUT.exonH + 2);
  tick.setAttribute("stroke", "#1f2937");
  tick.setAttribute("stroke-width", "0.8");
  svg.appendChild(tick);
  addText(svg, x, labelY + yShift, label, { anchor: "middle", size: 9, color: "#1f2937", bold: true });
}
function drawTxAxis(svg, exons, pxPerNt, widthPx, txLen) {
  // Bottom axis: tick marks at sensible round numbers along transcript length.
  const y = LAYOUT.trackY + LAYOUT.exonH + 8;
  // Choose ~5 tick positions in tx coordinates
  const nTicks = 5;
  const txTicks = [];
  for (let i = 0; i <= nTicks; i++) {
    txTicks.push(Math.round((i / nTicks) * txLen));
  }
  // Draw a baseline (only across exon area)
  const baseLine = document.createElementNS(NS, "line");
  baseLine.setAttribute("x1", LAYOUT.marginLeft);
  baseLine.setAttribute("y1", y);
  baseLine.setAttribute("x2", widthPx - LAYOUT.marginRight);
  baseLine.setAttribute("y2", y);
  baseLine.setAttribute("stroke", COLORS.axis);
  baseLine.setAttribute("stroke-width", "0.5");
  svg.appendChild(baseLine);

  for (const tx of txTicks) {
    const px = txToPx(Math.max(1, tx), exons, pxPerNt);
    if (px == null) continue;
    const tick = document.createElementNS(NS, "line");
    tick.setAttribute("x1", px); tick.setAttribute("y1", y);
    tick.setAttribute("x2", px); tick.setAttribute("y2", y + 4);
    tick.setAttribute("stroke", COLORS.axis);
    tick.setAttribute("stroke-width", "0.5");
    svg.appendChild(tick);
    addText(svg, px, y + 14, tx.toLocaleString(),
            { anchor: "middle", size: 9, color: COLORS.muted });
  }
  // Axis label
  addText(svg, widthPx / 2, y + 28, `Transcript position (nt) · ${txLen.toLocaleString()} nt total`,
          { anchor: "middle", size: 9, color: COLORS.muted });
}
