// transcript_stack.js — render multiple isoforms as parallel horizontal tracks
// in genomic coordinates, so alternative splicing patterns are visible at a glance.
//
// Design:
//   - All tracks share one x-axis = the gene's genomic span (min tx_start,
//     max tx_end across the isoforms passed in).
//   - Each track = one isoform: intron line + exon rectangles (CDS taller,
//     UTR shorter) + label on the left.
//   - Row selection: click a track to select it (via onSelect callback);
//     the selected row gets a coloured left border + slightly darker background.
//   - Strand handling: if strand === "-" we flip x mapping so 5' is always on
//     the left. Mixed-strand isoforms in the same view are rare but handled.

const COLORS = {
  cds:  "#4a5568",
  utr:  "#cbd5e0",
  axis: "#9ca3af",
  text: "#1f2937",
  muted: "#6b7280",
  nmd:   "#ef8a62",
  selected: "#2563eb",
  gencodeOnly: "#f3f4f6",
};

const NS = "http://www.w3.org/2000/svg";
const CTS = ["AT", "DD", "FB", "MV"];

/**
 * @param {HTMLElement} container
 * @param {Array<Object>} isoforms — array of isoform records from gene shard
 * @param {Object} options
 *   @param {string} selectedId — the currently selected isoform id
 *   @param {Function} onSelect — called with isoform.id when a track is clicked
 *   @param {number} trackH — per-track pixel height (default 22)
 *   @param {number} labelWidth — pixel width reserved for the left label column
 *   @param {boolean} showRuler — draw a genomic-coordinate ruler at the top
 */
export function renderStack(container, isoforms, options = {}) {
  container.innerHTML = "";
  const opts = {
    trackH: options.trackH ?? 22,
    labelWidth: options.labelWidth ?? 190,
    showRuler: options.showRuler ?? true,
    selectedId: options.selectedId ?? null,
    onSelect: options.onSelect ?? (() => {}),
  };

  if (!isoforms || isoforms.length === 0) {
    container.innerHTML = '<p style="margin:1em;color:#6b7280;font-size:0.9em">No isoforms match the current filter.</p>';
    return;
  }

  // Compute gene span across all displayed isoforms
  const starts = isoforms.map(i => i.tx_start).filter(v => v != null);
  const ends   = isoforms.map(i => i.tx_end).filter(v => v != null);
  if (starts.length === 0) {
    container.innerHTML = '<p style="margin:1em;color:#6b7280;font-size:0.9em">No genomic coordinates.</p>';
    return;
  }
  const gStart = Math.min(...starts);
  const gEnd   = Math.max(...ends);
  const gSpan  = gEnd - gStart + 1;

  // Genomic left-to-right (no strand flip). Direction of transcription is
  // conveyed by chevrons on the intron lines — Isopair convention: `>` for +,
  // `<` for −. See transcript_stack.js drawChevrons().
  const strands = isoforms.map(i => i.strand).filter(Boolean);
  const dominantStrand = strands.length > 0 ? strands[0] : "+";

  const widthPx = Math.max(container.clientWidth || 900, 500);
  const rulerH  = opts.showRuler ? 22 : 0;
  const headerH = rulerH + 4;
  const heightPx = headerH + isoforms.length * opts.trackH + 8;

  const trackAreaX = opts.labelWidth;
  const trackW     = widthPx - trackAreaX - 12;

  const g2x = (g) => trackAreaX + ((g - gStart) / gSpan) * trackW;

  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("width", widthPx);
  svg.setAttribute("height", heightPx);
  svg.setAttribute("viewBox", `0 0 ${widthPx} ${heightPx}`);

  // Coordinate ruler at the top, with strand indicator on the far right
  if (opts.showRuler) drawRuler(svg, gStart, gEnd, trackAreaX, trackW, widthPx, dominantStrand);

  // Draw each isoform track
  isoforms.forEach((iso, rowIdx) => {
    const y = headerH + rowIdx * opts.trackH;
    drawTrack(svg, iso, rowIdx, y, opts.trackH, {
      g2x, gStart, gEnd,
      trackAreaX, trackW,
      labelWidth: opts.labelWidth,
      widthPx,
      selected: iso.id === opts.selectedId,
      onSelect: opts.onSelect,
    });
  });

  container.appendChild(svg);
}

function drawRuler(svg, gStart, gEnd, x0, trackW, widthPx, strand) {
  const y = 14;
  const line = document.createElementNS(NS, "line");
  line.setAttribute("x1", x0); line.setAttribute("y1", y);
  line.setAttribute("x2", x0 + trackW); line.setAttribute("y2", y);
  line.setAttribute("stroke", COLORS.axis); line.setAttribute("stroke-width", "0.5");
  svg.appendChild(line);
  const nTicks = 5;
  for (let i = 0; i <= nTicks; i++) {
    const frac = i / nTicks;
    const g = Math.round(gStart + frac * (gEnd - gStart));
    const px = x0 + frac * trackW;
    const tick = document.createElementNS(NS, "line");
    tick.setAttribute("x1", px); tick.setAttribute("y1", y);
    tick.setAttribute("x2", px); tick.setAttribute("y2", y + 4);
    tick.setAttribute("stroke", COLORS.axis); tick.setAttribute("stroke-width", "0.5");
    svg.appendChild(tick);
    addText(svg, px, y - 3, g.toLocaleString(),
            { anchor: "middle", size: 9, color: COLORS.muted, family: "mono" });
  }
  // Strand indicator on the left-label column
  const arrow = strand === "-" ? "◀" : "▶";
  addText(svg, x0 - 6, y + 4,
          `${strand} strand ${arrow}`,
          { anchor: "end", size: 10, color: COLORS.muted, bold: true });
}

function drawTrack(svg, iso, rowIdx, y, trackH, ctx) {
  const { g2x, labelWidth, widthPx, selected, onSelect } = ctx;
  const centerY = y + trackH / 2;
  const exonStarts = [].concat(iso.exon_starts || []);
  const exonEnds   = [].concat(iso.exon_ends || []);
  const nExons = exonStarts.length;
  const isGencodeOnly = iso.gencode_only === true;
  const isNmd = iso.nmd_responsive && CTS.some(ct => iso.nmd_responsive[ct]);

  // Row background: selected → light blue; gencode-only → light grey
  if (selected || isGencodeOnly) {
    const bg = document.createElementNS(NS, "rect");
    bg.setAttribute("x", 0); bg.setAttribute("y", y);
    bg.setAttribute("width", widthPx); bg.setAttribute("height", trackH);
    bg.setAttribute("fill", selected ? "#dbeafe" : COLORS.gencodeOnly);
    bg.setAttribute("opacity", selected ? "0.6" : "0.5");
    svg.appendChild(bg);
  }

  // Left-edge accent for selected row
  if (selected) {
    const bar = document.createElementNS(NS, "rect");
    bar.setAttribute("x", 0); bar.setAttribute("y", y);
    bar.setAttribute("width", 3); bar.setAttribute("height", trackH);
    bar.setAttribute("fill", COLORS.selected);
    svg.appendChild(bar);
  }

  // Left label: NMD dot + isoform ID (truncated); tooltip = full id
  const labelX = 8;
  if (isNmd) {
    const dot = document.createElementNS(NS, "circle");
    dot.setAttribute("cx", labelX + 3); dot.setAttribute("cy", centerY);
    dot.setAttribute("r", 3); dot.setAttribute("fill", COLORS.nmd);
    dot.appendChild(makeTitle("NMD-responsive in ≥1 cell type"));
    svg.appendChild(dot);
  }
  const shortId = truncate(iso.id, 22);
  const idText = addText(svg, labelX + (isNmd ? 12 : 4), centerY + 4, shortId,
                          { anchor: "start", size: 10, color: isGencodeOnly ? COLORS.muted : COLORS.text,
                            family: "mono" });
  idText.style.cursor = "pointer";
  idText.appendChild(makeTitle(iso.id + (isGencodeOnly ? " (GENCODE-annotated, not detected)" : "")));

  if (nExons === 0) {
    // Just a placeholder line
    const dash = document.createElementNS(NS, "line");
    dash.setAttribute("x1", ctx.trackAreaX); dash.setAttribute("y1", centerY);
    dash.setAttribute("x2", ctx.trackAreaX + ctx.trackW); dash.setAttribute("y2", centerY);
    dash.setAttribute("stroke", COLORS.axis); dash.setAttribute("stroke-width", "0.5");
    dash.setAttribute("stroke-dasharray", "2,3");
    svg.appendChild(dash);
    return;
  }

  // Intron line (from leftmost exon to rightmost exon)
  const exGmin = Math.min(...exonStarts);
  const exGmax = Math.max(...exonEnds);
  const introLine = document.createElementNS(NS, "line");
  const xL = g2x(exGmin);
  const xR = g2x(exGmax);
  introLine.setAttribute("x1", xL); introLine.setAttribute("y1", centerY);
  introLine.setAttribute("x2", xR); introLine.setAttribute("y2", centerY);
  introLine.setAttribute("stroke", isGencodeOnly ? COLORS.muted : "#4b5563");
  introLine.setAttribute("stroke-width", "0.8");
  svg.appendChild(introLine);

  // Chevrons on introns indicating strand (Isopair convention)
  // For each intron between consecutive exons, place `>` (or `<`) at intervals
  drawChevronsOnIntrons(svg, exonStarts, exonEnds, iso.strand, centerY, ctx);

  // Exons (CDS taller, UTR shorter). Colours slightly desaturated for gencode-only.
  const cdsCol = isGencodeOnly ? "#9ca3af" : COLORS.cds;
  const utrCol = isGencodeOnly ? "#e5e7eb" : COLORS.utr;
  // Defensive: coerce to null if not a plain number (guards against stale
  // shards where the value is an empty object due to the jsonlite NULL→{} bug).
  const cdsStart = (typeof iso.cds_start === "number") ? iso.cds_start : null;
  const cdsStop  = (typeof iso.cds_stop  === "number") ? iso.cds_stop  : null;
  const hasCDS = cdsStart !== null && cdsStop !== null;
  const exonH  = Math.min(trackH - 6, 14);
  const utrH   = Math.round(exonH * 0.55);

  for (let i = 0; i < nExons; i++) {
    const eS = exonStarts[i], eE = exonEnds[i];
    const draw = (gFrom, gTo, isCDS) => {
      const xa = Math.min(g2x(gFrom), g2x(gTo));
      const xb = Math.max(g2x(gFrom), g2x(gTo));
      const w  = Math.max(1.5, xb - xa);
      const h  = isCDS ? exonH : utrH;
      const r = document.createElementNS(NS, "rect");
      r.setAttribute("x", xa);
      r.setAttribute("y", centerY - h / 2);
      r.setAttribute("width", w);
      r.setAttribute("height", h);
      r.setAttribute("fill", isCDS ? cdsCol : utrCol);
      r.setAttribute("stroke", "#374151");
      r.setAttribute("stroke-width", "0.3");
      r.appendChild(makeTitle(
        `Exon ${i + 1}${isCDS ? " (CDS)" : " (UTR)"} · ${gFrom.toLocaleString()}–${gTo.toLocaleString()} · ${(gTo - gFrom + 1).toLocaleString()} nt`
      ));
      svg.appendChild(r);
    };
    if (hasCDS) {
      const segs = [
        { gFrom: eS, gTo: Math.min(eE, cdsStart - 1), isCDS: false },
        { gFrom: Math.max(eS, cdsStart), gTo: Math.min(eE, cdsStop), isCDS: true },
        { gFrom: Math.max(eS, cdsStop + 1), gTo: eE, isCDS: false },
      ].filter(s => s.gFrom <= s.gTo);
      segs.forEach(s => draw(s.gFrom, s.gTo, s.isCDS));
    } else {
      draw(eS, eE, false);
    }
  }

  // Click-to-select over the whole row
  const hit = document.createElementNS(NS, "rect");
  hit.setAttribute("x", 0); hit.setAttribute("y", y);
  hit.setAttribute("width", widthPx); hit.setAttribute("height", trackH);
  hit.setAttribute("fill", "transparent");
  hit.style.cursor = "pointer";
  hit.addEventListener("click", () => onSelect(iso.id));
  svg.appendChild(hit);
}

// Isopair-style intron chevrons.
// Convention:
//   - Auto-spacing = max(3% of total genomic range on the track area, 200 bp)
//   - Arm width in genomic coords = 30% of spacing
//   - Skip introns shorter than 3× arm width (chevron would overlap exons)
//   - Chevron: `>` for + strand, `<` for − strand
//   - Draw as SVG path with two arm segments meeting at the point
function drawChevronsOnIntrons(svg, exonStarts, exonEnds, strand, centerY, ctx) {
  const nExons = exonStarts.length;
  if (nExons < 2) return;
  // Genomic-coord spacing derived from the total displayed span (not just this iso)
  const totalGSpan = (ctx.gEnd - ctx.gStart + 1);
  const spacingG   = Math.max(Math.round(totalGSpan * 0.03), 200);
  const armWidthG  = Math.round(spacingG * 0.3);
  const armHalfY   = 3.2;  // pixels
  const dir = (strand === "-") ? -1 : +1;

  // Order exons by genomic start so intron boundaries are unambiguous
  const idx = [...Array(nExons).keys()].sort((a, b) => exonStarts[a] - exonStarts[b]);
  for (let k = 0; k < idx.length - 1; k++) {
    const iA = idx[k], iB = idx[k + 1];
    const intronStartG = exonEnds[iA] + 1;
    const intronEndG   = exonStarts[iB] - 1;
    const intronLenG   = intronEndG - intronStartG + 1;
    if (intronLenG < armWidthG * 3) continue;

    const nChev = Math.max(1, Math.floor(intronLenG / spacingG));
    const stepG = intronLenG / (nChev + 1);
    for (let c = 1; c <= nChev; c++) {
      const cxG = intronStartG + c * stepG;
      const cxPx = ctx.g2x(cxG);
      const armPx = Math.max(3, (armWidthG / totalGSpan) * ctx.trackW);
      // Build the two chevron arms as a path
      const tipX = cxPx;
      const baseX = cxPx - dir * armPx;
      const path = document.createElementNS(NS, "path");
      path.setAttribute("d",
        `M ${baseX} ${centerY - armHalfY} L ${tipX} ${centerY} L ${baseX} ${centerY + armHalfY}`);
      path.setAttribute("fill", "none");
      path.setAttribute("stroke", "#4b5563");
      path.setAttribute("stroke-width", "1");
      path.setAttribute("stroke-linecap", "round");
      path.setAttribute("stroke-linejoin", "round");
      svg.appendChild(path);
    }
  }
}

// ── helpers ──
function addText(svg, x, y, str, { anchor = "start", size = 10, color = "#1f2937", family = "sans", bold = false } = {}) {
  const t = document.createElementNS(NS, "text");
  t.setAttribute("x", x);
  t.setAttribute("y", y);
  t.setAttribute("text-anchor", anchor);
  t.setAttribute("font-size", size);
  t.setAttribute("fill", color);
  t.setAttribute("font-family",
    family === "mono"
      ? "ui-monospace, SFMono-Regular, Menlo, monospace"
      : "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif");
  if (bold) t.setAttribute("font-weight", "600");
  t.textContent = str;
  svg.appendChild(t);
  return t;
}
function makeTitle(s) {
  const t = document.createElementNS(NS, "title");
  t.textContent = s;
  return t;
}
function truncate(s, n) {
  return s.length > n ? s.slice(0, n - 1) + "…" : s;
}

// Sort helper — importable from app.js so both the in-panel view and modal
// use the same ordering.
export function sortIsoforms(list) {
  return [...list].sort((a, b) => {
    const aNmd = a.nmd_responsive && CTS.some(ct => a.nmd_responsive[ct]);
    const bNmd = b.nmd_responsive && CTS.some(ct => b.nmd_responsive[ct]);
    if (aNmd !== bNmd) return aNmd ? -1 : 1;
    const aGc = a.gencode_only === true;
    const bGc = b.gencode_only === true;
    if (aGc !== bGc) return aGc ? 1 : -1;
    const am = Math.max(...CTS.map(ct => a.cpm?.[ct] ?? 0));
    const bm = Math.max(...CTS.map(ct => b.cpm?.[ct] ?? 0));
    if (bm !== am) return bm - am;
    return (b.n_exons ?? 0) - (a.n_exons ?? 0);
  });
}
