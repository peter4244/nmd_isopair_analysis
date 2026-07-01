// app.js — main controller for the NMD lung atlas.
import { renderTranscript } from "./transcript_viz.js";
import { renderStack, sortIsoforms } from "./transcript_stack.js";

const IN_PANEL_STACK_CAP = 8;

const CTS = ["AT", "DD", "FB", "MV"];
const CT_FULL = { AT: "Alveolar type 2", DD: "Large airway epithelial",
                   FB: "Fibroblast", MV: "Microvascular endothelial" };
const CT_SHORT = { AT: "AT2", DD: "LAE", FB: "Fib", MV: "MVE" };
const CT_COLOR = { AT: "#4c72b0", DD: "#dd8452", FB: "#55a868", MV: "#c44e52" };
const NMD_COLOR = "#ef8a62";   // canonical NMD palette
const PLOT_FONT = { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif", size: 11, color: "#1f2937" };

const state = {
  genesIndex: [],
  quantiles: null,
  manifest: null,
  currentGene: null,    // gene shard JSON
  currentIso: null,     // selected isoform record within currentGene
  filter: {
    minCpm: 1,          // default threshold — matches the mashr-DIE inclusion floor
    showAnnotated: false,
  },
};

// ── boot ──
init();

async function init() {
  try {
    // Fetch manifest first so we can use its data_version to bust caches on
    // the other JSON payloads (gene index + shards).
    const mf = await fetch("data/manifest.json", { cache: "no-store" }).then(r => r.json());
    const v = encodeURIComponent(mf.data_version || Date.now());
    state.manifest = mf;
    const [idx, qs] = await Promise.all([
      fetch(`data/genes_index.json?v=${v}`).then(r => r.json()),
      fetch(`data/quantiles.json?v=${v}`).then(r => r.json()),
    ]);
    state.genesIndex = idx;
    state.quantiles  = qs;
    document.getElementById("manifest-info").textContent =
      `${mf.n_genes.toLocaleString()} genes · ${mf.n_isoforms.toLocaleString()} isoforms · v${mf.data_version}`;
    wireSearch();
    wireWelcomeLinks();
    wireNav();
    wireFilterBar();
    wireStackModal();
  } catch (err) {
    document.getElementById("manifest-info").textContent = "data load failed";
    console.error(err);
  }
}

// ── search ──
function wireSearch() {
  const input = document.getElementById("gene-search");
  const results = document.getElementById("search-results");
  let timer = null;

  input.addEventListener("input", () => {
    clearTimeout(timer);
    timer = setTimeout(() => runSearch(input.value.trim()), 80);
  });
  input.addEventListener("keydown", (e) => {
    const items = [...results.querySelectorAll("li")];
    const idx = items.findIndex(li => li.classList.contains("highlight"));
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      const next = e.key === "ArrowDown" ? Math.min(idx + 1, items.length - 1) : Math.max(idx - 1, 0);
      items.forEach(li => li.classList.remove("highlight"));
      if (items[next]) items[next].classList.add("highlight");
    } else if (e.key === "Enter") {
      const sel = items.find(li => li.classList.contains("highlight")) || items[0];
      if (sel) selectGene(sel.dataset.geneId);
    } else if (e.key === "Escape") {
      results.innerHTML = "";
    }
  });
  document.addEventListener("click", (e) => {
    if (!e.target.closest(".search-section")) results.innerHTML = "";
  });
}

function runSearch(q) {
  const results = document.getElementById("search-results");
  if (!q) { results.innerHTML = ""; return; }
  const lc = q.toLowerCase();
  const hits = state.genesIndex
    .filter(g => (g.hgnc_symbol && g.hgnc_symbol.toLowerCase().includes(lc)) ||
                 g.gene_id.toLowerCase().includes(lc))
    .slice(0, 20);
  results.innerHTML = "";
  for (const g of hits) {
    const li = document.createElement("li");
    li.dataset.geneId = g.gene_id;
    li.innerHTML = `<span class="sym">${escapeHtml(g.hgnc_symbol || "(no symbol)")}</span>` +
                   `<span class="gid">${escapeHtml(g.gene_id)}</span>` +
                   (g.any_nmd_iso ? '<span class="badge">NMD</span>' : '') +
                   `<span class="gid"> · ${g.n_isoforms} isoform${g.n_isoforms > 1 ? "s" : ""}</span>`;
    li.addEventListener("click", () => selectGene(g.gene_id));
    results.appendChild(li);
  }
}

function wireWelcomeLinks() {
  document.querySelectorAll(".welcome a[data-gene]").forEach(a => {
    a.addEventListener("click", (e) => { e.preventDefault(); selectGene(a.dataset.gene); });
  });
}

// ── nav ──
function wireNav() {
  const goHome = (e) => { if (e) e.preventDefault(); resetToHome(); };
  document.getElementById("nav-home").addEventListener("click", goHome);
  document.getElementById("nav-home-title").addEventListener("click", goHome);
}
function resetToHome() {
  state.currentGene = null;
  state.currentIso  = null;
  // Reset filter UI so revisiting a gene from Home doesn't inherit the previous state
  state.filter.showAnnotated = false;
  const chk = document.getElementById("show-gencode-all");
  if (chk) chk.checked = false;
  document.getElementById("gene-panel").classList.add("hidden");
  document.getElementById("welcome").classList.remove("hidden");
  document.getElementById("gene-search").value = "";
  document.getElementById("search-results").innerHTML = "";
  // Also close the modal if it happens to be open
  const modal = document.getElementById("stack-modal");
  if (modal && !modal.classList.contains("hidden")) {
    modal.classList.add("hidden");
    document.body.style.overflow = "";
  }
  window.scrollTo({ top: 0, behavior: "smooth" });
}

// ── filter bar ──
function wireFilterBar() {
  const slider = document.getElementById("cpm-filter");
  const label  = document.getElementById("cpm-filter-value");
  const chk    = document.getElementById("show-gencode-all");
  slider.addEventListener("input", () => {
    state.filter.minCpm = parseFloat(slider.value);
    label.textContent = formatCpmThreshold(state.filter.minCpm);
    if (state.currentGene) renderGenePanel();
  });
  chk.addEventListener("change", () => {
    state.filter.showAnnotated = chk.checked;
    if (state.currentGene) renderGenePanel();
  });
  // Initialise label
  label.textContent = formatCpmThreshold(state.filter.minCpm);
}
function formatCpmThreshold(v) {
  return v === 0 ? "0 (all)" : v.toFixed(1).replace(/\.0$/, "");
}

// ── stacked-view modal ──
function wireStackModal() {
  const btn = document.getElementById("stack-expand");
  const closeBtn = document.getElementById("modal-close");
  const modal = document.getElementById("stack-modal");
  btn.addEventListener("click", openStackModal);
  closeBtn.addEventListener("click", closeStackModal);
  modal.addEventListener("click", (e) => {
    if (e.target === modal) closeStackModal();
  });
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && !modal.classList.contains("hidden")) closeStackModal();
  });
}
function openStackModal() {
  if (!state.currentGene) return;
  const modal = document.getElementById("stack-modal");
  modal.classList.remove("hidden");
  renderStackInto("modal-stack", getFilteredIsoforms(), { trackH: 28 });
  document.body.style.overflow = "hidden";
}
function closeStackModal() {
  document.getElementById("stack-modal").classList.add("hidden");
  document.body.style.overflow = "";
}

// ── gene selection ──
async function selectGene(geneId) {
  try {
    const v = encodeURIComponent(state.manifest?.data_version || Date.now());
    const r = await fetch(`data/genes/${encodeURIComponent(geneId)}.json?v=${v}`);
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    const shard = await r.json();
    state.currentGene = shard;
    state.currentIso  = pickInitialIso(shard);
    document.getElementById("welcome").classList.add("hidden");
    document.getElementById("gene-panel").classList.remove("hidden");
    document.getElementById("search-results").innerHTML = "";
    document.getElementById("gene-search").value = shard.hgnc_symbol || shard.gene_id;
    renderGenePanel();
  } catch (err) {
    alert(`Could not load gene ${geneId}: ${err.message}`);
    console.error(err);
  }
}

function pickInitialIso(shard) {
  // Prefer a DETECTED isoform (skip GENCODE-only rows), highest max-CPM first.
  const detected = shard.isoforms.filter(i => !i.gencode_only);
  const pool = detected.length > 0 ? detected : shard.isoforms;
  const ranked = [...pool].sort((a, b) => {
    const aMax = Math.max(...CTS.map(ct => (a.cpm?.[ct] ?? 0)));
    const bMax = Math.max(...CTS.map(ct => (b.cpm?.[ct] ?? 0)));
    return bMax - aMax;
  });
  return ranked[0];
}

function renderGenePanel() {
  const { currentGene, currentIso } = state;
  document.getElementById("gene-symbol").textContent = currentGene.hgnc_symbol || "(no symbol)";
  document.getElementById("gene-id").textContent = currentGene.gene_id;

  renderIsoformTable();
  renderSelectedIso();
  renderStackedView();
  // If the modal is open, keep it in sync with the current filter
  const modal = document.getElementById("stack-modal");
  if (modal && !modal.classList.contains("hidden")) {
    renderStackInto("modal-stack", getFilteredIsoforms(), { trackH: 28 });
  }
}

function renderIsoformTable() {
  const tbl = document.getElementById("isoform-table");
  const shown = getFilteredIsoforms();
  const detected = state.currentGene.isoforms.filter(i => !i.gencode_only);
  const annotatedOnly = state.currentGene.isoforms.filter(i => i.gencode_only === true);
  const shownDetected  = shown.filter(i => !i.gencode_only);
  const shownAnnotated = shown.filter(i => i.gencode_only === true);

  // Ensure the currently selected isoform is still in view (if filtered out, pick a new default)
  if (!shown.find(i => i.id === state.currentIso?.id)) {
    state.currentIso = shown[0] || null;
    if (state.currentIso) {
      renderSelectedIso();
    } else {
      // Nothing to show — clear the detail panes so they don't stay stale
      document.getElementById("iso-id").innerHTML = "—";
      document.getElementById("iso-meta").innerHTML =
        '<div class="hint">No isoforms match the current filter. Slide Min. CPM lower or check "Also show annotated but not detected."</div>';
      document.getElementById("transcript-viz").innerHTML = "";
      document.getElementById("cpm-chart").innerHTML = "";
      document.getElementById("logfc-chart").innerHTML = "";
    }
  }

  // Count badge
  const badge = document.getElementById("iso-count-badge");
  const totalDetected = detected.length;
  const totalAnnotated = annotatedOnly.length;
  badge.textContent =
    `${shownDetected.length} of ${totalDetected} detected` +
    (state.filter.showAnnotated && totalAnnotated ? ` + ${shownAnnotated.length} of ${totalAnnotated} annotated` : "");

  const header = `<thead><tr>
    <th>Isoform</th><th>Exons</th><th class="numeric">Max CPM</th>
    <th>NMD+ CTs</th>
    </tr></thead>`;
  const rows = shown.map(iso => {
    const nmdCts = iso.nmd_responsive
      ? (CTS.filter(ct => iso.nmd_responsive[ct]).join(",") || "—")
      : "—";
    const mx = maxCpm(iso);
    const sel = (iso.id === state.currentIso?.id) ? "selected" : "";
    const gc = iso.gencode_only ? "gencode-only" : "";
    const isNmd = nmdCts !== "—";
    return `<tr class="${sel} ${gc}" data-iso="${escapeHtml(iso.id)}">
      <td class="iso-id">${escapeHtml(iso.id)}</td>
      <td class="numeric">${iso.n_exons ?? "—"}</td>
      <td class="numeric">${mx > 0 ? mx.toFixed(2) : "—"}</td>
      <td class="${isNmd ? "nmd-true" : "nmd-false"}">${nmdCts}</td>
    </tr>`;
  }).join("");
  tbl.innerHTML = header + `<tbody>${rows}</tbody>`;
  tbl.querySelectorAll("tr[data-iso]").forEach(tr => {
    tr.addEventListener("click", () => {
      const iso = state.currentGene.isoforms.find(i => i.id === tr.dataset.iso);
      if (iso) { state.currentIso = iso; renderGenePanel(); }
    });
  });
}
function maxCpm(iso) { return Math.max(...CTS.map(ct => iso.cpm?.[ct] ?? 0)); }

// Single source of truth for what the filter shows — used by the table AND the
// stacked view + modal so they stay in sync.
function getFilteredIsoforms() {
  if (!state.currentGene) return [];
  const { minCpm, showAnnotated } = state.filter;
  const all = state.currentGene.isoforms;
  const detected = all.filter(i => !i.gencode_only)
                       .filter(i => maxCpm(i) >= minCpm);
  const annotated = showAnnotated
    ? all.filter(i => i.gencode_only === true)
    : [];
  return sortIsoforms([...detected, ...annotated]);
}

function renderStackedView() {
  const all = getFilteredIsoforms();
  const inPanel = all.slice(0, IN_PANEL_STACK_CAP);

  const badge = document.getElementById("stack-count-badge");
  const capNote = document.getElementById("stack-cap-note");
  const btn = document.getElementById("stack-expand");
  badge.textContent = `${all.length} matching`;
  if (all.length > IN_PANEL_STACK_CAP) {
    capNote.textContent = `Top ${IN_PANEL_STACK_CAP} shown — click Compare all for the rest.`;
    btn.style.display = "";
  } else {
    capNote.textContent = "";
    btn.style.display = all.length > 6 ? "" : "none"; // still let user pop-out for scale
  }

  renderStackInto("transcript-stack", inPanel, { trackH: 22 });
}

function renderStackInto(elId, isoforms, opts = {}) {
  const el = document.getElementById(elId);
  renderStack(el, isoforms, {
    selectedId: state.currentIso?.id,
    trackH: opts.trackH ?? 22,
    onSelect: (id) => {
      const iso = state.currentGene.isoforms.find(i => i.id === id);
      if (!iso) return;
      state.currentIso = iso;
      // Re-render everything that reflects selection
      renderGenePanel();
      // If modal is open, keep it in sync too
      const modal = document.getElementById("stack-modal");
      if (!modal.classList.contains("hidden")) {
        renderStackInto("modal-stack", getFilteredIsoforms(), { trackH: 28 });
      }
    },
  });
}

function renderSelectedIso() {
  const iso = state.currentIso;
  if (!iso) return;
  const isGencodeOnly = iso.gencode_only === true;

  // Header with GENCODE-only badge
  const idEl = document.getElementById("iso-id");
  idEl.innerHTML = escapeHtml(iso.id) +
    (isGencodeOnly ? ' <span class="badge-muted">GENCODE-annotated, not detected in our data</span>' : '');

  const meta = document.getElementById("iso-meta");
  const dash = (v) => (v === null || v === undefined || v === "") ? "—" : String(v);
  meta.innerHTML = [
    ["Chromosome",       dash(iso.chr)],
    ["Strand",           dash(iso.strand)],
    ["Length (bp)",
      iso.tx_start && iso.tx_end
        ? sumExonLengths(iso).toLocaleString() + " nt"
        : "—"],
    ["Exons",            dash(iso.n_exons)],
    ["GENCODE biotype",  dash(iso.gencode_biotype)],
    ["GENCODE tags",     dash(iso.tags)],
  ].map(([k, v]) => `<div><span class="label">${k}:</span> <span class="value">${escapeHtml(v)}</span></div>`).join("");

  // CDS provenance banner (below the meta grid) — always visible; source-labelled + coloured
  const cdsSrc = iso.cds_source || "none";
  const cdsDetail = iso.cds_source_detail || "";
  const cdsMap = {
    gencode_v49:            { label: "GENCODE v49",                   klass: "cds-src-gencode" },
    gencode_nf_placeholder: { label: "GENCODE placeholder (cds_start_NF / mRNA_start_NF)", klass: "cds-src-gencode-nf" },
    ref_aug_projection:     { label: "Reference-AUG projection",       klass: "cds-src-refaug" },
    td2:                    { label: "TD2 / SQANTI3 (caution)",        klass: "cds-src-td2" },
    none:                   { label: "No CDS annotated",               klass: "cds-src-none" },
  };
  const cdsInfo = cdsMap[cdsSrc] || cdsMap.none;
  const cdsPanel = document.getElementById("cds-provenance");
  cdsPanel.className = `cds-provenance ${cdsInfo.klass}`;
  cdsPanel.innerHTML =
    `<span class="cds-src-label">CDS source: ${escapeHtml(cdsInfo.label)}</span>` +
    (cdsDetail ? `<span class="cds-src-detail">${escapeHtml(cdsDetail)}</span>` : "");

  renderTranscript(document.getElementById("transcript-viz"), iso);
  if (isGencodeOnly) {
    renderNotDetectedPlaceholder("cpm-chart", "Not detected — no expression data for this transcript in our long-read dataset.");
    renderNotDetectedPlaceholder("logfc-chart", "Not detected — no NMD-response measurement for this transcript.");
  } else {
    renderCpmChart(iso);
    renderLogfcChart(iso);
  }
}
function sumExonLengths(iso) {
  const s = [].concat(iso.exon_starts || []);
  const e = [].concat(iso.exon_ends   || []);
  let total = 0;
  for (let i = 0; i < s.length; i++) total += e[i] - s[i] + 1;
  return total;
}
function renderNotDetectedPlaceholder(id, msg) {
  const el = document.getElementById(id);
  el.innerHTML = `<div class="placeholder-panel">${escapeHtml(msg)}</div>`;
}

// ── CPM chart with per-CT decile rank ──
// Each x-axis label includes the cell-type short name AND the decile rank
// underneath, so the rank info doesn't fight with the bar-value labels.
function renderCpmChart(iso) {
  const vals  = CTS.map(ct => iso.cpm[ct] ?? 0);
  const ranks = CTS.map(ct => deciles(state.quantiles[ct], iso.cpm[ct]));
  const tickLabels = CTS.map((ct, i) => {
    const r = ranks[i] != null ? `D${ranks[i]}` : "—";
    return `<b>${CT_SHORT[ct]}</b><br><span style="font-size:9px;color:#6b7280">${r}</span>`;
  });
  const colors = CTS.map(ct => CT_COLOR[ct]);
  Plotly.newPlot("cpm-chart", [{
    x: CTS, y: vals, type: "bar",
    marker: { color: colors },
    customdata: CTS.map((ct, i) => [CT_FULL[ct], ranks[i]]),
    hovertemplate: "<b>%{customdata[0]}</b><br>CPM: %{y:.2f}<br>Decile: D%{customdata[1]}<extra></extra>",
  }], {
    margin: { t: 14, r: 10, b: 58, l: 64 },
    xaxis: {
      tickvals: CTS, ticktext: tickLabels,
      tickfont: PLOT_FONT,
    },
    yaxis: {
      title: { text: "Expression (CPM)", standoff: 10, font: { ...PLOT_FONT, size: 12 } },
      tickfont: PLOT_FONT,
      rangemode: "tozero",
    },
    height: 230,
    showlegend: false,
    font: PLOT_FONT,
  }, { displayModeBar: false, responsive: true });
}

function deciles(qs, v) {
  if (qs == null || v == null || v <= 0) return null;
  const probs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
  // qs is an object { "0%":..., "10%":..., ... } from R's quantile()
  const arr = probs.map(p => qs[`${(p * 100).toFixed(0)}%`] ?? null);
  for (let i = 1; i < arr.length; i++) {
    if (v <= arr[i]) return i;  // decile 1..10
  }
  return 10;
}

// ── log2FC chart with NMD-responsive flag ──
// NMD-responsive cell types coloured orange and marked with a star inside the
// x-axis label; the bar value is shown on hover only (no double-labeling above).
function renderLogfcChart(iso) {
  const vals = CTS.map(ct => iso.logfc[ct] ?? 0);
  const responsive = CTS.map(ct => iso.nmd_responsive[ct] === true);
  const colors = responsive.map(r => r ? NMD_COLOR : "#9ca3af");
  const tickLabels = CTS.map((ct, i) => {
    const star = responsive[i] ? '<span style="color:#ef8a62">★</span>' : '';
    return `<b>${CT_SHORT[ct]}</b>${star ? "<br>" + star : ""}`;
  });
  Plotly.newPlot("logfc-chart", [{
    x: CTS, y: vals, type: "bar",
    marker: { color: colors },
    customdata: CTS.map((ct, i) => [CT_FULL[ct], responsive[i] ? "yes" : "no"]),
    hovertemplate: "<b>%{customdata[0]}</b><br>log₂FC: %{y:.2f}<br>NMD-responsive: %{customdata[1]}<extra></extra>",
  }], {
    margin: { t: 14, r: 10, b: 58, l: 64 },
    xaxis: {
      tickvals: CTS, ticktext: tickLabels,
      tickfont: PLOT_FONT,
    },
    yaxis: {
      title: { text: "log₂ fold-change", standoff: 10, font: { ...PLOT_FONT, size: 12 } },
      tickfont: PLOT_FONT,
      zeroline: true, zerolinecolor: "#9ca3af",
    },
    height: 230,
    showlegend: false,
    font: PLOT_FONT,
  }, { displayModeBar: false, responsive: true });
}

// ── tiny helpers ──
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}
