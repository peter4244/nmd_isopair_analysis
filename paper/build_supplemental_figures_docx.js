#!/usr/bin/env node
/*
 * Build a Word document containing SF25–SF43 for the NMD manuscript
 * Supplemental Figures (§4 and §5). Format matches the existing Yul-authored
 * SF1–SF23 document (Supplemental Figures NMD.pdf):
 *
 *   - "Manuscript Supplemental Figures" bold + underlined banner at the top
 *     of the first page.
 *   - Figures flow continuously (no forced page break per SF; Word will page
 *     them naturally). A small vertical spacer sits between one figure block
 *     and the next.
 *   - Each SF: centered figure image, then a left-aligned caption paragraph
 *     that begins with **SF## | Title.** inline (bold, ending at the period
 *     that closes the title), followed by the rest of the caption prose.
 *     Multi-paragraph captions preserve their paragraph breaks.
 *
 * Output: ~/claude_projects/nmd/paper/nmd_supplemental_figures_sf24_sf42.docx
 * (Co-located with the manuscript, generator, and other paper artifacts —
 * NOT Desktop; see feedback_default_output_location memory rule.)
 */

const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, ImageRun,
  AlignmentType, PageOrientation, UnderlineType,
} = require('docx');

const SFROOT = '/Users/petecastaldi/claude_projects/nmd/figures/SupplementalFigures';
const OUTPUT = '/Users/petecastaldi/claude_projects/nmd/paper/nmd_supplemental_figures_sf24_sf42.docx';

// Manifest: manuscript order → dir + png filename + legend filename
const SFS = [
  ['SF25_SpliceEventCategories',        'figure_sf25_splice_event_categories.png',       'figure_sf25_splice_event_categories_legend.md'],
  ['SF26_IsoformsPerGene',              'figure_sf26_isoforms_per_gene.png',             'figure_sf26_isoforms_per_gene_legend.md'],
  ['SF27_ReferenceShare',               'figure_sf27_reference_share.png',               'figure_sf27_reference_share_legend.md'],
  ['SF28_TranscriptLengthByRole',       'figure_sf28_transcript_length_by_role.png',     'figure_sf28_transcript_length_by_role_legend.md'],
  ['SF29_PairAnalysisFlowchart',        'figure_s_pair_analysis_flowchart.png',          'figure_s_pair_analysis_flowchart_legend.md'],
  ['SF30_GainDirectionByEvent',         'figure_sf30_gain_direction_by_event.png',       'figure_sf30_gain_direction_by_event_legend.md'],
  ['SF31_PTCDistanceDoseResponse',      'figure_sf31_ptc_distance_dose_response.png',    'figure_sf31_ptc_distance_dose_response_legend.md'],
  ['SF32_NMDEffectByEJCCount',          'figure_sf32_nmd_effect_by_ejc_count.png',       'figure_sf32_nmd_effect_by_ejc_count_legend.md'],
  ['SF33_CdsAnd3UTR_GENCODE',           'figure_sf33.png',                               'figure_sf33_legend.md'],
  ['SF34_TD2Bias_broad',                'figure_sf34.png',                               'figure_sf34_legend.md'],
  ['SF35_TD2Bias_occult',               'figure_sf35.png',                               'figure_sf35_legend.md'],
  ['SF36_CdsAnd3UTR_refAUG',            'figure_sf36.png',                               'figure_sf36_legend.md'],
  ['SF37_ShapAcrossWindows',            'figure_sf37_shap_across_windows.png',           'figure_sf37_shap_across_windows_legend.md'],
  ['SF38_StopCodonUsage',               'figure_s_stop_codon_usage.png',                 'figure_s_stop_codon_usage_legend.md'],
  ['SF39_AttentionDistribution',        'figure_s_attention_distribution.png',           'figure_s_attention_distribution_legend.md'],
  ['SF40_PTCSubclassBranchSHAP',        'figure_s_branch_shap_by_subclass.png',          'figure_s_branch_shap_by_subclass_legend.md'],
  ['SF41_PTCSubclassPerformance',       'figure_s_performance_by_subclass.png',          'figure_s_performance_by_subclass_legend.md'],
  ['SF42_GCcontentStopWindow',          'figure_sf42_gc_content_stop_window.png',        'figure_sf42_gc_content_stop_window_legend.md'],
  ['SF43_ModelComparison',              'figure_s_model_comparison.png',                 'figure_s_model_comparison_legend.md'],
];

// ── page geometry ─────────────────────────────────────────────────────
// US Letter portrait, 1" margins → 6.5" content width
const CONTENT_WIDTH_IN = 6.5;

// PNG dimensions via IHDR chunk (avoid extra image-size dep).
function pngSize(filepath) {
  const buf = fs.readFileSync(filepath);
  return { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) };
}

// Split a legend into { header, bodyParagraphs[] }.
// The legend begins with **SF## | Title.** followed by body prose. The body
// may itself be one or more paragraphs (blank-line separated).
function parseLegend(md) {
  const trimmed = md.trim();
  const m = trimmed.match(/^\*\*(.+?)\*\*\s*(.*)$/s);
  if (!m) return { header: null, bodyParagraphs: [trimmed] };
  const header = m[1].trim();
  const body   = m[2].trim();
  const bodyParagraphs = body.split(/\n\s*\n/).map(p => p.replace(/\n/g, ' ').trim());
  return { header, bodyParagraphs };
}

// Convert an inline-markdown string into a run of TextRuns preserving
// **bold**, *italic*, and `code` markers. Enough for our legends.
function markdownRuns(text) {
  const out = [];
  let i = 0;
  while (i < text.length) {
    if (text[i] === '\\' && i + 1 < text.length) {
      // Backslash escape: emit the next character literally (e.g. \* -> *),
      // so significance markers like *p / **p / ***p render as plain asterisks.
      out.push(new TextRun(text[i + 1]));
      i += 2;
    } else if (text.startsWith('**', i)) {
      const end = text.indexOf('**', i + 2);
      if (end === -1) { out.push(new TextRun(text.slice(i))); break; }
      out.push(new TextRun({ text: text.slice(i + 2, end), bold: true }));
      i = end + 2;
    } else if (text[i] === '*') {
      const end = text.indexOf('*', i + 1);
      if (end === -1) { out.push(new TextRun(text.slice(i))); break; }
      out.push(new TextRun({ text: text.slice(i + 1, end), italics: true }));
      i = end + 1;
    } else if (text[i] === '`') {
      const end = text.indexOf('`', i + 1);
      if (end === -1) { out.push(new TextRun(text.slice(i))); break; }
      out.push(new TextRun({ text: text.slice(i + 1, end), font: 'Courier New' }));
      i = end + 1;
    } else {
      const nextEsc    = text.indexOf('\\', i);
      const nextBold   = text.indexOf('**', i);
      const nextItalic = text.indexOf('*',  i);
      const nextCode   = text.indexOf('`',  i);
      const candidates = [nextEsc, nextBold, nextItalic, nextCode].filter(x => x !== -1);
      const next = candidates.length ? Math.min(...candidates) : text.length;
      out.push(new TextRun(text.slice(i, next)));
      i = next;
    }
  }
  return out;
}

// Build the flowing content.
const children = [];

// Top banner (Yul-style): bold + underlined "Manuscript Supplemental Figures".
children.push(new Paragraph({
  spacing: { after: 240 },
  children: [new TextRun({
    text: 'Manuscript Supplemental Figures',
    bold: true,
    underline: { type: UnderlineType.SINGLE },
  })],
}));

// Per-SF content
for (let idx = 0; idx < SFS.length; idx++) {
  const [dir, pngName, legendName] = SFS[idx];
  const pngPath    = path.join(SFROOT, dir, pngName);
  const legendPath = path.join(SFROOT, dir, legendName);

  const { w, h } = pngSize(pngPath);
  const aspect = w / h;
  const widthPx  = Math.round(CONTENT_WIDTH_IN * 96); // docx-js takes pixels @ 96 DPI
  const heightPx = Math.round(widthPx / aspect);

  const legendMd = fs.readFileSync(legendPath, 'utf8');
  const { header, bodyParagraphs } = parseLegend(legendMd);

  // Centered figure image.
  children.push(new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 120, after: 60 },
    children: [new ImageRun({
      type: 'png',
      data: fs.readFileSync(pngPath),
      transformation: { width: widthPx, height: heightPx },
      altText: { title: header || dir, description: header || dir, name: dir },
    })],
  }));

  // Caption paragraph: **SF## | Title.** inline, then first body paragraph
  // in the SAME paragraph — matches Yul's format where the header runs
  // inline with the description.
  const headerText = header
    ? (header.replace(/\.\.+$/, '.') + (header.endsWith('.') ? '' : '.'))
    : '';
  const firstBody = bodyParagraphs[0] || '';
  const firstRuns = [
    new TextRun({ text: headerText + ' ', bold: true }),
    ...markdownRuns(firstBody),
  ];
  children.push(new Paragraph({
    spacing: { after: 120 },
    children: firstRuns,
  }));

  // Additional body paragraphs (multi-paragraph captions).
  for (let j = 1; j < bodyParagraphs.length; j++) {
    children.push(new Paragraph({
      spacing: { after: 120 },
      children: markdownRuns(bodyParagraphs[j]),
    }));
  }

  // Vertical spacer between SF blocks — matches Yul's small gap.
  if (idx < SFS.length - 1) {
    children.push(new Paragraph({ spacing: { after: 240 }, children: [] }));
  }
}

// Assemble the document.
const doc = new Document({
  styles: {
    default: {
      document: { run: { font: 'Arial', size: 22 } }, // 11 pt body
    },
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840, orientation: PageOrientation.PORTRAIT },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children,
  }],
});

Packer.toBuffer(doc).then(buffer => {
  fs.writeFileSync(OUTPUT, buffer);
  console.log(`wrote ${OUTPUT}`);
  console.log(`size: ${(buffer.length / 1024).toFixed(1)} KB`);
});
