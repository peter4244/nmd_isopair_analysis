# NMD Lung Atlas — static-site implementation

Public, interactive browser for the 4-cell-type NMD long-read dataset
accompanying *Leshem et al.* 2026. **Architecture:** pure static site
(HTML/CSS/vanilla JS) with per-gene JSON shards — no server, no
database. Hosts on Cloudflare Pages / Netlify / GitHub Pages for $0
recurring.

## Layout

```
code/nmd_atlas/
├── README.md                  ← this file
├── export_atlas_data.R        ← R script: cache + mashr → JSON
└── public/                    ← what gets deployed
    ├── index.html
    ├── style.css
    ├── app.js                 ← main controller
    ├── transcript_viz.js      ← SVG transcript-structure component
    └── data/
        ├── genes_index.json   ← gene→symbol lookup (small, loaded upfront)
        ├── quantiles.json     ← per-CT CPM decile breakpoints
        ├── manifest.json
        └── genes/             ← one JSON shard per gene (lazy-loaded)
            └── <gene_id>.json
```

## Build the data

```bash
cd code/nmd_atlas

# Small test scope (10 featured genes; fast iteration on UI)
Rscript export_atlas_data.R --test

# Full scope (~25k genes, ~80k isoforms; ~minutes)
Rscript export_atlas_data.R
```

Outputs land under `public/data/`. The R script reads:

- `nmd_orf_model_v5_4ct/tmp/ref_orf_ptc_cache_with_nmd_2026.5.12.rds`
- `results/.../structures.rds`, `cds.rds`
- `isocall_dge/mashr/nmd_mashr_die_{at,dd,fb,mv}_2026.3.10.csv`
- `isocall_dge/mashr/mashr_isoform_lfsr_2026.3.10.csv`

Scope = isoforms with both a long-read-resolved structure AND any
expression/mashr signal (~80k isoforms across ~25k genes).

## Run locally

```bash
cd code/nmd_atlas/public
python3 -m http.server 8080
# open http://localhost:8080/
```

No build step; edit `app.js`, `style.css`, or `transcript_viz.js` and
refresh. Search any gene by symbol or Ensembl ID, or click a featured
gene link.

## Deploy

Any static host. Recommended: **Cloudflare Pages** (free, generous
limits, fast CDN, no maintenance for the manuscript-citation horizon).

```bash
# One-time setup with Wrangler CLI
npm install -g wrangler
wrangler pages project create nmd-lung-atlas
wrangler pages deploy public --project-name nmd-lung-atlas
```

Alternative one-line deploys:
- **Netlify**: `npx netlify-cli deploy --dir=public --prod`
- **GitHub Pages**: push `public/` to `gh-pages` branch
- **Vercel**: `npx vercel deploy public --prod`

## Cost projection

| Item | Cost |
|---|---|
| Hosting (Cloudflare Pages free tier) | $0 / mo |
| Custom domain (optional) | $10–15 / yr |
| Bandwidth at 500 hits/mo × ~2 MB/visitor | <2 GB/mo (well within all free tiers) |

Full data footprint at 25k-gene scope: ~200–500 MB on disk (one JSON shard per gene, ~20 KB average, plus a ~2 MB `genes_index.json` loaded on init). Each visitor fetches only the shards for the genes they look up, so per-visitor bandwidth is dominated by the index + a handful of shards (~3 MB uncompressed, ~500 KB gzipped).

## Page layout

- **Search**: substring match on gene symbol + Ensembl ID, with NMD-isoform badges.
- **Isoform table** (left): one row per isoform in the selected gene
  with max-CPM and the cell-type list of NMD-responsive calls.
- **Selected isoform** (right):
  - **Transcript structure** — SVG track with exons (CDS shaded dark,
    UTR shaded light), introns as lines with strand arrows, ATG +
    STOP annotated, hover tooltips on each segment.
  - **Expression (CPM)** — per-CT bar chart with each cell type's
    decile rank (D1–D10) annotated above each bar.
  - **NMD response** — per-CT log₂FC bar chart, NMD-responsive cell
    types coloured orange and starred (lfsr < 0.05 + positive
    posterior mean).

## Extending

- **More cell types** — add to the `CTS` constant in `export_atlas_data.R`,
  `app.js`, and `style.css`.
- **Better search** — drop in `fuse.js` (CDN) for fuzzy matching.
- **Comparison view** — small multiples of multiple isoforms within a
  gene. Hook into `state.currentGene.isoforms` in `app.js`.
- **Cross-gene queries** — `manifest.json` + `genes_index.json` are
  small; a future "compare across genes" page can load all shards on
  demand.

## Data versioning

`manifest.json` carries `data_version` and `generated_at`. Bump
`data_version` whenever you regenerate the data; the front-end footer
displays it so users can confirm which dataset they're viewing.
