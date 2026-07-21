# Environment manifest — versions that produced the manuscript figures

Captured 2026-07-20 from Pete's laptop, the machine that rendered the tracked
figure PNGs. **This is a record, not a lockfile.** Full `renv::init()` is deferred to
the curated tree (see consolidation plan Phase 6.5 rationale below).

## Why a manifest now, `renv` later
Running `renv::init()` against the current repo would snapshot the dependencies of the
~200 files the prune is about to delete — locking in exactly the cruft being removed.
The irreplaceable information is *which versions produced the published figures*, so we
capture that now and generate the real lockfile from the curated tree at snapshot time.

## Corroboration
The tracked figure PNGs embed `Matplotlib version3.11.0`, matching the matplotlib below.
This is what makes the Phase 5.2 **PNG byte-compare** gate viable: the renderer on this
machine is the one that produced the committed figures.

## R
```
R 4.5.2 (2025-10-31)   |   Bioconductor 3.22
```
| Package | Version | | Package | Version |
|---|---|---|---|---|
| data.table | 1.18.0 | | msigdbr | 25.1.1 |
| dplyr | 1.2.1 | | htmltools | 0.5.9 |
| ggplot2 | 4.0.1 | | readxl | 1.4.5 |
| edgeR | 4.8.2 | | knitr | 1.51 |
| tidyr | 1.3.1 | | jsonlite | 2.0.0 |
| tibble | 3.3.1 | | fgsea | 1.36.2 |
| patchwork | 1.3.2 | | Biostrings | 2.78.0 |
| scales | 1.4.0 | | AnnotationDbi | 1.72.0 |
| Isopair | 0.99.4 | | pathview | 1.50.0 |
| reshape2 | 1.4.5 | | GenomicRanges | 1.62.1 |
| org.Hs.eg.db | 3.22.0 | | rtracklayer | 1.70.1 |
| mashr | 0.2.79 | | matrixStats | 1.5.0 |
| limma | 3.66.0 | | tidyverse | 2.0.0 |
| ashr | 2.2.63 | | DT | 0.34.0 |

**`topGO` — NOT INSTALLED locally.** Used by `code/upstream/productive_response.Rmd`
(Yul runs it on Channing against `/udd/reyle/Rlibs`). Consequence: that chunk cannot be
re-run or verified on this laptop. Either install it here or record Channing's version
from Yul — otherwise the §3 GO enrichment is unreproducible outside her environment.

## Python
```
Python 3.14.4
```
| Package | Version |
|---|---|
| matplotlib | 3.11.0 |
| numpy | 2.4.6 |
| pandas | 3.0.3 |
| scipy | 1.18.0 |
| pillow (PIL) | 12.2.0 |
| seaborn | 0.13.2 |

## External tools (pinned in Methods text; versions not yet verified against runs)
nf-core/rnaseq v3.14.0 · Nextflow 24.04.4 · PacBio Isocall · SQANTI3 · TransDecoder2 (TD2)
· minimap2. **Action (plan 6.5):** confirm each against the actual run logs, and cache the
`pathview` KEGG `hsa04141` template in-repo — it is currently fetched from the network at
runtime and will drift.

## Not captured here
Yul's Channing environment (`/udd/reyle/Rlibs`), which is what actually produced the
§1–3 results. Her scripts set `.libPaths()` to it. **Ask Yul for `sessionInfo()` from the
run that generated the manuscript numbers** — the versions above are Pete-side and may
differ from hers.
