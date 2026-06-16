# isocall_dge/limma/

Limma-voom transcript-level differential expression results from the isocall pipeline.

## Layout

- `yul/` — **Yul's canonical limma run, dated 2026.3.3.** This is the official limma DGE for the manuscript (project ownership model: Yul owns limma + mashr DGE; Pete owns Isopair + the deep-learning model).
  Files (gitignored): `nmd_dge_{at,dd,ddali,do,fb,mv}_2026.3.3.csv`

## Removed 2026-06-16

`pete/` (Pete's parallel exploratory limma run, dated 2026.3.1) was deleted on 2026-06-16 as part of the confusion-cleanup pass. Yul's `yul/` is the canonical version; the parallel `pete/` directory was a methodological side-by-side that risked being mistaken for canonical. If reproduction of that run is needed in future it's regenerable from the count matrix + sample metadata.

## Why gitignored

CSVs are 22–41 MB each. They're regenerable from the count matrix + sample metadata, so they live on disk but stay out of git history (matches the existing pattern for `shortread_dge/*.csv`).
