# Data inputs needed to run Yul's upstream scripts (D1 inventory)

Extracted 2026-07-20 from the 22 imported scripts under `code/upstream/` that read a
`nmd_fig_data/` bundle. **None of these files are in either repo.** Full reproducibility
(the project target) requires each to be obtainable. Grouped by provenance + likely
deposit route.

## A. Yul's derived intermediates — the `nmd_fig_data/` bundle (need from Yul)
The DGEList / count objects her figure + analysis scripts read. These are her
**outputs**, not raw data, so GEO won't hold them — they need a Zenodo *data* record or
to ship with the code (if small enough).

| File | Read by | Notes |
|---|---|---|
| `dge_gene_unfiltered_2026.1.2.rds` | correlation / Fig 1A | SR gene DGEList |
| `dge_shortread_gene_2026.3.2.rds` + `_filtered_` | DGE, productive | SR gene |
| `dge_isoform_longread_2026.3.3.rds` + `_filtered_` | most figures | LR isoform DGEList (core) |
| `salmon.merged.gene_counts_length_scaled.rds` + `.FULL.rds` | Fig 1A | nf-core salmon output |
| `tan_tx_mashr_model.rds` | Tan reanalysis | fitted mashr object |
| `gmap_ENSGv115_2025.08.12.rds`, `gmap_txlevel_ENSGv115_...rds` | annotation | gene/tx maps (ENSEMBL v115) |

## B. Randell isoform-discovery pipeline outputs (`randell:` — separate deposit)
SQANTI3 / isocall products on Channing. Route: study data deposit + Randell attribution,
not Yul's Zenodo.
- `nmd_lungcells_corrected.cds.gff3`, `nmd_lungcells_corrected.fasta`,
  `nmd_lungcells_filtered.gtf`, `nmd_isocall.isoforms.gtf`

## C. Example-gene annotation (small; may ship with code)
- `srsf.gtf` / `srsf.cds.gff3` / `srsf.fasta`, `sr_simple_cds.gtf`,
  `sr_simple_subset.gtf`, `sub.isoforms.gtf`

## D. Already in Pete's repo (no action)
The per-CT + shared mashr CSVs (`shortread_dge/mashr/`, `isocall_dge/mashr/`) — the
scripts that need these read the same `2026.3.10` CSVs already tracked here.

---
**Open question for Yul (O-7):** what is the intended deposit for group A? A Zenodo data
record, or small enough to ship in-repo? And do the `dge_*.rds` de-identify cleanly
(no donor PHI beyond the `DD###X` sample codes already discussed)?
