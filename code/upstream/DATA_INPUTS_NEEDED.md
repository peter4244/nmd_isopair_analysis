# Data inputs needed to run Yul's upstream scripts (D1 inventory)

> **⚠ SCOPE CHANGED 2026-07-20 (plan D-5).** Pete set the standard: **raw reads (GEO
> GSE329233) + complete code chain; interim data files are NOT deposited.** Therefore
> **groups A, B, and C below are NO LONGER DELIVERABLES** — they are interim products.
> This file is retained as a *runtime* inventory (what you need on hand to actually
> execute these scripts), not as a deposit checklist. Group B's "no identified home"
> problem is moot — do not resurrect it.
>
> **The obligation this creates instead:** every step that *produces* these files must
> have its code in the repo. The open gaps are **M2** (nf-core launcher) and **M3**
> (isocall / SQANTI3 config) — see the plan.

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

## B. Isoform-discovery pipeline outputs (`randell:`) — **NO IDENTIFIED HOME**
SQANTI3 / isocall products under the Randell lab's Channing tree. Not Yul's outputs, so
not hers to supply.
- `nmd_lungcells_corrected.cds.gff3`, `nmd_lungcells_corrected.fasta`,
  `nmd_lungcells_filtered.gtf`, `nmd_isocall.isoforms.gtf`

**Open — do not record this as settled.** The manuscript's Data Availability lists only
**GEO GSE329233** + the Isopair GitHub link. GSE329233 is an expression record and will
not hold GTF/GFF3/FASTA pipeline products, so these files currently have no deposit
target — the same gap as the group-A bundle. Resolve with Scott Randell / Yul; do not
assume a separate "Randell deposit" exists.

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
