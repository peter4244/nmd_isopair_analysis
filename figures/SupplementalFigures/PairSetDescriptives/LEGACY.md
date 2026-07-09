# LEGACY — superseded 2026-07-08

This directory holds the original 3-panel PairSetDescriptives supplemental
figure (Panels A/B/C = isoforms-per-gene / reference-share / transcript-length).

The manuscript's SF numbering (2026-07-08) splits these three panels into
three independent supplemental figures:

| Original panel | Now |
|---|---|
| Panel A (isoforms per gene, n = 3,009) | `SF25_IsoformsPerGene/` |
| Panel B (reference share of gene expression) | `SF26_ReferenceShare/` |
| Panel C (transcript length by pair role) | `SF27_TranscriptLengthByRole/` |

Each of those three dirs is a standalone ggplot-mimic matplotlib panel with
its own legend and data — they are the canonical source going forward.

**This directory is retained for provenance only.** The `data_export.R` script
here is still the canonical source of the TSVs the SF25/26/27 scripts copy from
(`isoforms_per_gene.tsv`, `ref_expression_fraction.tsv`, `tx_length_by_role_long.tsv`,
`descriptives_summary.tsv`). To regenerate any of those TSVs, run
`data_export.R` here and copy the outputs into the three SF dirs.

Do **not** use `figure_s_pairset_descriptives.py` for the paper — the split
SF25/SF26/SF27 scripts are the shipping renders.
