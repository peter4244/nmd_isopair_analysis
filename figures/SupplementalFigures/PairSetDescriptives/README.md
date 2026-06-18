# Supplemental Figure — Pair-set descriptives (pop_BC)

Companion to §4¶1 + §4¶2 — descriptive characterisation of the 3,009-gene Isopair pair set (pop_BC).

## What it shows

Three panels:
- **(A)** Isoforms per gene, across 3,009 pair-genes. **Median 7** — verified match to §4¶1.
- **(B)** Reference-isoform share of parent-gene non-NMD expression. Median **30.8%**, with 37.6% of references ≥ 50%.
- **(C)** Transcript-length distribution per pair role (NMD comparator / Reference / Control comparator). Medians **3,049 / 2,893 / 2,762 nt** — verified exact match to §4¶2.

## Numbers — verified match to the revised §4¶1 (2026-06-18)

| Claim | Manuscript (revised) | Computed | Verdict |
|---|---|---|---|
| Pair genes | 3,009 | 3,009 | ✓ |
| Median isoforms/gene | 7 | 7 | ✓ |
| Ref. share — median | **31%** | 30.8% | ✓ |
| Ref. share — % above 50% | **38%** | 37.6% | ✓ |
| NMD comparator median length | 3,049 nt | 3,049 nt | ✓ |
| Reference median length | 2,893 nt | 2,893 nt | ✓ |
| Control comparator median length | 2,762 nt | 2,762 nt | ✓ |
| Length test | "p < 0.001" | KW p = 1.5×10⁻⁷; all pairwise p ≤ 0.034 | ✓ direction |

The §4¶1 sentence as confirmed by Pete on 2026-06-18 reads:

> "These genes had a median of 7 isoforms each (SFx - Isoform Count in Isoform Pair Sets), and the median reference isoform accounted for 31% of its parent gene expression, with 38% of the reference isoforms accounting for >50% of parent gene expression."

This SF (Panels A + B) is the canonical source for those two numbers. The "SFx" placeholder in the §4¶1 sentence resolves to this Supplemental Figure (Panel A specifically). Panel B documents the second descriptive (reference share); Panel C provides the §4¶2 transcript-length comparison.

## Files

| File | What |
|---|---|
| `data_export.R` | Builds Panels A/B/C TSVs from local Isopair pipeline objects (`profiles_c2/c4_allsamples.rds`, `structures.rds`, `expression_data.rds`, `nmd_classification.rds`, `dmso_samples.rds`); also writes `descriptives_summary.tsv` |
| `figure_s_pairset_descriptives.py` | 3-panel matplotlib figure (validator-clean) |
| `figure_s_pairset_descriptives_legend.md` | Manuscript-style legend |
| `data/isoforms_per_gene.tsv` | Per-gene isoform count |
| `data/ref_expression_fraction.tsv` | Per-gene reference fraction |
| `data/tx_length_by_role_long.tsv` | Long-form tx-length data (9,027 rows) |
| `data/pairwise_tx_length.tsv` | Wilcoxon contrasts for Panel C |
| `data/descriptives_summary.tsv` | Single-row summary used by Methods §1 |

## Regenerating

```bash
Rscript data_export.R
python3 figure_s_pairset_descriptives.py
```

## Cross-references

- **§4¶1** — primary text cites median 7 isoforms/gene + the 70%/75% reference numbers; verified match for n_pair_genes (3,009) and median 7.
- **§4¶2** — primary text cites the three transcript-length medians; verified exact match.
- **Figure 3 Panel A** — pair construction schematic; uses the same 3,009 numbers.
- **Methods §1 (Isoform pairs analysis)** — `paper/methods_updates_2026-06-18.md` already cites these numbers; once the reference-share discrepancy is resolved, that draft is the canonical home for the §4¶1 descriptive text.
