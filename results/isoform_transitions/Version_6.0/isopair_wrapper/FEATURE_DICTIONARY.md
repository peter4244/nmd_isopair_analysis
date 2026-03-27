# Feature Dictionary — NMD Prediction Model

## Canonical Feature Definitions

Each feature is computed from two CDS sources: **TD2** (TransDecoder2-predicted)
and **Reference** (dominant non-NMD isoform's CDS traced through the target).
The TD2 version uses the R variable name; the Reference version uses the
`ref_` prefix.

| R variable | Display name | Category | Type | Definition |
|---|---|---|---|---|
| `downstream_ejc` | Downstream EJC count | Stop codon position | Continuous (0–5) | Number of exon-exon junctions downstream of the stop codon. Truncated at 5. Source: `ptc.rds` (TD2) or `ref_cds_features_all.rds` (ref). |
| `atg_density` | 5'UTR ATG density | Start codon context | Continuous | ATG codons per 100 bp of 5'UTR. Source: `Isopair::scan5UtrFeatures()`. |
| `atg_count` | 5'UTR ATG count | Start codon context | Count | Total ATG codons in the 5'UTR. |
| `atg_strong_kozak` | Strong Kozak ATG count | Start codon context | Count | 5'UTR ATGs matching strong Kozak consensus (favorable nucleotides at positions -3 and +4). |
| `atg_orphan_count` | Orphan ATG count | Start codon context | Count | 5'UTR ATGs not associated with any predicted ORF (no downstream in-frame stop codon within the 5'UTR). |
| `uorf_count_overlapping` | Overlapping uORF count | Upstream ORF | Count | uORFs whose reading frame extends past the annotated CDS start codon. |
| `uorf_longest_nt` | Longest uORF (nt) | Upstream ORF | Continuous | Length in nucleotides of the longest uORF in the 5'UTR. |
| `uorf_count_inframe` | In-frame uORF count | Upstream ORF | Count | uORFs in the same reading frame as the main CDS. |
| `uorf_count_outframe` | Out-of-frame uORF count | Upstream ORF | Count | uORFs in a different reading frame from the main CDS. |
| `utr5_orf_coverage` | 5'UTR ORF coverage | Start codon context | Continuous (0–100%) | Percentage of the 5'UTR spanned by predicted ORFs (union of all uORF intervals). |
| `stop_density` | 5'UTR stop density | Start codon context | Continuous | Stop codons (all three frames) per 100 bp of 5'UTR. |
| `log_utr3_length` | 3'UTR length (log bp) | 3'UTR | Continuous | Log-transformed (log1p) exonic basepairs from the stop codon to the 3' transcript end. |

## Reference-CDS Prefix Convention

For the reference CDS version, all variable names are prefixed with `ref_`:
`ref_downstream_ejc`, `ref_atg_density`, `ref_atg_count`, etc.

## Source Scripts

| Script | Features produced | Output file |
|---|---|---|
| `05l_unified_model.R` | TD2 features (12) | `unified_model.rds` |
| `05t_ref_cds_features.R` | Reference CDS features (12) | `ref_cds_features_all.rds` |
| `05v_model_comparison.R` | Combined model (24 features) | `model_comparison.rds` |
| `05k_utr5_all_isoforms.R` | Raw 5'UTR scan (all isoforms) | `utr5_features_all.rds` |

## Feature Label Mapping (R code)

```r
feature_labels <- c(
  downstream_ejc           = "Downstream EJC count",
  uorf_count_overlapping   = "Overlapping uORF count",
  uorf_longest_nt          = "Longest uORF (nt)",
  uorf_count_inframe       = "In-frame uORF count",
  uorf_count_outframe      = "Out-of-frame uORF count",
  atg_density              = "5'UTR ATG density (per 100 bp)",
  atg_count                = "5'UTR ATG count",
  atg_strong_kozak         = "Strong Kozak ATG count",
  atg_orphan_count         = "Orphan ATG count",
  utr5_orf_coverage        = "5'UTR ORF coverage (%)",
  stop_density             = "5'UTR stop density (per 100 bp)",
  log_utr3_length          = "3'UTR length (log bp)"
)
# Reference CDS versions: paste0("ref_", names(feature_labels))
```
