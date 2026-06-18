# Supplemental Figure — Attention Distribution: NMD vs Control

Companion to the §5 manuscript sentence:

> "Attention analysis demonstrated that, while most of the attention was
> focused on ORF0 (the most likely CDS), attention was more broadly
> distributed for NMD than Control isoforms (SFx, panel A and B)."

## What it shows

- **Panel A** (`fig5b_attention_by_rank`) — boxplots of attention weight
  per ORF rank, NMD vs Control. Rank-0 dominance + visible NMD tail at
  higher ranks.
- **Panel B** (`fig5a_entropy`) — density of per-isoform Shannon entropy
  over the 5-ORF attention vector. NMD shifted right of Control = more
  distributed attention.

## Status

**DRAFT.** Both panels are embedded as PNGs copied from the legacy
`deep_nmd_model/figures/` directory (`fig5a_entropy.png`,
`fig5b_attention_by_rank.png`). The underlying attention exports
(`04_interpret_attention.py` → `uorf_attention_predictions.tsv`) don't
depend on the stop-codon bug fix and are therefore canonical as-is.
When the cluster is back, the report's §5 attention chunks can be
re-rendered natively against the latest model outputs without changing
the directional finding.

## Files

| File | What |
|---|---|
| `figure_s_attention_distribution.py` | 2-panel composite, side-by-side embed |
| `figure_s_attention_distribution_legend.md` | Manuscript-style legend |
| `data/legacy_fig5b_attention_by_rank.png` | Panel A source (2400×1500) |
| `data/legacy_fig5a_entropy.png` | Panel B source (2400×1500) |

## Regenerating

```bash
python3 figure_s_attention_distribution.py
```

## Cross-references

- Figure 5 Panel G (uORF attention at n=1,166) — a downstream
  interpretation of the same attention vectors, focused on the
  uORF-mediated NMD mechanism for PTC− retained isoforms.
- Legacy v5 Rmd §5 (attention analysis): `orf_model_report_v5.Rmd`
  Sections 5 / 9.9 — entropy quantification, rank distribution, and
  per-CDS-source attention comparisons.
