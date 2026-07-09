# Supplemental Figure — Stop Codon Usage: NMD vs Control

Companion to the §5 manuscript sentence:

> "We confirmed that, when comparing NMD to non-NMD isoforms, UGA is
> more commonly used (SFx)."

## What it shows

Within-class stop-codon frequencies on the deep-learning model's chr-1/3/5/7
paralog-free test set (rank-0 ORFs, 5-run weighted): **UGA 52.0% NMD vs
48.3% Control** (enriched), UAA 27.1% vs 30.3% (depleted), UAG 20.9% vs
21.3% (unchanged).

## Status

**DRAFT placeholder.** The rendered image is extracted from base64 in
`results/.../deep_nmd_model/orf_model_report_v5.html` — the v5 Rmd's
`fig4a_stop_codon` chunk, post-bug-fix (`scripts/patch_stop_codon.py`,
2026-04-30, fixed off-by-one in the ORFik stop-codon column).

The legacy rendered `fig3c_stop_identity.png` in the same directory is
**STALE** — generated from pre-bug-fix data, shows two codons (TAA + TGA)
with frequencies that contradict the post-fix numbers. Do not use it.

**Subtitle removed.** The legacy render's ggplot subtitle (`Chi-sq = …, p = …
… pairwise Fisher exact tests …`) was wider than the panel and clipped at
the right edge. We applied two fixes:

1. **Rmd source (`NMD_orf_model_v5_4ct/orf_model_report_v5.Rmd`)** —
   dropped the `subtitle = …` argument from the `fig4a-stop-codon-freq`
   chunk's `labs()` call so future re-renders (when the cluster is back)
   no longer produce that text.
2. **Local placeholder** — `figure_s_stop_codon_usage.py` crops the top
   ~8.5% of the embedded legacy PNG to remove the subtitle band on the
   draft figure.

When the cluster comes back, the canonical re-render is the
`fig4a-stop-codon-freq` chunk in `orf_model_report_v5.Rmd` (lines
~1108–1175). Output: `fig4a_stop_codon.{pdf,png}` at 7 × 6 in. The
`SUBTITLE_CROP_FRAC` constant in this folder's `.py` can be set to 0
once we re-render from the updated Rmd source — but the current 8.5%
crop is consistent with the corrected source.

## Files

| File | What |
|---|---|
| `figure_s_stop_codon_usage.py` | Embeds the post-fix legacy PNG as the figure body (matches the Panel E / F placeholder pattern in `figure5_dl_model/`) |
| `figure_s_stop_codon_usage_legend.md` | Manuscript-style legend |
| `data/legacy_fig4a_stop_codon.png` | Extracted base64 from the v5 HTML (1344 × 1152) |

## Regenerating

```bash
python3 figure_s_stop_codon_usage.py
```

## Cross-references

- Figure 5 Panel F (signed SHAP × input at the stop codon) — same
  directional finding: G at the variable stop position carries positive
  signed attribution.
- v5 Rmd Section 4.1 narrative: `orf_model_report_v5.Rmd` lines
  ~1099–1205 covers the frequency + SHAP analyses end-to-end.
- Bug-fix audit: `BUGFIX_STOP_CODON_2026-03-31.md` (in the model repo).
