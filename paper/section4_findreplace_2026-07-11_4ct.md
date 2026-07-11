# Manuscript find/replace pairs — 4-CT re-scope (2026-07-11)

**Apply to the Google Doc** (source of truth): https://docs.google.com/document/d/1Tz6coXnDwpGZaV1Jl11fmN_LVX1R8YQPf2y_7wzBLKs/edit
FIND strings taken from the `NMD manuscript 2026.2.5.pdf` snapshot — **verify each against the live Doc** before applying (the Doc may have drifted).
This supersedes and replaces the never-applied 6-CT `section4_findreplace_2026-07-10_referencefloor.md`. The whole `all_samples` analysis is now the four manuscript cell types (AT, DD, FB, MV) at the source — there is no 6-CT / "all sequenced libraries" state anywhere.
Provenance for every number: `isopair_wrapper/REFERENCE_FLOOR_NUMBERS_DELTA_4CT.md` + the (all-passing) verifiers `reproducibility/verify_pass7_new_rmd.R` (37/37) and `verify_cross_check_new_rmd_vs_figures.R` (57/57). All independently re-derived.

## §4 body prose

1. FIND: `This resulted in 3,009 genes from which the isoform`
   REPL: `This resulted in 1,548 genes from which the isoform`

2. FIND: `31% of its parent gene expression, with 38% of the reference isoforms accounting for >50%`
   REPL: `64% of its parent gene expression, with 68% of the reference isoforms accounting for >50%`
   (Headline correction — the pre-floor 31% was the displacement artifact; the 4-CT floored references are genuinely dominant. SF26: median 64.3%, ≥50% share 67.5%.)

3. FIND: `median length 2893 and 3049`
   REPL: `median length 2896 and 3018`
   FIND: `slightly shorter (median 2762 nucleotides`
   REPL: `slightly shorter (median 2742 nucleotides`
   (Reference-vs-NMD remains "similar"; Control still shorter than both. No wording change needed.)

4. FIND: `yielding 190 NMD and Control pairs in which every`
   REPL: `yielding 130 NMD and Control pairs in which every`

5. FIND: `we observed an 18-fold enrichment of PTCs in NMD susceptible isoforms (37.9% vs 2.1% PTC rate, p < 10⁻²⁰)`
   REPL: `we observed a 24-fold enrichment of PTCs in NMD susceptible isoforms (36.9% vs 1.5% PTC rate, p < 10⁻¹³)`
   (Fisher OR = 37.06, p = 1.58×10⁻¹⁴.)

6. ⚠ **PROSE JUDGMENT (Pete decide):** FIND: `peaked at 4-5 downstream EJCs from the PTC (SF31`
   REPL (proposed): `was broadly similar across 1–7+ downstream EJCs (SF31`
   (4-CT SF31 shows no 4–5 peak — the per-bin median trend is flat/noisy, nominal max at bin 6. The "peaked at 4-5" claim no longer holds. Reword or drop.)

7. FIND: `frameshift events were the leading mechanism (55%), followed by in-frame stop events (33%) and 3' UTR splicing (12%`
   REPL: `frameshift events were the leading mechanism (52%), followed by in-frame stop events (38%) and 3' UTR splicing (10%`
   (48 attributed events: 25 frameshift / 18 in-frame stop / 5 3′UTR splice.)

8. FIND: `Skipped exon was the most common event type (44% of PTC-causing events vs 14% in Controls, Fisher p < 10⁻⁷), with A5SS also significantly enriched (13% vs 4.5%, p < 10⁻²)`
   REPL: `Skipped exon was the most common event type (48% of PTC-causing events vs 14% in Controls, Fisher p < 10⁻⁶), with A5SS also significantly enriched (13% vs 4%, p < 0.05)`
   (SE 47.9% vs 14.4%, p = 6.9×10⁻⁷; A5SS 12.5% vs 4.1%, p = 2.8×10⁻².)

9. FIND: `yielded 1166 NMD and Control pairs. In 50% of these cases`
   REPL: `yielded 819 NMD and Control pairs. In 52% of these cases`
   (TD2 == reference AUG in 429/819 = 52.4%.)

10. FIND: `much stronger in the 492 novel NMD+ isoforms`
    REPL: `much stronger in the 348 novel NMD+ isoforms`

11. FIND: `99% (487/492) of the TD2-called CDS were`
    REPL: `99% (345/348) of the TD2-called CDS were`

12. FIND: `reference AUG was stronger in 78% of cases (384/492, SF34)`
    REPL: `reference AUG was stronger in 82% of cases (286/348, SF34)`

13. FIND: `we identified PTCs in 90% of the NMD+ isoforms (1050/1166)`
    REPL: `we identified PTCs in 92% of the NMD+ isoforms (756/819)`

14. FIND: `both the GENCODE restricted isoform set (n=190) and the more expansive pair set`
    REPL: `both the GENCODE restricted isoform set (n=130) and the more expansive pair set`

15. FIND: `median length 1290/800/948 for NMD+/PTC+, NMD+/PTC-, and Control`
    REPL: `median length 1324/722/922 for NMD+/PTC+, NMD+/PTC-, and Control`
    (3′UTR length via first non-premature stop, broad ref-AUG scope; SF35 Panel C.)

## Figure 3 legend

16. FIND: `all three isoforms are GENCODE annotated (n=190).`
    REPL: `all three isoforms are GENCODE annotated (n=130).`

17. FIND: `SE is twice as prevalent in NMD pairs (44.2 vs 21.2%, p ≈ 10⁻⁸¹)`
    REPL: `SE is more than twice as prevalent in NMD pairs (51.0 vs 21.6%, p ≈ 10⁻⁶⁵)`

18. FIND: `Direct PTC rate 37.9% NMD (72/190) vs 2.1% Control (4/190)`
    REPL: `Direct PTC rate 36.9% NMD (48/130) vs 1.5% Control (2/130)`

19. FIND: `n = 69 attributed events from 72 PTC+ pairs) vs all events in 190 Controls (light blue; n = 447 total events)`
    REPL: `n = 48 attributed events from 48 PTC+ pairs) vs all events in 130 Controls (light blue; n = 292 total events)`

20. FIND: `SE accounts for 43.5% of PTC-causing events vs 14.1% of Control events (Fisher p = 8×10⁻⁸); A5SS is enriched (13.0% vs 4.5%, p = 9×10⁻³); Alt TES is depleted (5.8% vs 26.8%, p = 3×10⁻⁵)`
    REPL: `SE accounts for 47.9% of PTC-causing events vs 14.4% of Control events (Fisher p = 7×10⁻⁷); A5SS is enriched (12.5% vs 4.1%, p = 3×10⁻²); Alt TES is depleted (4.2% vs 26.0%, p = 3×10⁻⁴)`

21. FIND: `69 attributed pairs split 38 frameshift (coral) / 23 in-frame stop (blue) / 8 3′UTR splice (teal) = 55% / 33% / 12%`
    REPL: `48 attributed pairs split 25 frameshift (coral) / 18 in-frame stop (blue) / 5 3′UTR splice (teal) = 52% / 38% / 10%`

## Figure 4 legend

22. FIND: `all in GENCODE (n = 72 NMD+/PTC+, 118 NMD+/PTC−, 190 Control)`
    REPL: `all in GENCODE (n = 48 NMD+/PTC+, 82 NMD+/PTC−, 130 Control)`

23. FIND: `Of n = 1,166 NMD comparator pairs in which the reference AUG`
    REPL: `Of n = 819 NMD comparator pairs in which the reference AUG`

24. FIND: `1,050 (90%) had a`
    REPL: `756 (92%) had a`

25. FIND: `116 had no downstream EJC (NMD+/PTC−); the 1,166 ref-AUG-traceable Control pairs`
    REPL: `63 had no downstream EJC (NMD+/PTC−); the 819 ref-AUG-traceable Control pairs`

## §5 — conditional (NOT found in the 2026.2.5 PDF; check the live Doc)

26. ⚠ IF the manuscript contains a sentence stating the START/ATG window is "roughly three times more important" in NMD+/PTC− than NMD+/PTC+ isoforms (per SF39):
    FIND (approx): `roughly three times more important in the NMD+/PTC-`
    REPL: `roughly twice as important in the NMD+/PTC-`
    (SF39 ATG-branch share ratio PTC−:PTC+ = 1.87× under 4-CT; still "roughly twice".)
    NB: the §5 "STOP site sequence was nearly three times as important" sentence is a DIFFERENT claim (STOP vs ATG branch importance, model-global) — **do NOT change it**; the model was not retrained.

27. **SF42 (predictor comparison) is no longer frozen** — it was regenerated on the 4-CT cohort (test set n = 415: 195 NMD+/PTC+, 16 NMD+/PTC−, 204 Control; test Spearman NMDetective-B 0.79, NMDEP 0.79, our model 0.77). IF §4/§5 prose cites any of these predictor-comparison numbers, update them to match; the 2026.2.5 PDF does not appear to cite them, so likely no change is needed. Verify against the live Doc.

## Methods — "Isoform Pairs Analysis" (documents the floor + the 4-CT cell-type scope)

M-1. **Document the floor + 4-CT scope** — FIND: `This design ensures that when we compare NMD-related splicing transitions to Control splicing transitions, these comparisons share an identical reference isoform.`
   REPL: `This design ensures that when we compare NMD-related splicing transitions to Control splicing transitions, these comparisons share an identical reference isoform. To ensure the reference represents the gene's dominant transcript rather than a minor non-NMD isoform, we retained only genes in which the selected reference isoform accounted for at least 25% of the gene's total isoform expression in DMSO (mean across the four cell types); genes whose highest-expressed strictly-non-NMD isoform fell below this threshold were excluded, yielding 1,548 gene-matched triplets.`
   (M-1 wording: "all sequenced libraries" → "the four cell types" — the isoform universe, normalization, reference selection, floor, and pair construction are all AT/DD/FB/MV.)

M-2. FIND: `all GENCODE-annotated coding transcripts (n = 190)`
    REPL: `all GENCODE-annotated coding transcripts (n = 130)`

M-3. FIND: `so the NMD and Control arms share the same gene set (n = 1,166)`
    REPL: `so the NMD and Control arms share the same gene set (n = 819)`

M-4. FIND: `the reference-anchored analysis reveals (n = 492)`
    REPL: `the reference-anchored analysis reveals (n = 348)`

M-5. The Isopair **vignette** does NOT need the floor — it documents the generic `generatePairsExpression` package function; the ≥25% floor is a project-specific post-selection filter in the wrapper `02_build_profiles_mashr.R`. The manuscript Methods (M-1) is the reader-facing home for it.

**Pre-existing Methods discrepancy (NOT re-scope-related; Pete's call):** the Methods say non-NMD = "adj.P.Val > 0.50", but the code uses > 0.30 (per ONBOARDING §6 / map claim 4.3). Flagging in passing; out of scope for this pass.

## Not changed (verified re-scope-independent)
Abstract (no §4 numbers), §1–§3, §5 model performance (AUC / AUPRC — model not retrained), branch-importance percentages, all model-global SHAP/attention/GC figures (SF36/37/38/41).
