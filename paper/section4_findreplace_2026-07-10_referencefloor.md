# Manuscript find/replace pairs — 25% reference-share floor (2026-07-10)

**Apply to the Google Doc** (source of truth): https://docs.google.com/document/d/1Tz6coXnDwpGZaV1Jl11fmN_LVX1R8YQPf2y_7wzBLKs/edit
FIND strings taken from the `NMD manuscript 2026.2.5.pdf` snapshot — **verify each against the live Doc** before applying (the Doc may have drifted).
Provenance for every number: `isopair_wrapper/REFERENCE_FLOOR_NUMBERS_DELTA.md` + `REFERENCE_FLOOR_VERIFICATION.md`. All independently verified.

## §4 body prose

1. FIND: `This resulted in 3,009 genes from which the isoform`
   REPL: `This resulted in 1,585 genes from which the isoform`

2. FIND: `31% of its parent gene expression, with 38% of the reference isoforms accounting for >50%`
   REPL: `67% of its parent gene expression, with 71% of the reference isoforms accounting for >50%`
   (This is the headline correction — the pre-floor 31% was the displacement artifact; post-floor references are genuinely dominant.)

3. FIND: `median length 2893 and 3049`
   REPL: `median length 2958 and 2991`
   FIND: `slightly shorter (median 2762 nucleotides`
   REPL: `slightly shorter (median 2808 nucleotides`
   (Reference-vs-NMD is now n.s. (p=0.058) but both remain "similar"; Control still shorter than both, p<0.001. No wording change needed.)

4. FIND: `yielding 190 NMD and Control pairs in which every`
   REPL: `yielding 136 NMD and Control pairs in which every`

5. FIND: `we observed an 18-fold enrichment of PTCs in NMD susceptible isoforms (37.9% vs 2.1% PTC rate, p < 10⁻²⁰)`
   REPL: `we observed a 54-fold enrichment of PTCs in NMD susceptible isoforms (39.7% vs 0.7% PTC rate, p < 10⁻¹⁷)`

6. ⚠ **PROSE JUDGMENT (Pete decide):** FIND: `peaked at 4-5 downstream EJCs from the PTC (SF31`
   REPL (proposed): `was broadly similar across 1–7+ downstream EJCs (SF31`
   (Post-floor SF31 no longer shows a 4–5 peak — the median trend is flat/noisy, nominal max at 6. The "peaked at 4-5" claim no longer holds. Reword or drop.)

7. FIND: `frameshift events were the leading mechanism (55%), followed by in-frame stop events (33%) and 3' UTR splicing (12%`
   REPL: `frameshift events were the leading mechanism (57%), followed by in-frame stop events (32%) and 3' UTR splicing (11%`

8. FIND: `Skipped exon was the most common event type (44% of PTC-causing events vs 14% in Controls, Fisher p < 10⁻⁷), with A5SS also significantly enriched (13% vs 4.5%, p < 10⁻²)`
   REPL: `Skipped exon was the most common event type (47% of PTC-causing events vs 15% in Controls, Fisher p < 10⁻⁶), with A5SS also significantly enriched (15% vs 4%, p < 10⁻²)`

9. FIND: `yielded 1166 NMD and Control pairs. In 50% of these cases`
   REPL: `yielded 888 NMD and Control pairs. In 52% of these cases`

10. FIND: `much stronger in the 492 novel NMD+ isoforms`
    REPL: `much stronger in the 380 novel NMD+ isoforms`

11. FIND: `99% (487/492) of the TD2-called CDS were`
    REPL: `99% (375/380) of the TD2-called CDS were`

12. FIND: `reference AUG was stronger in 78% of cases (384/492, SF34)`
    REPL: `reference AUG was stronger in 82% of cases (310/380, SF34)`

13. FIND: `we identified PTCs in 90% of the NMD+ isoforms (1050/1166)`
    REPL: `we identified PTCs in 92% of the NMD+ isoforms (818/888)`

14. FIND: `both the GENCODE restricted isoform set (n=190) and the more expansive pair set`
    REPL: `both the GENCODE restricted isoform set (n=136) and the more expansive pair set`

15. FIND: `median length 1290/800/948 for NMD+/PTC+, NMD+/PTC-, and Control`
    REPL: `median length 1311/834/948 for NMD+/PTC+, NMD+/PTC-, and Control`

## Figure 3 legend

16. FIND: `all three isoforms are GENCODE annotated (n=190).`
    REPL: `all three isoforms are GENCODE annotated (n=136).`

17. FIND: `SE is twice as prevalent in NMD pairs (44.2 vs 21.2%, p ≈ 10⁻⁸¹)`
    REPL: `SE is twice as prevalent in NMD pairs (51.8 vs 22.7%, p ≈ 10⁻⁶⁵)`

18. FIND: `Direct PTC rate 37.9% NMD (72/190) vs 2.1% Control (4/190)`
    REPL: `Direct PTC rate 39.7% NMD (54/136) vs 0.7% Control (1/136)`

19. FIND: `n = 69 attributed events from 72 PTC+ pairs) vs all events in 190 Controls (light blue; n = 447 total events)`
    REPL: `n = 53 attributed events from 54 PTC+ pairs) vs all events in 136 Controls (light blue; n = 297 total events)`

20. FIND: `SE accounts for 43.5% of PTC-causing events vs 14.1% of Control events (Fisher p = 8×10⁻⁸); A5SS is enriched (13.0% vs 4.5%, p = 9×10⁻³); Alt TES is depleted (5.8% vs 26.8%, p = 3×10⁻⁵)`
    REPL: `SE accounts for 47.2% of PTC-causing events vs 14.8% of Control events (Fisher p = 6×10⁻⁷); A5SS is enriched (15.1% vs 3.7%, p = 3×10⁻³); Alt TES is depleted (3.8% vs 26.6%, p = 7×10⁻⁵)`

21. FIND: `69 attributed pairs split 38 frameshift (coral) / 23 in-frame stop (blue) / 8 3′UTR splice (teal) = 55% / 33% / 12%`
    REPL: `53 attributed pairs split 30 frameshift (coral) / 17 in-frame stop (blue) / 6 3′UTR splice (teal) = 57% / 32% / 11%`

## Figure 4 legend

22. FIND: `all in GENCODE (n = 72 NMD+/PTC+, 118 NMD+/PTC−, 190 Control)`
    REPL: `all in GENCODE (n = 54 NMD+/PTC+, 82 NMD+/PTC−, 136 Control)`

23. FIND: `Of n = 1,166 NMD comparator pairs in which the reference AUG`
    REPL: `Of n = 888 NMD comparator pairs in which the reference AUG`

24. FIND: `1,050 (90%) had a`
    REPL: `818 (92%) had a`

25. FIND: `116 had no downstream EJC (NMD+/PTC−); the 1,166 ref-AUG-traceable Control pairs`
    REPL: `70 had no downstream EJC (NMD+/PTC−); the 888 ref-AUG-traceable Control pairs`

## §5 — conditional (NOT found in the 2026.2.5 PDF; check the live Doc)

26. ⚠ IF the manuscript contains a sentence stating the START/ATG window is "roughly three times more important" in NMD+/PTC− than NMD+/PTC+ isoforms (per SF39):
    FIND (approx): `roughly three times more important in the NMD+/PTC-`
    REPL: `roughly twice as important in the NMD+/PTC-`
    (SF39 ATG-branch ratio moved 2.21× → 1.98× post-floor.)
    NB: the §5 "STOP site sequence was nearly three times as important" sentence is a DIFFERENT claim (STOP vs ATG branch importance, model-global) — **do NOT change it**; the model was not retrained.

## Methods — "Isoform Pairs Analysis" (documents the floor; the Step-5 gap)

M-1. **Document the floor** — FIND: `This design ensures that when we compare NMD-related splicing transitions to Control splicing transitions, these comparisons share an identical reference isoform.`
   REPL: `This design ensures that when we compare NMD-related splicing transitions to Control splicing transitions, these comparisons share an identical reference isoform. To ensure the reference represents the gene's dominant transcript rather than a minor non-NMD isoform, we retained only genes in which the selected reference isoform accounted for at least 25% of the gene's total isoform expression in DMSO (mean across all sequenced libraries); genes whose highest-expressed strictly-non-NMD isoform fell below this threshold were excluded, yielding 1,585 gene-matched triplets.`

M-2. FIND: `all GENCODE-annotated coding transcripts (n = 190)`
    REPL: `all GENCODE-annotated coding transcripts (n = 136)`

M-3. FIND: `so the NMD and Control arms share the same gene set (n = 1,166)`
    REPL: `so the NMD and Control arms share the same gene set (n = 888)`

M-4. FIND: `the reference-anchored analysis reveals (n = 492)`
    REPL: `the reference-anchored analysis reveals (n = 380)`

M-5. The Isopair **vignette** does NOT need the floor — it documents the generic `generatePairsExpression` package function; the ≥25% floor is a project-specific post-selection filter in the wrapper `02_build_profiles_mashr.R` (already documented in code comments + `results_to_code_map.md` M11). The manuscript Methods (M-1) is the reader-facing home for it.

**Pre-existing Methods discrepancy (NOT floor-related; Pete's call):** the Methods say non-NMD = "adj.P.Val > 0.50", but the code uses > 0.30 (per ONBOARDING §6 / map claim 4.3). Flagging in passing; out of scope for this floor pass.

## Not changed (verified floor-independent)
Abstract (no §4 numbers), §1–§3, §5 model performance (AUC 0.93 / AUPRC 0.83), branch-importance percentages, all model-global SHAP/attention/GC figures (SF36/37/38/41), SF42 (frozen June comparison).
