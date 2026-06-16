# Manuscript audit — current Google Doc vs. canonical figures

**Source of truth:** Google Doc at `https://docs.google.com/document/d/1Tz6coXnDwpGZaV1Jl11fmN_LVX1R8YQPf2y_7wzBLKs/edit` (fetched 2026-06-15).
**Sections audited:** §3 (Transcriptional Output Lost and Productive Compensation Index) and §4 (Attributing NMD Susceptibility to Specific Splicing Events). §5 (deep-learning model) out of scope.
**Canonical figure scopes (current):**
- Figure 3 D/E/F + Figure 4 Panels A/B: all-3-ENST + coding-CDS, n = 190/190, own-GENCODE-stop classifier, headline 37.9%/2.1%/18×.
- Figure 4 Panels C/D: ENST-reference + ref-AUG-traceable, 1:1 gene-matched re-intersected, n = 1,050 PTC+ / 113 PTC− / 1,166 Control.
- TD2BiasEvidence supplement: n = 1,166 broad + n = 492 occult-PTC subset; ORF ratio 4.69× / 7.73×; Kozak p = 1.2×10⁻³⁸ / 2.6×10⁻³⁸; position 50% same / 43% downstream broad, 99.0% downstream occult-PTC.
- CDSand3UTR_GENCODEonly supplement: n = 72/118/190; CDS medians 760/834/968; 3′UTR translation-based 869/384/517; 3′UTR non-PTC-stop 524/384/517.

## Bottom line

**§4 is largely up to date.** Most stale-number issues from the previous audit (which had been reading a stale Rmd export, not the Google Doc) have already been fixed by Pete in the live document. The numerical claims, panel references, and supplement references all reconcile with the current canonical figures.

**Three items worth attention:**
1. The "novel NMD+ isoforms" framing in the §4 TD2-bias paragraph is technically not quite right (the supplement scope isn't restricted to novel comparators — see Item 1 below).
2. The 3′UTR-outside-the-gene-matched-set paragraph cites numbers (n=1,904 NMD / 22,335 non-NMD; 605 vs 919 nt 3′UTR medians) that aren't in any current figure or supplement deliverable — needs source verification.
3. Standard placeholder cleanup: "Figure X" → "Figure 3" / "Figure 4"; "SFx" → specific supplement number. (These are manuscript-stage placeholders; flagged for completeness, not as drift.)

§3 (PCI) prose was not in scope for the figure-side updates we did — it's about a separate analytical thread (productive compensation, GPR180 example, pathway enrichment) and does not reference Figures 3/4 numbers.

---

## §3 audit

§3 covers the Productive Compensation Index (PCI) story — % output lost across cell types, PCI definition, GPR180 example, pathway enrichment for PCI genes. **None of the §3 prose references Figures 3 or 4** or the canonical scopes from this audit. The figure references in §3 are to Figure 2A (output lost), Fig X (panel B / panel C / panel D / panel E / panel F) and various SFx supplements that I did not check in this round. No numerical drift identified against current Figures 3/4.

## §4 audit (paragraph by paragraph)

### ¶1 — Isopair construction

Quotes 3,009 genes / 7 isoforms / 70% / 75% / 3067 / 3107 / 2906 nt. **All upstream Isopair-pipeline numbers**, not in the canonical figure reference list. Trust as-is; verify against current Rmd output (`05_final_report_mashr.Rmd` sections covering pair construction) if you want to be fully rigorous. ✓ structurally consistent with the analysis.

### ¶2 — Sequence similarity + splicing event prevalence

Mostly qualitative. References "Fig X, B" and "Fig X, C" for sequence similarity and splice-event prevalence (these correspond to Figure 3 Panels B and C in the current layout). ✓ panel references correct.

### ¶3 — PTC analysis at all-3-ENST + coding-CDS scope

**All numbers check out:**
- "n = 301 NMD and 301 Control pairs" ✓ (all-3-ENST gene-matched)
- "190 NMD pairs and 190 Control pairs" ✓ (after coding-CDS filter + re-intersection)
- "18-fold enrichment of PTCs in NMD-susceptible isoforms (37.9% vs 2.1% PTC rate, p < 10⁻²⁰)" ✓
- Mechanism split: "frameshift (55%), in-frame stop (33%), 3'UTR splicing (12%)" ✓ matches 38/23/8 = 55.1% / 33.3% / 11.6%
- Event-type enrichment: "Skipped exon ... 44% of PTC-causing events vs 14% in Controls, Fisher p < 10⁻⁷" ✓ matches 43.5% / 14.1% / p=8×10⁻⁸
- "A5SS also significantly enriched (13% vs 4.5%, p < 10⁻²)" ✓ matches 13.0% / 4.5% / p=9×10⁻³
- "Fig 3, E and F" ✓ correct panel references

### ¶4 — NMD+/PTC− mechanism (the long paragraph)

**Mostly correct, one minor framing issue.**

Correct claims:
- "compared CDS length, 5' UTR and longest uORF length, and 3' UTR length" ✓ matches Section A + CDSand3UTR supplement
- "NMD+/PTC- isoforms had longer 5' UTRs and longer uORFs than both other groups (Figure X Panel A and B)" ✓ Figure 4 Panels A/B
- "no significant difference in CDS or 3'UTR length between the NMD+/PTC+ and NMD+/PTC- isoforms (SFx)" ✓ CDSand3UTR supplement
- "n=1166" ✓ Section C universe
- "50% ... TD2-called CDS used the same AUG as the reference" ✓ 50% same in supplement broad scope
- "492 novel NMD+ isoforms where the ORF defined by the reference AUG included a PTC" ✓ occult-PTC subset
- "99% (487/492)" ✓
- "Kozak score was stronger for the reference AUG in 78% of cases (384/492, SFx)" ✓
- "1050/1166" ✓ Section C PTC+ rate

**Item 1 — minor framing issue, optional fix:**

> "We then focused on the novel NMD+ isoforms (CDS called from TD2 predictions) where the reference isoform AUG start codon was also present, and we compared the TD2-called CDS to the ORF defined by the reference AUG (n=1166)."

The **n=1,166 Section C universe in TD2BiasEvidence is NOT restricted to novel comparators** — it includes ENST comparators where the comparator has its own GENCODE-annotated CDS. Calling the 1,166 "novel NMD+ isoforms" is technically inaccurate: ~80% of NMD c2 comparators in this scope are novel, but ~20% are ENST.

Three options:
- **(A) Drop "novel"** — recast as "We then focused on the NMD+ isoforms in the broader population (n=1,166) where the reference isoform's AUG was present in the comparator's exonic sequence ..."
- **(B) Restrict the analysis to novel comparators** — would require re-running the supplement at the novel-only scope (drops 1,166 → ~924). Not recommended; loses generality.
- **(C) Leave as-is** — the difference between "novel + ENST" and "novel only" doesn't change the bias direction or any of the quoted numbers materially, and "TD2-called CDS" is the right framing for the comparison that's actually shown.

Recommendation: **(A)**, drop "novel" — the cleanest fix. The supplement's universe is "all NMD comparator isoforms in the Section C scope," not just novel ones.

### ¶5 — 3'UTR aside (GENCODE-only outside-the-gene-matched-set)

**Item 2 — numbers not in any current deliverable, needs source verification:**

> "we analysed only isoforms with GENCODE-annotated CDS (n=1,904 NMD isoforms and 22,335 non-NMD), and we observed that the 3' UTRs were shorter in the NMD susceptible gene set (605 [299-1,119] nt in NMD versus 919 [360-1980] in non-NMD, p<0.001, SFx)."

This paragraph cites numbers from a separate analysis (NMD-vs-non-NMD 3'UTR at the GENCODE-CDS-only scope, outside the gene-matched Isopair set). Neither the n's (1,904 / 22,335) nor the medians (605 / 919) appear in any current Figure 3 / Figure 4 / supplement deliverable. The analysis these numbers came from was last run in the Rmd pipeline (probably `05_final_report_mashr.Rmd` section on "3'UTR length distribution outside the gene-matched set" or similar) and may or may not still be current at the new manuscript scopes.

**Action:** confirm whether this analysis still exists in current form. If it does, the numbers are likely fine (the GENCODE-only scope is methodology-independent of the gene-matched analysis). If the Rmd that produces it has drifted, the numbers need rechecking.

## Placeholder cleanup (manuscript-stage, not drift)

These are typical pre-submission placeholders that need final filling, listed here for completeness:

- "Fig X" → "Figure 3" or "Figure 4" depending on context. Current mapping:
  - "Fig X, A" / "Fig X, B" / "Fig X, C" / "Fig X, D" → **Figure 3** Panels A / B / C / D
  - "Fig 3, E and F" → already correctly numbered, ✓
  - "Figure X Panel A and B" → **Figure 4** Panels A and B
  - "Figure X, Panels C and D" → **Figure 4** Panels C and D
- "SFx" → specific supplement letter:
  - Splice-event prevalence aside ("Gain Direction by Event Type") → existing supplement
  - PTC distance dose response → existing supplement
  - NMD Effect Size by Number of EJCs → existing supplement
  - CDS / 3'UTR no-significant-difference aside → **CDSand3UTR_GENCODEonly** supplement
  - TD2 bias references (4 mentions across ¶4) → **TD2BiasEvidence** supplement (single supplement; current 2×3 layout)
  - 3'UTR outside-gene-matched aside (¶5 SFx) → existing supplement (verify)
  - Isopair Splice Events / Isoform Count / Expression Levels / Transcript Length / NMD pairs flowchart → existing supplements

## Summary table

| Item | Status | Action |
|---|---|---|
| §3 PCI prose | Out of scope (no Fig 3/4 refs) | None |
| §4 ¶1 Isopair construction | ✓ correct | Verify upstream Rmd produces 3,067 / 3,107 / 2,906 |
| §4 ¶2 splice prevalence | ✓ correct | None |
| §4 ¶3 PTC analysis (n=190) | ✓ all numbers correct | None |
| §4 ¶4 NMD+/PTC− mechanism | ✓ mostly correct | **Item 1**: drop "novel" framing in the TD2-bias sentence |
| §4 ¶5 3'UTR outside-set aside | ⚠ unverified numbers | **Item 2**: confirm n=1,904 / 22,335 / 605 / 919 still produced by current Rmd |
| "Fig X" placeholders | Standard pre-submission | Fill in final numbering |
| "SFx" placeholders | Standard pre-submission | Map to final supplement letters |

## Provenance

- This audit fetched the Google Doc via WebFetch on 2026-06-15.
- The figures + supplements it audits against are the current state in `figures/multipanel/figure{3,4}_*` and `figures/SupplementalFigures/{CDSand3UTR_GENCODEonly, TD2BiasEvidence}/`.
- Supersedes and replaces the prior multiple find/replace docs (`section3_findreplace_2026-06-15.md`, `section3_findreplace_2026-06-15_v2.md`, `section4_findreplace_2026-06-13.md`, `section4_findreplace_2026-06-15.md`) which were deleted as part of this consolidation.
- The earlier "audit" that incorrectly read `paper_outputs/manuscript_working.md` (an old Rmd export) and produced a long list of false-positive "stale numbers" was based on a stale local file, NOT the live Google Doc. Those findings do NOT apply to the current manuscript.
