# Manuscript-Wide TD2 Bias Audit

**Status:** in progress — opened 2026-06-14 during §4 / Figure 4 deep-dive.

## The general issue

Figure 3 of the manuscript establishes that **TD2 (TransDecoder2) systematically prefers PTC-avoiding ORFs on novel PacBio transcripts.** Where the true biological CDS would produce a downstream EJC > 50 nt past the stop (i.e., a PTC), TD2 frequently selects an alternative, non-canonical downstream ATG that produces a non-PTC ORF instead.

This means any manuscript analysis that uses TD2-derived CDS information (CDS boundaries, stop positions, 5′UTR boundaries, 3′UTR boundaries, PTC calls) on the comparator side of an NMD-vs-Control pair carries a directional bias:
- "5′UTR" computed from TD2's chosen start is systematically *longer* than the true biological 5′UTR for PTC-containing comparators.
- "3′UTR" computed from TD2's chosen stop is systematically *shorter* than the true PTC-containing 3′UTR (TD2 prefers non-PTC stops).
- "Number of uORFs / longest uORF in 5′UTR" computed from a TD2-defined 5′UTR includes ORFs that lie in the true protein-coding region.
- "PTC+ vs PTC−" classification based on TD2-derived stop positions is biased *toward* PTC− calls.

The bias source is documented (Figure 3, §3 of the manuscript). The implication — that downstream analyses using TD2 CDS information inherit this bias — has not yet been propagated through the manuscript and the analysis pipeline.

This audit tracks which analyses are affected and what remediation each needs.

## What §4 / Figure 4 work uncovered

The §4 analysis revealed three concrete instances of this issue:

1. **`utr5_features_all.rds`** — currently scans the 5′UTR using TD2-defined comparator CDS boundary. For PTC-containing comparators (and any comparator where TD2 chose a downstream start), the scanned region inflates the 5′UTR and likely double-counts true-CDS ORFs as "uORFs". §4 plan: add a ref-AUG-projected scan (`utr5_features_refaug.rds`) as primary, retain TD2-scan as sensitivity. **Mitigation in progress.**

2. **`original_ptc` field (Isopair `comp_has_ptc`)** — computed from each comparator's TD2-called CDS via the 50-nt rule. Inherits TD2 bias. §4 plan: do NOT use `original_ptc` in mechanism classification; classification is 100% ref-AUG-derived. **Mitigation in `mechanism_class.R`.** Documented in `figures/multipanel/figure4_ptcneg_and_model/RATIONALE.md` §3, §7.

3. **3′UTR length on PTC+ pairs** — has a separate measurement-bias issue independent of TD2 (`tx_len − PTC_stop` includes coding sequence between PTC and original stop). For PTC+, §4 does not compare 3′UTR length as an independent finding. Methods text required. **Documented in RATIONALE.md §4.2, §7.**

## Audit checklist — candidates to check

These analyses or RDS files are TD2-dependent or potentially TD2-affected. For each, the audit must determine: (a) where TD2 enters, (b) bias direction, (c) whether published claims are affected, (d) recommended remediation.

### Direct TD2-derived inputs (canonical CDS source)

- [ ] `data_mashr/cds.rds` (comparator-side CDS coordinates from SQANTI3 + TD2 for novel isoforms)
- [ ] `data_mashr/analysis_cache/ptc.rds` (downstream of `cds.rds`)
- [ ] `data_mashr/analysis_cache/ptc_c2_allsamples.rds` (per-pair PTC summary, see §4 finding)
- [ ] `data_mashr/analysis_cache/utr5_features_all.rds` (5′UTR scanned using TD2 5′UTR boundary)
- [ ] `data_mashr/analysis_cache/ref_cds_features_all.rds` (per-isoform CDS features; mix of GENCODE + TD2)
- [ ] `data_mashr/analysis_cache/orfik_scan.rds` (ORFik uORF scan)
- [ ] `data_mashr/analysis_cache/orf_landscape.rds`, `orf_landscape_per_orf.rds`

### Analyses likely affected — to triage

- [ ] **§2 NMD response analyses** — do any quantify "NMD via PTC" using TD2 CDS calls? If yes, are claims directionally affected?
- [ ] **§3 PTC accounting / Output Lost analyses** — these are Figure 3 territory. The TD2 bias claim is itself the §3 finding; analyses that USE that conclusion should be fine. But analyses that quantify "X% of NMD-via-PTC" using TD2 CDS may need scope restrictions.
- [ ] **§4 Isopair-attribution analyses** — covered by the new §4 framework; pre-registered in RATIONALE.md.
- [ ] **§5 deep-learning model (Figure 5)** — model input features include CDS-derived features. To what extent does TD2 bias enter the training inputs? Does the model learn around the bias, or amplify it? Could the model's "PTC prediction" performance be inflated by the bias on the training set?
  - The model is trained on labels derived from limma-voom-mashr-driven NMD-responsive classification (independent of TD2) but uses *sequence features* including CDS-derived ones. TD2 bias enters via the structural-feature branch.
  - Specifically check: 5′UTR / 3′UTR length features, ORF length features, EJC distance features.
- [ ] **Any pathway / functional enrichment analyses** that condition on PTC status — those use TD2's PTC call.
- [ ] **Splice event attribution analyses** (Figure 3 Panel F) — these condition on PTC+ pairs. TD2-bias check.
- [ ] **Per-cell-type DGE / DIE analyses** — independent of TD2 CDS, should be safe.

### Public-facing prose claims to check

- [ ] Methods text: every place CDS / PTC / UTR length is computed — does the text clearly identify which calls are TD2 vs GENCODE-derived?
- [ ] Results text: are claims framed in a TD2-bias-aware manner where applicable?
- [ ] Figure legends: per-panel transparency about which CDS source was used.

## Remediation patterns

Recommended remediation per affected analysis:
- **Scope restriction** — restrict to ENST-only references (GENCODE-CDS-on-reference) where the reference-side input is the dominant determinant. Used for §4.
- **Ref-AUG anchoring** — recompute the analysis using ref-AUG-projected positions on the comparator. Used for §4 5′UTR scan. Requires writing parallel analyses.
- **Documentation only** — for analyses where TD2 bias is acknowledged and works against the claim (i.e., bias is conservative for the finding). Document explicitly in methods.
- **Removal** — for claims that cannot survive a bias correction, remove from the manuscript or relabel as exploratory.

## Triage status (to fill in as audit proceeds)

| Analysis / file | TD2 enters | Bias direction | Affects published claim? | Remediation | Status |
|---|---|---|---|---|---|
| `utr5_features_all.rds` 5′UTR scan | comparator's TD2 CDS boundary | Inflates 5′UTR length / uORF count for effectively_ptc (TD2 picks downstream non-PTC start); shrinks 5′UTR for no_downstream_ejc (TD2 picks upstream uORF-ATG). Either direction mis-characterizes the true uORF burden. | Yes — §4 mechanism candidate | **Done**. Ref-AUG-projected sibling scan saved to `utr5_features_refaug.rds`. Both scans reported; ref-AUG-anchored is primary for §4 per RATIONALE.md §4.1. Original `utr5_features_all.rds` preserved untouched. Confirmed magnitude: effectively_ptc NMD median longest uORF dropped 198 → 42 nt under ref-AUG; no_downstream_ejc NMD median rose 234 → 570 nt. | Done (Task #47, 2026-06-15) |
| `original_ptc` / `ptc.rds` | comparator's TD2 CDS | Biased toward PTC− calls on novel | Yes — §4 classification | §4 uses ref-AUG-derived classifier only | Done (mechanism_class.R) |
| 3′UTR on PTC+ | `tx_len − PTC_stop` (definitional, separate from TD2) | Inflates 3′UTR on PTC+ | Yes — affects PTC+ vs PTC− 3′UTR comparison | Methods caveat + PTC+ comparison demoted | Pre-registered (RATIONALE.md §4.2) |
| `ref_cds_features_all.rds` | mixed (GENCODE + TD2 for novel refs) | varies by row | Maybe — used in §5 model features | TBD | Open |
| §5 DL model inputs | CDS-derived features (5′UTR length, ORF length, etc.) | varies | Unknown — needs investigation | TBD | Open |
| §3 Figure 3 panels | self-aware (TD2 bias IS the topic) | n/a | No — bias is the finding | None | Done (Figure 3 verified) |
| Splice event attribution | conditions on PTC+ | conditions inherit TD2 bias | Maybe | TBD | Open |
| Per-CT DGE / DIE | independent of TD2 CDS | n/a | No | None | OK |
| Methods text TD2 disclosure | text-level | n/a | Possibly missing transparency | Methods audit | Open |

## Audit trail

- 2026-06-14 — TD2 bias on `original_ptc` confirmed and addressed in §4 via ref-AUG-derived classification (RATIONALE.md). General audit document opened to track the broader manuscript-wide implications.
- 2026-06-15 — Considered a comprehensive reference-set re-anchoring scheme (top-ENST-per-gene from re-built pair structure, summed mean CPM across 4 CTs). Computed candidate counts (16,718 genes overall; 9,186 with non-NMD filter; ~11.6% reference change vs OLD allsamples when matching non-NMD filter). Decision: defer. For §4 (and any other analyses requiring TD2-bias mitigation), use the simpler "filter existing pop_BC to ENST-reference" approach. Maintains Isopair pair-building convention, preserves Stage-2 gene-matching 1:1, minimal Fig 3 disruption (panels D/E/F only). Tracked as Task #50 for possible later revisit.
- 2026-06-15 — Built `utr5_features_refaug.rds` (Task #47). Empirical TD2-vs-ref-AUG comparison on 4,052 pop_traceable comparator isoforms confirms the TD2 bias is large and directionally inconsistent: effectively_ptc isoforms had longest 5′UTR ORF inflated ~5× by TD2 (median 198 → 42 nt under ref-AUG); no_downstream_ejc isoforms had it shrunk by TD2 (median 234 → 570 nt under ref-AUG). Both errors push apparent uORF burden in opposite directions, supporting the methodological choice that the ref-AUG-anchored scan is the right primary measure for §4. truncated_no_ejc cases are mostly insensitive to boundary choice (~95% identical). Manuscript-wide audit triage table updated to mark this analysis "Done"; §5 DL-model input audit still open.
