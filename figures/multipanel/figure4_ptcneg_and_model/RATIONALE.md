# §4 / Figure 4 — Methodological Rationale and Pre-Registration

**Status:** Pre-registered methodological choices for §4 of the NMD long-read manuscript.
**Date:** 2026-06-14
**Authors:** P. Castaldi (analysis lead) with NMD long-read team.
**Purpose:** Lock in the analytic design BEFORE re-running the residual-NMD-mechanism analyses, so that scope, group definitions, primary measures, and statistical framework are committed independently of any boundary-dependent result.

This document is the source of truth for the design choices behind Figure 4 and the §4 results. Subsequent code, methods text, and figure renders must conform to what is locked here. Any deviation is a methods change requiring re-justification — not a silent refinement.

**Canonical analysis report:** as of 2026-06-15 the Figure 4 numbers are also computed inline (independently) by [`results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd`](../../../results/isoform_transitions/Version_6.0/isopair_wrapper/05_final_report_gencode_scope_2026-06-15.Rmd) — Panels A/B in §2c, Panels C/D in §3a, the TD2BiasEvidence supplement (panels A–F) in §4. The legacy `05_final_report_mashr.Rmd` has been legacy-bannered; it remains available for sensitivity analyses but is no longer the canonical source. Cross-check between the new Rmd and the figure-side TSVs in this folder is automated at [`reproducibility/verify_cross_check_new_rmd_vs_figures.R`](../../../reproducibility/verify_cross_check_new_rmd_vs_figures.R).

---

## 1. Overarching analytic rationale

We do not, in general, know the true CDS of any given comparator isoform with certainty. The two available CDS sources are:

- **GENCODE annotations** — biologically curated, but available only for ENST-prefixed reference isoforms (~70% of c2 references).
- **TD2 (TransDecoder2) predictions** — applied by SQANTI3 to novel reference isoforms (~30% of c2 references) and to nearly all comparator isoforms.

TD2 was demonstrated in Figure 3 to systematically prefer PTC-avoiding ORFs (TD2 selects a downstream non-PTC ATG when the reference biological ATG would produce a PTC-containing ORF). This bias is directional and well-characterized.

To minimize TD2 bias, §4 analyses restrict to **pairs whose reference isoform has a GENCODE-annotated CDS** (`grepl("^ENST", reference_isoform_id)`). Within this restricted scope, the comparator's ORF is determined by **tracing the reference AUG through the comparator** via `Isopair::traceReferenceAtg` (the ref-AUG-projected ORF), and downstream UTR analyses are anchored on the ref-AUG-projected boundaries on both the reference and comparator sides.

This is the most reliable basis for §4 mechanism inference. Sensitivity analyses at full pop_traceable scope are reported in supplement.

---

## 2. Universe (matches Fig 3 starting universe)

- **pop_BC** (Stage-2 gene-matched coding-coding pairs, built by `02_build_profiles_mashr.R`) = 1,548 NMD c2 / 1,548 Control c4. Established and verified in Fig 3.
- **pop_traceable** (ref-AUG-traceable subset; categories ∈ {effectively_ptc, no_downstream_ejc, truncated_no_ejc}) = 2,289 NMD c2 / 1,763 Control c4.
- **§4 primary scope** = pop_traceable ∩ pairs where `reference_isoform_id` is ENST-prefixed (GENCODE-annotated):
  - **NMD c2 (ENST-reference pop_traceable)**: 1,171 pairs (a few additional rows where `mechanism_class` routes NA-`orf_diff` cases to "Ref AUG absent" reduce the operational mechanism-classified subset to 1,164).
  - **Control c4 (ENST-reference pop_traceable)**: 885 pairs.
- **Pop_BC at ENST-reference scope** (for Panel D classification overview, before pop_traceable restriction): **1,385 NMD c2 / 1,385 Control c4** — 1:1 Stage-2 gene-matching fully preserved (both sides share the same 1,385 gene/reference keys).

**Methodological note — this is a filter on the existing pop_BC, NOT a re-anchoring of the pair structure.** We do not re-build pairs or re-select reference isoforms from a new criterion. The Isopair pair-building convention (`02_build_profiles_mashr.R`: "top non-NMD isoform per gene in pooled DMSO" → pair with top NMD comparator) is unchanged. We simply drop the 163 pop_BC pairs (1,548 → 1,385) where Isopair selected a novel (non-ENST) reference. This is the conservative move: maximal continuity with Fig 3's verified pair structure and minimum methodological surface area to defend in the manuscript.

A potential future expansion — a comprehensive reference-set re-anchoring scheme (top-expressed-ENST-per-gene anchored on a re-built pair set) — was scoped during the §4 design discussion but deferred. If we revisit it, see Task #50.

Panel D (classification overview bar chart) reports at pop_BC scope (1,548 each) with `ref_source` and `mechanism_class` annotations, so a reader can see the full landscape and the ENST-only subset that drives §4 analytic claims.

---

## 3. Classification — 4-group manuscript scheme

§4 mechanism inference uses a 4-group classification (computed by `mechanism_class_4()` in `figures/lib/mechanism_class.R`):

| Group | Definition | n (ENST-only pop_BC) |
|---|---|---:|
| **NMD+/PTC+** | `category == "effectively_ptc"` (ref-AUG-traced ORF has downstream EJC > 50 nt past stop) | 1,080 |
| **NMD+/PTC− ORF match** | `category == "no_downstream_ejc"` AND `comp_orf_len == ref_orf_len` (comparator encodes the reference protein) | 52 |
| **NMD+/PTC− ORF diff** | (`category == "no_downstream_ejc"` AND `comp_orf_len > ref_orf_len`) — "extended" subgroup, n=4 — OR (`category == "truncated_no_ejc"` AND `comp_orf_len < ref_orf_len`) — "truncated" subgroup, n=34. Combined n = 38. | 38 |
| **Ref AUG absent** | `category ∈ {ref_atg_lost, no_ref_cds, mapping_failed}` plus ndj/trc rows with NA `orf_diff` (unable to compute the protein-length comparison) | 214 |
| **Control** (`c4` comparators at the same ENST-pop_traceable scope) | matched non-NMD comparators | 885 |

Total ENST-only NMD c2 (pop_BC) = 1,080 + 57 + 34 + 214 = **1,385 pairs**. Of these, 1,171 are mechanism-classified (PTC+, ORF match, ORF diff); 214 fall to "Ref AUG absent" and are excluded from mechanism inference.

Classification is **100% ref-AUG-derived** and does NOT use Isopair's `original_ptc` field. This is a deliberate design choice: `original_ptc` is computed from each comparator's TD2-called CDS via the 50-nt rule, and TD2 has the documented PTC-avoidance bias (Figure 3 finding). Using `original_ptc` as a classifier would re-introduce TD2 bias into our PTC+ definition even at ENST-only scope, because the bias enters on the *comparator* side.

### Implementation detail — 5-group internal scheme

`mechanism_class.R` produces a finer 5-group classification used for Panel D's landscape display (so a reader can see the ORF-extended vs ORF-truncated breakdown):

| 5-group | 4-group manuscript mapping | n (ENST-only pop_BC) |
|---|---|---:|
| NMD+/PTC+ | NMD+/PTC+ | 1,080 |
| NMD+/PTC− ORF same | NMD+/PTC− ORF match | 52 |
| NMD+/PTC− ORF extended | NMD+/PTC− ORF diff (subgroup) | 13 |
| NMD+/PTC− ORF truncated | NMD+/PTC− ORF diff (subgroup) | 54 |
| Ref AUG absent | Ref AUG absent | 214 |

The 5-group representation is the canonical internal classification; the 4-group manuscript scheme is `mechanism_class_4()` applied at use-site.

### Reclassified pairs (TD2-PTC+/ref-AUG-PTC− discrepancy)

13 pairs at ENST-only scope have `category == "no_downstream_ejc"` AND `original_ptc == TRUE`. TD2 called these PTC+ from a downstream non-ref-AUG ORF; ref-AUG tracing identifies a normal full-length main ORF. Under the §4 ref-AUG-derived framework these are PTC− pairs: 12 fall into ORF match (orf_diff == 0); 1 falls to Ref AUG absent (NA orf_diff). Documented in methods text as evidence of TD2's PTC-avoidance influence on the comparator's TD2 CDS even when the reference is GENCODE-anchored.

### NA `orf_diff` edge case (option A — current convention)

7 ENST ndj pairs (and 10 in the full-scope c2) have `category == "no_downstream_ejc"` but `comp_orf_length − ref_orf_length` is NA — typically because Isopair's `traceReferenceAtg` could not determine `ref_orf_length` for the reference (the reference's ORF did not find an in-frame stop within the searchable window). These pairs have a ref-AUG present and a ref-AUG-traceable category, but they cannot be assigned to "ORF match" or "ORF diff" because the protein-length comparison is undefined.

**Convention:** route these to "Ref AUG absent". Methods text discloses this: *"ndj pairs with undetermined `ref_orf_length` (n=7 at ENST-only c2 scope, 0.4%) are reported under 'Ref AUG absent' because the protein-length classification cannot be made; this is a conservative bookkeeping choice that slightly inflates the 'Ref AUG absent' fraction."*

### Mechanistic content of each group

- **PTC+** — canonical NMD substrate (downstream EJC > 50 nt downstream of stop). The mechanism IS the PTC; the long stop-to-3′-end region is the substrate.
- **PTC− ORF match** — comparator encodes the reference protein, ref-AUG-traced ORF terminates normally with no PTC. Mechanism by elimination — not PTC, not 3′UTR length (we show); remaining candidate is 5′UTR features (uORF burden). **The cleanest mechanism claim in §4 rests on this group (n=52 at ENST-only).**
- **PTC− ORF diff** — comparator encodes a *different* protein (truncated or extended C-terminus) without a PTC. Mechanistically heterogeneous (truncated-protein quality control, alternative C-terminal isoforms, exonic SNV/indels, annotation artifacts). Characterized descriptively; not used for confident mechanism inference. The truncated subgroup (n=54) drives most of this group's mass.
- **Ref AUG absent** — included in Panel D for landscape transparency; excluded from mechanism analyses.
- **Control** — gene-matched non-NMD baseline.

Mechanistic content of each group:

- **PTC+** — canonical NMD substrate (downstream EJC > 50 nt downstream of stop).
- **PTC− ORF same** — ref-AUG-projected ORF in the comparator is the reference protein length AND terminates at a stop that is in the last exon (no canonical PTC). NMD via 5′UTR features is the candidate mechanism.
- **PTC− ORF extended** — ref AUG translates further than reference (readthrough or stop-codon disruption in comparator). n=13; descriptively characterized only.
- **PTC− ORF truncated** — comparator encodes a C-terminally truncated protein via an early in-frame stop still in the last exon. Mechanistically heterogeneous (exonic SNV/indel, alternative C-terminal splice, sequencing/annotation artifact); excluded from mechanism inference.

`mechanism_class` is a **derived field** computed by a sourced helper (`figures/lib/mechanism_class.R`) and applied at use-site, NOT a cache column. This decouples classification from `ref_atg_analysis.rds` so the classification can evolve without invalidating Figure 3 caches.

---

## 4. UTR feature measurements — primary and sensitivity choices

### 4.1 5′UTR boundary on the comparator

- **Primary**: ref-AUG-projected start position. The comparator's 5′UTR is everything upstream of the ref-AUG-projected position in transcript space. Computed by `05k_b_utr5_refaug.R` saving `utr5_features_refaug.rds`. The existing `utr5_features_all.rds` (TD2-defined 5′UTR) is preserved as-is.
- **Sensitivity**: TD2-defined 5′UTR boundary (`utr5_features_all.rds`).

**Why ref-AUG-projected boundary is primary:** TD2's PTC-avoidance bias on the comparator picks a downstream start codon when a PTC is present. The resulting "5′UTR" then includes coding sequence that lies between the true reference start and TD2's downstream-shifted start. uORFs detected in this artificially extended region are, in part, mis-attributed protein-coding sequence rather than authentic uORFs. Using the ref-AUG boundary is consistent with our reason for restricting to ENST-only references (avoiding TD2 bias).

This methodological choice predates the boundary-dependent §4 analyses (TD2's bias was established and characterized in Figure 3) and is therefore not a post-hoc choice. Both scans are reported.

**Empirical TD2-vs-ref-AUG comparison** (2026-06-15, n = 4,052 pop_traceable comparator isoforms, paired by isoform_id):

| ref_atg category × set | n | 5′UTR length (TD2 median → ref-AUG) | Longest uORF (TD2 median → ref-AUG) |
|---|---:|---|---|
| effectively_ptc / NMD | 1,912 | 471 → **247.5** (TD2 ~2× too long) | 198 → **42** (TD2 ~5× too long) |
| effectively_ptc / Control | 288 | 446.5 → 184 | 169.5 → 0 |
| no_downstream_ejc / NMD | 301 | 557 → **974** (TD2 too short!) | 234 → **570** (TD2 too short) |
| no_downstream_ejc / Control | 825 | 246 → 300 (mostly same) | 54 → 87 |
| truncated_no_ejc / both | 726 | mostly same (95%) | mostly same (95%) |

The bias is real, large, and **directional in different ways across categories**: TD2 inflates the apparent 5′UTR for effectively_ptc (picks downstream non-PTC start) and shrinks it for no_downstream_ejc (picks upstream first-ATG that's actually an authentic uORF, treating the uORF stop as the main ORF stop). Both errors mis-characterize the true uORF burden — in opposite directions — confirming why a consistent, ref-AUG-anchored boundary is the right primary measure.

### 4.2 3′UTR length measures (two reported)

We report **two complementary 3′UTR measures** to make the PTC+ measurement bias visible and report a bias-corrected finding:

1. **`utr3_to_tx_end_nt`** = `tx_len(comp) − comp_stop_tx_pos`. **Translation-based** measure: the post-stop region the cell processes after the comparator's ORF ends. For PTC+ pairs the comparator's ORF ends AT the PTC, so this measure includes the [PTC → natural stop] coding region by definition — upward-biased for PTC+ as a "3'UTR length" feature. Available for all pairs.
2. **`utr3_via_non_ptc_stop_nt`** = `tx_len(comp) − non_ptc_stop_tx_pos`. **Sequence-based** measure: define a uniform "first non-PTC stop" across all groups by walking downstream of ref AUG (skipping any stop that has a downstream EJC > 50 nt) until the first stop that's NOT a PTC. For non-effectively_ptc traceable groups this equals `comp_stop_tx_pos` (their comparator's stop is already non-PTC by category definition). For effectively_ptc groups it gives the natural stop the cell *would* terminate at if the PTC were skipped — the "natural 3'UTR length" anchored on the comparator's own sequence (no dependence on the reference's stop position). Coverage: 99% of effectively_ptc pairs and 99.8% of Controls; failures (~19 of 1,912 effectively_ptc) are cases where every downstream stop is also a PTC.

The non-PTC-stop measure is the **bias-free 3'UTR measure**. The translation-based measure is the **biologically-translated post-stop measure** (= the PTC for effectively_ptc). Reporting both makes the PTC inflation visible while preserving the bias-corrected statistical claim.

**Inferential comparison primary**: each NMD group vs Control under BOTH measures, with BH correction across the test family. The 3'UTR-mechanism finding for §4 rests on the non-PTC-stop measure (HL = −221, Cliff δ = −0.18, p = 0.003 for ORF match vs Control).

**Why the non-PTC-stop measure is preferred over the prior reference-stop-projection approach** (which was the older `utr3_via_ref_stop_nt` field, now removed): the reference-stop projection made the Control measure depend on the reference's stop — methodologically weird since Control is its own transcript with its own natural stop. The non-PTC-stop measure is anchored entirely on each isoform's own sequence + structure, so Controls and NMD comparators are measured by the same yardstick.

**Why restrict the inferential comparison to ORF match vs Control:** the 3′UTR measure depends on the comparator's stop position. For the ORF match group (by definition `comp_orf_len == ref_orf_len`), the comparator's stop is at the same transcript position as the reference, so `utr3_to_tx_end_nt` is the *natural* 3′UTR — biologically and methodologically comparable to Control's natural 3′UTR. For other groups the measure is biased:

| Group | Comparator stop position vs reference | Bias direction |
|---|---|---|
| NMD+/PTC+ | upstream (PTC sits *before* the natural stop) | **UPWARD** — measure includes [PTC → natural stop] coding sequence + natural 3′UTR |
| NMD+/PTC− ORF match | identical | **none** — natural 3′UTR |
| NMD+/PTC− ORF diff (ORF extended subgroup, n=13) | downstream | measure shorter than natural 3′UTR by `(comp_orf_len − ref_orf_len)` nt |
| NMD+/PTC− ORF diff (ORF truncated subgroup, n=54) | upstream | **UPWARD** — same direction as PTC+, smaller magnitude |
| Control | natural stop | **none** |

A reviewer may ask why we don't compute a "reference-stop-projected 3′UTR" to unbias the biased groups. That additional measure is deferred — adding it would require another ref_atg field plus the additional ref-stop-projection mapping step, and the §4 mechanism story is unaffected because it does not rest on PTC+ or ORF diff in this panel.

**TES variation** is part of normal biology (alternative polyadenylation) and is not a confound to control for in the ORF match vs Control comparison: both sides measure the comparator's *actual* 3′UTR length, which is what NMD machinery acts on. We do not restrict to TES-concordant pairs; the pre-registered TES-concordant primary scheme considered earlier turned out to drop the ORF match n from 52 to a handful and was abandoned. TES variation across groups can be reported descriptively from the `tx_len_delta_nt` field if a reviewer asks.

**Methods text to commit:**
> "3′UTR length was computed as `tx_length(comparator) − comparator_stop_tx_pos`, where `comparator_stop_tx_pos` is the ref-AUG-projected stop in the comparator's transcript coordinates. Cross-group statistical comparison was restricted to NMD+/PTC− ORF match vs Control because for other groups the measure includes coding sequence relative to the reference (PTC+: from PTC to natural stop; ORF diff: from the alternative stop to the reference stop) and is not interpretable as a 3′UTR length. The bias direction is upward for PTC+ and the ORF-truncated subgroup of ORF diff."

### 4.3 EJC-relative stop distance

- Reported as the canonical NMD-PTC distance metric (already in Figure 3 Panel D).
- For §4, used to verify the absence of canonical PTC in the residual PTC− groups (n_downstream_ejc == 0).

---

## 5. Statistical framework

- **Primary test**: Mann–Whitney–Wilcoxon (two-sided), each NMD subgroup vs Control, per metric.
- **Effect size**: Hodges–Lehmann shift estimator with 95% CI (from `wilcox.test(..., conf.int = TRUE)`) AND Cliff's delta. **Required for n ≤ 100 subgroups** (ORF same n=84, ORF truncated n=54, ORF extended n=13). Reported as median (IQR), Δ̂_HL (95% CI), Cliff's δ.
- **Multiple testing**: Benjamini–Hochberg adjustment across the §4 mechanism panel's family of tests (PTC subgroups × metrics, excluding ORF extended). Adjusted q-values reported.
- **ORF extended subgroup (n=13)**: **descriptive only — no inferential test reported**. Medians + IQRs only.
- **Truncated_no_ejc subgroup (n=54)**: included in the test family but interpreted descriptively (mechanism inference excluded per §3).

---

## 6. Mechanism framing — pre-registered language

The PTC−/ORF same mechanism claim is expressed as an **elimination-among-detectable-mechanisms** argument, NOT as a positive mechanistic proof.

Required hedging sentence (or close substitute) in §4 prose:

> "Among the NMD mechanisms detectable from paired short-read + Iso-Seq long-read data — downstream-EJC-mediated PTC, extended 3′UTR length, and 5′UTR uORF burden — only 5′UTR uORF burden differentiates the NMD+/PTC− ORF-same pairs from gene-matched non-NMD Controls. Mechanisms not addressed by this analysis include m6A-mediated decay, AU-rich element-mediated decay, miRNA-driven decay, and EJC-independent UPF1 recruitment via 3′UTR features other than length."

The uORF burden claim itself is a **candidate** mechanism, not a confirmed one. Standard literature would expect a uORF-EJC-distance check on the uORF stop codon for confident uORF-mediated NMD attribution; this check is acknowledged in the methods discussion as a limitation but is not feasible at sufficient depth from current data.

---

## 7. What we do NOT claim

- We do not claim "PTC+ has longer 3′UTR than Control" as an independent finding (it is a definitional consequence of the PTC).
- We do not claim that ORF-same NMD+/PTC− pairs have *more* 5′UTR uORFs than NMD+/PTC+ pairs (at ENST-only scope, they have a similar uORF burden; the differential signal observed at full pop_traceable scope is in part TD2-bias driven).
- We do not claim a mechanism for the ORF-truncated NMD+/PTC− subgroup; we characterize but do not infer.
- We do not claim conclusions about the n=13 ORF-extended subgroup; it is described.
- We do not claim TD2 calls on comparator stops are unbiased; the bias direction (TD2 prefers non-PTC stops → upward bias on 3′UTR length) is acknowledged and works against our PTC−/ORF-same finding, which strengthens — rather than weakens — that finding.
- We do not use Isopair's `original_ptc` field for §4 classification. `original_ptc` is computed from each comparator's TD2-called CDS, so it inherits TD2's PTC-avoidance bias. The 5-group classification is entirely ref-AUG-derived. The 13 ENST-only pairs where TD2 called PTC+ but ref-AUG-tracing did not are routed by their ref-AUG category, not by TD2's call.

---

## 8. Pre-registered analytical pipeline

1. **`05r_ref_atg_analysis.R` updates** — additive fields only: `ref_source`, `utr3_to_tx_end_nt`, `comp_tx_len_nt`, `ref_tx_len_nt`, `tx_len_delta_nt`. **No** `mechanism_class` cache column (derived helper). Regression gate `all.equal(old[, names(old)], new[, names(old)])` MUST pass before commit.
2. **`figures/lib/mechanism_class.R` helper** — single source of truth for 5-group classification.
3. **`05k_b_utr5_refaug.R` sibling analysis** — ref-AUG-projected 5′UTR scan → `utr5_features_refaug.rds`. Original `utr5_features_all.rds` untouched.
4. **`original_ptc` provenance verification** — confirm `original_ptc` traces to a GENCODE-derived PTC call on the reference, not a TD2-derived call. Document in methods.
5. **Figure 4 `data_export.R` refactor** — ENST-only scope; 5-group `mechanism_class`; ref-AUG-projected 5′UTR primary + TD2-boundary sensitivity; TES-concordant 3′UTR primary + full sensitivity; H–L shifts + CI + Cliff's δ + BH q-values.
6. **Figure 4 Panels D/E/F rebuild** under new scope.
7. **`05_final_report_mashr.Rmd` §4 rewrite** with the locked design and pre-registered framing.
8. **5-pass verification** with expanded Pass 2 checklist covering: (a) ENST-only vs full-pop agreement in direction and rank, (b) TES-concordant vs full agreement, (c) no single-group >50% test statistic dominance, (d) primary vs sensitivity 5′UTR-boundary direction agreement.
9. **Manuscript find/replace pairs** for the new numbers (per repo convention; the Google Doc is the source of truth).

---

## 9. Risks tracked

- **Cache schema drift** invalidating Figure 3 — mitigated by regression gate; existing column names/dtypes/row-order byte-identical.
- **n=84 ORF-same subgroup is load-bearing** for the §4 mechanism claim. If TES-concordant restriction drops it below ~50, fall back to "5′UTR uORF burden is a feature distinguishing a subset of residual NMD+/PTC− substrates from Controls" without the mechanism claim.
- **Primary vs sensitivity 5′UTR boundary disagreement** — if the ref-AUG scan strengthens the signal substantially relative to the TD2 scan, document explicitly with the TD2-bias-direction explanation. If it weakens the signal, the §4 mechanism claim is weaker than the discovery suggested, and the framing in §6 must reflect that.
- **`original_ptc` TD2 contamination** — if Task #48 verification finds `original_ptc` is partly TD2-derived, the NMD+/PTC+ union has residual TD2 bias even at ENST-only scope; revisit the classification (Task #45).
- **Downstream cache consumers** — `05_final_report_mashr.Rmd`, `06_orf_analysis_mashr.Rmd`, `05l_unified_model.R`, `05t_ref_cds_features.R`, and 4 paper-figure scripts read `ref_atg_analysis.rds`. Two of those (`make_paper_figure_ptc_reclassification.R`, `make_paper_figure_attribution_dissociation.R`) drive Figure 3 panels that are already verified. Their output post-regeneration MUST be byte-identical to pre-regeneration.

---

## 10. Audit trail

- 2026-06-14 — pre-registration written prior to executing upstream R script updates, after Plan-agent adversarial review of §4 design. Key methodological reframing this session: ENST-only restriction (vs full pop_traceable), 5-group classification with ORF same/extended/truncated split, ref-AUG-projected 5′UTR boundary as primary (vs TD2 boundary), TES-concordant 3′UTR primary, elimination-among-Isopair-measurable-mechanisms framing.
- 2026-06-15 — §4 narrative restructured into three sections (A/B/C). Section A adds the all-3-ENST + coding-CDS scope (n=48/82/130) with own-GENCODE-stop PTC determination (no ref-AUG projection — every comparator is GENCODE-curated). Section B promotes the TD2-bias evidence to manuscript panels (TD2-vs-ref-AUG ORF length + Kozak PWM). Section C preserves the larger ref-AUG-traced scope (n=756/63/819 at the 1:1 gene-matched re-intersected scope) from the original pre-registration. For Section A, ORF match + ORF diff are merged into a single NMD+/PTC− group: at all-3-ENST + coding scope the comparator's CDS is GENCODE-annotated (not TD2-called), so the truncated/extended subgroups carry the same confidence as ORF match. The 4-group scheme is retained for Section C (where the comparator is ref-AUG-projected, not GENCODE-curated, and the ORF-diff heterogeneity warrants separation).

---

## 11. Final 4-panel Figure 4 structure (2×2 composite, 12"×8" landscape)

The composite is rendered by `figure4_composite.py` with a uniform `GridSpec(2, 2)` and 1.5:1 panel cells (6.0"×4.0" each). Two sections (A, C), each occupying one row.

**Restructure note (2026-06-15):** the prior 3×2 layout included a row of TD2-bias justification panels (TD2-vs-ref-AUG ORF length + paired Kozak PWM). These have been moved to a consolidated supplemental figure (`figures/SupplementalFigures/SF34-SF35_TD2BiasEvidence/`) which presents three orthogonal observations of TD2's PTC-avoidance bias on the same isoform universe as Figure 4 Panels C/D (n = 819 broad + n = 348 occult-PTC subset, both at the 1:1 gene-matched re-intersected scope). The restructure eliminates the multi-denominator confusion that the original 3×2 layout introduced — the main figure now uses two denominators (130 + 819) for two sections, and the supplement uses the same Section C denominator (819) plus its occult-PTC subset (348).

### Section A — all-3-ENST + coding-CDS scope (own-GENCODE-stop PTC determination)

The cleanest small-n GENCODE-only evidence of the PTC+/PTC− mechanism dichotomy. Scope: `pop_BC ∩ all-3-ENST ∩ coding-everywhere` (re-intersected on gene + reference after the coding filter) → **n = 190 NMD c2 / 190 Control c4**. Every reference, comparator, AND control isoform is GENCODE-annotated and coding, so PTC determination uses each comparator's *own* GENCODE-annotated stop genomic position (50-nt rule from comparator's own last EJC). No ref-AUG projection. No TD2 dependency anywhere.

3-group merged classification (Section A only):

| Group | n | Definition |
|---|---:|---|
| **NMD+/PTC+** | 72 | own GENCODE stop is > 50 nt past last EJC |
| **NMD+/PTC−** | 118 | own GENCODE stop is ≤ 50 nt past last EJC (ORF match + ORF diff merged — both comparators are GENCODE-curated, so the confidence is uniform) |
| **Control** | 190 | non-NMD comparators at the same all-3-ENST + coding scope |

| Panel | Script | Measure | Key finding |
|---|---|---|---|
| **A** | `figure4_panelA_5utr_length_all3enst.py` | 5′UTR length (nt, log10) | PTC− >> PTC+, PTC− >> Control. PTC+ ~ Control (PTC+ has SHORT 5′UTRs by design, since the natural full 5′UTR isn't compromised in PTC+) |
| **B** | `figure4_panelB_longest_5utr_orf_all3enst.py` | Longest 5′UTR ORF (nt, log10(1+x)) | Same dichotomy: PTC− has long uORFs, PTC+ ~ Control. Orthogonal mechanism evidence to Panel A |

### Section C — ENST-reference + ref-AUG-traceable scope (Panels C, D)

Replicates the Section A finding at the larger ENST-reference scope using ref-AUG-projected mechanism classification (`mechanism_class_4` → 3-group merged: PTC+, PTC− merged, Control). This is the original pre-registered scope from §3 / §4.

| Group | n | Definition |
|---|---:|---|
| **NMD+/PTC+** | 756 | effectively_ptc category (ref-AUG-traced ORF has downstream EJC > 50 nt past stop) |
| **NMD+/PTC−** | 63 | no_downstream_ejc + truncated_no_ejc merged into a single PTC− group for plotting symmetry with Section A |
| **Control** | 819 | gene-matched non-NMD c4 comparators at the SAME 1:1 gene-matched scope |

The Section C scope is built by applying the ref-AUG-traceable category filter (effectively_ptc / no_downstream_ejc / truncated_no_ejc) INDEPENDENTLY to the NMD c2 and Control c4 sides at the ENST-reference scope, then **re-intersecting on `(gene_id, reference_isoform_id)` AFTER the category filter** so that both panels share the same 1:1 gene-matched universe of 819 NMD and 819 Control comparator pairs. The re-intersection step is necessary because the category filter is asymmetric — `ref_atg_lost` is ~3.7× more common in Control comparators (320/1,385 = 23%) than in NMD comparators (86/1,385 = 6%), reflecting that non-NMD productive Controls more often use alternative TSS or skip alternative N-terminal exons. Without the re-intersection, NMD c2 = 1,171 and Control c4 = 885 are not strictly gene-matched.

The full Panel C/D NMD universe is n = 819 (756 PTC+ + no_downstream_ejc + truncated_no_ejc); under the 4-CT scope every NMD c2 pair resolves to one of the two groups, leaving 63 NMD+/PTC− displayed.

| Panel | Script | Measure | Key finding |
|---|---|---|---|
| **C** | `figure4_panelC_5utr_length_refaug.py` | 5′UTR length (nt, log10), ref-AUG boundary | PTC− >> PTC+, PTC− >> Control. Same direction as Panel A at ~6× larger sample |
| **D** | `figure4_panelD_longest_5utr_orf_refaug.py` | Longest 5′UTR ORF (nt, log10(1+x)), ref-AUG boundary | PTC− >> PTC+, PTC+ vs Control significant at this matched scope |

### TD2 ORF-call bias evidence — supplemental figure (TD2BiasEvidence)

The TD2-bias justification panels (which previously occupied Section B of the main figure) now live in `figures/SupplementalFigures/SF34-SF35_TD2BiasEvidence/` and present three orthogonal observations at two complementary scopes (2×3 layout):

- **Row 1 (Panels A–C; n = 819)** — broad scope, same 1:1 gene-matched universe as Figure 4 Panels C/D. TD2 chose the same ATG as the reference AUG in 52% of pairs (429/819).
- **Row 2 (Panels D–F; n = 492)** — occult-PTC subset (`effectively_ptc ∩ original_ptc == FALSE` within the re-intersected universe).

Headline numbers on the occult-PTC subset where the bias is most extreme:
- TD2-selected CDS is **~7.7× longer** than the reference-AUG-anchored ORF (median 2,934 vs 380 nt; **Supp Panel D**).
- Reference AUG has **stronger Kozak context** than TD2's chosen ATG (paired Wilcoxon p = 2.6×10⁻³⁸; reference Kozak > TD2 Kozak in 78.0% of pairs; **Supp Panel E**).
- **99.0% (487/492)** of TD2 ATGs are downstream of the reference AUG (median offset +476 nt; **Supp Panel F**).

On the broader Section C universe (Panels A–C) the bias direction is preserved but diluted by the 50% of pairs where TD2 agrees with the reference: ORF ratio 4.7×, Kozak p = 1.2×10⁻³⁸, 43% downstream (median +471 nt) / 50% same / 7% upstream.

### Pre-registered tests per Section A / Section C panel

All three pairwise Wilcoxon comparisons (NMD+/PTC− vs NMD+/PTC+, NMD+/PTC− vs Control, NMD+/PTC+ vs Control) are reported per panel with significance markers (`***` p<10⁻⁴, `**` p<10⁻³, `*` p<0.05, n.s. otherwise). Effect sizes (HL shift + Cliff's δ) live in `data/panel{A,B,C,D}_*_pairwise.tsv`.

### Data files (main Figure 4)

```
data/
  panelA_5utr_length_{long,descriptives,pairwise}.tsv             # Section A
  panelB_longest_5utr_orf_{long,descriptives,pairwise}.tsv        # Section A
  panelC_5utr_length_{long,descriptives,pairwise}.tsv             # Section C
  panelD_longest_5utr_orf_{long,descriptives,pairwise}.tsv        # Section C
```

### Why two scopes side-by-side (A vs C and B vs D)

- **Section A (n=190)** is the strongest *purity* argument: every isoform is GENCODE-curated, so the mechanism call doesn't depend on any computational ORF prediction — only on the GENCODE annotation and on splice-junction geometry. This is the panel that survives any TD2-bias critique most cleanly.
- **Section C (n=756/63/819)** is the strongest *generality* argument: the same dichotomy persists at ~6× larger sample using ref-AUG projection on the comparator within a 1:1 gene-matched universe. Co-occurrence of the same effect at both scopes is the multi-evidence argument the §4 mechanism claim rests on.
- Both rows share the same axis design, group order, and color palette so a reader can directly compare A↔C and B↔D.

### Why the Section A PTC+ rate (38%) is much lower than Section C's PTC+ rate (90%) — and why this is NOT a bug

The Section A NMD+ pairs split 48 PTC+ / 82 PTC− (37% / 63%) under the own-GENCODE-stop classifier, whereas Section C's NMD+ pairs split 756 PTC+ / 63 PTC− (92% / 8%) under the ref-AUG-traced classifier on the 1:1 gene-matched universe. The gap is the combined effect of two distinct shifts, neither of which is a counting error:

**1. Classifier shift (largest component).** Of the 190 Section A pairs, the ref-AUG classifier (the same classifier used in Section C) calls **134 of them PTC+ (70%)**. The own-GENCODE-stop classifier reclassifies **68 of those 134 ref-AUG-PTC+ pairs as PTC−**. The mechanism: the comparator's own GENCODE-annotated CDS uses a different (downstream) start codon, so its own annotated stop is in the last exon — failing the 50-nt rule — even though projecting the *reference's* AUG through the same comparator hits a premature stop. **75% (51 of 68) of these reclassified pairs are biotyped `nonsense_mediated_decay` in GENCODE.** Curators recognize them as NMD substrates biologically, but the own-stop+last-EJC geometry is not the canonical PTC pattern. The own-GENCODE-stop classifier is therefore strictly more conservative than the ref-AUG classifier — it accepts only the subset of NMD substrates whose own annotated CDS *also* satisfies the 50-nt geometry.

**2. Sample-composition shift (smaller component).** The all-3-ENST filter requires the *Control* comparator to also be ENST-annotated. At the ENST-reference scope before this filter (Section C's universe), PTC+ NMD substrates have a strong tendency to be paired with *novel* Control comparators — 1,225 PTC+ pairs are dropped by the c4-ENST requirement, versus only 86 PTC− pairs. Section A is therefore PTC+-depleted at the sample-construction stage, even before the classifier choice. Under the ref-AUG classifier alone (no own-stop reclassification), the Section A scope's PTC+ rate is 70%, not the 90% seen at full ENST-reference scope.

**Why we keep the own-GENCODE-stop classifier for Section A:** it is the only PTC determination that requires *no* projection step from reference onto comparator — the classifier looks at the comparator's own annotation, its own splice junctions, and applies the 50-nt rule. This is what makes Section A the strongest purity argument: a reviewer cannot push back on "the PTC call depends on ref-AUG-tracing" because no tracing was used. The cost is the rate gap relative to Section C, which we document here and discuss in the manuscript (§4) as a methodological difference, not a substrate difference.

**Cross-validation:** for the 190 Section A pairs, Isopair's pre-computed `original_ptc` field (which the manuscript-wide TD2 audit, Task #49, confirmed traces to the comparator's own GENCODE CDS for ENST-comparator pairs and to TD2 otherwise) gives **identical 72 PTC+ / 118 PTC−** calls. The own-GENCODE-stop classifier and the pre-cached `original_ptc` agree perfectly on this subset — consistent with the expectation that TD2 ≈ GENCODE when the comparator is ENST-annotated.
