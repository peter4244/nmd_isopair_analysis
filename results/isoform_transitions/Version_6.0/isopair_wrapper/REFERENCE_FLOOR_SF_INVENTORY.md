# SF24–SF42 floor-impact inventory (2026-07-10)

Full audit of every supplemental figure I own (SF24–42) for reference-floor dependence,
data producer, and current state. Built after discovering I'd initially missed SF25/27/28/29.

## Status legend
- **DONE** — regenerated on floored data, visually/data-verified, committed.
- **TODO** — floor-affected, still shows pre-floor data.
- **N/A** — not floor-affected (static, or model-global output; model is NOT retrained).
- **FLAG** — separate deliverable; floor-dependence needs a Pete decision.

| SF | Title | Floor-dep? | Producer | Data state | Status |
|---|---|---|---|---|---|
| SF24 | SpliceEventCategories | No (static glossary) | illustrative coords from Rmd glossary | n/a | **N/A** |
| SF25 | IsoformsPerGene | **Yes** (pop_BC 3,009→1,585) | `code/` cluster (isoform_transition_splicing_analysis.Rmd / filter_major_isoforms.R) — **trace needed** | 07-08 stale | **TODO** |
| SF26 | ReferenceShare | Yes | rebuilt `data_export.R` (this work) | fresh | **DONE** |
| SF27 | TranscriptLengthByRole | **Yes** (ref/NMD/ctrl lengths on pop_BC) | `code/` cluster — **trace needed** | 07-08 stale | **TODO** |
| SF28 | PairAnalysisFlowchart | **Yes** (cohort cascade) | `build_flowchart.R` — **hardcoded** 3,009/190/1,166 in DOT | inline | **TODO** (swap → 1,585/136/888) |
| SF29 | GainDirectionByEvent | **Yes** (pop_BC pairs) | **deprecated legacy `05_final_report_mashr.Rmd`** — provenance problem | 07-08 stale | **TODO** (needs floored producer) |
| SF30 | PTCDistanceDoseResponse | Yes | SF30 `data_export.R` | fresh | **DONE** |
| SF31 | NMDEffectByEJCCount | Yes | SF30 `data_export.R` (writes SF31 too) | fresh | **DONE** |
| SF32 | CdsAnd3UTR_GENCODE | Yes | SF32 `data_export.R` → CDSand3UTR_GENCODEonly | fresh | **DONE** |
| SF33 | TD2Bias_broad | Yes | SF33 `data_export.R` → TD2BiasEvidence | fresh | **DONE** |
| SF34 | TD2Bias_occult | Yes | SF33 `data_export.R` | fresh | **DONE** |
| SF35 | CdsAnd3UTR_refAUG | Yes | SF32 `data_export.R` | fresh | **DONE** |
| SF36 | ShapAcrossWindows | No (model-global SHAP profile, test set) | model `results_4ct/shap_profile_*` | 07-08 (fine) | **N/A** |
| SF37 | StopCodonUsage | No (model test-set NMD calls, n=2,268; "class"=NMD/non-NMD, not PTC subclass) | model `results_4ct/stop_codon_freq_by_class_sf37` | model-repo | **N/A** |
| SF38 | AttentionDistribution | No (model test-set attention) | model `results_4ct/uorf_attention_predictions` | model-repo | **N/A** |
| SF39 | PTCSubclassBranchSHAP | Yes | SF39 `data_export.R` | fresh | **DONE** |
| SF40 | PTCSubclassPerformance | Yes | SF39 `data_export.R` (PERF_OUT fixed) | fresh | **DONE** |
| SF41 | GCcontentStopWindow | No (model-global GC across windows) | model `results_4ct/gc_content_across_*` | 07-08 (fine) | **N/A** |
| SF42 | ModelComparison | **Maybe** (points colored by PTC subclass) | separate `code/nmd_predictor_comparison/*_2026.6.20.tsv` | 06-20 | **FLAG** |

## Summary
- **DONE (9):** SF26, 30, 31, 32, 33, 34, 35, 39, 40.
- **TODO (4):** SF25, SF27, SF28, SF29 — floor-affected, still pre-floor.
- **N/A (5):** SF24, SF36, SF37, SF38, SF41 — static or model-global (model not retrained; the floor changes neither the model's test-set SHAP/attention/GC/stop-codon outputs nor a glossary).
- **FLAG (1):** SF42 — depends on the separate `nmd_predictor_comparison` deliverable; PTC-subclass coloring *may* shift under the floor. Pete decision: is that comparison re-run against floored subclasses, or is it a fixed June analysis?

## TODO detail / risks
- **SF28** — easiest: hardcoded DOT numbers (3,009/190/1,166 → 1,585/136/888) plus `N_C2<-nrow(profiles_c2)` auto-updates from floored profiles. Also check the "674/1,166 dropped" and "111/111" figures.
- **SF25 / SF27** — producer is an older `code/` script, not the isopair_wrapper pipeline; must confirm it reads floored caches or give each a rebuilt `data_export.R` (as done for SF26).
- **SF29** — its data comes from the **deprecated** legacy Rmd; needs a proper floored producer (do NOT re-bless the legacy Rmd).
