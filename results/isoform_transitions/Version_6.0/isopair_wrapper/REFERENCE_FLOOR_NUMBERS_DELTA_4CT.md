# 4-CT re-scope number delta — SSOT for R3/R4/R5/R6

Captured 2026-07-11 from the 4-CT rebuild (decision **b**: `01` restricts to AT/DD/FB/MV at
source; `02`+floor / `03b --force` / `05r` / `05k` / `05k_b` re-run on the 4-CT `expr_mat`; then
`code/nmd_predictor_comparison/` 01→04 re-run). **Supersedes the 6-CT `REFERENCE_FLOOR_NUMBERS_DELTA.md`
entirely.** Three columns: pre-floor (original 6-CT), 6-CT-floor (2026-07-10, now void), **4-CT (current)**.

| Quantity | pre-floor (6-CT) | 6-CT floor (void) | **4-CT (current)** | Source |
|---|---|---|---|---|
| Filtered isoform universe (`expr_mat`) | 123,790 | 123,790 | **95,623** | 01 |
| samples | 36 | 36 | **26** (13 DMSO + 13 Smg1i) | 01 |
| pop_BC genes (all_samples C2) | 3,009 | 1,585 | **1,548** | 02 |
| per-CT C2 AT/DD/FB/MV | 2,710/2,907/2,583/2,756 | 1,445/1,487/1,358/1,487 | **1,432/1,476/1,340/1,466** | 02 |
| C4 floor drop (all_samples) | — | — | **8,134 → 5,001 (dropped 3,133)** | 02 |
| n=190 gencode_all3 (each arm) | 190 | 136 | **130** | figure5 export |
| n=1,166 ref-AUG (each arm) | 1,166 | 888 | **819** | figure5 n1166 export |
| occult-PTC scope | 492 | 380 | **348** | SF33 export |
| ref_atg NMD C2 pairs | — | — | **1,410** | 05r |
| ref_atg Control C4 pairs | — | — | **1,365** | 05r |
| effectively_ptc (ref-AUG) | — | — | **1,180 (83.7%)** — orig_ptc T=624 / F=556 | 05r |
| Total PTC-mediated | — | — | **1,210 / 1,410 (85.8%)** | 05r |

## Predictor comparison (SF42) — 4-CT, DATESTAMP 2026.7.11
| Quantity | old (2026.6.20, pre-floor) | **4-CT (2026.7.11)** | Source |
|---|---|---|---|
| Cohort rows | 2,332 (1,166+1,166) | **1,638 (819 NMD + 819 Control)** | 01 |
| — NMD+/PTC+ / NMD+/PTC- | 1,050 / 116 | **756 / 63** | 01 |
| our_model_prob coverage | ~95% (114 dropped) | **1,570 / 1,638 = 95.85% (68 dropped: 21/9/38)** | 01 (all.x merge) |
| head-to-head (all 3 score) | — | **1,570** | 04 |
| head-to-head:test | 561? | **415** (195 PTC+ / 16 PTC- / 204 Control) | 04 |
| pooled Spearman our / ndB / nmdep | — | **0.7700 / 0.8042 / 0.8041** | 04 |
| test Spearman our / ndB / nmdep | — | **0.7666 / 0.7906 / 0.7930** | 04 |

**SF42 legend n (frozen prose to update):** old 561/255/30/276 → **new test-split 415 = 195/16/204**
(verify exact scope against `figure_s_model_comparison.py` at R3).

## Figure 3 (headline PTC enrichment) — 4-CT
| Quantity | pre-floor | 6-CT floor (void) | **4-CT** |
|---|---|---|---|
| gencode_all3 NMD PTC+ / PTC- / Control | 72/118/190 | 54/82/136 | **48 / 82 / 130** |
| NMD PTC rate | 37.9% | 39.7% (54/136) | **36.9% (48/130)** |
| Control PTC rate | 2.1% (4/190) | 0.7% (1/136) | **1.5% (2/130)** |
| PTC fold-enrichment | 18× | 54× | **24.0×** |
| Fisher OR | 28.2 | 88.1 | **37.06** (95% CI 9.3–323) |
| Fisher p | 1.88e-20 | 5.15e-18 | **1.58e-14** |
| Panel E: PTC-attributed events (NMD) / Control baseline | — | — | **48 / 292** |

## TBD — fill in during R3 as each figure/Rmd regenerates
n=1,166 PTC rate (of 819), Kozak paired-Wilcoxon p (broad + occult of 348), SF39 branch-SHAP subclass n,
transcript-length medians (SF27), reference-share median, %ENST reference / %novel NMD, SE prevalence,
PTC-cause mech split, A5SS attribution, TD2==ref-AUG %.
