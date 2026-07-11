# Reference-floor number delta (pre → post) — SSOT for prose/verifier/map/manuscript updates

Captured 2026-07-10 from the floored rebuild (02 + 03b --force + 05r + 05k_b) and the
regenerated `data_export.R` outputs. Every downstream prose/assertion/verifier value must move
pre → post per this table. "TBD" = to be read from the rendered Rmd / remaining panel exports.

| Quantity | Pre-floor | Post-floor | Source |
|---|---|---|---|
| pop_BC genes (all_samples) | 3,009 | **1,585** | 02 / 03b gene-matching |
| per-CT C2 AT/DD/FB/MV | 2,710/2,907/2,583/2,756 | 1,445/1,487/1,358/1,487 | 02 |
| n=190 gencode_all3 (each arm) | 190 | **136** | figure3/figure5 export |
| — NMD+/PTC+ | 72 | **54** | figure5 export |
| — NMD+/PTC− | 118 | **82** | figure5 export |
| — Control | 190 | **136** | figure5 export |
| NMD PTC rate | 37.9% | **39.7%** (54/136) | figure3 export |
| Control PTC rate | 2.1% (4/190) | **0.7%** (1/136) | figure3 export |
| PTC fold-enrichment | 18× | **54×** | figure3 export |
| Fisher OR | 28.2 | **88.1** | computed |
| Fisher p | 1.88e-20 | **5.15e-18** | computed |
| n=1,166 refaug (each arm) | 1,166 | **888** | figure5 n1166 export |
| — NMD+/PTC+ | 1,050 | **818** | figure4 export |
| — NMD+/PTC− | 116 | **70** | figure4 export |
| n=1,166 PTC rate | 90.1% (1050/1166) | **92.1%** (818/888) | figure4 export |
| occult-PTC scope | 492 | **380** | SF33 export |
| Kozak broad paired-Wilcoxon p | 1.21e-38 | **2.92e-39** | SF33 export |
| Kozak occult paired-Wilcoxon p | 2.6e-38 | **5.17e-38** | SF34 export |
| SF39 branch-SHAP subclass n (PTC+/PTC−/Ctrl) | 1,016/95/1,107 | **794/60/851** | SF39 |
| SF39 ATG-branch share PTC+ / PTC− | 9.0% / 18.1% | **9.0% / 17.7%** | SF39 |
| **§5 ATG-branch ratio PTC−:PTC+** ("roughly 3×") | 2.21× | **1.98×** (~2×) → soften prose to "~twice" | SF39 |
| **§4 NMD-effect peak vs downstream-EJC count** ("peaked at 4-5") | 4–5 EJCs | **flat/noisy; nominal peak at 6** (bin5 dips to 2.03) → reconsider "peaked at 4-5" wording | SF31 |
| Reference-share median | 31% | **TBD (≥25% by construction; expect ~50%+)** | Rmd §1 render |
| Transcript-length medians (ref/NMD/ctrl) | 2893/3049/2762 | **2958/2991/2808** (SF27); NMD-vs-ref now n.s. (p=0.058) but still "similar"; Control still shorter p<10⁻⁵ | SF27 |
| Median isoforms per pop_BC gene | 7 | **7** (unchanged) | SF25 |
| %ENST reference / %novel NMD | 69.7 / 76.6 | TBD | Rmd render |
| SE prevalence NMD vs Ctrl | 44.2 vs 21.2% | TBD | figure3 Panel C |
| PTC-cause mech split (fs/ifs/3′utr) | 55/33/12% | TBD | figure3 Panel F |
| A5SS attribution enrichment | 13.0 vs 4.5%, p=9e-3 | **TBD — watch (n=54 PTC+)** | figure3 Panel E |
| TD2 == ref AUG | 49.6% (578/1166) | TBD | Rmd §4 |
| TD2 downstream (occult) | 99.0% (487/492) | TBD (of 380) | Rmd §4c |
| Kozak ref>TD2 (occult) | 78.0% (384/492) | TBD (of 380) | SF34 |

**Interpretation:** the floor sharpened the headline PTC enrichment (18×→54×, OR 28→88) by removing
displaced/low-expression references. p only marginally less extreme (10⁻²⁰→10⁻¹⁸) from smaller n.
