**SF39 | Branch-level Shapley attributions of the deep-learning model, stratified by isoform PTC subclass.** Shapley attributions of the model's NMD prediction across its three input branches — the start-codon window, the stop-codon window, and a set of structural features — measured at the branch-aggregation point of the trained network (the same attribution method as in Figure 5 Panel C). Attributions are stratified across three isoform subclasses defined in Figure 4 Panels C and D: NMD susceptible with a premature termination codon (NMD+/PTC+, n = 1,016), NMD susceptible without a premature termination codon (NMD+/PTC−, n = 95), and matched non-NMD Controls (n = 1,107). Attributions are computed for every isoform in the subgroup.

Each panel matches the style of Figure 5 Panel C: bars are ordered structural / stop / start; the value inside each bar is the mean absolute Shapley attribution at that branch; the percentage above each bar is that branch's share of the within-subclass total attribution. The vertical scale is identical across the three panels so heights can be compared directly.

(A) NMD susceptible with a premature termination codon (NMD+/PTC+). Structural features dominate (62 % of total branch attribution); the stop-codon window contributes 29 %; the start-codon window contributes 9 %. This pattern mirrors the pooled attribution shown in Figure 5 Panel C — when a premature termination codon is present, the structural and stop-codon features carry the prediction.

(B) NMD susceptible without a premature termination codon (NMD+/PTC−). The structural share drops to 53 % (no premature termination codon means no excess downstream-junction signal); the stop share is essentially unchanged at 29 %; the start-codon-window share doubles to 18 %. The model relies more heavily on start-codon context to call NMD in these isoforms.

(C) Non-NMD Controls. The proportional profile resembles NMD+/PTC+ (62 % / 27 % / 11 %) but at lower absolute magnitude across all three branches.

The start-codon-window share in NMD+/PTC− is ~2-fold higher than in NMD+/PTC+ (18 % vs 9 %), quantifying the shift toward start-codon information when a premature termination codon is not present.

Abbreviations. NMD, nonsense-mediated decay; PTC, premature termination codon.
