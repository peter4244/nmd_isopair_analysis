# Reference-share floor — scientific rationale (internal)

## Problem
Isopair's reference isoform is selected as rank-1 by DMSO mean CPM within the **strict
non-NMD** candidate pool (`nmd_class[[ct]]$non_nmd`; adj.P.Val > 0.30 in all 4 CTs). The gene's
true dominant transcript is often a **"neither"-gap** isoform — not NMD-responsive, but also
failing the strict non-NMD bar in ≥1 CT — which is filtered out of the candidate pool *before*
ranking. Consequently the selected reference is NOT the gene's dominant transcript for ~47% of
pop_BC genes (1,413/3,009), and is frequently a minor (often novel, often <1 CPM) isoform.

The `generatePairsExpression` function is faithful — it correctly ranks whatever pool it is
handed. The displacement is a consequence of the upstream pool restriction in the wrapper, not a
bug in the package. It is not an error, but an unexpected consequence of the filtering.

## Fingerprint in the current manuscript
§4¶1 reports the reference "accounted for a median 31% of parent gene expression" — this *is* the
displacement (median reference share ≈ 30.8%). If the reference were the gene dominant, the
median would be ~64%. The legacy report (`05_final_report_mashr.Rmd`) reported ~70% because it
measured a *different* quantity: `max_cpm/total_cpm` (the true dominant's share), then asserted
"the dominant non-NMD isoform serves as the reference" — an assumption false for ~47% of genes.

## Decision (Pete, 2026-07-10)
**Family A + 25% floor + all-isoform denominator.** Keep the strict-non-NMD selection (we want a
*truly non-NMD* reference — this is the more important property than being the absolute
highest-expressed isoform), and additionally require the selected reference to account for ≥25%
of total gene DMSO expression. Genes whose strict-non-NMD reference falls below 25% are **dropped**
(not re-anchored to a "neither"-gap isoform, which would not be truly non-NMD). The 25% floor at
the median of the GENCODE-reference share distribution removes ~92% of novel/displaced references
while retaining ~70% of GENCODE references, with modest loss to the headline subsets.

Why 25% and not 50%: 50% halves even the GENCODE-referenced headline subsets (n=190→~100,
n=1,166→~640) because the GENCODE-ref median share sits at 50%, risking the secondary Fig 3E/F
event-attribution stats. 25% (n=190→~134, n=1,166→~861) fixes essentially the same displacement
with roughly half the power loss.

Why all-isoform denominator: "% of overall gene expression" reads most naturally as a share of
the gene's total output. Caveat: `expr_mat` is already the 5%-condition-stratified filtered matrix
(`01_prepare_data_mashr.R:129-165`), so the denominator is "all isoforms passing the 5% filter,"
not literally every isoform. Methods must state this precisely.

## Exposure (verified, all-samples DMSO basis)
- pop_BC displaced references: 1,413/3,009 (47.0%); of displaced, 94.7% displaced by a "neither"-gap isoform.
- Novel references: median share 7.9%, 92% dropped at 25% floor.
- GENCODE references: median share 50.7%, ~70% retained at 25% floor.
- Headline subsets have ZERO novel references (n=190 all-GENCODE; n=1,166 ref-AUG-traceable).
