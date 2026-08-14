# Figure 3D — should it use the reference AUG-traceable universe?

**For Yul.** Drafted 2026-08-11 at Pete's request. Nothing here is committed to the manuscript and
nothing in the shipped panel has been changed. This is a description of a finding, the evidence
behind it, the parts that are not settled, and what you would need to decide.

---

## The short version

Figure 3D currently draws 130 NMD and 130 Control pairs. Redrawn on the reference AUG-traceable
universe — about 819 pairs per arm — the same measurement separates far more strongly, and the
picture matches what the text already claims.

| | current panel | reference AUG-traceable |
|---|---:|---:|
| pairs per arm | 130 | 819 |
| NMD median distance | −57 nt | **+421 nt** |
| Control median distance | −151 nt | −105 nt |
| NMD with stop ≥ 50 nt upstream of the last EJC | 36.9% | **92.3%** |
| Control, same | 1.5% | 16.0% |
| fold enrichment | 24× | 5.8× |

`panelD_819_mockup.png` in this directory is the redrawn version. `figure3_panelD_stop_codon_distance.png`
one level up is the current one. Same quantity, same sign convention, same colors, so they can be
compared directly.

**The trade is real and it runs in both directions.** The separation becomes much more visible and
rests on roughly six times more data, but the fold enrichment falls from 24× to 5.8×, because the
Control rate rises from 1.5% to 16%. Whether that reads as a stronger result depends on whether the
sentence leads with the percentages or with the fold.

---

## Why the two scopes differ

Both panels measure the same thing: the distance from a comparator isoform's stop codon to its last
exon-exon junction. Positive means the stop sits upstream of that junction — a premature termination
codon under the 50-nt rule.

What differs is which isoforms are eligible, and that turns on how the coding sequence is called.

**The current panel avoids the problem by shrinking.** It uses only pairs where all three isoforms
(reference, NMD comparator, Control comparator) are GENCODE-annotated, so each one carries its own
annotated CDS and no CDS prediction is needed. That is a clean design and it is why the panel is
small: 130 pairs.

**The expanded scope solves it instead by tracing.** It requires only that the reference is
annotated, then projects the reference start codon onto each comparator's spliced sequence to define
the comparator's reading frame. This is the "reference AUG tracing" the manuscript already describes
at line 74 of the Methods.

The reason this matters is that the alternative — using TD2's predicted CDS for novel isoforms —
is biased, and the size of the bias is measurable:

> Among NMD comparators, TD2 called **no** premature stop for 756 isoforms. Reference AUG tracing
> finds a premature stop in **556 of them**. The opposite error — TD2 calling a premature stop where
> tracing finds none — happens **30** times.

An 18-to-1 asymmetry, in the direction of hiding exactly the isoforms the figure is about. This is
the same bias the paper reports in SF34–35 and quantifies as 355 occult-PTC isoforms; the 556 here
is a broader count over a different scope, not a contradiction.

---

## What is NOT settled

**The cohort is 819, and the published number is 833.** The manuscript reports the reference
AUG-traceable set as n=833 (768 PTC+ and 65 PTC−). Reconstructing that cohort from the analysis
cache gives **819** (756 PTC+ and 63 PTC−). The PTC rate agrees to a tenth of a point — 92.3% against
the published 92.2% — but 14 pairs are unaccounted for.

I stopped there rather than adding filters until the number matched, because that procedure
manufactures agreement rather than finding it. The likely explanation is a vintage difference: the
`ref_atg_analysis.rds` used here is dated 2026-07-11 and the published figure may come from a later
build. **Before any of this reaches the manuscript, the exact cohort should come from whatever code
built the published 833**, not from the reconstruction described below.

Worth knowing: there are at least three reference AUG-traceable populations in circulation, and they
are not interchangeable.

| population | size | where it is defined |
|---|---:|---|
| `pop_traceable` | 2,289 NMD / 1,763 Control | `figure4_ptcneg_and_model/RATIONALE.md` §2 |
| §4 Section C | 1,050 PTC+ / 113 PTC− / 1,166 Control | `figure4_ptcneg_and_model/data_export.R` line 9 |
| the published 833 | 768 PTC+ / 65 PTC− | Supplemental Methods, "reference AUG-traceable scope" |
| this reconstruction | 756 PTC+ / 63 PTC− | `build819.R` in this directory |

**The Control tail is real.** The expanded Control curve has visible density beyond the threshold
that the 130-pair version does not. That is not an artifact — it is what a 16% rate looks like. It
should be explained rather than cropped: these are non-NMD comparators whose traced stop does sit
upstream of a junction, and the interesting question is why they are not degraded.

---

## How to reproduce what is in this directory

Everything runs locally; no cluster access needed. Files are in this directory.

1. **`ejc.R`** — computes the last exon-exon junction position in transcript coordinates for every
   isoform, from `structures.rds`. It validates itself against the 260 rows whose last-EJC position
   is already shipped in the current panel's data file, and reports **260/260 exact agreement**. If
   that check ever fails, nothing downstream should be trusted.

2. **`build819.R`** — builds the cohort and computes the distances. The cohort filter is: drop
   pairs where reference AUG tracing failed (`mapping_failed`, `ref_atg_lost`), keep pairs whose
   reference is GENCODE-annotated, then re-intersect the NMD and Control arms on (gene, reference)
   so both arms cover the same genes. Writes `panelD_819.tsv`.

3. **`panelD_819.py`** — draws the mockup. Mirrors the shipped panel's styling deliberately.

```bash
cd /Users/petecastaldi/claude_projects/nmd/figures/multipanel/figure3_isopair_and_ptc/scope_investigation
Rscript ejc.R          # expect: 260/260 exact agreement
Rscript build819.R     # expect: 819 per arm, NMD 92.3% / Control 16.0%
python3 panelD_819.py
```

**Sign convention, since it is easy to invert:** `distance = last_ejc_tx_pos − stop_tx_pos`, so
*positive means the stop codon is upstream of the last EJC*. This is the convention the shipped
panel uses; it was confirmed against its data file rather than assumed.

---

## What would change in the manuscript

The sentence carrying Figure 3D is at line 72 of `NMD manuscript 2026.7.17.md`:

> "…we observed a 24-fold enrichment of PTCs in NMD susceptible isoforms (36.4% vs 1.5% PTC rate,
> p < 10⁻¹³) and a strong enrichment of stop codons upstream of the terminal EJC (**Fig 3D**)."

Every number in it is scope-specific. On the expanded universe the percentages become roughly
92% and 16%, and the fold becomes about 6×. The framing sentence just before it — "We first limited
our analysis to the GENCODE reference set" — is what makes the current panel correct as it stands,
so the panel is not wrong; it is answering a deliberately conservative question.

Note that line 74 **already** reports the expanded result: "we identified PTCs in 92% of the NMD+
isoforms (768/833)". So the paper contains both scopes already. The question is only which one
Figure 3D should show.

Three ways this could go, and this is your call and Pete's rather than something to settle from the
data:

- **Leave 3D alone.** The conservative panel is defensible and the strong result is already in the
  text two paragraphs later.
- **Replace 3D with the expanded version.** The figure then matches the claim that the section
  builds toward, on six times more data — at the cost of the 24-fold headline.
- **Show both.** The GENCODE-restricted panel as the conservative result and the expanded one beside
  it. Costs space in Figure 3, and makes the scope distinction visible rather than buried in Methods.

---

## Questions worth answering before deciding

1. **Where do the missing 14 pairs go?** Which script built the published 833, and does re-running
   it today still give 833?
2. **Does the p-value survive?** The 24-fold claim carries p < 10⁻¹³ on 130 pairs. The expanded
   comparison has not been tested; with 819 per arm it will be extreme, but it should be computed
   rather than assumed.
3. **What are the 16% of Controls with an upstream stop?** If they are concentrated in one
   structural category, that is worth a sentence — and possibly a supplemental panel.
4. **Should the fold or the percentages lead?** 24× on a 1.5% baseline and 5.8× on a 16% baseline
   are both true; they answer different questions.
5. **Does Figure 4 need to stay consistent?** Figure 4 already works at a reference AUG-traceable
   scope. If 3D moves there too, the two figures become directly comparable, which may be an
   argument in itself.

---

## Provenance

Every number above was measured from the analysis cache at
`results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/`, not carried over from any
document. The comparison figures are `panelD_819_mockup.png` (this directory) and
`figure3_panelD_stop_codon_distance.png` (one level up). The mockup has **not** been through the
figure validator gate and is not a deliverable; it exists to show what the scope change looks like.
