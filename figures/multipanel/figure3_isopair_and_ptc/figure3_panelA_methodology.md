# Figure 3 Panel A — methodology

**Title:** Isoform pair construction
**Render script:** `figure3_panelA_pair_concept.py`
**Source (R original):** `nmd/results/isoform_transitions/Version_6.0/isopair_wrapper/paper_figures/scripts/make_pair_concept_figure.R`
**Scope:** N/A (conceptual schematic, no data)

## Headline claim

For each gene with ≥1 NMD susceptible isoform, Isopair constructs two paired contrasts that share the same reference isoform: a NMD pair (reference vs. NMD comparator) and a Control pair (reference vs. next-best non-NMD comparator). Splice events distinguishing each comparator from the reference are detected by Isopair and validated through reference reconstruction.

## Source data

None — conceptual schematic. The gene model coordinates are fixed in the script:

- Reference: 6 exons (200–350, 600–800, 1000–1200, 1400–1550, 1800–2000, 2200–2500)
- NMD comparator: 5 exons (exon 4 at 1400–1550 is skipped)
- Control comparator: 6 exons but exon 3 is A3SS-lengthened (940–1200 instead of 1000–1200)

## Population

N/A — schematic does not represent any specific subset of the data.

## Computation

Pure ggplot/matplotlib drawing. No data computation.

## Caveats / limitations

- The depicted skipped-exon and A3SS events are illustrative, chosen to match the most common event types observed in Panel C (skipped exon is the strongest NMD-enriched event; A3SS is also enriched).
- Panel sizing tuned for the standalone view; composite-layout scaling may need title or annotation adjustments.

## Cross-references

- `figures/lib/principles.md` — figure-making principles (two-font, validator-clean, etc.)
- Isopair package source: `~/claude_projects/Isopair/R/pair-generation.R`, `primary-orf.R`

## Composite cross-reference

This panel appears in the Figure 3 composite as the labeled cell — see `figure3_composite_methodology.md` for the layout and per-cell mapping. The composite embeds the panel's pre-rendered PNG (`figure3_panelA.png`); regenerating the composite does NOT re-run this panel. Re-render this panel script first if its data or rendering has changed.
