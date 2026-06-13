# §4 prose updates — Google Doc find/replace pairs

**Date:** 2026-06-13 (PM revision, supersedes the AM version)
**Source of new numbers:** canonical Isopair analysis chain at the **ref-AUG-traceable scope** (per Pete's policy: TD2 CDS annotations are unreliable; PTC analyses limit to isoform pairs where reference-AUG tracing can be performed). Same source of truth as Figure 3 Panels D/E/F.

**Why this replaces the AM version:** the previous find/replace doc was computed at the canonical 2,496-pair scope (Stage 2 + coding-coding + gene_id-only re-match) with both TD2-PTC+ and ref-AUG-recovered combined. Pete's clarified policy moves us off that mixed scope onto a clean single scope where ref-AUG-traced ORFs are the canonical PTC determination.

## Numbers in the §4 prose that need to update

| Old prose | New value | Where in §4 |
|---|---|---|
| "3,674 genes" | **3,099 eligible genes** (3,009 yielding NMD pair sets) | Paragraph 1 |
| "44%" PTC rate (NMD) | **83.5%** (ref-AUG-traceable scope) | Paragraph 3 |
| "5%" PTC rate (Control) | **16.3%** (ref-AUG-traceable scope) | Paragraph 3 |
| "15-fold enrichment" | **~5.1-fold** (still highly significant; OR can be recomputed) | Paragraph 3 |
| "56%" NMD+/PTC− | superseded — no longer the right way to slice the data | Paragraph 4 |
| "88%" of NMD+/PTC− contain ref AUG | superseded — the new scope is defined BY ref-AUG-traceability | Paragraph 4 |
| "72%" of those ORFs contain early PTCs | superseded — see above | Paragraph 4 |
| "85%" combined PTC explanation | **83.5%** (now the direct PTC rate among ref-AUG-traceable NMD pairs) | Paragraph 4 |
| "TD2 CDS" / "TD2's CDS caller" | **"GENCODE+TD2 CDS"** or **"initial CDS"** — and the "anti-PTC bias" critique is specifically about TD2 (novel isoforms); GENCODE-annotated CDS isn't biased the same way | Throughout §4 |

## Suggested find/replace pairs (Google Doc-friendly)

Apply in order.

### Pair 1 — paragraph 1 gene count

FIND:
> This resulted in 3,674 genes from which the isoform pairs were drawn.

REPLACE WITH:
> This resulted in 3,099 eligible genes — 3,009 of which yielded gene-matched NMD-vs-Control isoform pair sets.

### Pair 2 — paragraph 3 PTC rate (replace the whole sentence; the new framing is now ref-AUG-traceable)

FIND:
> we observed a strong 15-fold enrichment of PTCs in NMD susceptible isoforms (44% vs 5% PTC rate, p < 0.001)

REPLACE WITH:
> among the 2,289 NMD comparator isoforms where reference-AUG tracing could be performed, 83.5% (1,912) have a premature termination codon (PTC), compared with 16.3% (288/1,763) of Control comparators — a 5.1-fold enrichment (Fisher's exact, p < 10⁻³⁰⁰).

### Pair 3 — paragraph 3 PTC-distance dose-response (renders Figure 3 Panel D)

FIND:
> Looking at all stop codons, we observed a large enrichment of stop codons far upstream of the terminal exon junction in NMD isoforms

REPLACE WITH:
> Looking at the reference-AUG-traced stop codons, we observed a large enrichment of stop codons far upstream of the terminal exon junction in NMD isoforms

(The "looking at all stop codons" framing is no longer accurate because we're now using a specific stop-codon source — the ref-AUG-traced one.)

### Pair 4 — paragraph 4 prose — REWRITE rather than find/replace

The current §4 paragraph 4 describes a sequential rescue logic (TD2 → ref-AUG reveals occult PTCs → 85% combined). Under the new framing the analytical population IS the ref-AUG-traceable subset by design, so this whole paragraph should be restructured. Proposed replacement text:

> Because TD2's CDS calls for novel isoforms are systematically biased against premature stop codons, we limit our main PTC analyses to the **2,289 NMD comparator isoforms where reference-AUG tracing can be performed** (i.e., the reference isoform has an annotated AUG, and that AUG is exonic in the comparator transcript). For each such comparator, the stop codon is taken from the reference-AUG-traced ORF, sidestepping TD2's anti-PTC bias for novel isoforms. Under this scope, **83.5%** of NMD comparator isoforms carry a PTC. The 18% of NMD pairs excluded by the ref-AUG-traceability filter (ref-AUG lost or mapping failed) represent a category we cannot resolve from the available data; their distribution by alternative mechanisms is discussed in §4 paragraph 5.

(Pete to edit/refine.)

### Pair 5 — terminology consistency

FIND (replace-all):
> the predicted CDS

REPLACE WITH:
> the initial GENCODE+TD2 CDS

FIND (replace-all):
> TD2 program

REPLACE WITH:
> TransDecoder2 (TD2)

FIND (replace-all in §4):
> TD2's strong prioritization of ORF length

REPLACE WITH:
> TD2's strong prioritization of ORF length (applies to novel isoforms; GENCODE-annotated CDS for known isoforms is not affected by this bias)

## Sanity check on the new numbers (audit trail)

Computed by `figures/multipanel/figure3_isopair_and_ptc/data_export.R` on 2026-06-13 PM. Population structure:

| Layer | NMD n | Control n | Notes |
|---|---|---|---|
| 1. Stage 2 gene-matched (pop_BC) | 3,009 | 3,009 | (Panels B, C use this) |
| 2. ref-AUG-traceable (pop_traceable) | 2,289 | 1,763 | Categories: effectively_ptc + no_downstream_ejc + truncated_no_ejc (Panel D) |
| 3. ref-AUG PTC+ (pop_ptc_plus) | 1,912 | 288 | Category: effectively_ptc (Panels E, F) |
| PTC rate (NMD): | **83.5%** | | 1,912 / 2,289 |
| PTC rate (Control): | | **16.3%** | 288 / 1,763 |

## Related figure / methodology files (also updated)

- `figures/multipanel/figure3_isopair_and_ptc/data_export.R` — new population logic
- `figures/multipanel/figure3_isopair_and_ptc/panel_e_compute.R` — attribution combined across two sources at the new scope
- `figures/multipanel/figure3_isopair_and_ptc/figure3_panel{A,B,C,D,E,F}_methodology.md` — will be rewritten under Task #22 (in progress)
- The bar-chart supplement (`figure3_supp_ptc_rate_by_cds.{py,pdf,png}`) is now SUPERSEDED — the same story is now in Panel D as the headline; the supplement can be deleted (or kept showing TD2 vs combined for didactic purposes; Pete to decide)
- Upstream Rmd updates (`05_final_report_mashr.Rmd`, `compute_ptc_rates_row()`, `goal2-ptc-mechanisms`) — Task #23, deferred

## Outstanding decisions for Pete to make on §4 prose

1. **Should we keep the "TD2 vs ref-AUG rescue" narrative or remove it?** The new scope makes ref-AUG-traceable the canonical population, so the "TD2 missed PTCs but we found them" rescue narrative is less central. Could be moved to a supplement / methods note.
2. **The 5'UTR / uORF analysis in paragraph 5** — population for that should probably also be the ref-AUG-traceable subset. If so, the subgroup counts there will change.
3. **The 3'UTR-length analysis in paragraph 6** — already on a different population (GENCODE-annotated CDS only), so unaffected.
