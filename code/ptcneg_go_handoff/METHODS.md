# METHODS — NMD+/PTC- (reference-AUG) gene-set hand-off for GO enrichment

Analysis: `ptcneg_go_handoff.Rmd` (this directory). Produces the NMD+/PTC- target
gene set and a matched background universe for GO over-representation (topGO).

## 1. Scope and rationale

We identify **NMD-responsive isoforms that lack a premature termination codon
(PTC)** and the genes that carry them, then export them with a matched gene
universe for GO enrichment. Two deliberate choices:

- **PTC is determined per isoform from the reference AUG, not from TD2/SQANTI.**
  TransDecoder/SQANTI predicted CDS is biased against PTC-containing ORFs
  (CLAUDE.md §"SQANTI/TD2 CDS bias"); its `predicted_NMD` field is therefore not
  used. PTC status comes from the reference-AUG projection (§3).
- **Per-isoform, not pairwise.** Unlike the Isopair NMD-vs-control pair analysis,
  PTC here is a property of a single isoform relative to its gene's reference AUG.

## 2. Inputs

| Input | Path | Role |
|---|---|---|
| Reference-AUG ORF/PTC cache | `nmd_orf_model_v5_4ct/tmp/ref_orf_ptc_cache_with_nmd_2026.5.12.rds` (`$trace_aug`) | PTC determination only |
| Isoform mashr DIE (4 CT) | `isocall_dge/mashr/nmd_mashr_die_{at,dd,fb,mv}_2026.3.10.csv` | NMD status (canonical) |
| Gene annotation | `org.Hs.eg.db` (ENSEMBL → SYMBOL, ENTREZID) | identifier mapping |

The cache is built by `nmd_orf_model_v5_4ct/tmp/prep_dominant_ref_orf_ptc_2026.5.12.R`
(reference-AUG trace) and `…/augment_nmd_die_2026.5.12.R` (DIE merge); the trace
engine is `Isopair::traceReferenceAtg()` (`Isopair/R/ptc-attribution.R`). `$trace_aug`
has exactly one row per isoform (key `comparator_isoform_id`).

## 3. Reference-AUG PTC determination

Implemented in `Isopair::traceReferenceAtg()`; invoked with `ejc_threshold = 50`
and `resolve_alt_start = FALSE`.

1. **Reference selection.** For each gene the **reference isoform** is its
   dominant DMSO-expressed isoform, scored as the mean across the 4 cell types
   (AT, DD, FB, MV) of the per-CT mean DMSO CPM (equal weight per CT), restricted
   to isoforms present in GENCODE v49 with an annotated CDS. The reference AUG is
   that isoform's annotated CDS start codon (strand-aware genomic coordinate).
   Genes whose dominant isoform is not a GENCODE-CDS transcript are excluded.
2. **AUG compatibility** (`ref_atg_exonic_in_comp`). TRUE iff all three bases of
   the reference ATG codon fall within exonic sequence of the target isoform
   (`.isCodonExonic`). TRUE = the reference start is retained/traceable; FALSE =
   at least one base is intronic/absent (`category = ref_atg_lost`).
3. **ORF trace.** The reference AUG genomic coordinate is mapped into target
   transcript coordinates, the base is confirmed to be `ATG`, and the sequence is
   walked codon-by-codon to the first in-frame stop (`comp_stop_tx_pos`). This
   uses the isoform's own spliced sequence only — **no TD2/SQANTI CDS**.
4. **EJC count** (`n_downstream_ejc`). Number of exon-exon junctions located
   > 50 nt 3′ of the stop codon (canonical 50-nt rule).
5. **Category.** Among AUG-compatible, successfully traced isoforms:
   - `effectively_ptc` — `n_downstream_ejc > 0` ⇒ **PTC+** (predicted NMD substrate).
   - `truncated_no_ejc` — `n_downstream_ejc == 0` and ORF shorter than the reference ⇒ **PTC-**.
   - `no_downstream_ejc` — `n_downstream_ejc == 0`, ORF not shorter ⇒ **PTC-**.
   - `mapping_failed` — AUG exonic but trace unresolved (ORF/EJC = NA); excluded.
   - `ref_atg_lost` — AUG not exonic; excluded by the compatibility filter.

**PTC- (no PTC)** is defined as `n_downstream_ejc == 0`
(`truncated_no_ejc` ∪ `no_downstream_ejc`). `effectively_ptc` is PTC+ and excluded.

## 4. NMD status (canonical, direction-aware)

NMD-responsiveness is taken from the per-cell-type mashr DIE CSVs, using the
project-canonical flag `nmd_responsive == TRUE` ⇔ **mashr `lfsr < 0.05` AND
posterior mean > 0** (CLAUDE.md §"Data file conventions"): isoforms significantly
**up-regulated** under SMG1 inhibition, the direction expected of NMD substrates.
An isoform is **NMD+** if `nmd_responsive == TRUE` in **≥ 1** of AT, DD, FB, MV
(`n_ct_nmd ≥ 1`).

The cache's own `nmd_responsive_any` column is **not** used: it is defined as
`adj.P.Val < 0.05` in ≥ 1 cell type with **no direction filter**
(`augment_nmd_die_2026.5.12.R`), which admits down-regulated isoforms (63% of the
would-be target) and is not an NMD-substrate criterion. Substituting it inflates
the target ~2.7× and contaminates it with wrong-direction genes; this report
avoids that by sourcing NMD status from the canonical mashr flag.

## 5. Target and universe

- **Target (NMD+/PTC-)** = AUG-compatible ∧ PTC- ∧ NMD+. Result: 467 isoforms /
  362 genes (`no_downstream_ejc` 248 iso / `truncated_no_ejc` 219 iso).
- **Background universe** = genes with ≥ 1 reference-AUG-compatible isoform
  (10,978 genes) — the population eligible to be called NMD+/PTC-, so ORA is not
  diluted by genes that could never enter the target. Independent of NMD status.
- Gene IDs are Ensembl (version stripped). Symbols and Entrez IDs from
  `org.Hs.eg.db` (≈99.8% coverage); Entrez is the topGO keytype.

## 6. Deliverables (this directory, date-stamped)

| File | Contents |
|---|---|
| `nmd_ptcneg_refAUG_target_genes_*.csv` | 362 target genes; IDs + supporting-isoform detail (GO target) |
| `nmd_ptcneg_refAUG_background_genes_*.csv` | 10,978 universe genes; IDs + `is_target` flag (topGO `allGenes`) |
| `nmd_ptcneg_refAUG_target_symbols_*.txt` | target HGNC symbols (Enrichr / g:Profiler) |
| `nmd_ptcneg_refAUG_target_isoforms_*.csv` | 467 target isoforms; per-isoform ref-AUG trace |

A ready-to-run topGO snippet is in `ptcneg_go_handoff.Rmd` §5.

## 7. Reproducibility notes & limitations

- The cache and its generating scripts live under a **gitignored `tmp/`**, and
  `Isopair` is a separate repo — a fresh checkout of the `nmd` repo alone cannot
  regenerate the cache. Provenance paths are recorded above; regenerate from the
  `NMD_orf_model_v5_4ct` repo when needed.
- "Reference" is the **dominant expressed** GENCODE-CDS isoform (4-CT equal-weight
  mean), not the MANE/canonical isoform; a gene dominant in one CT but absent in
  others can still anchor the reference AUG.
- The PTC- fraction reported among NMD+ AUG-compatible isoforms (6.3%) is a
  per-isoform transcriptome-wide quantity and is **not** the manuscript's "~15%",
  which is the PTC-negative share of gene-matched NMD *pairs* in Isopair — a
  different population and denominator.

## 8. Verification status

Subjected to the 5-step report verification (factual accuracy, result
correctness, documentation accuracy, reproducibility, METHODS). Step 2 caught an
NMD-direction defect in an earlier version (use of the cache's direction-agnostic
flag); corrected here by adopting the canonical mashr definition. The in-report
`verify` chunk additionally asserts internal set consistency at render time.
