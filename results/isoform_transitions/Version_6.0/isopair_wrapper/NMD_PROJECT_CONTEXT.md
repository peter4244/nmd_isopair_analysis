# NMD Isoform Project Context — For Sequence-Based Modeling

## Project Overview

This project studies **Nonsense-Mediated mRNA Decay (NMD)** in human lung cell lines using PacBio IsoSeq long-read RNA-seq. NMD is a post-transcriptional surveillance mechanism that degrades mRNAs containing premature termination codons (PTCs). We inhibited the NMD pathway with SMG1 inhibitor (SMG1i) across 6 lung cell types (AT2, DD, DD_ALI, DO, FB, MV) and compared to DMSO controls.

**Key experimental design:**
- 36 samples: 6 cell types × 2 treatments (DMSO, SMG1i) × 3 donors
- PacBio IsoSeq count matrix: ~110K isoforms × 36 samples
- NMD classification via mashr multi-condition differential isoform expression
- NMD-sensitive: mashr posterior logFC > 0, lfsr < 0.05 (union across 5 main cell types)
- Non-NMD: adj.P.Val > 0.50 (intersection across 5 main cell types, conservative)

## Key Findings from Structural Analysis (completed)

### The CDS Misprediction Problem
- TransDecoder2 (TD2) predicts CDS boundaries for each isoform
- For NMD isoforms, TD2 systematically selects the **wrong reading frame** — choosing a longer alternative ORF that avoids the PTC
- This means TD2's downstream EJC count (the primary PTC indicator) is incorrect for ~40% of NMD isoforms
- The gene's **dominant non-NMD isoform** provides the biologically correct reading frame

### Reference-CDS Approach
- For each gene, the dominant non-NMD isoform (highest DMSO CPM) serves as the CDS anchor
- The reference ATG is traced through each target isoform's transcript
- If exonic: walk codons to the first stop, count downstream EJCs
- This correctly identifies PTCs that TD2 missed
- 69.7% of "PTC-negative" NMD isoforms are effectively PTC+ from the reference reading frame

### Combined Prediction Model (AUC = 0.94)
- 24 features: 12 from TD2 CDS + 12 from reference CDS
- Each set includes: downstream EJC count, 5'UTR ATG density/count/Kozak, uORF counts, ORF coverage, stop density, 3'UTR length
- Elastic net (alpha=0.5), holdout chr 1,3,5,7 (~26%), paralog-free test set
- **Key result**: Reference-CDS PTC alone (1 feature) achieves AUC = 0.896, nearly matching TD2's full 11-feature model (0.917)
- TD2's 5'UTR features compensate for incorrect PTC calls, not independent biology
- The structural model achieves AUC = 0.94 but has a ceiling — sequence-level features may capture what structure misses

### SHAP Decomposition (NMD subpopulations)
- k=4 clusters among holdout NMD isoforms
- Clusters separate by which CDS source drives the prediction
- Some NMD isoforms have no PTC from either CDS source — these may require sequence-based approaches

### What Structure Can't Explain
- ~15-20% of NMD pairs remain unexplained by PTC mechanisms
- These have no downstream EJCs from either TD2 or reference CDS
- Possible mechanisms: long 3'UTR NMD (faux UTR model), uORF-mediated NMD, sequence-specific features
- A CNN model on sequence alone achieved AUC ~0.73 (see below)

## Data Locations and Formats

### Transcript Sequences
- **SQANTI corrected FASTA**: `/proj/regeps/regep00/studies/ExternalCellLines/data/longread/mrna/Randell_Lung_Cells_2025/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta`
  - ~230K transcript sequences
  - Headers are isoform IDs (ENST or ENSG.novel format)
  - Full-length corrected transcript sequences

### Key Data Files (on local machine, may need transfer)
All under `results/isoform_transitions/Version_6.0/isopair_wrapper/data_mashr/`:

- `expression_data.rds` — CPM matrix (isoform × sample), TMM-normalized
- `sample_metadata.rds` — sample_id, ct, donor, treatment
- `nmd_classification.rds` — list: `$all_samples$nmd` and `$all_samples$non_nmd` (isoform ID vectors)
- `cds.rds` — isoform_id, coding_status, cds_start, cds_stop, strand
- `structures.rds` — isoform_id, gene_id, chr, strand, exon_starts, exon_ends
- `gene_map.rds` — isoform_id → gene_id mapping
- `ptc.rds` — isoform_id, n_downstream_ejcs, has_ptc, ptc_distance, orf_length

Under `data_mashr/analysis_cache/`:
- `unified_model.rds` — TD2 feature matrix + elastic net models (~60K isoforms)
- `ref_cds_features_all.rds` — reference-CDS features for all isoforms
- `model_comparison.rds` — combined model with train/test splits, all model variants
- `utr5_features_all.rds` — 5'UTR features from Isopair::scan5UtrFeatures
- `orfik_scan.rds` — comprehensive ORF scan (2.4M ORFs, 61K transcripts)

### Isoform ID Formats
- ENST: `ENST00000665867.2` — keep version (matches count matrix)
- Novel: `ENSG00000196878.16.novel26` — do NOT strip
- Gene IDs: `ENSG00000196878.16` — versioned
- Novel detection: `grepl("\\.novel\\d+$", ...)`

### NMD Classification Numbers
- NMD-sensitive (all_samples): ~9,274 coding isoforms
- Non-NMD (all_samples): ~52,423 coding isoforms
- Total classified coding: ~61,697

## Prior CNN Work

A preliminary CNN model was trained on transcript sequences:
- Input: one-hot encoded transcript sequences (padded/truncated)
- Architecture: 1D CNN with multiple filter sizes
- AUC ~0.73 on sequence alone
- The 0.16 AUC gap between sequence (0.73) and the structural unified model (0.89) suggests explicit ORF-level reasoning is needed
- The combined structural model reaches 0.94 — can sequence close or exceed this?

Relevant files:
- `data_mashr/analysis_cache/cnn_data.tsv` — labeled sequences for CNN training
- `data_mashr/analysis_cache/cnn_nmd_model.pt` — trained PyTorch model
- `data_mashr/analysis_cache/cnn_nmd_results.tsv` — CNN predictions

## Modeling Recommendations for Sequence-Based Approach

### Population
Use the same matched population as the structural model (42,742 isoforms with complete features from both CDS sources) for fair comparison. Same holdout chromosomes (1, 3, 5, 7) and paralog removal.

### Promising Directions
1. **Hybrid model**: Sequence features + structural features. Can sequence capture what PTC-based features miss?
2. **Attention-based models**: Transformer or attention over transcript sequence may learn reading frame and splice junction positions
3. **ORF-aware sequence encoding**: Encode the transcript with CDS boundaries marked (from both TD2 and reference CDS) — gives the model explicit reading frame information
4. **Focus on the hard cases**: The ~15-20% of NMD isoforms unexplained by PTC mechanisms are the most interesting target for sequence models

### Important Considerations
- **Data leakage**: Paralogs with >80% protein identity across train/test chromosomes must be removed from test set
- **Class imbalance**: ~16% NMD prevalence — use AUC, not accuracy
- **Holdout discipline**: chr 1,3,5,7 NEVER used for any training or feature selection
- **Baseline to beat**: AUC = 0.94 (combined structural model), 0.73 (prior CNN)

## User Preferences (from prior sessions)

- Share reasoning before acting; flag uncertainty explicitly
- Verify against actual data before presenting claims
- Never hardcode computed values — use dynamic references
- Reproducible research: trace from source data
- Don't run full dataset without direction
- Plan → review → implement → check workflow
- Channing cluster: linux12h (many jobs, 12h limit), linux01/r8 (no limit, 20 concurrent)
