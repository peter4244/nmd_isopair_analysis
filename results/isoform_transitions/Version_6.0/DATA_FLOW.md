# Data Flow Documentation - Version 6.0

**Purpose:** Track the flow of data through the isoform transitions analysis pipeline, from raw inputs to final analysis-ready datasets.

**Date:** 2026-02-20

---

## **PRIMARY DATA SOURCES**

### **1. Isoform Expression Data**

**File:** `/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds`

**Type:** DGEList object (edgeR)

**Contents:**
- **Counts matrix:** 1,746,281 isoforms × 38 samples
- **Sample metadata:** Cell types (at2, dd, ddali, do, fb, mv), treatments (DMSO, Smg1i)
- **Isoform metadata:**
  - `txid`: Transcript IDs (ENST*, PB.*)
  - `gene_id_ens115_sqanti`: Gene IDs (ENSG*, novelGene_*, fusions)
  - `biotype`: Transcript biotype annotations

**Design:** ~treatment + ct + ct*treatment

**Gene Categories:**
- 87,079 GENCODE genes (ENSG[0-9]+)
- 141,069 PacBio Novel genes (novelGene_*)
- 3,717 Fusion/Chimeric genes (ENSG_ENSG)
- **Total:** 231,865 gene IDs

### **2. GENCODE Annotations (GTF/GFF3)**

**File:** `/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz`

**Format:** GFF3 (compressed)

**Contents:**
- **Exon coordinates** for GENCODE transcripts (ENST*)
- **CDS annotations** for coding transcripts
- Version: GENCODE v49

**Used for:**
- Extracting exon structures for ENST isoforms
- CDS start/stop coordinates
- Transcript biotype information

### **3. SQANTI/PacBio Annotations (GFF3)**

**File:** `/Users/petecastaldi/claude_projects/nmd/reference_files/sqanti3_corrected.cds.gff3`

**Format:** GFF3

**Contents:**
- **Exon coordinates** for PacBio isoforms (PB.*)
- **CDS predictions** (ORF annotations) for novel isoforms
- SQANTI3 structural classifications

**Used for:**
- Extracting exon structures for PB.* isoforms
- ORF predictions for novel transcripts
- Novel isoform structural information

---

## **PIPELINE ARCHITECTURE**

Scripts are organized into `scripts/core/` (generic, reusable) and `scripts/nmd/` (NMD study-specific).

### **Phase 1: Data Preparation (nmd/01, core/02-05)**

Builds analysis-ready data structures from raw inputs (unfiltered).

### **Phase 2: Filtering (nmd/06)**

Applies expression filtering (filterByExpr) and gene category filtering to match DGE analysis.

### **Phase 3: Classification and Pair Generation (nmd/07)**

Classifies isoforms as NMD-sensitive or non-NMD, generates comparison pairs (C1-C4) across 7 sample sets, deduplicates.

### **Phase 4: Event Detection (core/08) + Reconstruction Validation (core/09)**

Comprehensive splicing event detection on filtered data, followed by reconstruction validation.

### **Phase 5: Analysis (nmd/09-14)**

Performs statistical analysis on splicing choice profiles.

---

## **DATA FLOW: PHASE 1 (Preparation)**

### **Script nmd/01: Prepare Expression Data (Memory Optimized)**

**Input:**
- `dge_isoform_nofilter_2026.2.7.rds` (DGEList)

**Process:**
1. Load DGEList (1,746,281 isoforms)
2. Calculate CPM using edgeR normalization
3. **Filter to "major isoforms"** (≥5% expression in any sample)
   - Process gene-by-gene to avoid memory issues
   - For each gene, calculate: `proportion = isoform_CPM / gene_total_CPM`
   - Keep isoforms with `max_proportion >= 0.05` across all samples
4. Identify dominant isoform per gene (highest mean expression)

**Output:**
- `data/expression_data.rds` - CPM matrix in tidy format (502,994 isoforms)
- `data/dominant_isoforms.rds` - Dominant isoform per gene (192,156 genes)
- `data/sample_metadata.rds` - Sample annotations (38 samples)
- `data/filtering_stats.rds` - Per-gene filtering statistics

**Key Decisions:**
- **5% threshold:** Very permissive - keeps isoforms ≥5% in ANY sample
- **Gene-by-gene processing:** Memory efficient for large datasets
- **Major isoforms:** 502,994 / 1,746,281 = 28.8% kept

**Notes:**
- This is an isoform-level filter, NOT a gene-level expression filter
- All 231,865 genes retained in filtering_stats (even if 0 major isoforms)
- Focuses analysis on well-expressed isoforms

---

### **Script core/02: Extract Isoform Structures**

**Input:**
- `data/expression_data.rds` (502,994 isoforms)
- `gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz`
- `sqanti3_corrected.cds.gff3`

**Process:**
1. Extract list of isoform IDs from expression_data
2. Parse GENCODE GFF3:
   - Filter to isoforms in expression_data
   - Extract exon features (transcript_id, exon_number, start, end, strand)
   - Handle version numbers (match versionless IDs)
3. Parse SQANTI/PacBio GFF3:
   - Extract exon features for PB.* isoforms
   - Same structure as GENCODE
4. Combine into unified structure table

**Output:**
- `data/isoform_structures.rds` (15MB)
  - Columns: `isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end`
  - 150,781 GENCODE isoforms (818,103 exons)
  - 352,213 PacBio isoforms (1,894,853 exons)
  - **Total:** 502,994 isoform structures (100% coverage)

**Key Decisions:**
- Strip version numbers for matching (ENST00000000233.10 → ENST00000000233)
- Keep original versioned IDs in output for traceability
- Use rtracklayer for robust GFF parsing

**Notes:**
- Mean 5.4 exons/isoform
- Range: 1-361 exons per isoform

---

### **Script core/03: Build Union Exons**

**Input:**
- `data/isoform_structures.rds`

**Process:**
1. For each gene (89,208 genes with structures):
   - Extract all exons from all isoforms
   - Collect all unique exon boundaries and classify as START, END, or BOTH
   - Split into atomic segments using boundary-type-aware logic:
     - START boundary: split before → [..., b-1] then [b, ...]
     - END boundary: split after → [..., b] then [b+1, ...]
     - BOTH boundary: three segments → [..., b-1], [b, b], [b+1, ...]
   - Filter to segments covered by at least one exon
   - Number segments sequentially per gene
2. Create isoform-to-union-exon mapping:
   - Uses containment check: atomic UE fully within isoform exon (UE_start >= exon_start AND UE_end <= exon_end)

**Output:**
- `data/union_exons.rds`
  - Atomic union exons (more fine-grained than merged approach) from **89,208 genes**
  - Columns: `gene_id, union_exon_id, exon_index, seqnames, start, end, strand`
  - Union exon models for 89,208 genes
- `data/isoform_union_mapping.rds`
  - Columns: `isoform_id, gene_id, union_exon_id, present`
  - Isoform-union exon presence/absence
  - Binary mapping of which union exons are present in each isoform

**Key Concept:**
- **Atomic union exons:** Non-overlapping segments split at every unique exon boundary
- Fundamental property: every isoform exon decomposes into a set of complete atomic UEs via strict containment
- Used by event detection pipeline for isoform comparison; reconstruction uses direct event coordinates

**Notes:**
- 0 overlapping union exons (validated in tests)
- Processed in batches (1000 genes) with progress tracking

---

### **Script core/04: Extract CDS Annotations**

**Input:**
- `gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz`
- `sqanti3_corrected.cds.gff3`
- `data/isoform_structures.rds`

**Process:**
1. Parse GENCODE GFF3:
   - Extract CDS features for each ENST transcript
   - Calculate: `cds_start = min(CDS_start)`, `cds_stop = max(CDS_stop)`
   - Mark as coding if CDS present
2. Parse SQANTI GFF3:
   - Extract CDS predictions for PB.* isoforms
   - Same structure as GENCODE
3. Combine and classify:
   - `coding`: Has CDS coordinates
   - `non_coding`: No CDS or explicitly non-coding
   - `unknown`: Ambiguous cases
4. Calculate ORF length for coding isoforms

**Output:**
- `data/isoform_cds_metadata.rds` (2.8MB)
  - **235,022 coding isoforms** (46.7%) with CDS coordinates
  - **267,972 unknown** (53.3%) - no CDS in GFF (likely non-coding)
  - Mean ORF length: 1,256 bp (median: 996 bp, range: 3-107,820 bp)
  - Columns: `isoform_id, coding_status, cds_start, cds_stop, orf_length`
  - CDS annotations for all 502,994 isoforms

**Coding Status Breakdown:**
- Coding: 235,022 (46.7%)
  - **GENCODE:** 54,952 (36.4% of 150,781 GENCODE isoforms)
    - Have CDS start/stop coordinates from GENCODE GFF3
  - **PacBio:** 180,070 (51.1% of 352,213 PacBio isoforms)
    - Have CDS coordinates from SQANTI3 ORF predictions
    - SQANTI3 lacks start_codon/stop_codon features, inferred using strand-aware logic
  - ORF length calculated from sum of CDS feature lengths
- Non-coding: 0 (0.0%)
  - No explicit non-coding annotations found
- Unknown: 267,972 (53.3%)
  - No CDS features in GFF
  - Includes: true non-coding RNAs, novel isoforms without ORF predictions
  - Cannot distinguish UTR vs CDS for these isoforms

**Key Decisions:**
- GENCODE: Use annotated CDS features from GFF3
- PacBio: Use SQANTI3 ORF predictions from sqanti3_corrected.cds.gff3 (GTF format, 2.0 GB)
- Unknown status when no CDS features found
- **ID handling:**
  - GENCODE: Strip version numbers for matching (ENST00000415666.1 → ENST00000415666)
  - PacBio: Keep full ID intact - format is `PB.gene.isoform` (e.g., PB.2144.377 where .377 is isoform ID, NOT version)

**Coordinate System (Strand-Aware):**
- `cds_start` = **minimum genomic coordinate** of CDS (works for both strands)
- `cds_stop` = **maximum genomic coordinate** of CDS (works for both strands)
- This allows strand-independent range checking in downstream scripts
- Rationale: Genomic coordinates always increase left→right, but on minus strand, biological 5'→3' runs right→left
- Example (minus strand):
  - Start codon at genomic position 40,856,538 (biological 5' end)
  - Stop codon at genomic position 40,826,538 (biological 3' end)
  - Stored as: cds_start=40,826,538 (min), cds_stop=40,856,540 (max)
  - Allows simple range check: `exon_start >= cds_start && exon_end <= cds_stop` → CDS

**Notes:**
- Coordinates in genomic space (not transcript space)
- Required for distinguishing UTR vs CDS changes across both strands
- 36% of GENCODE isoforms are coding, 51% of PacBio isoforms are coding
- ORF length calculated from sum of CDS feature lengths (per-exon CDS coordinates, strand-independent)

---

### **Script core/05: Annotate Region Types**

**Input:**
- `data/union_exons.rds`
- `data/isoform_union_mapping.rds`
- `data/isoform_cds_metadata.rds`

**Process:**
1. For each (isoform, union_exon) pair where isoform contains union_exon:
2. Lookup isoform coding status
3. If non-coding or unknown: `region_type = "non_coding"` or `"unknown"`
4. If coding, apply region type logic:
   ```
   if (exon_end < cds_start):
       region_type = "5'UTR"
   elif (exon_start > cds_stop):
       region_type = "3'UTR"
   elif (exon contains BOTH cds_start AND cds_stop):
       region_type = "contains_orf_start_stop"
   elif (exon contains cds_start):
       region_type = "contains_orf_start"
   elif (exon contains cds_stop):
       region_type = "contains_orf_stop"
   elif (exon_start >= cds_start AND exon_end <= cds_stop):
       region_type = "CDS"
   else:
       region_type = "unknown"
   ```

**Output:**
- `data/isoform_union_exons_annotated.rds` (20MB)
  - **2,712,956 (isoform, union_exon) mappings** annotated
  - Columns: `isoform_id, gene_id, union_exon_id, exon_index, seqnames, start, end, strand, present, region_type`
  - Region type annotations for all isoform-union exon pairs

**Region Type Distribution (across all mappings):**
- **CDS: 1,350,947 (49.8%)** - correctly identified coding exons across both strands
- non_coding: 564,548 (20.8%) - isoforms without CDS annotations
- contains_orf_start: 235,954 (8.7%) - exons spanning start codon boundary
- contains_orf_stop: 192,992 (7.1%) - exons spanning stop codon boundary
- 5'UTR: 185,592 (6.8%) - before CDS
- 3'UTR: 182,828 (6.7%) - after CDS

**Mappings Per Union Exon:**
- Mean: 4.7 isoforms per exon (median: 2)
- Range: 1 to 502 isoforms
- 42.3% are isoform-specific (used by only 1 isoform)
- 34.2% shared by 2-5 isoforms
- 23.5% shared by 6+ isoforms (constitutive/highly shared)

**Region Type Combinations Per Exon:**
- **73.3% of exons have single region type** (419,231 exons - consistent across all isoforms)
  - [non_coding]: 202,910 exons (35.5%) - only in non-coding isoforms
  - **[CDS]: 128,269 exons (22.4%)** - constitutive coding exons (always CDS)
  - [contains_orf_start]: 27,082 exons (4.7%) - always span start codon
  - [3'UTR]: 21,031 exons (3.7%) - always 3'UTR
  - [5'UTR]: 20,536 exons (3.6%) - always 5'UTR
  - [contains_orf_stop]: 19,403 exons (3.4%) - always span stop codon

- **26.7% of exons have multiple region types** (152,723 exons - functionally different across isoforms)
  - [CDS, non_coding]: 23,191 exons (4.1%) - coding in some isoforms, non-coding in others
  - **[3'UTR, CDS]: 17,272 exons (3.0%)** - CDS vs 3'UTR (**NMD-relevant: premature termination!**)
  - **[5'UTR, CDS]: 16,963 exons (3.0%)** - CDS vs 5'UTR (alternative start sites)
  - [CDS, contains_orf_start]: 9,835 exons (1.7%)
  - [CDS, contains_orf_stop]: 9,768 exons (1.7%)
  - [contains_orf_start, non_coding]: 7,266 exons (1.3%)
  - [5'UTR, contains_orf_start]: 5,754 exons (1.0%)
  - [5'UTR, CDS, contains_orf_start]: 5,227 exons (0.9%)
  - [3'UTR, contains_orf_stop]: 5,199 exons (0.9%)
  - [3'UTR, CDS, contains_orf_stop]: 5,153 exons (0.9%)
  - Plus 10+ more combinations (each <0.7%)

**Key Concepts:**

*Region type is per (isoform, union_exon) pair:*
- Same union exon can have different region types in different isoforms
- 26.7% of exons show functional switching between isoforms
- These represent biologically meaningful isoform differences

*Critical NMD-relevant patterns:*
- **[3'UTR, CDS]: 17,272 exons (3.0%)** - The most important pattern for NMD
  - Exons that are CDS in some isoforms but 3'UTR in others
  - Indicates premature termination codons (PTCs) in some isoforms
  - Classic NMD trigger: coding exon becomes 3'UTR due to upstream stop codon
- **[5'UTR, CDS]: 16,963 exons (3.0%)** - Alternative translation start sites
  - May affect NMD through uORFs or altered CDS boundaries
- **[CDS, non_coding]: 23,191 exons (4.1%)** - Entire isoforms lose coding capacity

*Strand-aware coordinate system:*
- Uses genomic coordinate ranges (cds_start=min, cds_stop=max)
- Works correctly for both plus and minus strand genes
- Before fix: minus strand CDS exons misclassified as 5'UTR (strand bug)
- After fix: CDS mappings increased from 24.9% to 49.8% (doubled!)

**Validation:**
- ✓ 100% of mappings successfully annotated (2,712,956 / 2,712,956)
- ✓ Strand-aware coordinate system validated on minus strand genes
- ✓ CDS percentage doubled after strand fix (24.9% → 49.8%)
- ✓ Region type distributions match biological expectations

---

### **Script nmd/06: Filter to Analysis-Ready Subset**

**Input:**
- `dge_isoform_nofilter_2026.2.7.rds` (original DGEList)
- `data/*.rds` (7 unfiltered files from scripts 01-05)

**Process:**

**Step 1: Apply filterByExpr**
- Build design matrix: `~treatment + ct + ct*treatment` (38 samples, 12 coefficients)
- Apply `filterByExpr(dge, design, min.count=5, min.total.count=10, min.prop=0)`
- **Matches long-read DGE analysis exactly**
- Result: 300,209 / 1,746,281 isoforms pass (17.2%)

**Step 2: Exclude Fusion/Chimeric Genes**
- Categorize gene IDs:
  - **KEEP:** GENCODE (^ENSG[0-9]+$) + PacBio Novel (^novelGene_)
  - **EXCLUDE:** Fusion/Chimeric (e.g., ENSG_ENSG)
- After filterByExpr: 1,485 fusion genes present
- Final result: 296,680 isoforms from 59,386 genes

**Step 3: Filter All Downstream Data**
- Create isoform_ids_to_keep (296,680) and gene_ids_to_keep (59,386) lists
- Filter all 7 data files consistently:
  1. expression_data → expression_data_filtered.rds (37MB)
  2. dominant_isoforms → dominant_isoforms_filtered.rds (1.1MB)
  3. isoform_structures → isoform_structures_filtered.rds (6.2MB)
  4. union_exons → union_exons_filtered.rds (3.0MB, 223,300 exons from 26,553 genes)
  5. isoform_union_mapping → isoform_union_mapping_filtered.rds (4.1MB)
  6. isoform_cds_metadata → isoform_cds_metadata_filtered.rds (1.2MB)
  7. isoform_union_exons_annotated → isoform_union_exons_annotated_filtered.rds (4.4MB, 283,722 mappings)

**Output:**
- `data/*_filtered.rds` (7 files)
- `results/filtering_report_script06.txt` (detailed statistics)

**Filtering Summary:**

| Metric | Original | After Script 01 | After Script 06 | % of Original | % of Script 01 |
|--------|----------|-----------------|-----------------|---------------|----------------|
| **Isoforms** | 1,746,281 | 502,994 (28.8%) | 296,680 | 17.0% | 59.0% |
| **Genes** | 231,865 | 192,156 (82.9%) | 59,386 | 25.6% | 30.9% |
| **GENCODE Genes** | 87,079 | - | 29,768 | 34.2% | - |
| **Novel Genes** | 141,069 | - | 29,618 | 21.0% | - |
| **Fusion Genes** | 3,717 | - | 0 | 0% | - |
| **Union Exons** | 571,954 | - | 223,300 | 39.1% | - |
| **Annotated Mappings** | 2,712,956 | - | 283,722 | 10.5% | - |

**Gene Categories (Filtered):**
- GENCODE: 29,768 genes (50.1%)
- PacBio Novel: 29,618 genes (49.9%)
- Balanced representation

**Key Decisions:**
- **Expression filtering:** Matches DGE analysis exactly (reproducible)
- **Gene category filtering:** Excludes fusions, keeps novel genes
- **Consistent filtering:** All downstream data filtered with same criteria
- **Ordering change:** Filtering BEFORE event detection (more efficient)

**Rationale:**
- Focus on adequately expressed isoforms (filterByExpr validated by DGE)
- Exclude fusion artifacts (ENSG_ENSG chimeras)
- Keep PacBio novel genes (genuine discoveries, not artifacts)
- Ensures all downstream analyses use identical gene/isoform sets
- **Efficiency:** Run intensive event detection on 41% fewer union exons

---

## **DATA FLOW: PHASE 2 (Filtering)**

*See Script nmd/06 above*

---

## **DATA FLOW: PHASE 3 (Classification + Pair Generation)**

### **Script nmd/07: Classify Isoforms and Generate Pairs**

**Input:**
- `data/expression_data_filtered.rds` — tidy CPM
- `data/sample_metadata.rds` — sample annotations
- `data/isoform_structures_filtered.rds` — for validation
- 6 DE CSVs: `longread_dge/nmd_dge_{ct}_2026.1.18.csv`

**Process:**
1. Classify isoforms as NMD-sensitive (`adj.P.Val < 0.05 AND logFC > 0`) or non-NMD (`adj.P.Val > 0.95`)
2. Apply 5% expression filter (treatment-agnostic for C1/C2, DMSO-only for C3/C4)
3. Calculate dominance proportions (NMD in Smg1i samples, non-NMD in DMSO samples)
4. Generate comparison pairs for 4 comparison types x 7 runs (all_samples + 6 cell types)
5. Deduplicate pairs across all 28 sets

**Four Comparison Types:**

| | Dominant | Comparator | Pairs/gene |
|---|---|---|---|
| **C1** | Dominant non-NMD (DMSO) | Dominant NMD-sensitive (Smg1i) | 1 |
| **C2** | Top non-NMD by CPM (DMSO) | Top NMD-sensitive by CPM (Smg1i) | 1 |
| **C3** | Dominant non-NMD (DMSO) | All other non-NMD | Multiple |
| **C4** | Dominant non-NMD (DMSO) | Next-best non-NMD by CPM | 1 |

**Output:**
- `comparisons/{C1,C2,C3,C4}/{run}/classification.rds` — per-run classification
- `comparisons/{C1,C2,C3,C4}/{run}/pairs.tsv` — per-run pairs
- `comparisons/deduplicated/all_pairs.tsv` — unique pairs across all 28 sets

**CLI Flags:**
- `--test N`: Limit to first N genes

---

## **DATA FLOW: PHASE 4 (Event Detection + Reconstruction Validation)**

### **Script core/08: Extract Splicing Choice Profiles**

**Input:**
- `data/isoform_union_exons_annotated_filtered.rds`
- `data/dominant_isoforms_filtered.rds`
- `data/union_exons_filtered.rds`
- `data/isoform_structures_filtered.rds`

**Process:**
1. For each gene:
   - Identify dominant isoform (from dominant_isoforms, or from `--pairs-file`)
   - Get dominant's union exon pattern
2. For each non-dominant isoform in gene:
   - Compare to dominant isoform via hierarchical event detection
   - Extract comprehensive splicing events:
     - **TSS/TES changes:** First/last exon boundaries
     - **Exon inclusion:** Missing/added union exons (by region type)
     - **Boundary modifications:** A5SS/A3SS detection
     - **Splicing dysfunction:** IR/Partial_IR detection
     - **Complexity metrics:** Number of exons, junctions, length
3. Optionally reconstruct and verify each pair on-the-fly (`--reconstruction_check`)

**Output:**
- `data/splicing_choice_profiles.rds`
- `data/splicing_choice_profiles_intermediate.rds` (after batch 1, for development)

**CLI Flags:**
- `--test N`: Limit to first N genes
- `--reconstruction_check`: On-the-fly reconstruction verification per pair. Adds `reconstruction_status` and `reconstruction_reason` columns to profiles. Prints PASS/FAIL/ERROR summary at end. Eliminates need to run Script core/09 separately.
- `--pairs-file <path>`: Specify explicit contrast pairs (TSV with columns: `gene_id`, `dominant_isoform_id`, `comparator_isoform_id`). When used, only the specified pairs are processed instead of auto-generating all non-dominant vs dominant pairs. Skips loading `dominant_isoforms_filtered.rds`.
- `--output <path>`: Custom output path for splicing profiles RDS (default: `data/splicing_choice_profiles.rds`).

**Key Features:**
- **Comprehensive event detection:** A5SS, A3SS, SE, Missing_Internal, IR, Partial_IR, IR_diff
- **Hierarchical detection order:** IR → boundary shifts → exon skipping → terminal events
- **On-the-fly reconstruction check:** Validates detection accuracy without separate core/09 run
- **Explicit pairs mode:** Enables targeted re-analysis of specific isoform pairs
- **Runs on filtered data:** More efficient, focuses on expressed isoforms

---

### **Script core/09: Validate Reconstruction (Standalone)**

**Purpose:** Standalone batch reconstruction validation. Useful for validating profiles generated without `--reconstruction_check`, or for re-verification after code changes.

**Input:**
- `data/splicing_choice_profiles.rds` (events per pair)
- `data/isoform_structures_filtered.rds`
- `data/union_exons_filtered.rds`

**Process:**
1. For each isoform pair with detected events:
   - Extract comparator exons
   - Apply events to reconstruct the dominant isoform
   - Compare reconstructed exons to original dominant exons
2. Classification:
   - **PASS**: All exons match (same count, same coordinates within tolerance)
   - **FAIL**: Exon mismatch (count or coordinate differences)
   - **ERROR**: Reconstruction could not complete (e.g., no comparator exons)

**Output:**
- `data/reconstruction_verification.rds`
- `logs/reconstruction_verification.log`

**Validation Results:**
- Curated test suite (synthetic + real failures): 126/126 (100%)
- GENCODE test data: 4258/4274 = 99.6%

**Shared Libraries Used:**
- `scripts/core/event_detection_functions.R` -- Hierarchical event detection, detection thresholds (TSS_TOLERANCE, TES_TOLERANCE)
- `scripts/core/reconstruction_functions.R` -- Reconstruction (`reconstruct_dominant_v2`) and verification (`verify_transcript`)

---

## **DATA FLOW: PHASE 5 (Analysis)**

*Scripts nmd/09-14*

**Will use:** `data/splicing_choice_profiles.rds` or `comparisons/deduplicated/all_splicing_profiles.rds`

**Analysis scripts:**
- `nmd/09_analyze_complexity_relationship.R`: Complexity analysis
- `nmd/10_analyze_cooccurrence.R`: Co-occurrence analysis
- `nmd/11_analyze_spatial_patterns.R`: Spatial pattern analysis
- `nmd/12_analyze_functional_context.R`: Functional context analysis
- `nmd/13_analyze_patterns.R`: Pattern classification
- `nmd/14_generate_report.Rmd`: Final report

**Outputs:** `results/*.tsv`, `figures/*.pdf`

*See ANALYSIS_PLAN.md for detailed analysis specifications*

---

## **KEY DATA CHARACTERISTICS**

### **Final Analysis-Ready Dataset (After Script 06):**

**Isoforms:** 296,680
- 150,781 GENCODE (ENST*)
- 145,899 PacBio (PB.*)

**Genes:** 59,386
- 29,768 GENCODE genes
- 29,618 PacBio novel genes
- 0 fusion genes

**Splicing Profiles:** (to be calculated post-filtering)
- Non-dominant vs dominant comparisons
- Covering genes with 2+ major isoforms

**Samples:** 38
- 6 cell types: at2, dd, ddali, do, fb, mv
- 2 treatments: DMSO, Smg1i
- Paired design

**Design Matrix:**
- Formula: ~treatment + ct + ct*treatment
- Interaction term tests cell-type-specific treatment effects

---

## **FILTERING PHILOSOPHY**

### **Two-Stage Filtering:**

**Script 01 (Isoform-Level):**
- **Purpose:** Reduce noise, focus on well-expressed isoforms
- **Threshold:** ≥5% within gene in at least one sample
- **Scope:** Isoform abundance relative to gene
- **Result:** 502,994 "major isoforms"

**Script 06 (Gene-Level):**
- **Purpose:** Match DGE analysis, exclude artifacts
- **Threshold:** filterByExpr (expression across samples)
- **Scope:** Absolute gene/isoform expression levels
- **Result:** 296,680 adequately expressed isoforms

### **Why Two Filters?**

1. **Memory efficiency:** Script 01 reduces dataset size early
2. **Biological focus:** Script 01 removes low-abundance noise within genes
3. **Statistical validity:** Script 06 ensures adequate expression for DGE
4. **Reproducibility:** Script 06 matches published DGE analysis exactly

### **Excluded Categories:**

❌ **Low-abundance isoforms:** <5% within gene (Script 01)
❌ **Lowly expressed genes:** Fail filterByExpr (Script 06)
❌ **Fusion/chimeric genes:** ENSG_ENSG artifacts (Script 06)

✓ **Included Categories:**

✅ **GENCODE isoforms:** Well-annotated reference transcripts
✅ **PacBio novel isoforms:** Novel discoveries from long-read data
✅ **Adequately expressed:** Pass filterByExpr across samples
✅ **Major isoforms:** ≥5% abundance within gene

---

## **FILE ORGANIZATION**

```
/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/

├── data/
│   ├── expression_data.rds                              [nmd/01 - unfiltered]
│   ├── dominant_isoforms.rds                            [nmd/01 - unfiltered]
│   ├── sample_metadata.rds                              [nmd/01]
│   ├── filtering_stats.rds                              [nmd/01]
│   ├── isoform_structures.rds                           [core/02 - unfiltered]
│   ├── union_exons.rds                                  [core/03 - unfiltered]
│   ├── isoform_union_mapping.rds                        [core/03 - unfiltered]
│   ├── isoform_cds_metadata.rds                         [core/04 - unfiltered]
│   ├── isoform_union_exons_annotated.rds                [core/05 - unfiltered]
│   ├── expression_data_filtered.rds                     [nmd/06 - FILTERED]
│   ├── dominant_isoforms_filtered.rds                   [nmd/06 - FILTERED]
│   ├── isoform_structures_filtered.rds                  [nmd/06 - FILTERED]
│   ├── union_exons_filtered.rds                         [nmd/06 - FILTERED]
│   ├── isoform_union_mapping_filtered.rds               [nmd/06 - FILTERED]
│   ├── isoform_cds_metadata_filtered.rds                [nmd/06 - FILTERED]
│   ├── isoform_union_exons_annotated_filtered.rds       [nmd/06 - FILTERED]
│   ├── splicing_choice_profiles.rds                     [core/08]
│   └── reconstruction_verification.rds                  [core/09]
│
├── comparisons/                                         [nmd/07 + core/08]
│   ├── {C1,C2,C3,C4}/{run}/                            [Per-comparison outputs]
│   │   ├── classification.rds
│   │   └── pairs.tsv
│   └── deduplicated/
│       ├── all_pairs.tsv                                [Unique pairs across all sets]
│       └── all_splicing_profiles.rds                    [core/08 on deduplicated pairs]
│
├── results/
│   └── filtering_report_script06.txt                    [Script 06 report]
│
├── scripts/
│   ├── core/                                            [Generic, reusable]
│   │   ├── 02_extract_isoform_structures.R              [Data prep]
│   │   ├── 03_build_union_exons.R                       [Data prep]
│   │   ├── 04_extract_cds_annotations.R                 [Data prep]
│   │   ├── 05_annotate_region_types.R                   [Data prep]
│   │   ├── 08_extract_splicing_profiles.R               [Event detection]
│   │   ├── 09_validate_reconstruction.R                 [Reconstruction validation]
│   │   ├── event_detection_functions.R                  [Shared library]
│   │   ├── reconstruction_functions.R                   [Shared library]
│   │   └── visualization_functions.R                    [Shared library]
│   ├── nmd/                                             [NMD study-specific]
│   │   ├── 01_prepare_expression_data.R                 [Data prep]
│   │   ├── 06_filter_to_analysis_subset.R               [Filtering]
│   │   ├── 07_classify_and_pair.R                       [Classification + pairing]
│   │   └── 09-14_*.R                                    [Downstream analysis]
│   ├── tests/                                           [Validation suite]
│   ├── dev/                                             [Development/diagnostic]
│   └── archive/                                         [Archived earlier versions]
│
└── logs/
    └── *.log                                             [Execution logs]
```

---

## **QUALITY CONTROL CHECKPOINTS**

### **nmd/01:**
✓ 100% of genes have filtering stats recorded
✓ Mean dominant proportion: 91.3%

### **core/02:**
✓ 100% coverage: all 502,994 major isoforms have structures
✓ 0 missing isoforms

### **core/03:**
✓ 0 overlapping union exons (validated in tests)
✓ All union exons non-overlapping within gene

### **core/04:**
✓ Coding status assigned for all isoforms
✓ CDS coordinates validated

### **core/05:**
✓ Region type distribution as expected (5'UTR → CDS → 3'UTR)
✓ No unexpected "unknown" region types

### **nmd/06:**
✓ Exact match: 300,209 isoforms from filterByExpr (reproducible)
✓ Gene categories validated (GENCODE, Novel, Fusion)
✓ All 7 data files filtered consistently

### **core/08:**
✓ 102,928 profiles created (full dataset, legacy mode)
✓ All profiles have complexity category assigned
✓ Mean 5.0 differences per profile
✓ Supports explicit pairs mode (--pairs-file) for comparison-based analysis

### **core/09:**
✓ Curated test suite (synthetic + real failures): 126/126 (100%)
✓ GENCODE test data: 4258/4274 = 99.6% (0 FAILs, 16 ERRORs)
✓ All ERRORs are "No comparator exons" edge cases
✓ All event types (including IR) use direct-coordinate reconstruction

---

## **NOTES AND DECISIONS LOG**

### **2026-02-14:**
- Implemented scripts 01-05 (data preparation) + 07 (splicing profiles)
- Script 01 v2 created due to OOM error on full dataset
- Fixed version number matching in script 02
- Fixed function masking (S4Vectors vs dplyr) in scripts

### **2026-02-15:**
- Added script 06 (filtering step)
- Clarified gene categorization (novelGene vs fusion genes)
- Verified filterByExpr reproducibility (300,209 isoforms)
- Documented filtering philosophy and rationale
- Created DATA_FLOW.md to track data transformations

### **2026-02-20:**
- Rewrote Script 03 for atomic union exon construction (BOTH boundary fix)
- Added shared libraries: scripts/core/event_detection_functions.R, scripts/core/reconstruction_functions.R
- Added reconstruction validation (now core/09)
- Created SPLICE_FUNCTION_CATALOG.md and EVENT_DETECTION_FLOW.md

---

## **FUTURE WORK**

### **Immediate:**
- Run script 07 to generate splicing choice profiles
- Implement analysis scripts 09-13

### **Potential Enhancements:**
- Cell-type-specific dominant isoforms (currently overall dominant)
- Sequence-based frameshift analysis
- NMD prediction (PTC identification)
- Protein domain impact analysis

---

**END OF DATA_FLOW.md**
