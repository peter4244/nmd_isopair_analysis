# Isoform Choice Analysis: Splicing Decisions When Cells Shift from Dominant Isoforms

**Version:** 6.0
**Date:** 2026-02-09 (Updated: 2026-02-20)
**Analysis Framework:** Dominant isoform-centric splicing choice characterization

---

## **IMPLEMENTATION STATUS** (Updated 2026-02-20)

### **Phase 1: Data Preparation** ✅ **COMPLETE**

| Script | Status | Output | Notes |
|--------|--------|--------|-------|
| 01 | ✅ Complete | `dge_isoform_nofilter_2026.2.7.rds` | Expression data loaded |
| 02 | ✅ Complete | `isoform_structures.rds` | GENCODE + SQANTI structures extracted |
| 03 | ✅ Complete | `union_exons.rds`, `isoform_union_mapping.rds` | Atomic union exon models built for all genes |
| 04 | ✅ Complete | `isoform_cds_metadata.rds` | CDS annotations extracted |
| 05 | ✅ Complete | `isoform_union_exons_annotated.rds` | Region types assigned (5'UTR/CDS/3'UTR) |
| 06 | ✅ Complete | `*_filtered.rds` | Analysis subset filtered (18,754 dominant, 24,787 non-dominant) |

### **Phase 2: Event Detection & Reconstruction Validation** ✅ **COMPLETE**

| Script | Status | Output | Notes |
|--------|--------|--------|-------|
| 07 | ✅ Complete | `splicing_choice_profiles.rds` | ~24,645 profiles from ~24,000 genes |
| 08 | ✅ Complete | Reconstruction verification | Reconstructs dominant from comparator + events, verifies match |

**Event Detection Implementation**:
- ✅ Alt_TSS / Alt_TES detection (terminal boundary events)
- ✅ A5SS / A3SS detection (<100bp boundary differences)
- ✅ Partial_IR_5 / Partial_IR_3 detection (>=100bp shared boundary differences, split by side)
- ✅ IR detection (exon spanning >=2 other exons)
- ✅ IR_diff_5 / IR_diff_3 / IR_diff_5_3 detection (intron retention with boundary mismatch)
- ✅ SE detection (skipped exon with strict flanking check)
- ✅ Missing_Internal detection (absent exon, flanking condition not met)

**Terminal Boundary Rules** (Critical Implementation):
- ✅ Asymmetric terminal exon handling with `second_event` emission
- ✅ Dual-boundary decomposition for internal exons with both boundaries differing
- ✅ LOSS-before-GAIN ordering in reconstruction (`arrange(direction != "LOSS", event_type)`)

**Validation** ✅ **PASSING**:
- Curated test suite (synthetic + real failures): **126/126 tests passing** (100%) — both GTF-built and production UEs
- Real data (GENCODE): **4258/4274 = 99.6%** -- 0 FAILs + 16 ERRORs (all "No comparator exons")
- Coverage: Both strands, all event types, IR subtypes, edge cases, multi-event scenarios
- Validation framework: `scripts/tests/run_tests.R`, `scripts/tests/test_data/`

### **Phase 3: Analysis Execution** 🔄 **PENDING**

| Script | Status | Purpose |
|--------|--------|---------|
| 09 | ⏳ Pending | Complexity relationship analysis |
| 10 | ⏳ Pending | Co-occurrence analysis |
| 11 | ⏳ Pending | Spatial organization |
| 12 | ⏳ Pending | Functional context |
| 13 | ⏳ Pending | Pattern classification + report generation |

### **Key Achievements**

1. **Complete atomic union exon framework**: All 23,483 genes with atomic union exon models
2. **Validated event detection**: 126/126 curated tests passing (100% accuracy, synthetic + real failures)
3. **Reconstruction validation**: All event types reconstruct via direct coordinates (no UE lookups), verified round-trip correctness
4. **Comprehensive event type coverage**: Alt_TSS, Alt_TES, A5SS, A3SS, Partial_IR_5, Partial_IR_3, IR, IR_diff_5, IR_diff_3, IR_diff_5_3, SE, Missing_Internal
5. **Shared function libraries**: `event_detection_functions.R`, `reconstruction_functions.R`, `create_union_exons_and_junctions.R`
6. **Robust terminal detection**: TSS/TES detection uses both coordinate tolerance and exon overlap check
7. **Production-ready**: Processing full dataset (~24,645 profiles)

### **Documentation**

- ✅ **ANALYSIS_PLAN.md**: Complete analysis framework (this document)
- ✅ **METHODS.md**: Detailed algorithmic documentation with pseudocode
- ✅ **SPLICE_FUNCTION_CATALOG.md**: Catalog of all splicing event detection functions
- ✅ **EVENT_DETECTION_FLOW.md**: End-to-end event detection workflow documentation
- ✅ **Testing framework**: Synthetic validation (44 tests) + reconstruction verification
- ✅ **Validation**: 44/44 synthetic tests, 99.6% real data pass rate

### **Next Steps**

1. Re-run 1000-gene subset pipeline to verify fixes at scale
2. Run full NMD dataset with updated reconstruction logic
3. Begin Script 09: Analyze complexity vs event count relationship
4. Use Script 09 results to inform stratification strategy for Scripts 10-13
5. Continue with downstream statistical analyses (co-occurrence, spatial, functional)

---

## **SECTION 1: BIOLOGICAL FRAMEWORK & MOTIVATION**

### **Core Premise**
The dominant isoform in a gene plays an important functional role - this is why it's dominant. When cells shift away from the dominant isoform to produce alternative (non-dominant) isoforms, they are making deliberate assembly choices.

### **Central Research Question**
**What is the cell trying to achieve when it shifts away from the dominant isoform?**

### **Two-Stage Investigation**

**Stage 1 (Current Analysis):** Catalog the assembly choices
- What combinations of TSS, TES, and internal splicing choices create the alternative isoforms?
- How do these choices co-occur?
- Are there spatial/topological patterns to these choices?
- **Focus:** Understanding HOW cells make isoform assembly decisions

**Stage 2 (Future Analysis):** Functional impact
- What sequence regions are added or removed?
- How do these changes affect protein domains, NMD targeting, regulatory elements?
- What functional outcomes result from these assembly choices?
- **Focus:** Understanding WHY cells make these decisions

### **Current Scope**
This analysis focuses exclusively on **Stage 1**: characterizing the splicing choice repertoire and patterns before assessing functional consequences.

---

## **SECTION 2: CONCEPTUAL FRAMEWORK - HIERARCHICAL SPLICING CHOICES**

### **Overview**
When assembling an isoform different from the dominant, the cell deploys choices across multiple hierarchical levels. We organize these into a biologically interpretable framework.

### **The Cell's Splicing Choice Repertoire**

#### **Level 1: Transcript Boundaries**
Decisions about where transcription starts and ends.

**A. TSS (Transcription Start Site) Choice**

**Primary classification: Region usage**
- **Same region:** TSS in the same exon as dominant
  - *Then classify modification:*
    - **Same position:** Identical TSS
    - **Internal modification:** TSS further into the exon → shorter 5' end
    - **Extended modification:** TSS upstream within the exon → longer 5' end

- **Different region:** TSS in a completely different exon than dominant
  - *Biological interpretation:* Alternative promoter usage, different first exon
  - *Characterization:*
    - Which exon contains the alternative TSS?
    - How many exons upstream/downstream from dominant TSS?
    - UTR vs CDS impact

**B. TES (Transcription End Site) Choice**

**Primary classification: Region usage**
- **Same region:** TES in the same exon as dominant
  - *Then classify modification:*
    - **Same position:** Identical TES
    - **Internal modification:** TES earlier in the exon → shorter 3' end
    - **Extended modification:** TES later within the exon → longer 3' end

- **Different region:** TES in a completely different exon than dominant
  - *Biological interpretation:* Alternative polyadenylation signal, different terminal exon
  - *Characterization:*
    - Which exon contains the alternative TES?
    - How many exons upstream/downstream from dominant TES?
    - UTR vs CDS impact

**C. TSS/TES Characterization Framework**

**Track three dimensions:**
1. **Outside boundary:** extended / internal / same (the terminus itself)
2. **Inside boundary:** extended / internal / same (the junction with next/previous exon)
3. **Overlap:** yes / no (do terminal exons share sequence?)

**Functional summary for overlap cases:**
- **UTR:** more / less / same
- **CDS:** more / less / same
- **ORF affected:** yes / no (TSS past start codon, or TES before stop codon)

**D. UTR vs CDS Classification**

To understand functional implications, we must distinguish:
- **5' UTR changes:** Affect translation efficiency, regulatory elements (uORFs, IRES, etc.)
- **CDS (Coding Sequence) changes:** Directly affect protein sequence and domains
- **3' UTR changes:** Affect mRNA stability, localization, miRNA binding sites

**Data requirements:**
- CDS start/stop coordinates for each isoform
- **GENCODE annotations:** Available for ENST transcripts
- **SQANTI predictions:** ORF predictions for novel PacBio isoforms
- **Union exon annotation:** Classify each union exon as 5'UTR / CDS / 3'UTR / contains_orf_start / contains_orf_stop / non_coding / unknown

#### **Level 2: Internal Exon Choices** (Exon-centric view)

**A. Exon Inclusion Decisions:** *"Do I want this exon?"*

Fundamental yes/no decisions about exon presence:
- **SE (Skipped Exon):** Single exon excluded relative to dominant
- **MXE (Mutually Exclusive Exons):** Choose exon A OR exon B (but not both)
- **Runs of MISSING:** Multiple consecutive exons excluded
- **Exon addition:** Exons present in non-dominant but absent in dominant

*Biological interpretation:* These represent discrete choices about including/excluding functional units (exons).

**B. Exon Boundary Modifications:** *"Do I want to change the boundaries of this exon?"*

Fine-tuning exon length and composition:
- **A5SS (Alternative 5' Splice Site):** Adjust where the exon ends (donor site)
- **A3SS (Alternative 3' Splice Site):** Adjust where the exon starts (acceptor site)

*Biological interpretation:* These represent choices to include/exclude partial exon sequences, affecting the exact protein sequence encoded.

**With UTR/CDS classification:**
- A5SS/A3SS in 5' UTR → affects regulatory sequences
- A5SS/A3SS in CDS → affects protein sequence (may maintain/disrupt frame)
- A5SS/A3SS in 3' UTR → affects regulatory elements

**C. Splicing Dysfunction:** *"Have things gone haywire?"*

Failures in normal splicing or regulated NMD-targeting:
- **IR (Intron Retention):** Complete failure to remove an intron
- **Partial_IR (Partial Intron Retention):** Partial retention of intronic sequence

*Biological interpretation:* May represent splicing errors, quality control issues, or deliberate NMD-targeting mechanisms.

**With UTR/CDS classification:**
- IR in CDS → likely triggers NMD (introduces PTCs)
- IR in UTR → may escape NMD, affect regulation

#### **Level 3: Spatial Topology of Boundary Modifications**

When both A5SS and A3SS events occur in the same non-dominant isoform, their spatial relationship may reveal coordination:

**Face-to-Face (F2F):** A3SS and A5SS on opposite sides of the SAME intron
```
Exon A [====] ----intron---- [====] Exon B
              ^A3SS    ^A5SS
           (acceptor)(donor)
```
*Interpretation:* Coordinated regulation of both boundaries of a single intron

**Back-to-Back (B2B):** Alternative splice sites on sequential exons
```
Example: A5SS on consecutive exons
Exon A [====]^donor1 --intron-- [====]^donor2 Exon B
```
*Interpretation:* Clustering or cascading of splice site variation

**Distributed:** Multiple exons separating A5SS and A3SS events
*Interpretation:* Independent or spatially unrelated boundary modifications

---

## **SECTION 3: ANALYSIS UNIT & DATA SOURCES**

### **Unit of Analysis: Non-Dominant vs Dominant Comparisons**

**Core approach:** For each gene, identify the dominant isoform (highest overall expression), then characterize each non-dominant isoform relative to this dominant reference.

**Scope:** Overall dominant isoform (across all cell types)

### **Data Sources**

#### **1. Primary Data: DGEList Object**
- **File:** `/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds`
- **Contains:**
  - Isoform counts (all isoforms, unfiltered)
  - Sample metadata (phenotype data in `$samples` slot)
  - Normalization factors
- **Use:**
  - Calculate isoform expression levels
  - Identify dominant isoforms per gene
  - Extract sample information (cell types, treatments, subjects)

#### **2. CDS/ORF Annotations**
**GENCODE (ENST transcripts):**
- **File:** `/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz`
- **Contains:** CDS features for annotated transcripts

**SQANTI (PacBio isoforms):**
- **File:** `/Users/petecastaldi/claude_projects/nmd/reference_files/sqanti3_corrected.cds.gff3`
- **Contains:** ORF predictions for novel isoforms

**Use:** Annotate union exons as 5'UTR / CDS / 3'UTR; determine ORF impact

#### **3. Isoform Structure (Exon Coordinates)**
- **Source:** GENCODE GFF3 (for ENST) and SQANTI/PacBio GFF (for PacBio isoforms)
- **Extract:** Exon coordinates for each isoform
- **Use:** Build union exon structure, identify exon inclusion/exclusion patterns

### **Data Preparation: Building Analysis-Ready Structures**

**Note:** Version 6.0 builds all structures from scratch using original source data only. No dependencies on Version 5.0 intermediate files.

#### **Step 1: Extract Isoform Structures from GFF Files**
- Parse GENCODE GFF3 → exon coordinates for ENST transcripts
- Parse PacBio GFF → exon coordinates for PacBio isoforms
- Create isoform structure table:
  ```
  isoform_id | gene_id | seqnames | strand | exon_number | exon_start | exon_end
  ```

#### **Step 2: Build Union Exon Structure per Gene**
- For each gene, merge overlapping exons across all isoforms
- Create union exons (maximal exon set)
- Assign union exon IDs and indices (5' → 3')
- Create mapping: (isoform_id, union_exon_id) → presence/absence

#### **Step 3: Extract CDS Annotations**
- Parse GENCODE GFF3 → CDS start/stop for ENST transcripts
- Parse SQANTI GFF3 → CDS start/stop for PacBio isoforms
- Create isoform CDS metadata:
  ```
  isoform_id | coding_status | cds_start | cds_stop | orf_length
  ```
  Where `coding_status` = "coding" / "non-coding" / "unknown"

#### **Step 4: Annotate Union Exons with Region Type**
- For each (isoform, union_exon) pair:
  - If isoform is non-coding → `region_type = "non_coding"`
  - If isoform is coding:
    - Determine overlap with CDS region
    - Assign: "5'UTR" / "CDS" / "3'UTR" / "contains_orf_start" / "contains_orf_stop" / "unknown"

**Region type assignment logic:**
```
if (exon_end < cds_start):
    region_type = "5'UTR"

elif (exon_start > cds_stop):
    region_type = "3'UTR"

elif (exon contains cds_start):
    region_type = "contains_orf_start"

elif (exon contains cds_stop):
    region_type = "contains_orf_stop"

elif (exon_start >= cds_start AND exon_end <= cds_stop):
    region_type = "CDS"

else:
    region_type = "unknown"
```

**Key clarification:** `region_type` is a property of the **(isoform, union_exon) pair**, not the union exon alone, because the same union exon plays different functional roles in different isoforms.

#### **Step 5: Identify Dominant Isoforms**
- Calculate isoform expression levels (CPM or TPM from DGEList)
- For each gene, identify dominant isoform:
  - Highest overall expression across all samples
  - Or cell-type-specific dominance (future refinement)

**Output:** Analysis-ready data structure with:
- Union exon structures per gene
- Isoform-union exon mapping with region types
- Dominant isoform identification
- All built fresh from original source data

---

## **SECTION 4: SPLICING CHOICE PROFILE STRUCTURE**

### **Goal**
For each non-dominant isoform, create a comprehensive profile characterizing all differences from the dominant isoform.

### **Profile Components**

#### **1. Transcript Boundaries**

**TSS (Transcription Start Site):**
```
tss_overlap: yes / no
tss_outside_boundary: extended / internal / same
tss_inside_boundary: extended / internal / same
tss_utr_impact: more / less / same
tss_cds_impact: more / less / same
tss_orf_affected: yes / no
```

**TES (Transcription End Site):**
```
tes_overlap: yes / no
tes_outside_boundary: extended / internal / same
tes_inside_boundary: extended / internal / same
tes_utr_impact: more / less / same
tes_cds_impact: more / less / same
tes_orf_affected: yes / no
```

#### **2. Exon Inclusion Decisions**

```
n_exons_removed: count (in dominant, missing in non-dominant)
n_exons_added: count (missing in dominant, present in non-dominant)
removal_pattern: "none" / "single" / "clustered" / "distributed"
addition_pattern: "none" / "single" / "clustered" / "distributed"

# By region type:
exons_removed_5utr: count
exons_removed_cds: count
exons_removed_3utr: count
exons_removed_orf_start: count
exons_removed_orf_stop: count
exons_added_5utr: count
exons_added_cds: count
exons_added_3utr: count

# Specific patterns:
involves_mxe: yes / no
consecutive_missing: max_run_length
```

#### **3. Exon Boundary Modifications**

```
n_a5ss: count
n_a3ss: count

# By region type:
a5ss_in_5utr: count
a5ss_in_cds: count
a5ss_in_3utr: count
a5ss_in_orf_start: count
a5ss_in_orf_stop: count
a3ss_in_5utr: count
a3ss_in_cds: count
a3ss_in_3utr: count
a3ss_in_orf_start: count
a3ss_in_orf_stop: count

# Spatial relationships (if both A5SS and A3SS present):
a5ss_a3ss_topology: "face_to_face" / "back_to_back" / "distributed" / "NA"
a5ss_a3ss_distance_exons: count / NA
```

#### **4. Splicing Dysfunction**

```
n_ir: count
n_partial_ir: count

# By region type:
ir_in_cds: count (likely NMD targets)
ir_in_utr: count
partial_ir_in_cds: count
partial_ir_in_utr: count

likely_nmd_target: yes / no / uncertain
```

#### **5. Spatial Summary**

```
total_changes: count (sum of all non-constitutive categories)
change_distribution: "5prime_clustered" / "3prime_clustered" / "distributed"
distance_tss_to_first_change: bp / NA
distance_last_change_to_tes: bp / NA
```

#### **6. Complexity Metrics**

```
complexity_n_exons: count (number of exons in non-dominant)
complexity_n_junctions: count
complexity_length: bp (transcript length)
complexity_n_union_exons_in_gene: count (gene-level complexity)
```

---

## **SECTION 5: RESEARCH QUESTIONS**

### **Category 1: Co-occurrence Analysis**
*Do specific splicing choices occur together more than expected?*

#### **1A. Within-Hierarchy Co-occurrence**

**Transcript boundaries:**
- TSS + TES: Do boundary changes at both ends co-occur?
- What proportion of TSS changes also have TES changes (and vice versa)?
- **What proportion of TSS and TES changes result in CDS changes (affect ORF)?**

**Exon inclusion:**
- Do removed exons cluster spatially (consecutive MISSING)?
- Distribution: single vs runs of MISSING exons

**Boundary modifications:**
- **Heterotypic:** A5SS + A3SS co-occurrence (different splice site types)
- **Homotypic:** A5SS + A5SS co-occurrence? A3SS + A3SS co-occurrence? (same type)

#### **1B. Cross-Hierarchy Co-occurrence**

**Boundaries × Internal events:**
- TSS changes × exon removal/addition
- TES changes × exon removal/addition
- TSS/TES changes × A5SS/A3SS events

**Exon inclusion × Boundary modifications:**
- Are MISSING exons flanked by A5SS/A3SS events?
- (Tests: exon skipping implemented via splice site choice)

**Dysfunction × Other changes:**
- IR/Partial_IR × other splicing changes
- (Tests: NMD-targeted isoforms have multiple alterations)

---

### **Category 2: Spatial Organization**
*Where do splicing changes occur along the transcript?*

#### **2A. Positional Bias** (all internal event types)

**For each event type (A5SS, A3SS, SE, MXE, IR, Partial_IR):**
- Enriched at 5' end? 3' end? Or uniformly distributed?
- Metrics: absolute position from TSS, relative position (0-1), exon index
- Statistical test: observed vs uniform distribution (KS test)

#### **2B. Spatial Topology** (multi-event patterns)

**When A5SS + A3SS co-occur:**
- **Face-to-face:** Both on same intron (coordinated boundary regulation)
- **Back-to-back:** On sequential exons (clustered variation)
- **Distributed:** Separated by multiple exons
- Enrichment test: observed vs expected frequencies

**When A5SS + A5SS co-occur (or A3SS + A3SS):**
- Do they preferentially affect sequential exons (back-to-back)?
- (Tests: hotspot regions for splice site variation)

#### **2C. Proximity Analysis**

**Distance measurements:**
- TSS to first internal change
- Last internal change to TES
- A5SS/A3SS to nearest MISSING exon
- Between co-occurring events of same/different types

---

### **Category 3: Functional Context**
*What regions and consequences are affected?*

#### **3A. Regional Distribution (UTR vs CDS)**

**For all event types, quantify:**

**Within coding isoforms:**
- % in 5' UTR
- % in CDS
- % in 3' UTR
- % in exons containing ORF start
- % in exons containing ORF stop

**Across all isoforms:**
- % in non-coding isoforms
- % with unknown coding status or region

**Comparisons:**
- Are certain event types enriched in specific regions?
- **Exon boundary modification susceptibility:**
  - **Are exons containing ORF start/stop more or less susceptible to A5SS/A3SS (exon boundary modifications)?**
  - Statistical test: Compare A5SS/A3SS frequency in ORF boundary exons vs regular CDS exons

**Note:** This specifically tests **exon boundary** modifications (A5SS/A3SS), not CDS/ORF boundary changes, in exons that happen to contain ORF start/stop codons.

#### **3B. ORF Impact**

**TSS changes:**
- % that move past start codon (ORF affected)
- Results: non-coding isoforms, alternative start codons

**TES changes:**
- % that truncate before stop codon (ORF affected)
- Results: incomplete ORFs, likely NMD targets

**CDS-altering events:**
- Exon removals in CDS: count and frequency
- IR/Partial_IR in CDS: count and frequency
- **Note:** Detailed frameshift analysis (in-frame vs frameshift) and sequence-based NMD prediction (PTC identification) deferred to future work

---

### **Category 4: Pattern Classification**
*What combinations of splicing choices are common?*

#### **4A. Isoform Complexity vs Splicing Event Count** (Foundational)

**Goal:** Establish baseline relationship before analyzing patterns

**Complexity metrics to test:**
1. Number of internal exons
2. Number of junctions
3. Transcript length (bp)
4. Number of union exons in gene

**Event count:** Total splicing changes distinguishing non-dominant from dominant

**Key analysis (for each metric):**
- Plot: Complexity vs number of events
- Fit models: Linear, polynomial, log-linear
- Report: R², relationship type, significance
- **Select primary metric:** Highest R² for stratification in 4B-4D

**Interpretation:**
- Establishes null expectation for event counts given complexity
- Informs whether to use stratification or regression adjustment

#### **4B. Pattern Analysis (Complexity-Controlled)**

**Within complexity strata or adjusted for complexity:**

**Frequency of profile types:**
- "Boundary-only" (TSS/TES change, internal unchanged)
- "Inclusion-only" (exon addition/removal, boundaries same)
- "Modification-only" (A5SS/A3SS, no exon changes)
- "Dysfunction" (IR/Partial_IR present)
- "Combined" (multiple categories)

**Analysis approach:**
- **Stratification:** Group by complexity bins (e.g., 2-5 exons, 6-10 exons, 11-20 exons, 20+ exons)
- **Or regression adjustment:** Control for complexity as covariate
- Compare pattern frequencies within strata

#### **4C. Most Common Combinations** (Complexity-Controlled)

**Top patterns by frequency (within complexity strata):**
- Example: In simple isoforms (2-5 exons): "Internal TSS + 1 exon removed"
- Example: In complex isoforms (20+ exons): "Internal TSS + 3 exons removed + 2 A3SS events"

**Questions:**
- Do different complexity levels favor different pattern types?
- Are certain patterns only feasible in complex isoforms?
- Does pattern diversity increase with complexity?

#### **4D. Complexity-Independent Pattern Features**

**After controlling for complexity:**
- Are certain event combinations still enriched (beyond what complexity predicts)?
- Hierarchical clustering of patterns (adjusted for complexity)
- Identification of "canonical" patterns that transcend complexity

---

### **Category 5: Statistical Framework**
*Rigorous testing approach for all questions*

#### **5A. Co-occurrence Testing (Category 1)**

**Approach 1: Simple Contingency Tables** (Quick overview)

**Unit of analysis:**
- Gene-level or isoform-level (specify per question)

**Statistical test:**
- **2×2 contingency table:**
  ```
                Feature B present | Feature B absent
  Feature A present    n_both           n_A_only
  Feature A absent     n_B_only         n_neither
  ```
- **Test:** Fisher's exact test (small n) or Chi-square test (large n)
- **Effect size:** Odds ratio (OR) = (n_both × n_neither) / (n_A_only × n_B_only)
- **Report:** OR, 95% CI, p-value, FDR q-value

**Interpretation:** Crude association, unadjusted for isoform structure

---

**Approach 2: Regression-Based Framework** (Structure-controlled)

**Rationale:** Simple contingency tables may show "co-occurrence" due to:
- Both events being more likely in complex isoforms (more opportunities)
- Both events favoring same union exons (hotspot effects)

**Step 1: Baseline Event Rate Modeling**

For each event type (A3SS, A5SS, SE, etc.), fit separately:

```
Logistic regression with LASSO/ridge regularization:

event_occurs ~ union_exon_id (as factor) + n_exons_in_isoform
```

**Where:**
- `event_occurs`: binary (1 if event type occurs in this isoform, 0 otherwise)
- `union_exon_id`: categorical variable for each union exon (high-dimensional)
- `n_exons_in_isoform`: isoform complexity metric
- **Regularization:** LASSO or ridge penalty to handle high-dimensional union exon factors
- **Cross-validation:** Select optimal penalty parameter

**Output:**
- Baseline probability of each event type given union exon composition
- Identifies "hotspot" union exons for each event type

**Step 2: Co-occurrence Testing (Controlling for Structure)**

Test if Event A and Event B co-occur beyond structural expectations:

```
Logistic regression with clustered robust SEs:

has_event_B ~ has_event_A + union_exon_id + n_exons
            + cluster(gene_id)  # robust standard errors
```

**Where:**
- `has_event_A`: binary predictor (Event A present)
- `union_exon_id`: controls for exon composition
- `n_exons`: controls for isoform complexity
- `cluster(gene_id)`: robust SEs accounting for multiple isoforms per gene

**Test:** Is coefficient for `has_event_A` significantly ≠ 0?

**Effect size:** Adjusted odds ratio for `has_event_A`

**Step 3: Comparative Reporting**

**For each event pair (A, B), report:**
- **Crude OR** (from Approach 1): Unadjusted association
- **Adjusted OR** (from Approach 2): Structure-controlled association
- **Interpretation:**
  - Crude OR ≈ Adjusted OR → minimal confounding by structure
  - Crude OR > Adjusted OR → confounded by shared complexity/composition
  - Adjusted OR significant → true biological co-occurrence

**Report table:**
```
Event_A | Event_B | Crude_OR | Adj_OR | 95%CI | P_crude | P_adj | FDR_q
A3SS    | A5SS    | 2.5      | 1.8    | 1.3-2.5 | 0.001 | 0.02  | 0.05
```

---

**Both approaches complement each other:**
- **Approach 1:** Fast, interpretable, generates hypotheses
- **Approach 2:** Rigorous, controls confounders, confirms true associations

---

#### **5B. Spatial Analysis Testing (Category 2)**

**2A. Positional Bias:**
- **Test:** Kolmogorov-Smirnov test
- **Null:** Uniform distribution along transcript (0-1 relative position)
- **Alternative:** Two-sided (enriched at either end)
- **Effect size:** Median position, mean position, D-statistic
- **Report:** KS D-statistic, p-value, median position

**2B. Topology Enrichment:**
- **Null model construction:**
  - **Simulation-based:** Randomly place events along transcript, calculate expected F2F/B2B/distributed frequencies
  - **Analytical:** Expected frequency = f(gene structure, n_events)
    - F2F: probability both events on same intron
    - B2B: probability events on adjacent exons
    - Distributed: 1 - F2F - B2B
- **Test:** Chi-square goodness-of-fit (observed vs expected frequencies)
- **Effect size:** Enrichment ratio = observed / expected
- **Report:** Observed counts, expected counts, enrichment ratio, p-value

**2C. Proximity Analysis:**
- **Test:** Permutation test
  - **Null:** Randomly permute event positions while preserving gene structure
  - **Test statistic:** Mean distance between co-occurring events
  - **Permutations:** 10,000
- **Alternative:** Wilcoxon rank-sum test (observed vs simulated distances)
- **Effect size:** Median distance, mean distance
- **Report:** Observed median, expected median (from permutations), p-value

#### **5C. Regional Distribution Testing (Category 3)**

**3A. Regional Enrichment:**
- **Null model:** Events distributed proportional to region size
  - Expected frequency in CDS = (CDS length) / (total transcript length)
  - Expected frequency in 5'UTR = (5'UTR length) / (total transcript length)
  - Expected frequency in 3'UTR = (3'UTR length) / (total transcript length)
- **Test:** Multinomial goodness-of-fit test or Chi-square test
- **Per-region test:** Binomial test (observed vs expected in focal region)
- **Effect size:** Enrichment ratio = (observed %) / (expected %)
- **Report:** Observed %, expected %, enrichment ratio, p-value

**ORF boundary exon susceptibility (A5SS/A3SS):**
- **Test:** Two-sample proportion test
  - Group 1: ORF boundary exons (contains start or stop)
  - Group 2: Regular CDS exons (no ORF boundary)
  - Proportion: % with A5SS or A3SS event
- **Alternative:** Fisher's exact test (2×2 table)
- **Effect size:** Risk ratio = (prop in boundary exons) / (prop in regular CDS)
- **Report:** Proportions, risk ratio, 95% CI, p-value

#### **5D. Pattern Classification (Category 4)**

**4A. Complexity-Event Relationship:**
- **Model fitting:**
  - Linear: events = β₀ + β₁ × complexity
  - Polynomial: events = β₀ + β₁ × complexity + β₂ × complexity²
  - Log-linear: events = β₀ + β₁ × log(complexity)
- **Model selection:** AIC, BIC, adjusted R²
- **Diagnostics:** Residual plots, normality tests
- **Report:** Best-fit model, R², coefficients, p-values

**4B-4D. Complexity-Controlled Pattern Analysis:**
- **Stratification approach:**
  - Define complexity bins (e.g., quartiles of primary complexity metric)
  - Test pattern frequencies within each bin
  - Compare pattern frequencies across bins (Chi-square or Fisher's exact)
- **Regression approach:**
  - Logistic regression: pattern_type ~ complexity + other_covariates
  - Test if patterns remain significant after controlling for complexity
- **Report:** Pattern frequencies by bin, OR adjusted for complexity, p-values

#### **5E. Multiple Testing Correction**

**Strategy:**
- **Within-category correction:** Apply FDR correction (Benjamini-Hochberg) within each category
  - Category 1: All co-occurrence tests
  - Category 2: All spatial tests
  - Category 3: All regional tests
- **Global correction:** Optional additional FDR correction across all tests
- **Threshold:** FDR q < 0.05 for significance
- **Report:** Both raw p-values and FDR q-values

#### **5F. Effect Sizes and Reporting Standards**

**Always report:**
- **Sample sizes:** n per group, n per stratum
- **Effect sizes:** OR, enrichment ratio, risk ratio (not just p-values)
- **Confidence intervals:** 95% CI for effect sizes
- **Power considerations:** Note when tests are underpowered (small n)

**Handling sparse data:**
- **Minimum n threshold:** Require n ≥ 10 per cell for Chi-square, otherwise use Fisher's exact
- **Flag low-power tests:** Note when n < 30 per group
- **Pooling strategy:** Consider pooling rare event types or complexity bins if sparse

---

## **SECTION 6: IMPLEMENTATION PLAN**

### **Overview**

Practical workflow for implementing the analysis framework defined in Sections 1-5.

**Version 6.0 Strategy:**
- Build all structures fresh from original source data
- No dependencies on Version 5.0 intermediate files
- Can reference Version 5.0 code for implementation patterns

**Estimated total time:** 3-5 hours (depending on dataset size and computational resources)

---

### **PHASE 1: Data Preparation** (60-90 min)

**Goal:** Build analysis-ready structures from original source data

#### **Step 1.1: Load and Prepare DGEList**

**Input:** `/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds`

**Script:** `01_prepare_dge_data.R`

**Process:**
1. Load DGEList object
2. Extract sample metadata from `$samples` slot
3. Calculate expression levels (CPM or TPM)
4. Identify dominant isoform per gene:
   - Group by gene_id
   - Rank isoforms by mean expression across samples
   - Flag dominant (highest expression)
5. Filter to relevant samples/conditions if needed (e.g., DMSO only)

**Output:**
- `expression_data.rds` (isoform × sample matrix with CPM/TPM)
- `dominant_isoforms.rds` (gene_id, dominant_isoform_id, mean_expression)
- `sample_metadata.rds` (sample info from DGEList)

---

#### **Step 1.2: Extract Isoform Structures from GFF Files**

**Input:**
- GENCODE: `/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz`
- SQANTI: `/Users/petecastaldi/claude_projects/nmd/reference_files/sqanti3_corrected.cds.gff3`
- PacBio GFF (if needed): Path TBD

**Script:** `02_extract_isoform_structures.R`

**Process:**
1. Parse GENCODE GFF3:
   - Extract exon features for transcripts
   - Store: transcript_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end
   - Filter to transcripts present in DGEList

2. Parse SQANTI/PacBio GFF:
   - Extract exon features for PacBio isoforms
   - Same structure as GENCODE

3. Combine into unified isoform structure table

**Output:** `isoform_structures.rds`
```
Columns: isoform_id, gene_id, seqnames, strand, exon_number, exon_start, exon_end
```

---

#### **Step 1.3: Build Union Exon Structure per Gene**

**Input:** `isoform_structures.rds`

**Script:** `03_build_union_exons.R`

**Process:**
1. For each gene:
   - Extract all exons from all isoforms
   - Build **atomic** union exons (not merged). The algorithm collects all unique exon boundaries, classifies each as START, END, or BOTH, and splits into atomic segments:
     - START boundary → split before: [..., b-1] then [b, ...]
     - END boundary → split after: [..., b] then [b+1, ...]
     - BOTH boundary → three segments: [..., b-1], [b, b], [b+1, ...]
   - This produces the smallest segments that are either fully included or fully excluded in any isoform. Each atomic segment can be claimed by any isoform exon via strict containment (UE_start >= exon_start AND UE_end <= exon_end).
   - Create union exons with unique IDs
   - Order union exons 5' → 3' (respect strand)
   - Assign union exon indices

2. Create isoform-union exon mapping:
   - For each isoform, determine which union exons are present
   - Binary matrix: (isoform_id, union_exon_id) → 1 (present) or 0 (absent)

**Output:**
- `union_exons.rds` (gene_id, union_exon_id, exon_index, seqnames, start, end, strand)
- `isoform_union_mapping.rds` (isoform_id, union_exon_id, present)

---

#### **Step 1.4: Extract CDS Annotations**

**Input:**
- GENCODE GFF3
- SQANTI GFF3

**Script:** `04_extract_cds_annotations.R`

**Process:**
1. Parse GENCODE GFF3:
   - Extract CDS features for each transcript
   - Store: transcript_id, cds_start, cds_stop

2. Parse SQANTI GFF3:
   - Extract CDS predictions for PacBio isoforms
   - Store: isoform_id, cds_start, cds_stop, coding_potential

3. Determine coding status:
   - If has CDS coordinates → "coding"
   - If explicitly non-coding or no CDS → "non_coding"
   - If ambiguous → "unknown"

4. Calculate ORF length for coding isoforms

**Output:** `isoform_cds_metadata.rds`
```
Columns: isoform_id, coding_status, cds_start, cds_stop, orf_length
```

---

#### **Step 1.5: Annotate Union Exons with Region Type**

**Input:**
- `union_exons.rds`
- `isoform_union_mapping.rds`
- `isoform_cds_metadata.rds`

**Script:** `05_annotate_region_types.R`

**Process:**

For each (isoform, union_exon) pair where isoform contains union_exon:

1. Lookup isoform coding status
2. If non-coding or unknown → assign "non_coding" or "unknown"
3. If coding:
   - Get union exon coordinates (start, end)
   - Get CDS coordinates (cds_start, cds_stop)
   - Apply region type logic:
   ```
   if (exon_end < cds_start):
       region_type = "5'UTR"

   elif (exon_start > cds_stop):
       region_type = "3'UTR"

   elif (exon_start <= cds_start AND exon_end >= cds_start AND exon_end < cds_stop):
       region_type = "contains_orf_start"

   elif (exon_start > cds_start AND exon_start <= cds_stop AND exon_end >= cds_stop):
       region_type = "contains_orf_stop"

   elif (exon_start >= cds_start AND exon_end <= cds_stop):
       region_type = "CDS"

   else:
       region_type = "unknown"
   ```

**Output:** `isoform_union_exons_annotated.rds`
```
Columns: isoform_id, gene_id, union_exon_id, exon_index, seqnames, start, end, strand,
         present (1/0), region_type
```

**Validation:**
- Check distribution of region_types
- Verify coding isoforms have expected UTR/CDS pattern (5'UTR → CDS → 3'UTR)
- Flag unexpected "unknown" assignments

---

### **PHASE 2: Dominant vs Non-Dominant Comparison** (60-90 min)

**Goal:** Extract splicing choice profiles for each non-dominant isoform

#### **Step 2.1: Extract Splicing Choice Profiles**

**Input:**
- `isoform_union_exons_annotated.rds`
- `dominant_isoforms.rds`
- `isoform_cds_metadata.rds`

**Script:** `07_extract_splicing_profiles.R`

**Process:**

For each gene:
1. Identify dominant isoform
2. Get dominant's union exon pattern (which exons present, region types)
3. For each non-dominant isoform:
   - Get non-dominant's union exon pattern
   - **Compare patterns to extract:**

**A. TSS/TES Changes:**
- First exon comparison (TSS):
  - Same exon or different exon? (overlap yes/no)
  - If same exon:
    - Outside boundary: compare TSS positions (extended/internal/same)
    - Inside boundary: compare 3' end of first exon (extended/internal/same)
    - UTR impact, CDS impact (compare region_types)
    - ORF affected: does TSS move past start codon?
  - If different exon:
    - Different region (alternative promoter)

- Last exon comparison (TES):
  - Same logic as TSS

**B. Exon Inclusion:**
- Identify MISSING exons (present in dominant, absent in non-dominant)
- Identify ADDED exons (absent in dominant, present in non-dominant)
- Count by region_type
- Spatial pattern: consecutive runs vs distributed

**C. Boundary Modifications:**
- Identify union exons with different boundaries in dominant vs non-dominant
  - This requires comparing exact exon coordinates, not just presence/absence
  - A5SS: same 3' boundary, different 5' boundary
  - A3SS: same 5' boundary, different 3' boundary
- Count by region_type
- If both A5SS and A3SS: calculate topology (F2F, B2B, distributed)

**D. Dysfunction:**
- Identify IR, Partial_IR (requires junction-level analysis or intron retention detection)
- Count by region_type

**E. Spatial Metrics:**
- Distance from TSS to first internal change
- Distance from last internal change to TES
- Relative positions (0-1 along transcript)

**F. Complexity Metrics:**
- n_exons, n_junctions, transcript_length for non-dominant isoform

**Output:** `splicing_choice_profiles.rds`

**Structure:**
```
Columns:
gene_id | dominant_isoform_id | nondominant_isoform_id |

# TSS
tss_overlap | tss_outside_boundary | tss_inside_boundary |
tss_utr_impact | tss_cds_impact | tss_orf_affected |

# TES
tes_overlap | tes_outside_boundary | tes_inside_boundary |
tes_utr_impact | tes_cds_impact | tes_orf_affected |

# Exon inclusion
n_exons_removed | n_exons_added |
n_exons_removed_5utr | n_exons_removed_cds | n_exons_removed_3utr | ...
removal_pattern | addition_pattern |

# Boundary modifications
n_a5ss | n_a3ss |
n_a5ss_5utr | n_a5ss_cds | n_a5ss_3utr | ...
a5ss_a3ss_topology | a5ss_a3ss_distance_exons |

# Dysfunction
n_ir | n_partial_ir |
ir_in_cds | ir_in_utr | ...

# Spatial
total_changes | change_distribution |
distance_tss_to_first_change | distance_last_change_to_tes |

# Complexity
complexity_n_exons | complexity_n_junctions | complexity_length |
complexity_n_union_exons_in_gene
```

---

### **PHASE 3: Analysis Execution** (90-120 min)

**Input:** `splicing_choice_profiles.rds`

#### **Step 3.1: Category 4 First - Complexity Analysis** (15 min)

**Why first:** Informs stratification strategy for subsequent analyses

**Script:** `09_analyze_complexity_relationship.R`

**Process:**
1. For each complexity metric (n_exons, n_junctions, length):
   - Scatter plot: complexity vs total_changes
   - Fit models: linear, polynomial, log-linear
   - Select best model (R², AIC)

2. Select primary complexity metric (highest R²)

3. Define complexity bins for stratification:
   - Quartiles or custom bins based on distribution

**Outputs:**
- `results/complexity_relationship_results.rds`
- `figures/complexity_vs_events.pdf`
- `results/complexity_bins_definition.rds`

---

#### **Step 3.2: Category 1 - Co-occurrence Analysis** (45 min)

**Script:** `10_analyze_cooccurrence.R`

**Process:**

**Approach 1: Simple contingency tables (20 min)**
- For each event pair of interest:
  - Build 2×2 table
  - Fisher's exact or Chi-square test
  - Calculate OR, 95% CI, p-value
  - Store results

- Apply FDR correction across all tests

**Approach 2: Regression-based (25 min)**
- **Step 1:** Fit baseline models (per event type):
  - Prepare data: (isoform_id, union_exon_id presence/absence, event_occurs)
  - Use glmnet for LASSO/ridge regression:
    ```
    event_occurs ~ union_exon_id (as factor) + n_exons
    ```
  - Cross-validation for penalty selection
  - Extract coefficients (identify hotspot union exons)

- **Step 2:** Co-occurrence models:
  - For each event pair (A, B):
    ```
    glm (logistic) with robust SEs:
    has_event_B ~ has_event_A + union_exon_id + n_exons
    ```
    - Use sandwich package for clustered SEs (cluster by gene_id)
    - Test coefficient for has_event_A
    - Extract adjusted OR, 95% CI, p-value

- **Step 3:** Compare crude vs adjusted ORs
  - Create comparison table
  - Identify confounded vs true associations

**Outputs:**
- `results/cooccurrence_crude_results.tsv`
- `results/cooccurrence_adjusted_results.tsv`
- `results/cooccurrence_comparison.tsv`
- `figures/cooccurrence_heatmap.pdf` (crude vs adjusted)

---

#### **Step 3.3: Category 2 - Spatial Organization** (20 min)

**Script:** `11_analyze_spatial_patterns.R`

**Process:**

**2A. Positional Bias (10 min):**
- For each event type:
  - Extract relative positions (0-1 along transcript)
  - KS test vs uniform distribution
  - Report median position, D-statistic, p-value

**2B. Topology Enrichment (5 min):**
- Filter to isoforms with A5SS + A3SS
- Classify topology (F2F, B2B, distributed)
- Simulate expected frequencies:
  - Random placement given gene structure
  - 10,000 simulations
- Chi-square test: observed vs expected
- Report enrichment ratios

**2C. Proximity (5 min):**
- Calculate pairwise distances between co-occurring events
- Permutation test:
  - Shuffle event positions within isoforms
  - Recalculate distances
  - 10,000 permutations
- Compare observed vs expected distributions

**Outputs:**
- `results/positional_bias_results.tsv`
- `results/topology_enrichment_results.tsv`
- `results/proximity_analysis_results.tsv`
- `figures/spatial_patterns.pdf`

---

#### **Step 3.4: Category 3 - Functional Context** (15 min)

**Script:** `12_analyze_functional_context.R`

**Process:**

**3A. Regional Distribution (10 min):**
- For each event type:
  - Count by region_type
  - Calculate expected % based on region sizes
    - Need total length of each region type across all isoforms
  - Binomial test for each region
  - Report enrichment ratios

- ORF boundary exon susceptibility:
  - Extract A5SS/A3SS events
  - Classify union exons: ORF boundary vs regular CDS
  - Two-proportion test
  - Report risk ratio

**3B. ORF Impact (5 min):**
- TSS changes: % with tss_orf_affected == TRUE
- TES changes: % with tes_orf_affected == TRUE
- Summary by region type

**Outputs:**
- `results/regional_distribution_results.tsv`
- `results/orf_impact_summary.tsv`
- `figures/regional_enrichment.pdf`

---

#### **Step 3.5: Category 4 - Pattern Classification** (10 min)

**Script:** `13_analyze_patterns_and_report.R`

**Process:**

**4B-4D: Pattern Analysis (using complexity bins from 3.1):**
- Classify each profile:
  - Boundary-only: TSS/TES changes, no internal
  - Inclusion-only: exon changes, no TSS/TES
  - Modification-only: A5SS/A3SS, no exon changes
  - Dysfunction: IR/Partial_IR present
  - Combined: multiple categories

- Within each complexity bin:
  - Count pattern frequencies
  - Identify top combinations

- Across bins:
  - Chi-square test for pattern × complexity association
  - Report patterns that differ by complexity

**Outputs:**
- `results/pattern_frequencies_by_complexity.tsv`
- `results/top_patterns.tsv`
- `figures/pattern_classification.pdf`

---

### **PHASE 4: Visualization and Reporting** (30-45 min)

**Script:** `13_analyze_patterns_and_report.Rmd` (integrated with Step 3.5)

**Process:**
1. Load all result files
2. Generate summary statistics across all categories
3. Create integrated figures:
   - Co-occurrence heatmaps (crude vs adjusted)
   - Spatial distribution plots
   - Regional enrichment barplots
   - Complexity vs events scatter with best-fit line
   - Pattern frequency charts by complexity
   - Topology enrichment barplots

4. Compile tables:
   - Top co-occurring event pairs (both approaches)
   - Top patterns by complexity
   - Regional enrichment summary
   - ORF impact summary

5. Render HTML report with narrative interpretation

**Output:** `report/isoform_choice_analysis_report.html`

---

### **PHASE 5: Output Organization**

**Directory structure:**
```
/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/Version_6.0/
├── ANALYSIS_PLAN.md                     # This document
├── data/                                 # Prepared data structures
│   ├── expression_data.rds
│   ├── dominant_isoforms.rds
│   ├── sample_metadata.rds
│   ├── isoform_structures.rds
│   ├── union_exons.rds
│   ├── isoform_union_mapping.rds
│   ├── isoform_cds_metadata.rds
│   ├── isoform_union_exons_annotated.rds
│   └── splicing_choice_profiles.rds
├── results/                              # Analysis results
│   ├── complexity_relationship_results.rds
│   ├── complexity_bins_definition.rds
│   ├── cooccurrence_crude_results.tsv
│   ├── cooccurrence_adjusted_results.tsv
│   ├── cooccurrence_comparison.tsv
│   ├── positional_bias_results.tsv
│   ├── topology_enrichment_results.tsv
│   ├── proximity_analysis_results.tsv
│   ├── regional_distribution_results.tsv
│   ├── orf_impact_summary.tsv
│   ├── pattern_frequencies_by_complexity.tsv
│   └── top_patterns.tsv
├── figures/                              # Visualizations
│   ├── complexity_vs_events.pdf
│   ├── cooccurrence_heatmap.pdf
│   ├── spatial_patterns.pdf
│   ├── regional_enrichment.pdf
│   └── pattern_classification.pdf
├── scripts/                              # Implementation scripts
│   ├── 01_prepare_dge_data_v2.R
│   ├── 02_extract_isoform_structures.R
│   ├── 03_build_union_exons.R
│   ├── 04_extract_cds_annotations.R
│   ├── 05_annotate_region_types.R
│   ├── 06_filter_to_analysis_subset.R
│   ├── 07_extract_splicing_profiles.R
│   ├── 08_validate_reconstruction.R
│   ├── 09_analyze_complexity_relationship.R
│   ├── 10_analyze_cooccurrence.R
│   ├── 11_analyze_spatial_patterns.R
│   ├── 12_analyze_functional_context.R
│   ├── 13_analyze_patterns_and_report.R
│   ├── event_detection_functions.R          # Shared library
│   ├── reconstruction_functions.R           # Shared library
│   └── create_union_exons_and_junctions.R   # Shared library
└── report/                               # Final output
    └── isoform_choice_analysis_report.html
```

---

### **Execution Order Summary**

**Critical path:**
1. **Phase 1 (Data Preparation)** → must complete sequentially:
   - Steps 1.1 → 1.2 → 1.3 (isoform structures → union exons)
   - Steps 1.4 → 1.5 (CDS annotations → region type annotation)

2. **Phase 2 (Profile Extraction)** → depends on Phase 1 completion

3. **Within Phase 3:**
   - **Step 3.1 (complexity) FIRST** → informs stratification
   - Steps 3.2-3.5 can run in parallel after 3.1 completes

4. **Phase 4 (Reporting)** → integrates all Phase 3 results

---

## **DATA SOURCES SUMMARY**

**Version 6.0 Original Source Data:**

1. **DGEList Object:**
   - Path: `/Users/petecastaldi/claude_projects/nmd/rds/dge_isoform_nofilter_2026.2.7.rds`
   - Contains: Counts, sample metadata, normalization factors

2. **GENCODE GFF3:**
   - Path: `/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz`
   - Contains: Exon and CDS coordinates for ENST transcripts

3. **SQANTI GFF3:**
   - Path: `/Users/petecastaldi/claude_projects/nmd/reference_files/sqanti3_corrected.cds.gff3`
   - Contains: ORF predictions for PacBio isoforms

**No dependencies on Version 5.0 intermediate files**

---

## **IMPLEMENTATION NOTES**

### **Key Differences from Version 5.0**

**Version 5.0 (Event Detection):**
- Focus: Pairwise isoform comparisons (reference vs comparison)
- Output: Event list with arbitrary reference/comparison designations
- Level 1: Pairwise exon assignments
- Level 2: Cumulative exon assignments across all comparisons

**Version 6.0 (Dominant-Centric Choice Analysis):**
- Focus: Non-dominant vs dominant comparisons (fixed reference point)
- Output: Splicing choice profiles for each non-dominant
- No Level 1/Level 2 distinction (direct comparison to dominant)
- Builds union exons fresh, not reusing Version 5.0 structures

### **Code Reuse Strategy**

**Can adapt from Version 5.0:**
- GFF parsing functions (exon extraction, CDS extraction)
- Union exon construction logic
- Spatial distance calculations
- Statistical testing frameworks

**Must write new for Version 6.0:**
- Dominant isoform identification
- TSS/TES overlap classification (region usage framework)
- Splicing choice profile extraction (Section 4 structure)
- Complexity-controlled co-occurrence testing (regression approach)

### **Future Enhancements**

**Stage 2 Extensions (after completing Stage 1):**
- Sequence-based frameshift analysis (requires FASTA)
- NMD prediction (PTC identification, >50-55bp rule)
- Protein domain impact analysis
- Functional annotation enrichment
- Cell-type-specific dominant isoforms

**Boundary Analysis:**
- Detailed A5SS/A3SS vs Partial_IR distance distributions
- Test if 100bp threshold is justified or arbitrary
- Continuous distance modeling

---

## **END OF ANALYSIS PLAN**

**Document Status:** Complete
**Ready for Implementation:** Yes
**Next Step:** Begin Phase 1 - Data Preparation
