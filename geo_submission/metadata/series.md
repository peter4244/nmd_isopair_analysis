# GEO SuperSeries — SERIES section

This file contains the text blocks for the SERIES section of the GEO
`seq_template.xlsx` spreadsheet. Paste each block into the matching
field in the template for the **SuperSeries** and for each **SubSeries**.

---

## SuperSeries (overarching study)

**Series title**
Long read RNA sequencing in primary lung cell types reveals principles on nonsense-mediated decay.

**Summary (abstract)**
Nonsense-mediated decay (NMD) is a highly conserved mechanism in which RNA
transcripts are targeted for degradation after ribosomal engagement. One of the
functions of NMD is to "proofread" transcripts in order to identify and degrade
those that contain premature termination codons (PTCs). However, NMD can be
triggered for reasons beyond PTCs, including long 3' UTRs and the presence of
upstream open reading frames. In general, NMD-triggering events lead to stalled
or delayed translation, and while some of the RNA sequence characteristics that
cause NMD are known, many remain unknown. Consequently, our ability to predict
the extent to which specific transcripts will be degraded by NMD both
qualitatively (will a transcript be degraded by NMD?) and quantitatively (to
what extent will it be degraded?) remains imperfect, despite the widespread use
of heuristics such as the 50-nucleotide rule.

In addition to its role as transcript proofreader, NMD seems to have been
integrated into other biological processes as a means of rapid regulation of
protein levels. For example, NMD plays a major role in cellular stress, where
it provides baseline suppression of stress response genes. When cells
experience stress, NMD is inhibited, resulting in a rapid increase in stress
response RNAs, providing clear evidence that cells have utilized NMD as a
regulatory mechanism.

We hypothesized that, in addition to degrading PTC-containing mRNAs, NMD
regulates multiple biological processes in primary lung cell types. To study
this hypothesis, we performed transcriptome-wide discovery of NMD-targeted
transcripts in five different primary lung cell types using short and
long-read RNA sequencing in normal culture conditions and in the presence of
NMD inhibition. We identified thousands of transcripts that are degraded by
NMD, indicating that NMD is not limited to rare, aberrant transcripts but
rather is widespread and affects highly expressed genes and transcripts. We
observed significant enrichment of NMD affecting genes involved in the
unfolded protein response and alternative splicing, as well as more cell
type-specific pathways and processes. Using deep learning, we trained a
highly accurate model to predict NMD from transcript sequence, identifying
novel RNA sequence features responsible for NMD, including sequence elements
at the stop codon that promote readthrough translation.

**Overall design**
Six primary lung cell populations (alveolar type 2 cells [AT2], differentiated
airway epithelial cells in submerged culture [DD], differentiated airway
epithelial cells at air-liquid interface [DD_ALI], undifferentiated donor
airway epithelial cells [DO], lung fibroblasts [FB], and lung microvascular
endothelial cells [MV]) from independent human donors were treated for 6 hours
at 37°C with 0.3 μM of the SMG1 inhibitor SMG1i or with an equivalent volume
of DMSO vehicle control in culture medium. Total RNA was profiled with two
orthogonal platforms: (1) paired-end short-read RNA sequencing on Illumina
NovaSeq 6000 and (2) long-read
PacBio IsoSeq (Kinnex/MAS-Seq) on a Revio instrument. The study is organized as
a SuperSeries comprising two SubSeries, one per platform. Each donor
contributes a matched DMSO / SMG1i pair within a given cell type.

**Contributor(s)**
Castaldi,Peter,J

**Supplementary file(s)** (SuperSeries-level)
- sample_metadata_longread.csv  (full per-sample phenotype table, long-read)
- sample_metadata_shortread.csv (full per-sample phenotype table, short-read)

---

## SubSeries A — Short-read RNA sequencing

**Series title**
Short-read RNA sequencing of primary lung cell types under SMG1 inhibition

**Summary (abstract)**
Paired-end Illumina RNA sequencing of six primary lung cell populations
treated with the SMG1 inhibitor SMG1i or DMSO vehicle. Gene-level quantification
was performed with the nf-core/rnaseq pipeline (STAR alignment to GRCh38,
Salmon quantification against GENCODE v49). See SuperSeries for overall design
and background.

**Overall design**
See SuperSeries. This SubSeries contains only the short-read Illumina samples
(paired-end FASTQs) and the Salmon gene-level count matrix.

**Contributor(s)**
Castaldi,Peter,J

---

## SubSeries B — Long-read RNA sequencing

**Series title**
PacBio IsoSeq long-read RNA sequencing of primary lung cell types under SMG1 inhibition

**Summary (abstract)**
PacBio HiFi IsoSeq (Kinnex/MAS-Seq) long-read RNA sequencing of six primary lung
cell populations treated with the SMG1 inhibitor SMG1i or DMSO vehicle. Reads
were aligned to GRCh38 with minimap2 (splice:hq preset). Joint isoform discovery
and transcript-level quantification were performed with isocall against
GENCODE v49, including both reference-annotated and novel isoforms.
See SuperSeries for overall design and background.

**Overall design**
See SuperSeries. This SubSeries contains only the long-read PacBio samples
(demultiplexed per-sample HiFi unaligned BAMs, one BAM per sample) and the
isocall gene- and transcript-level count matrices.

**Contributor(s)**
Castaldi,Peter,J

---

## TODO items before submission
- Add co-authors to the Contributor(s) line if more than one.
