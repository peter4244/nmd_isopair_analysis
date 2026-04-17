# GEO — PROTOCOLS & DATA PROCESSING sections

Paste each block into the matching field of the GEO `seq_template.xlsx`
spreadsheet (one spreadsheet per SubSeries).

---

## SHORT-READ SUBSERIES — PROTOCOLS

**Growth protocol** (TODO — confirm / expand with Randell lab wording)
Primary human lung cell populations were obtained from deceased donor lungs
and expanded in culture at the Randell lab (UNC). AT2 cells were cultured as
TODO. Differentiated airway epithelial cells were cultured either in submerged
conditions (DD) or at an air-liquid interface (DD_ALI). Undifferentiated donor
airway epithelial cells (DO) were maintained in expansion medium. Lung
fibroblasts (FB) and lung microvascular endothelial cells (MV) were cultured in
TODO medium. See [TODO: reference to published culture protocol].

**Treatment protocol**
Cultured cells were incubated with 0.3 μM of the SMG1 inhibitor SMG1i or an
equivalent volume of DMSO (vehicle control) in culture medium at 37°C for
6 hours prior to RNA harvest. Matched DMSO and SMG1i samples were collected
from the same donor within each cell type.

**Extracted molecule**
polyA RNA

**Extraction protocol**
After washing with PBS, cells were lysed in RLT Plus buffer and total RNA
was extracted with the RNeasy Plus Mini Kit (Qiagen) per the manufacturer's
instructions.

**Library construction protocol**
Illumina stranded mRNA libraries were prepared from polyA-selected RNA using
the NEBNext UltraExpress RNA Library Prep Kit (New England Biolabs).

**Library strategy**
RNA-Seq

---

## SHORT-READ SUBSERIES — DATA PROCESSING PIPELINE

(These lines go into the "data processing" rows of the GEO spreadsheet, one
line per row.)

1. RNA-seq data were processed through the nf-core/rnaseq pipeline (<https://nf-co.re/rnaseq/>) version 3.14.0 using Nextflow 24.04.4.
2. Quality control: read quality and trimming were performed with Trim Galore 0.6.7 (cutadapt 4.6) with default parameters. Alignment quality and other QC metrics were generated and evaluated with FastQC 0.12.1, samtools 1.17, and MultiQC. Samples had 81–244 million raw read pairs (median 121 million) and 96–98% alignment rate to GRCh38 (median 97%).
3. Alignment: stranded reads were aligned to the human GRCh38 reference genome using the STAR aligner 2.7.10a (two-pass mode) and the GENCODE 49 GTF with default parameters. Junctional reads were assigned according to default parameters.
4. Quantification: isoform counts were estimated with Salmon 1.10.1 with default parameters. Gene- and isoform-level counts were generated with the tximport Bioconductor package (bioconductor-tximeta 1.12.0, r-base 4.1.3).
5. Genome_build: GRCh38 primary assembly (GENCODE release 49).
6. Supplementary_files_format_and_content: `salmon.merged.gene_counts_length_scaled.tsv` is a tab-separated matrix of length-scaled Salmon gene-level counts, rows = Ensembl gene IDs (versioned), columns = samples; `sample_metadata_shortread.csv` gives per-sample phenotype (cell type, donor, treatment, library size).

---

## LONG-READ SUBSERIES — PROTOCOLS

**Growth protocol**
Same as short-read (see above).

**Treatment protocol**
Same as short-read (see above).

**Extracted molecule**
polyA RNA

**Extraction protocol**
After washing with PBS, cells were lysed in RLT Plus buffer and total RNA
was extracted with the RNeasy Plus Mini Kit (Qiagen) per the manufacturer's
instructions.

**Library construction protocol**
cDNA PacBio Kinnex libraries were constructed from isolated RNA according to
PacBio protocol 103-238-700 REV07 with the following modifications: 1.0X
SMRTbell cleanup beads were used instead of 0.9X to avoid enrichment of long
transcripts in order to achieve properly formed concatenated libraries. Prior
to concatenation, each of the eight Kinnex PCR reactions was pooled by equal
mass based on its Qubit concentration, instead of including a fixed volume
from each tube. The Kinnex FL RNA libraries were sequenced on a PacBio Revio
instrument with Revio SMRT cells and SPRQ chemistry to an average depth of
[fill in] FLNC reads.

**Library strategy**
OTHER (PacBio Iso-Seq / Kinnex long-read mRNA sequencing)

---

## LONG-READ SUBSERIES — DATA PROCESSING PIPELINE

1. PacBio HiFi reads were generated on the Revio instrument and demultiplexed per sample into IsoSeqX-barcoded unaligned BAM files.
2. Long reads were aligned to the human genome (`GRCh38.primary_assembly.genome.fa`) with minimap2 using the command `minimap2 -ax splice:hq --junc-bed known_junctions.bed --MD --cs=long -uf GRCh38_splice.mmi input.fastq.gz | samtools sort -o out.bam` [1].
3. Transcript-level quantification was performed on the transcriptome-aligned BAM files using the PacBio isocall pipeline (<https://github.com/PacificBiosciences/isocall>) with parameters `-j 20 --seq-tech pac-bio-hifi --filter-group no-filters --model-coverage`. Novel isoforms are represented in the output as `<ENSG_id>.novel<N>`; reference transcripts retain their `ENST…` IDs.
4. Per-sample transcript-level counts were compiled into a transcript × sample matrix; gene-level counts were obtained by summing across transcripts per gene.
5. Genome_build: GRCh38 primary assembly (GENCODE release 49).
6. Supplementary_files_format_and_content: `isocall_transcript_counts.tsv` is a tab-separated matrix of integer transcript-level read counts (rows = transcript IDs including `ENST…` reference transcripts and `ENSG….novel<N>` novel isoforms, columns = samples). `isocall_gene_counts.tsv` is the gene-level aggregation. `isocall_isoforms.gtf.gz` is the isocall-called isoform annotation (reference + novel). `sample_metadata_longread.csv` gives per-sample phenotype.

[1] Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34:3094–3100.

---

## TODO items

Still outstanding before submission:
- Fill culture medium / kit details for each cell type (short-read Growth protocol).
- Fill average FLNC read depth for the long-read libraries.

