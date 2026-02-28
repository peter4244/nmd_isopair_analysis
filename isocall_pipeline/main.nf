#!/usr/bin/env nextflow

/*
 * isocall Pipeline
 *
 * Joint RNA isoform calling across PacBio IsoSeq samples
 * using isocall v0.15.0
 *
 * Steps:
 *   1. PREP_ISOFORMS  — GTF → ref_seq.isoforms.gz
 *   2. FAIDX           — genome.fa → genome.fa.fai
 *   3. PROFILE          — aligned.bam → sample.gz (per-sample)
 *   4. MERGE            — all sample.gz → merged.gz
 *   5. CALL             — merged + isoforms + ref → GTF + counts
 */

nextflow.enable.dsl=2

// Module imports
include { PREP_ISOFORMS } from './modules/prep_isoforms'
include { FAIDX }          from './modules/faidx'
include { PROFILE }        from './modules/profile'
include { MERGE }          from './modules/merge'
include { CALL }           from './modules/call'

// Validate required params
if (!params.samplesheet)    { error "Missing required parameter: --samplesheet" }
if (!params.genome_fasta)   { error "Missing required parameter: --genome_fasta" }
if (!params.annotation_gtf) { error "Missing required parameter: --annotation_gtf" }

workflow {
    // Parse samplesheet → channel of [sample, bam]
    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam, checkIfExists: true)) }

    // Reference files
    genome_fasta   = file(params.genome_fasta, checkIfExists: true)
    annotation_gtf = file(params.annotation_gtf, checkIfExists: true)

    // Optional isocall config
    isocall_config = params.isocall_config
        ? file(params.isocall_config, checkIfExists: true)
        : file('NO_CONFIG')

    // Step 1 & 2: Reference prep (parallel)
    PREP_ISOFORMS(annotation_gtf)
    FAIDX(genome_fasta)

    // Step 3: Per-sample profiling (parallel fan-out)
    PROFILE(samples_ch)

    // Step 4: Merge all profiles
    MERGE(PROFILE.out.profile.collect())

    // Step 5: Joint isoform calling
    CALL(
        MERGE.out.merged,
        PREP_ISOFORMS.out.isoforms,
        genome_fasta,
        FAIDX.out.fai,
        isocall_config
    )
}
