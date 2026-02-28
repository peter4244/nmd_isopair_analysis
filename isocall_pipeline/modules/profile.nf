process PROFILE {
    tag "${sample}"
    publishDir "${params.outdir}/profiles", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    path "${sample}.gz", emit: profile

    script:
    """
    isocall profile \\
        --reads ${bam} \\
        --sample ${sample} \\
        --output ${sample}.gz
    """
}
