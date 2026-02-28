process MERGE {
    tag "merge"
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    path profiles

    output:
    path "merged.gz", emit: merged

    script:
    """
    isocall merge \\
        --profiles ${profiles} \\
        --output merged.gz
    """
}
