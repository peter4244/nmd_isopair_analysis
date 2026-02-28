process PREP_ISOFORMS {
    tag "prep_isoforms"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path gtf

    output:
    path "ref_seq.isoforms.gz", emit: isoforms

    script:
    // isocall requires gzip-compressed GTF input
    def gtf_gz = gtf.name.endsWith('.gz') ? gtf : "${gtf}.gz"
    def compress_cmd = gtf.name.endsWith('.gz') ? "" : "bgzip -c ${gtf} > ${gtf_gz}"
    """
    ${compress_cmd}
    isocall prep-isoforms \\
        --gtf ${gtf_gz} \\
        --output ref_seq.isoforms.gz
    """
}
