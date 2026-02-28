process CALL {
    tag "call"
    publishDir "${params.outdir}/call", mode: 'copy'

    input:
    path merged
    path isoforms
    path fasta
    path fai
    path config

    output:
    path "${params.output_prefix}.isoforms.gtf.gz", emit: gtf
    path "${params.output_prefix}.count_matrix.txt", emit: counts

    script:
    def config_arg = config.name != 'NO_CONFIG' ? "--config ${config}" : ''
    """
    isocall call \\
        --merged-profile ${merged} \\
        --known-isoforms ${isoforms} \\
        --reference ${fasta} \\
        --output-prefix ${params.output_prefix} \\
        --threads ${task.cpus} \\
        ${config_arg}
    """
}
