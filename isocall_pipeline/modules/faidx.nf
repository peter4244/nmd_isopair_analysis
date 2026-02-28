process FAIDX {
    tag "faidx"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai

    script:
    """
    samtools faidx ${fasta}
    """
}
