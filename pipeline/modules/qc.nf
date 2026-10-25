/*
 * Read-level chimera screen, as in the published QC job (2025-08-22).
 * Reads enter without dereplication or ;size= annotations, as in the published
 * job. The cluster-level filter that drove exclusion is pbaa's uchime score.
 */

process CHIMERA_QC {
    tag "${meta.id}"
    container 'quay.io/biocontainers/vsearch:2.18.0--h95f258a_0'
    cpus 8
    publishDir "${params.outdir}/chimera_qc", mode: 'copy', pattern: '*.{uchime.txt,log}'

    input:
    tuple val(meta), path(fastq)

    output:
    path "${meta.id}.uchime.txt"
    path "${meta.id}.vsearch.log"

    script:
    """
    vsearch --fastx_filter ${fastq} --fastq_qmax 93 --fastaout reads.fasta
    vsearch --uchime_denovo reads.fasta \\
        --nonchimeras /dev/null --chimeras ${meta.id}.chimeras.fasta \\
        --uchimeout ${meta.id}.uchime.txt --threads ${task.cpus} 2> ${meta.id}.vsearch.log
    """
}
