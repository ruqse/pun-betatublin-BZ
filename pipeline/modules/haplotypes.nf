/*
 * Cluster-based haplotypes. Parameters are those of the published run
 * (job scripts dated 2025-03-28; identity matrices 2025-04-02).
 */

process INDEX_FASTQ {
    tag "${meta.id}"
    container 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'

    input:
    tuple val(meta), path(fastq), path(guide)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), path("${meta.id}.fastq.gz.fai"), path("${meta.id}.fastq.gz.gzi"), path(guide), path("${guide}.fai"), emit: fastq

    script:
    // pbaa needs a bgzip-indexed FASTQ and an indexed guide.
    """
    gunzip -c ${fastq} | bgzip -i -I ${meta.id}.fastq.gz.gzi > ${meta.id}.fastq.gz
    samtools fqidx ${meta.id}.fastq.gz
    samtools faidx ${guide}
    """
}

process PBAA {
    tag "${meta.id}"
    container 'quay.io/biocontainers/pbaa:1.0.3--hdfd78af_0'
    publishDir "${params.outdir}/pbaa", mode: 'copy', pattern: '*{passed,failed}_cluster_sequences.fasta'

    input:
    tuple val(meta), path(fastq), path(fai), path(gzi), path(guide), path(guide_fai)

    output:
    tuple val(meta), path("${meta.id}.consensus_pbaa_passed_cluster_sequences.fasta"), emit: passed
    path "${meta.id}.consensus_pbaa_failed_cluster_sequences.fasta"
    path "${meta.id}.consensus_pbaa_read_info.txt"

    script:
    // --max-reads-per-guide = the normalised per-gene read count, so all reads are used.
    """
    echo ${fastq} > ${meta.id}.fofn
    pbaa cluster ${guide} ${meta.id}.fofn ${meta.id}.consensus_pbaa \\
        --min-read-qv 20 \\
        --max-reads-per-guide ${meta.max_reads} \\
        --min-cluster-read-count 5 \\
        --min-cluster-frequency 0.001 \\
        --max-uchime-score 1 \\
        --log-level INFO
    """
}

process MAP_TO_GUIDE {
    tag "${meta.id}"
    container 'quay.io/biocontainers/minimap2:2.26--he4a0461_2'

    input:
    tuple val(meta), path(passed), path(guide)

    output:
    tuple val(meta), path(passed), path("${meta.id}.guide.sam"), emit: sam

    script:
    """
    minimap2 -ax asm5 ${guide} ${passed} > ${meta.id}.guide.sam
    """
}

process ORIENT {
    tag "${meta.id}"
    container 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'
    publishDir "${params.outdir}/haplotypes", mode: 'copy'

    input:
    tuple val(meta), path(passed), path(sam)

    output:
    tuple val(meta), path("${meta.id}.consensus_oriented_sequences.fasta"), emit: fasta

    script:
    // Reverse-complement haplotypes that map to the guide's minus strand (SAM flag 16),
    // iterating records in coordinate order as the published job did.
    // pbaa names clusters sample-<barcodes>_guide-<scaffold>_cluster-N_...; downstream
    // files use Pun-bt-<gene>-<farm>_cluster-N_...
    """
    samtools sort -o sorted.bam ${sam}
    samtools view sorted.bam > records.sam
    orient_haplotypes.awk records.sam ${passed} \\
        | sed -E 's/^>sample-[^_]+_guide-[^_]+_/>Pun-bt-${meta.gene}-${meta.farm}_/' \\
        > ${meta.id}.consensus_oriented_sequences.fasta
    test "\$(grep -c '^>' ${meta.id}.consensus_oriented_sequences.fasta)" -eq "\$(grep -c '^>' ${passed})"
    """
}

process EXTRACT_CDS {
    tag "${meta.id}"
    // Published CDS came from run_build_cluster_cds.sh (2025-03-31): MAFFT
    // --adjustdirection of the reference CDS plus all haplotypes of one sample,
    // then the haplotype bases at reference (non-gap) columns.
    container 'quay.io/biocontainers/mafft:7.407--1'
    publishDir "${params.outdir}/haplotype_cds", mode: 'copy'

    input:
    tuple val(meta), path(passed), path(cds_ref)

    output:
    tuple val(meta), path("${meta.id}_CDS.fasta"), emit: cds

    script:
    """
    export TMPDIR=\$PWD   # MAFFT needs a writable temp dir inside the container
    cat ${cds_ref} ${passed} > input.fasta
    mafft --quiet --adjustdirection input.fasta > aligned.fasta
    extract_cds.awk aligned.fasta \
        | sed -E 's/^>(_R_)?sample-[^_]+_guide-[^_]+_/>\1Pun-bt-${meta.gene}-${meta.farm}_/' \
        > ${meta.id}_CDS.fasta
    test "\$(grep -c '^>' ${meta.id}_CDS.fasta)" -eq "\$(grep -c '^>' ${passed})"
    """
}

process ALIGN_HAPLOTYPES {
    tag "Pun-bt-${gene}"
    container 'docker://staphb/mafft:7.526'
    cpus 4
    publishDir "${params.outdir}/identity_matrices", mode: 'copy'

    input:
    tuple val(gene), path(resistant), path(susceptible)

    output:
    tuple val(gene), path("Pun-bt-${gene}_cluster_align.fasta"), emit: aln

    script:
    """
    export TMPDIR=\$PWD
    cat ${resistant} ${susceptible} > Pun-bt-${gene}_cluster.fasta
    mafft --thread ${task.cpus} --anysymbol --bl 62 --op 1.53 --ep 0.123 --reorder \\
        --retree 2 --treeout --maxiterate 2 --nuc Pun-bt-${gene}_cluster.fasta \\
        > Pun-bt-${gene}_cluster_align.fasta
    """
}

process IDENTITY_MATRIX {
    tag "Pun-bt-${gene}"
    container 'docker://staphb/clustalo:1.2.4'
    publishDir "${params.outdir}/identity_matrices", mode: 'copy'

    input:
    tuple val(gene), path(aln)

    output:
    path "${gene}_matrix.txt", emit: matrix

    script:
    """
    clustalo --infile=${aln} --percent-id --full --distmat-out=${gene}_matrix.txt > /dev/null
    """
}
