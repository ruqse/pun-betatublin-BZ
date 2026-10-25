/*
 * Read-based variant calling. Parameters are those of the published run
 * (job scripts dated 2025-03-20; versions read from the BAM @PG and VCF headers).
 */

process ALIGN {
    tag "${meta.id}"
    container 'quay.io/biocontainers/pbmm2:1.13.0--h9ee0642_0'
    cpus 8

    input:
    tuple val(meta), path(fastq)
    path reference

    output:
    tuple val(meta), path("${meta.id}.reads.sorted.bam"), emit: bam

    script:
    """
    pbmm2 align ${reference} ${fastq} ${meta.id}.reads.sorted.bam \\
        --preset CCS --sort --num-threads ${task.cpus}
    """
}

process CALL {
    tag "${meta.id}"
    // The published VCF headers read "freeBayes v1.3.2-dirty".
    container 'quay.io/biocontainers/freebayes:1.3.2--py38h40864fe_2'

    input:
    tuple val(meta), path(bam)
    path reference

    output:
    tuple val(meta), path("${meta.id}.reads.vcf"), emit: vcf

    script:
    """
    freebayes \\
        --fasta-reference ${reference} \\
        --pooled-continuous \\
        --min-coverage 10 \\
        --min-alternate-count 5 \\
        --min-alternate-fraction 0.05 \\
        --limit-coverage 1000 \\
        --ploidy 2 \\
        ${bam} > ${meta.id}.reads.vcf
    """
}

process FILTER {
    tag "${meta.id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/reads_vcf", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.reads.filtered.vcf.gz"), path("${meta.id}.reads.filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools filter -e 'QUAL < 50' ${vcf} -Oz -o ${meta.id}.reads.filtered.vcf.gz
    bcftools index --tbi ${meta.id}.reads.filtered.vcf.gz
    """
}

process PREPARE_VEP_GFF {
    container 'quay.io/biocontainers/htslib:1.20--h81da01d_0'

    input:
    path gff

    output:
    tuple path("annotation.sorted.gff3.gz"), path("annotation.sorted.gff3.gz.tbi"), emit: gff

    script:
    """
    (zcat -f ${gff} | grep '^#' ; zcat -f ${gff} | grep -v '^#' | sort -k1,1 -k4,4n -k5,5n -t\$'\\t') \\
        | bgzip -c > annotation.sorted.gff3.gz
    tabix -p gff annotation.sorted.gff3.gz
    """
}

process ANNOTATE {
    tag "${meta.id}"
    // The published annotation used the WormBase ParaSite VEP web tool (release 19,
    // Ensembl 111). This offline equivalent uses VEP 111 with the WBPS19 GFF3.
    container 'quay.io/biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0'
    publishDir "${params.outdir}/reads_annot", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(gff), path(gff_tbi)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.reads.txt"), emit: txt

    script:
    """
    vep --input_file ${vcf} --format vcf --gff ${gff} --fasta ${fasta} \\
        --tab --no_stats --force_overwrite --symbol --canonical --numbers --hgvs \\
        --fields "Uploaded_variation,Location,Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,REF_ALLELE,UPLOADED_ALLELE,DISTANCE,STRAND,FLAGS,SYMBOL_SOURCE,HGNC_ID,CANONICAL" \\
        --output_file ${meta.id}.reads.txt
    """
}
