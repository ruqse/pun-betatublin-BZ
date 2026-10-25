#!/usr/bin/env nextflow
/*
 * ============================================================================
 * Parascaris univalens beta-tubulin amplicon pipeline (PacBio HiFi)
 * ============================================================================
 *
 * Re-implements, with the same tools, versions and parameters, the SLURM jobs
 * that produced the published results:
 *   reads      : pbmm2 -> freebayes (pooled) -> bcftools QUAL >= 50 -> [VEP]
 *   haplotypes : pbaa cluster -> orient to guide -> rename -> CDS extraction
 *                -> MAFFT + Clustal Omega percent-identity matrix per gene
 *   qc         : vsearch --uchime_denovo on reads
 *
 * Input reads are the per-gene normalised HiFi reads deposited in ENA
 * (PRJEB102621); see bin/normalise.sh for how they were made.
 *
 * Usage:
 *   nextflow run pipeline/main.nf -profile uppmax \
 *       --reads_dir /path/to/fastq --project <slurm-account> --outdir results_pipeline
 * ============================================================================
 */

nextflow.enable.dsl = 2

include { ALIGN; CALL; FILTER; PREPARE_VEP_GFF; ANNOTATE } from './modules/reads.nf'
include { INDEX_FASTQ; PBAA; MAP_TO_GUIDE; ORIENT; EXTRACT_CDS; ALIGN_HAPLOTYPES; IDENTITY_MATRIX } from './modules/haplotypes.nf'
include { CHIMERA_QC } from './modules/qc.nf'

log.info """\
    P. univalens beta-tubulin pipeline
    ==================================
    samplesheet : ${params.samplesheet}
    reads_dir   : ${params.reads_dir}
    reference   : ${params.reference}
    outdir      : ${params.outdir}
    vep         : ${params.run_vep}
    """.stripIndent()

def guide_for(meta) { file("${params.ref_dir}/Pun-bt-${meta.gene}_genomic_REF.fasta", checkIfExists: true) }
def cds_ref_for(meta) { file("${params.ref_dir}/Pun-bt-${meta.gene}_CDS_REF.fasta", checkIfExists: true) }

workflow {
    if (!params.reads_dir) error "Set --reads_dir to the directory holding the HiFi FASTQ files"
    if (workflow.profile.contains('uppmax') && !params.project) error "uppmax profile: set --project <SLURM account>"

    samples = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, gene: row.gene, farm: row.farm,
                        max_reads: row.max_reads_per_guide as int]
            tuple(meta, file("${params.reads_dir}/${row.fastq}", checkIfExists: true))
        }

    reference = file(params.reference, checkIfExists: true)

    // Read-based variant calling
    ALIGN(samples, reference)
    CALL(ALIGN.out.bam, reference)
    FILTER(CALL.out.vcf)
    if (params.run_vep) {
        PREPARE_VEP_GFF(file(params.vep_gff, checkIfExists: true))
        ANNOTATE(FILTER.out.vcf, PREPARE_VEP_GFF.out.gff, file(params.vep_fasta, checkIfExists: true))
    }

    // Cluster-based haplotypes
    INDEX_FASTQ(samples.map { meta, fq -> tuple(meta, fq, guide_for(meta)) })
    PBAA(INDEX_FASTQ.out.fastq)
    MAP_TO_GUIDE(PBAA.out.passed.map { meta, fa -> tuple(meta, fa, guide_for(meta)) })
    ORIENT(MAP_TO_GUIDE.out.sam)
    EXTRACT_CDS(PBAA.out.passed.map { meta, fa -> tuple(meta, fa, cds_ref_for(meta)) })

    // One identity matrix per gene from Resistant + Susceptible haplotypes
    per_gene = ORIENT.out.fasta
        .map { meta, fa -> tuple(meta.gene, meta.farm, fa) }
        .groupTuple(size: 2)
        .map { gene, farms, fas ->
            def byFarm = [farms, fas].transpose().collectEntries { [(it[0]): it[1]] }
            tuple(gene, byFarm.Resistant, byFarm.Susceptible)
        }
    ALIGN_HAPLOTYPES(per_gene)
    IDENTITY_MATRIX(ALIGN_HAPLOTYPES.out.aln)

    // Read-level chimera screen
    if (params.run_chimera_qc) CHIMERA_QC(samples)
}
