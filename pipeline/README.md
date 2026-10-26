# Processing pipeline (Nextflow)

This pipeline goes from PacBio HiFi amplicon reads to the files in `../data/haplotypes/`, the read-based VCFs and the haplotype CDS. It uses the tools, versions and parameters of the jobs that produced the published results. Each process runs in a pinned container.

| Process | Tool (container version) | Key parameters | Output |
|---|---|---|---|
| ALIGN | pbmm2 1.13.0 | `--preset CCS --sort`; reference = the six β-tubulin scaffolds | BAM |
| CALL | freebayes 1.3.2 | `--pooled-continuous --min-coverage 10 --min-alternate-count 5 --min-alternate-fraction 0.05 --limit-coverage 1000 --ploidy 2` | VCF |
| FILTER | bcftools 1.20 | `-e 'QUAL < 50'` | `reads_vcf/` |
| ANNOTATE (optional, `--run_vep`) | Ensembl VEP 111 | WBPS19 GFF3 + genome | `reads_annot/` |
| INDEX_FASTQ | samtools 1.9 | bgzip + `fqidx` | indexed reads |
| PBAA | pbaa 1.0.3 | `--min-read-qv 20 --max-reads-per-guide <normalised count> --min-cluster-read-count 5 --min-cluster-frequency 0.001 --max-uchime-score 1`; guide = genomic region of the gene | `pbaa/` |
| MAP_TO_GUIDE, ORIENT | minimap2 2.26, samtools 1.9 | `-ax asm5`; minus-strand haplotypes reverse-complemented | `haplotypes/` |
| EXTRACT_CDS | MAFFT 7.407 | `--adjustdirection` alignment with the reference CDS; haplotype bases at reference columns | `haplotype_cds/` |
| ALIGN_HAPLOTYPES | MAFFT 7.526 | `--anysymbol --bl 62 --op 1.53 --ep 0.123 --reorder --retree 2 --maxiterate 2 --nuc` | `identity_matrices/` |
| IDENTITY_MATRIX | Clustal Omega 1.2.4 | `--percent-id --full` | `identity_matrices/{gene}_matrix.txt` |
| CHIMERA_QC | vsearch 2.18.0 | `--uchime_denovo` on reads | `chimera_qc/` |

Versions come from the published output headers (`@PG` in the BAMs; `##source` and `##bcftools_filterVersion` in the VCFs) and from the containers the published jobs used.

## Input reads

Input is the 14 per-gene normalised HiFi FASTQ files from ENA PRJEB102621, named as in `../data/metadata/samplesheet.csv`.

`bin/normalise.sh` documents how they were made from the demultiplexed reads. For each gene, the farm with more reads was down-sampled with `seqkit sample -s 123 -n <reads of the other farm>`, and the other farm was copied unchanged. Run on the demultiplexed files, the script reproduces 12 of the 14 deposited read sets exactly (identical read-name checksums). It does not reproduce Pun-bt-E Susceptible (105,376 reads) or Pun-bt-F Susceptible (63,881 reads). Those two were down-sampled earlier the same day by a command that was not recorded, and none of the SeqKit variants tested reproduces them. ENA holds those exact files, so every downstream step can still be reproduced from public data.

## Running

```bash
# Six-scaffold reference (pbmm2/freebayes), rebuilt from WormBase ParaSite
bash pipeline/bin/prepare_reference.sh data/reference

# Local (Singularity/Apptainer) or SLURM (UPPMAX profile)
nextflow run pipeline/main.nf --reads_dir /path/to/fastq
nextflow run pipeline/main.nf -profile uppmax --project <slurm-account> --reads_dir /path/to/fastq

# Quick end-to-end test on the smallest gene (Pun-bt-C, 13,295 reads per farm)
nextflow run pipeline/main.nf -profile uppmax,test --project <slurm-account> --reads_dir /path/to/fastq
```

Other profiles:
- `-profile docker` uses Docker instead of Singularity.
- `--run_vep --vep_gff <WBPS19 GFF3> --vep_fasta <WBPS19 genome>` adds VEP annotation.

## Checking a run against the published files

```bash
python3 pipeline/tests/compare_to_published.py results_pipeline
```

The script compares each output with the published files in `../data/` and prints PASS or FAIL:

| Check | What must match |
|---|---|
| Haplotypes (headers and sequences) | `../data/haplotypes/*_haplotypes.fasta` |
| Identity matrices | `../data/haplotypes/*_matrix.txt` |
| CDS | `../data/haplotype_cds/`, and their translations the curated protein alignments in `../data/cds_alignments/` |
| VCF records | `../data/variants/vcf/` |

It also reports how many reads vsearch flagged as chimeric.

## Verification

Re-runs on Pelle, checked with `tests/compare_to_published.py`:

| Output | Gene | Result |
|---|---|---|
| VCF records | Pun-bt-A | identical (329 Resistant, 130 Susceptible) |
| VCF records | Pun-bt-C | identical (397 Resistant, 125 Susceptible) |
| Haplotype sequences | Pun-bt-C | identical (15 Resistant, 36 Susceptible) |
| Identity matrix | Pun-bt-C | identical (51 × 51) |
| CDS | Pun-bt-C | identical (51 of 51) |
| vsearch read screen | all seven, both farms | 0 reads flagged |

freebayes 1.3.2 and pbaa 1.0.3 therefore reproduce the published files.

The published job scripts for the other genes run the same commands; only the gene, farm, FASTQ, reference and read cap differ. With those fields masked, all 14 read-based scripts are identical, and so are the pbaa scripts for A, B, C, F and G. The pbaa scripts for Pun-bt-D and Pun-bt-E were edited after their published outputs were made; the pipeline uses the gene-region guide for them, as for the other genes.

## Notes

- **Chimera screens.** The read-level vsearch screen runs on reads without abundance annotations, as the published job did. The cluster-level chimera filter that drove downstream exclusion is pbaa's uchime score.
- **VEP.** The published VEP tables came from the WormBase ParaSite VEP web tool. `ANNOTATE` is an offline equivalent and is not part of the published record.
