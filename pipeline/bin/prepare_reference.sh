#!/usr/bin/env bash
# =============================================================================
# prepare_reference.sh - Build the 6-scaffold mapping reference used by pbmm2/freebayes
# =============================================================================
# Usage: bash pipeline/bin/prepare_reference.sh [outdir]
#
# Downloads the P. univalens ASM225920v1 genome from WormBase ParaSite (WBPS17; the sequences of these six scaffolds were checked identical, by md5, to the reference used)
# and extracts the six scaffolds carrying the seven beta-tubulin genes
# (Pun-bt-E and Pun-bt-F share PgB04). Requires curl and samtools.
# Writes <outdir>/Betatublin_scaffolds_genomic.fa (+ .fai).
set -euo pipefail
OUTDIR="${1:-data/reference}"
URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/parascaris_univalens/PRJNA386823/parascaris_univalens.PRJNA386823.WBPS17.genomic.fa.gz"
SCAFFOLDS=(PgR003 PgR007 PgB04 PgB10 PgR045 PgE153)
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

curl -fsSL "$URL" | gunzip -c > "$TMP/genome.fa"
samtools faidx "$TMP/genome.fa"
samtools faidx "$TMP/genome.fa" "${SCAFFOLDS[@]}" > "$OUTDIR/Betatublin_scaffolds_genomic.fa"
samtools faidx "$OUTDIR/Betatublin_scaffolds_genomic.fa"
cut -f1,2 "$OUTDIR/Betatublin_scaffolds_genomic.fa.fai"
