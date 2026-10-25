#!/usr/bin/env bash
# =============================================================================
# normalise.sh - Per-gene read normalisation between the two farms
# =============================================================================
# Usage: bash pipeline/bin/normalise.sh <full_reads_dir> <out_dir> [samplesheet]
#
# For each gene, the farm with more reads is randomly down-sampled to the read
# count of the other farm with SeqKit (seed 123); the smaller farm is copied
# unchanged. These normalised reads are the ones deposited in ENA (PRJEB102621),
# so the pipeline normally starts from them. The full demultiplexed reads are
# not public; this script documents the step. It reproduces 12 of the 14
# deposited files exactly; Pun-bt-E and Pun-bt-F Susceptible were made by an
# unrecorded earlier command (see pipeline/README.md). Requires seqkit 2.4.0.
set -euo pipefail
IN="${1:?Usage: $0 <full_reads_dir> <out_dir> [samplesheet]}"
OUT="${2:?Usage: $0 <full_reads_dir> <out_dir> [samplesheet]}"
SHEET="${3:-data/metadata/samplesheet.csv}"
mkdir -p "$OUT"

count() { seqkit stats -T "$1" | awk 'NR == 2 {print $4}'; }

for gene in A B C D E F G; do
    r=$(awk -F, -v g="$gene" '$2 == g && $3 == "Resistant"   {print $4}' "$SHEET")
    s=$(awk -F, -v g="$gene" '$2 == g && $3 == "Susceptible" {print $4}' "$SHEET")
    nr=$(count "$IN/$r"); ns=$(count "$IN/$s")
    if (( nr <= ns )); then small="$r"; big="$s"; n="$nr"; else small="$s"; big="$r"; n="$ns"; fi
    cp "$IN/$small" "$OUT/$small"
    # One-pass -n (no -2): this exact call reproduces the deposited read sets.
    seqkit sample -s 123 -n "$n" "$IN/$big" -o "$OUT/$big"
    echo "Pun-bt-$gene: $small kept ($n reads); $big down-sampled to about $n reads"
done
