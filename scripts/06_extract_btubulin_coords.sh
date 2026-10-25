#!/usr/bin/env bash
# =============================================================================
# 06_extract_btubulin_coords.sh - scaffold coordinates of the seven beta-tubulin genes
# =============================================================================
# Usage: bash scripts/06_extract_btubulin_coords.sh <path/to/annotations.gff3>
#
# Reads the WormBase ParaSite gene IDs from data/metadata/gene_map.tsv, pulls
# each gene's scaffold, start, end and strand from the GFF3 "gene" features,
# and records host-scaffold length from the ##sequence-region directives.
# Assembly summary statistics are computed from the same directives.
#
# Annotation (P. univalens ASM225920v1, BioProject PRJNA386823, release WBPS19):
#   https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/parascaris_univalens/PRJNA386823/parascaris_univalens.PRJNA386823.WBPS19.annotations.gff3.gz
#   (gunzip before use)
#
# Outputs:
#   data/reference/btubulin_coords.tsv          pun_name, wormbase_gene_id, scaffold,
#                                               scaffold_length, start, end, strand
#   data/reference/assembly_stats.csv  n_scaffolds, total_mb, n50_mb, l50
# =============================================================================
set -euo pipefail
GFF="${1:?Usage: bash $0 <annotations.gff3>}"

GENE_MAP="data/metadata/gene_map.tsv"
COORDS="data/reference/btubulin_coords.tsv"
STATS="data/reference/assembly_stats.csv"

[[ -f "$GENE_MAP" ]] || { echo "Run from the repository root ($GENE_MAP not found)" >&2; exit 1; }
[[ -f "$GFF" ]] || { echo "GFF3 not found: $GFF" >&2; exit 1; }
mkdir -p "$(dirname "$STATS")"
tmp="$COORDS.tmp"
trap 'rm -f "$tmp"' EXIT

# --- Gene coordinates --------------------------------------------------------
# Labels come from gene_map.tsv (looked up by WormBase ID), never hardcoded.
{
  printf "pun_name\twormbase_gene_id\tscaffold\tscaffold_length\tstart\tend\tstrand\n"
  awk -F'\t' -v OFS='\t' '
    FNR == NR {
      if (FNR == 1) { for (i = 1; i <= NF; i++) col[$i] = i; next }
      label[$(col["wbps_gene_id"])] = $(col["label"]); n_wanted++; next
    }
    /^##sequence-region/ { split($0, f, " "); len[f[2]] = f[4] - f[3] + 1; next }
    $3 == "gene" && match($9, /ID=gene:[^;]+/) {
      id = substr($9, RSTART + 8, RLENGTH - 8)
      if (id in label) { n++; row[n] = label[id] OFS id OFS $1; scaf[n] = $1; rest[n] = $4 OFS $5 OFS $7 }
    }
    END {
      if (n != n_wanted) { print "Found " n " of " n_wanted " gene IDs in the GFF3" > "/dev/stderr"; exit 1 }
      for (i = 1; i <= n; i++) print row[i], len[scaf[i]], rest[i]
    }' "$GENE_MAP" "$GFF" | LC_ALL=C sort -t$'\t' -k3,3 -k5,5n
} > "$tmp"
mv "$tmp" "$COORDS"

# --- Assembly statistics -----------------------------------------------------
# N50 = length of the scaffold at which the length-sorted cumulative sum first
# reaches half the assembly; L50 = number of scaffolds needed to get there.
grep '^##sequence-region' "$GFF" | awk '{print $4 - $3 + 1}' | sort -nr | awk '
  { len[NR] = $1; total += $1 }
  END {
    for (i = 1; i <= NR; i++) { cum += len[i]; if (cum >= total / 2) break }
    print "n_scaffolds,total_mb,n50_mb,l50"
    printf "%d,%.6f,%.6f,%d\n", NR, total / 1e6, len[i] / 1e6, i
  }' > "$STATS"

echo "Wrote $COORDS"
column -t -s $'\t' "$COORDS"
echo "Wrote $STATS"
cat "$STATS"
