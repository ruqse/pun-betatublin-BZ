#!/usr/bin/env bash
# =============================================================================
# run_all.sh - Regenerate every figure and table of the paper from data/
# =============================================================================
# Usage (from the repository root, inside the envs/r-figures.yml environment):
#   bash scripts/run_all.sh
# Scripts 06_extract_btubulin_coords.sh and 07a_make_tx2gene.R are one-off
# helpers that need the WormBase ParaSite GFF3; their outputs are in data/.
set -euo pipefail
mkdir -p results/logs
[[ -f scripts/functions.R ]] || { echo "Run from the repository root" >&2; exit 1; }

for s in 01_haplotype_diversity 02_metahaplotype_clustering 03_intronic_variants \
         04_coding_variants 05_gst_pcoa 06_scaffold_map 07_expression_tbz \
         08_tables_metadata; do
    echo "== ${s}"
    Rscript "scripts/${s}.R" > "results/logs/${s}.log" 2>&1 \
        || { echo "FAILED: ${s} (see results/logs/${s}.log)" >&2; exit 1; }
done
Rscript scripts/qa/check_manuscript_numbers.R
