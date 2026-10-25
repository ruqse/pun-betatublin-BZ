#!/usr/bin/env Rscript
# ============================================================================
# Table 1 (genes and primers) and Table S1 (sequencing depth and normalisation)
# ============================================================================
#
# Purpose:
#   Table 1 lists the seven beta-tubulin genes with Roose et al. (2021) names,
#   WormBase ParaSite IDs, earlier names (Martin et al. 2021), amplicon primers
#   and sizes. Table S1 gives HiFi read counts per gene and farm before and
#   after per-gene normalisation (the farm with more reads was down-sampled to
#   the read count of the other farm; see pipeline/bin/normalise.sh).
#
# Inputs:
#   data/metadata/gene_map.tsv      gene names, IDs, primers, amplicon sizes
#   data/metadata/read_counts.tsv   seqkit stats of the demultiplexed and the
#                                   normalised (ENA) FASTQ files
#
# Outputs:
#   tables/Table1.csv               Table 1
#   tables/S1_Table.csv             Table S1
#   results/sequencing/read_summary.csv   run-level totals quoted in the Results
#
# Usage (from the repository root):
#   Rscript scripts/08_tables_metadata.R
# ============================================================================

source(file.path("scripts", "functions.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Table 1 - genes and primers     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gene_map <- load_gene_map()

table1 <- gene_map %>%
  dplyr::transmute(
    `Nomenclature (Roose et al., 2021)` = label,
    `WormBase ParaSite ID`              = wbps_gene_id,
    `Name in Martin et al. (2021)`      = tbb_name_martin2021,
    `Forward primer`                    = primer_F,
    `Reverse primer`                    = primer_R,
    `Amplicon size (kb)`                = amplicon_kb
  )
write_table(table1, "Table1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Table S1 - read depth           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

counts <- readr::read_tsv(file.path(DATA_DIR, "metadata", "read_counts.tsv"), show_col_types = FALSE)

s1 <- counts %>%
  dplyr::select(gene, farm, reads_demultiplexed, reads_normalised) %>%
  tidyr::pivot_wider(names_from = farm, values_from = c(reads_demultiplexed, reads_normalised)) %>%
  dplyr::left_join(gene_map %>% dplyr::select(gene, label), by = "gene") %>%
  dplyr::mutate(
    reads_normalised_lower_farm = pmin(reads_demultiplexed_Resistant, reads_demultiplexed_Susceptible),
    fold_reduction_higher_farm  = round(pmax(reads_demultiplexed_Resistant, reads_demultiplexed_Susceptible) /
                                          reads_normalised_lower_farm, 2)
  ) %>%
  dplyr::select(gene = label,
                reads_demultiplexed_Resistant, reads_demultiplexed_Susceptible,
                reads_normalised_Resistant, reads_normalised_Susceptible,
                reads_normalised_lower_farm, fold_reduction_higher_farm)
write_table(s1, "S1_Table")

# Down-sampling with seqkit sample -n is approximate, so the down-sampled farm
# can end a few reads short of the target; the shortfall is reported here.
stopifnot(all(s1$reads_normalised_Resistant <= s1$reads_demultiplexed_Resistant),
          all(s1$reads_normalised_Susceptible <= s1$reads_demultiplexed_Susceptible))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Run-level totals (Results)      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.create(file.path("results", "sequencing"), showWarnings = FALSE, recursive = TRUE)
read_summary <- data.frame(
  demultiplexed_reads      = sum(counts$reads_demultiplexed),
  per_barcode_min          = min(counts$reads_demultiplexed),
  per_barcode_max          = max(counts$reads_demultiplexed),
  per_barcode_mean         = mean(counts$reads_demultiplexed),
  per_barcode_sd           = stats::sd(counts$reads_demultiplexed),
  susceptible_over_resistant = mean(counts$reads_demultiplexed[counts$farm == "Susceptible"]) /
                               mean(counts$reads_demultiplexed[counts$farm == "Resistant"]),
  normalised_median        = stats::median(s1$reads_normalised_lower_farm),
  normalised_mean          = mean(s1$reads_normalised_lower_farm),
  normalised_sd            = stats::sd(s1$reads_normalised_lower_farm)
)
readr::write_csv(read_summary, file.path("results", "sequencing", "read_summary.csv"))
print(read_summary)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
