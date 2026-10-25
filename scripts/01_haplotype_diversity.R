#!/usr/bin/env Rscript
# ============================================================================
# Haplotype dominance per beta-tubulin gene and farm (n95 and top1)
# ============================================================================
#
# Purpose:
#   Summarises how concentrated the pbaa haplotype (cluster) pool is in each
#   gene x farm: n95, the number of haplotypes needed to reach 95% of reads,
#   and top1, the read frequency of the dominant haplotype.
#
# Inputs:
#   data/haplotypes/{A..G}_cluster_info.txt   pbaa cluster headers (farm, reads, UCHIME)
#   data/metadata/gene_map.tsv                gene labels
#
# Outputs:
#   figures/Fig1.tif (+ .png)                     Figure 1: n95 (A) and dominant haplotype frequency (B)
#   results/haplotype_diversity/dominance_summary.csv   gene, farm, n95, top1
#
# Usage (from the repository root):
#   Rscript scripts/01_haplotype_diversity.R
# ============================================================================

source(file.path("scripts", "functions.R"))

RESULTS_DIR <- file.path("results", "haplotype_diversity")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_LABELS <- with(load_gene_map(), stats::setNames(label, gene))
MIN_CLUSTER_FREQ <- 0.005

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Haplotype frequencies           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Chimeras removed first; tiny clusters (< 0.5% by the pbaa header frequency)
# are dropped within gene x farm, then frequencies are re-normalised by reads.
haplotypes <- dplyr::bind_rows(lapply(GENES, read_cluster_info)) %>%
  dplyr::group_by(gene, farm) %>%
  dplyr::filter(cluster_freq >= MIN_CLUSTER_FREQ) %>%
  dplyr::mutate(cluster_freq = read_count / sum(read_count)) %>%
  dplyr::ungroup()

dominance <- haplotypes %>%
  dplyr::group_by(gene, farm) %>%
  dplyr::arrange(gene, farm, dplyr::desc(cluster_freq)) %>%
  dplyr::mutate(cum = cumsum(cluster_freq), rank = dplyr::row_number()) %>%
  dplyr::summarise(n95 = min(rank[cum >= 0.95]),
                   top1 = max(cluster_freq),
                   .groups = "drop") %>%
  dplyr::mutate(gene = factor(GENE_LABELS[gene], levels = GENE_LABELS))

readr::write_csv(dominance, file.path(RESULTS_DIR, "dominance_summary.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Fig.1 - n95 and dominant haplotype   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dodge <- ggplot2::position_dodge(width = 0.7)

p_n95 <- ggplot2::ggplot(dominance, ggplot2::aes(gene, n95, fill = farm)) +
  ggplot2::geom_col(position = dodge, width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = n95), position = dodge,
                     vjust = -0.4, size = 3) +
  scale_fill_farm() +
  ggplot2::coord_cartesian(ylim = c(0, max(dominance$n95) * 1.12)) +
  ggplot2::labs(x = NULL, y = "Haplotypes to reach 95% of reads (n95)",
                fill = "Farm") +
  theme_pub()

p_top1 <- ggplot2::ggplot(dominance, ggplot2::aes(gene, top1, fill = farm)) +
  ggplot2::geom_col(position = dodge, width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = scales::percent(top1, accuracy = 1)),
                     position = dodge, vjust = -0.4, size = 3) +
  scale_fill_farm() +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::coord_cartesian(ylim = c(0, max(dominance$top1) * 1.12)) +
  ggplot2::labs(x = "Gene", y = "Dominant haplotype frequency", fill = "Farm") +
  theme_pub()

fig1 <- (p_n95 / p_top1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(legend.position = "bottom")

save_fig(fig1, "Fig1", width = 7.5, height = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
