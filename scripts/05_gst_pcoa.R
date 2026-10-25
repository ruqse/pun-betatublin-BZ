#!/usr/bin/env Rscript
# ============================================================================
# Between-farm differentiation per beta-tubulin gene: Nei's G_ST and PCoA
# ============================================================================
#
# Purpose:
#   Tests whether Farm-R vs Farm-S differentiation at the candidate genes
#   (Pun-bt-B, Pun-bt-C) exceeds that at the five non-candidate paralogs
#   (Pun-bt-A, D, E, F, G), which serve as an internal reference for
#   farm-level population structure.
#   (1) Per-gene Nei's G_ST from read-weighted haplotype (pbaa cluster)
#       frequencies.
#   (2) Classical MDS / PCoA of each gene's pairwise identity matrix,
#       shown for Pun-bt-B, C, F and G.
#
# Inputs:
#   data/haplotypes/{A..G}_cluster_info.txt   pbaa cluster headers (farm, reads, UCHIME)
#   data/haplotypes/{A..G}_matrix.txt         Clustal Omega percent-identity matrices
#   data/metadata/gene_map.tsv                gene labels
#
# Outputs:
#   figures/S6_Fig.tif (+ .png)   Figure S6: PCoA of Pun-bt-B, C, F, G haplotypes
#   figures/S7_Fig.tif (+ .png)   Figure S7: Nei's G_ST per gene
#   tables/S3_Table.csv           Table S3: per-gene G_ST, H_T, H_within, H_R, H_S, reads
#
# Usage (from the repository root):
#   Rscript scripts/05_gst_pcoa.R
# ============================================================================

source(file.path("scripts", "functions.R"))

GENE_LABELS     <- with(load_gene_map(), stats::setNames(label, gene))
CANDIDATE_GENES <- c("B", "C")
NON_CANDIDATES  <- setdiff(GENES, CANDIDATE_GENES)
PCOA_GENES      <- c("B", "C", "F", "G")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Functions                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Nei's G_ST = (H_T - H_within) / H_T. Each pbaa cluster is one haplotype
# (allele); read counts are allele-count weights within each farm pool.
# H_within is the read-weighted mean of the two within-farm heterozygosities.
compute_gst <- function(info) {
  wide <- info %>%
    dplyr::group_by(farm, name) %>%
    dplyr::summarise(reads = sum(read_count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = farm, values_from = reads, values_fill = 0)
  stopifnot(all(FARMS %in% colnames(wide)))

  het <- function(counts) {
    p <- counts / sum(counts)
    1 - sum(p^2)
  }
  reads_R <- sum(wide$Resistant)
  reads_S <- sum(wide$Susceptible)
  H_R <- het(wide$Resistant)
  H_S <- het(wide$Susceptible)
  H_within <- reads_R / (reads_R + reads_S) * H_R +
              reads_S / (reads_R + reads_S) * H_S
  H_T <- het(rowSums(wide[, FARMS]))

  # H_S is the Susceptible-farm heterozygosity (not Nei's mean within-
  # population H_S, which is H_within here).
  data.frame(n_clusters = nrow(wide), H_T = H_T, H_within = H_within,
             H_R = H_R, H_S = H_S, Gst = (H_T - H_within) / H_T,
             reads_R = reads_R, reads_S = reads_S)
}

# Distance = 100 - percent identity, projected to two axes with cmdscale.
# Variance explained uses the positive eigenvalues only.
compute_pcoa <- function(identity, info) {
  common <- intersect(rownames(identity), info$name)
  stopifnot(length(common) >= 3)
  mds <- cmdscale(as.dist(100 - identity[common, common]), k = 2, eig = TRUE)
  pos_eig <- mds$eig[mds$eig > 0]
  points <- data.frame(name = rownames(mds$points),
                       PC1 = mds$points[, 1], PC2 = mds$points[, 2]) %>%
    dplyr::left_join(info[, c("name", "farm", "read_count")], by = "name")
  list(points = points, var_exp = pos_eig[1:2] / sum(pos_eig) * 100)
}

# All panels share one size scale so patchwork collects a single legend.
plot_pcoa <- function(pcoa, label, size_limits) {
  ggplot2::ggplot(pcoa$points,
                  ggplot2::aes(x = PC1, y = PC2, fill = farm, size = read_count)) +
    # shape 21 with a black stroke keeps yellow points visible on white
    ggplot2::geom_point(shape = 21, colour = "black", stroke = 0.4, alpha = 0.85) +
    scale_fill_farm(name = "Farm") +
    ggplot2::scale_size_area(name = "Reads", max_size = 6, limits = size_limits,
                             breaks = c(1e3, 1e4, 1e5), labels = scales::comma) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(size = 4, shape = 21, colour = "black", stroke = 0.4))) +
    ggplot2::labs(title = label,
                  x = sprintf("PCo1 (%.1f%%)", pcoa$var_exp[1]),
                  y = sprintf("PCo2 (%.1f%%)", pcoa$var_exp[2])) +
    theme_pub(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Load haplotypes                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Chimeric clusters (positive UCHIME score) are dropped by read_cluster_info().
infos <- stats::setNames(lapply(GENES, read_cluster_info), GENES)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Table S3 - G_ST per gene        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gst_table <- lapply(GENES, function(g) {
  cbind(gene = GENE_LABELS[[g]], compute_gst(infos[[g]]))
}) %>%
  dplyr::bind_rows() %>%
  dplyr::select(gene, n_clusters, reads_R, reads_S, H_T, H_within, H_R, H_S, Gst) %>%
  dplyr::arrange(dplyr::desc(Gst))

write_table(gst_table, "S3_Table")
print(gst_table, digits = 3, row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.S7 - G_ST per gene          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gst_plot_df <- gst_table %>%
  dplyr::mutate(gene = factor(gene, levels = GENE_LABELS[GENES]),
                candidate = gene %in% GENE_LABELS[CANDIDATE_GENES])

# Solid line: highest non-candidate G_ST (primary comparator).
# Dashed line: mean of Pun-bt-F and Pun-bt-G.
max_non_candidate <- max(gst_plot_df$Gst[gst_plot_df$gene %in% GENE_LABELS[NON_CANDIDATES]])
fg_mean <- mean(gst_plot_df$Gst[gst_plot_df$gene %in% GENE_LABELS[c("F", "G")]])
gst_b <- gst_plot_df$Gst[gst_plot_df$gene == GENE_LABELS[["B"]]]
message(sprintf("Highest non-candidate G_ST = %.3f; F/G mean = %.3f; Pun-bt-B / highest non-candidate = %.2f",
                max_non_candidate, fg_mean, gst_b / max_non_candidate))

p_gst <- ggplot2::ggplot(gst_plot_df, ggplot2::aes(x = gene, y = Gst, fill = candidate)) +
  ggplot2::geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  ggplot2::geom_hline(yintercept = max_non_candidate, linetype = "solid",
                      colour = "grey40", linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = fg_mean, linetype = "dashed",
                      colour = "grey60", linewidth = 0.35) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = FARM_COLORS[["Resistant"]],
                                        `FALSE` = "#BFC6C6"), guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::labs(x = NULL, y = expression(G[st]~"(between-farm)")) +
  theme_pub() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

save_fig(p_gst, "S7_Fig", width = 7.5, height = 4.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.S6 - PCoA, B C F G          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

size_limits <- c(1, max(unlist(lapply(infos[PCOA_GENES], `[[`, "read_count"))))

pcoa_plots <- lapply(PCOA_GENES, function(g) {
  pcoa <- compute_pcoa(read_identity_matrix(g), infos[[g]])
  message(sprintf("%s PCo1 = %.1f%%, PCo2 = %.1f%%",
                  GENE_LABELS[[g]], pcoa$var_exp[1], pcoa$var_exp[2]))
  plot_pcoa(pcoa, GENE_LABELS[[g]], size_limits)
})

p_pcoa <- patchwork::wrap_plots(pcoa_plots, ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")

save_fig(p_pcoa, "S6_Fig", width = 7.5, height = 7)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
