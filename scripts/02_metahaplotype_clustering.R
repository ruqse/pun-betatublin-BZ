#!/usr/bin/env Rscript
# ============================================================================
# Meta-haplotype clustering and farm enrichment per beta-tubulin gene
# ============================================================================
#
# Purpose:
#   Groups each gene's pbaa haplotypes into meta-haplotypes (Group-1, ...) by
#   complete-linkage clustering of pairwise distances (100 - percent identity),
#   with the cut height chosen by a fixed rule over natural-break, Otsu and
#   silhouette candidates. Tests each group for read-weighted enrichment in
#   Farm-R versus Farm-S (2x2 Fisher's exact test, BH within gene).
#
# Inputs:
#   data/haplotypes/{A..G}_cluster_info.txt   pbaa cluster headers (farm, reads, UCHIME, diversity)
#   data/haplotypes/{A..G}_matrix.txt         Clustal Omega percent-identity matrices
#   data/metadata/gene_map.tsv                gene labels
#
# Outputs:
#   figures/Fig2.tif (+ .png)      Figure 2: PCA of Pun-bt-B (A) and Pun-bt-C (B) haplotypes
#   figures/Fig3.tif (+ .png)      Figure 3: forest plot of differentiated groups, all genes
#   figures/S1_Fig.tif (+ .png)    Figure S1: Pun-bt-B dendrogram with candidate cut heights
#   figures/S2_Fig.tif (+ .png)    Figure S2: Pun-bt-B group distribution by farm
#   figures/S3_Fig.tif (+ .png)    Figure S3: Pun-bt-B mean pi by group and farm
#   figures/S4_Fig.tif (+ .png)    Figure S4: Pun-bt-C group distribution by farm
#   figures/S5_Fig.tif (+ .png)    Figure S5: Pun-bt-C mean pi by group and farm
#   results/metahaplotype/threshold_selection.csv   candidate and chosen cut heights per gene
#   results/metahaplotype/group_assignments.csv     haplotype -> group
#   results/metahaplotype/group_enrichment.csv      Fisher test, fold change, class per group
#   results/metahaplotype/pca_variance.csv          PC1/PC2 percent variance per gene
#   results/metahaplotype/group_diversity.csv       mean pbaa diversity (pi) per group x farm
#
# Usage (from the repository root):
#   Rscript scripts/02_metahaplotype_clustering.R
# ============================================================================

source(file.path("scripts", "functions.R"))
check_packages(c("cluster", "RColorBrewer"))

RESULTS_DIR <- file.path("results", "metahaplotype")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_LABELS <- with(load_gene_map(), stats::setNames(label, gene))

SILHOUETTE_MIN     <- 0.30
GAP_RATIO_MIN      <- 3.0
CONSENSUS_MAX_DIFF <- 5
FC_THRESHOLD       <- 2
FDR_THRESHOLD      <- 0.01
EXCLUSIVE_CAP      <- 5  # log2 position of farm-exclusive groups in Fig. 3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Functions                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Candidate cut heights and the hierarchical decision rule (Methods).
select_threshold <- function(dist_mat, hc) {
  d <- dist_mat[lower.tri(dist_mat)]

  # Natural break: upper edge of the largest gap between sorted unique distances.
  sorted <- sort(unique(round(d, 2)))
  gaps <- diff(sorted)
  natural <- sorted[-1][which.max(gaps)]
  gap_ratio <- if (length(gaps) > 0 && mean(gaps) > 0) max(gaps) / mean(gaps) else 1

  # Otsu: histogram break maximising between-class variance.
  h <- hist(d, breaks = 100, plot = FALSE)
  prob <- h$counts / sum(h$counts)
  cum_prob <- cumsum(prob)
  cum_mean <- cumsum(prob * h$mids)
  total_mean <- sum(prob * h$mids)
  between_var <- (total_mean * cum_prob - cum_mean)^2 / (cum_prob * (1 - cum_prob))
  between_var[is.nan(between_var)] <- 0
  otsu <- h$breaks[which.max(between_var) + 1]

  # Silhouette is scored on the natural and Otsu heights plus the 90th and
  # 95th percentiles of pairwise distance; the best-scoring height is kept.
  candidates <- unique(c(natural, otsu, stats::quantile(d, c(0.9, 0.95))))
  sil <- vapply(candidates, function(cut_h) {
    groups <- stats::cutree(hc, h = cut_h)
    k <- length(unique(groups))
    if (k > 1 && k < length(groups) && k <= 50) {
      mean(cluster::silhouette(groups, stats::as.dist(dist_mat))[, "sil_width"])
    } else {
      -1
    }
  }, numeric(1))
  has_sil <- any(sil > -1)
  silhouette <- if (has_sil) candidates[which.max(sil)] else NA_real_
  best_sil <- if (has_sil) max(sil) else -1

  if (has_sil && best_sil > SILHOUETTE_MIN) {
    method <- "Silhouette"; chosen <- silhouette
  } else if (gap_ratio > GAP_RATIO_MIN) {
    method <- "Natural Breaks"; chosen <- natural
  } else if (abs(natural - otsu) < CONSENSUS_MAX_DIFF) {
    method <- "Consensus (Natural + Otsu)"; chosen <- mean(c(natural, otsu))
  } else {
    method <- "Otsu's Method"; chosen <- otsu
  }

  data.frame(natural = natural, otsu = otsu, silhouette = silhouette,
             best_silhouette_width = best_sil, gap_ratio = gap_ratio,
             chosen_method = method, chosen_height = chosen)
}

# Read-weighted 2x2 Fisher test per group (group reads vs all other reads,
# Farm-R vs Farm-S), BH within gene. Fold change is a ratio of within-farm
# proportions with a pseudocount = smallest non-zero proportion in the gene
# (capped at 0.001), so farm-exclusive groups get a finite value.
test_enrichment <- function(info) {
  farm_totals <- vapply(FARMS, function(f) sum(info$read_count[info$farm == f]), numeric(1))
  stopifnot(all(farm_totals > 0))

  wide <- info %>%
    dplyr::group_by(group, farm) %>%
    dplyr::summarise(reads = sum(read_count), .groups = "drop") %>%
    dplyr::mutate(proportion = reads / farm_totals[farm]) %>%
    tidyr::pivot_wider(names_from = farm, values_from = c(reads, proportion),
                       values_fill = 0)

  min_prop <- min(c(wide$proportion_Resistant[wide$proportion_Resistant > 0],
                    wide$proportion_Susceptible[wide$proportion_Susceptible > 0],
                    0.001))

  wide$p_value <- vapply(seq_len(nrow(wide)), function(i) {
    stats::fisher.test(matrix(c(wide$reads_Resistant[i],
                                farm_totals[["Resistant"]] - wide$reads_Resistant[i],
                                wide$reads_Susceptible[i],
                                farm_totals[["Susceptible"]] - wide$reads_Susceptible[i]),
                              nrow = 2))$p.value
  }, numeric(1))

  # Strict > and < as in the analysis run for the paper (Methods text says
  # >= and <=); no group lies exactly on a boundary.
  wide %>%
    dplyr::mutate(
      fold_change = (proportion_Resistant + min_prop) / (proportion_Susceptible + min_prop),
      fdr = stats::p.adjust(p_value, method = "BH"),
      resistant_only = proportion_Resistant > 0 & proportion_Susceptible == 0,
      susceptible_only = proportion_Susceptible > 0 & proportion_Resistant == 0,
      class = dplyr::case_when(
        resistant_only ~ "Resistant-exclusive",
        susceptible_only ~ "Susceptible-exclusive",
        fold_change > FC_THRESHOLD & fdr < FDR_THRESHOLD ~ "Enriched",
        fold_change < 1 / FC_THRESHOLD & fdr < FDR_THRESHOLD ~ "Depleted",
        TRUE ~ "Not significant"))
}

analyse_gene <- function(gene) {
  info <- read_cluster_info(gene)
  identity <- read_identity_matrix(gene)
  keep <- rownames(identity) %in% info$name
  identity <- identity[keep, keep]
  stopifnot(setequal(rownames(identity), info$name))

  dist_mat <- 100 - identity
  hc <- stats::hclust(stats::as.dist(dist_mat), method = "complete")
  threshold <- select_threshold(dist_mat, hc)
  groups <- stats::cutree(hc, h = threshold$chosen_height)
  info$group <- paste0("Group-", groups[match(info$name, rownames(dist_mat))])
  threshold$n_groups <- length(unique(groups))

  # PCA treats each haplotype's distance profile as a variable vector.
  pca <- stats::prcomp(as.data.frame(dist_mat), scale. = TRUE)
  var_pct <- pca$sdev^2 / sum(pca$sdev^2) * 100
  pca_points <- data.frame(name = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2]) %>%
    dplyr::left_join(info[, c("name", "farm", "group")], by = "name")

  list(info = info, hc = hc, threshold = threshold,
       enrichment = test_enrichment(info), pca = pca_points,
       pca_variance = data.frame(PC1_pct = var_pct[1], PC2_pct = var_pct[2]))
}

# Data frame keyed by gene label, from one list element of every gene.
collect <- function(results, element) {
  dplyr::bind_rows(lapply(names(results), function(g) {
    df <- results[[g]][[element]]
    df$gene <- NULL
    cbind(gene = GENE_LABELS[[g]], df)
  }))
}

group_levels <- function(x) {
  u <- unique(x)
  u[order(as.integer(sub("Group-", "", u)))]
}

# Segment coordinates of an hclust tree, leaves at height 0 in plotting order.
dendrogram_segments <- function(hc) {
  n <- length(hc$order)
  leaf_x <- match(seq_len(n), hc$order)
  node_x <- numeric(n - 1)
  child <- function(j) if (j < 0) c(leaf_x[-j], 0) else c(node_x[j], hc$height[j])
  segments <- vector("list", n - 1)
  for (i in seq_len(n - 1)) {
    l <- child(hc$merge[i, 1])
    r <- child(hc$merge[i, 2])
    node_x[i] <- (l[1] + r[1]) / 2
    segments[[i]] <- data.frame(x = c(l[1], r[1], l[1]), xend = c(l[1], r[1], r[1]),
                                y = c(l[2], r[2], hc$height[i]),
                                yend = rep(hc$height[i], 3))
  }
  do.call(rbind, segments)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Run all genes and tables        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

results <- stats::setNames(lapply(GENES, analyse_gene), GENES)

threshold_selection <- collect(results, "threshold")
print(threshold_selection)
readr::write_csv(threshold_selection, file.path(RESULTS_DIR, "threshold_selection.csv"))

group_assignments <- collect(results, "info") %>%
  dplyr::select(gene, haplotype = name, farm, read_count, group)
readr::write_csv(group_assignments, file.path(RESULTS_DIR, "group_assignments.csv"))

enrichment <- collect(results, "enrichment")
readr::write_csv(enrichment %>%
                   dplyr::select(gene, group, reads_resistant = reads_Resistant,
                                 reads_susceptible = reads_Susceptible,
                                 fold_change, p_value, fdr, class) %>%
                   dplyr::arrange(gene, as.integer(sub("Group-", "", group))),
                 file.path(RESULTS_DIR, "group_enrichment.csv"))

readr::write_csv(collect(results, "pca_variance"), file.path(RESULTS_DIR, "pca_variance.csv"))

group_diversity <- collect(results, "info") %>%
  dplyr::group_by(gene, group, farm) %>%
  dplyr::summarise(pi = mean(diversity), .groups = "drop") %>%
  dplyr::arrange(gene, as.integer(sub("Group-", "", group)), farm)
readr::write_csv(group_diversity, file.path(RESULTS_DIR, "group_diversity.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.2 - PCA of B and C          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_pca <- function(gene) {
  res <- results[[gene]]
  pts <- res$pca
  lev <- group_levels(pts$group)
  pts$group <- factor(pts$group, levels = lev)
  group_colors <- stats::setNames(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(lev)), lev)
  ggplot2::ggplot(pts, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(ggplot2::aes(color = group, shape = farm), size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_shape_manual(values = c(Resistant = 16, Susceptible = 17)) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", res$pca_variance$PC1_pct),
                  y = sprintf("PC2 (%.1f%%)", res$pca_variance$PC2_pct),
                  color = "Haplotype group", shape = "Farm") +
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1, title.position = "top"),
                    color = ggplot2::guide_legend(order = 2, ncol = 2, title.position = "top")) +
    theme_pub(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
}

fig2 <- (plot_pca("B") | plot_pca("C")) + patchwork::plot_annotation(tag_levels = "A")
save_fig(fig2, "Fig2", width = 7.5, height = 4.2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Fig.3 - Forest plot, all genes       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Rows are ordered by effect size within each gene panel.
forest <- enrichment %>%
  dplyr::filter(class != "Not significant") %>%
  dplyr::mutate(
    effect_size = dplyr::case_when(resistant_only ~ EXCLUSIVE_CAP,
                                   susceptible_only ~ -EXCLUSIVE_CAP,
                                   TRUE ~ log2(fold_change)),
    is_exclusive = resistant_only | susceptible_only,
    direction = ifelse(effect_size > 0, "Enriched", "Depleted"),
    gene = factor(gene, levels = GENE_LABELS),
    row_id = stats::reorder(paste(gene, group, sep = "__"), effect_size))
stopifnot(nrow(forest) > 0)

fig3 <- ggplot2::ggplot(forest, ggplot2::aes(effect_size, row_id)) +
  ggplot2::geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray70",
                      linewidth = 0.4) +
  ggplot2::geom_vline(xintercept = 0, color = "gray20", linewidth = 0.8) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect_size, yend = row_id),
                        color = "gray50", linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(fill = direction, shape = is_exclusive),
                      size = 3, stroke = 0.6, color = "gray20") +
  ggplot2::facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  ggplot2::scale_y_discrete(labels = function(x) sub(".*__", "", x)) +
  ggplot2::scale_fill_manual(
    values = c(Enriched = FARM_COLORS[["Resistant"]], Depleted = FARM_COLORS[["Susceptible"]]),
    labels = c(Enriched = "Enriched in\nresistant farm", Depleted = "Depleted in\nresistant farm"),
    name = "Association") +
  ggplot2::scale_shape_manual(
    values = c("FALSE" = 21, "TRUE" = 23),
    labels = c("FALSE" = "Present in\nboth farms", "TRUE" = "Exclusive to\none farm"),
    name = "Farm status") +
  ggplot2::scale_x_continuous(
    breaks = c(-5, -3, -2, -1, 0, 1, 2, 3, 5),
    labels = c("≤32×", "8×", "4×", "2×", "1",
               "2×", "4×", "8×", "≥32×"),
    limits = c(-5.5, 5.5)) +
  ggplot2::labs(x = "Fold change in resistant vs susceptible farm (log₂ scale)", y = NULL) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(title.position = "top",
                                 override.aes = list(shape = 21, size = 4)),
    shape = ggplot2::guide_legend(title.position = "top",
                                  override.aes = list(fill = "white", size = 4))) +
  theme_pub(base_size = 9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold"),
                 axis.text.x = ggplot2::element_text(size = 6.5),
                 legend.position = "bottom")

save_fig(fig3, "Fig3", width = 7.5, height = 7)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      S1 Fig - Pun-bt-B dendrogram         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hc_b <- results[["B"]]$hc
thr_b <- results[["B"]]$threshold
info_b <- results[["B"]]$info
k_b <- thr_b$n_groups

leaves <- data.frame(x = seq_along(hc_b$order), name = hc_b$labels[hc_b$order]) %>%
  dplyr::left_join(info_b[, c("name", "farm", "group")], by = "name") %>%
  dplyr::mutate(label = sub(paste0("^", GENE_LABELS[["B"]], "-"), "", name))

# Group boxes drawn as in stats::rect.hclust: top midway between the k-th and
# (k-1)-th highest merges.
boxes <- leaves %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(xmin = min(x) - 0.4, xmax = max(x) + 0.4, .groups = "drop") %>%
  dplyr::mutate(ymin = -0.03, ymax = mean(rev(hc_b$height)[(k_b - 1):k_b]))

cut_lines <- data.frame(
  method = c(sprintf("Silhouette (selected): %.2f", thr_b$silhouette),
             sprintf("Natural breaks: %.2f", thr_b$natural),
             sprintf("Otsu: %.2f", thr_b$otsu)),
  height = c(thr_b$silhouette, thr_b$natural, thr_b$otsu))
cut_lines$method <- factor(cut_lines$method, levels = cut_lines$method)

s1 <- ggplot2::ggplot() +
  ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                     fill = NA, color = "red", linewidth = 0.4) +
  ggplot2::geom_segment(data = dendrogram_segments(hc_b),
                        ggplot2::aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
  ggplot2::geom_hline(data = cut_lines,
                      ggplot2::aes(yintercept = height, color = method, linetype = method),
                      linewidth = 0.6) +
  ggplot2::geom_point(data = leaves, ggplot2::aes(x = x, y = 0, fill = farm),
                      shape = 21, size = 1.8, stroke = 0.3) +
  ggplot2::scale_color_manual(values = c("red", "blue", "darkgreen"), name = "Cut height") +
  ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Cut height") +
  scale_fill_farm(name = "Farm") +
  ggplot2::scale_x_continuous(breaks = leaves$x, labels = leaves$label,
                              expand = ggplot2::expansion(add = 0.6)) +
  ggplot2::labs(x = NULL, y = "Height (100 - % identity)") +
  theme_pub(base_size = 9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
                 legend.position = "top", legend.box = "horizontal")

save_fig(s1, "S1_Fig", width = 7.5, height = 6.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   S2-S5 Fig - group distribution and pi   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


plot_distribution <- function(gene) {
  dat <- results[[gene]]$info %>%
    dplyr::group_by(farm) %>%
    dplyr::mutate(farm_reads = sum(read_count)) %>%
    dplyr::group_by(group, farm) %>%
    dplyr::summarise(proportion = 100 * sum(read_count) / dplyr::first(farm_reads),
                     .groups = "drop") %>%
    tidyr::complete(group, farm, fill = list(proportion = 0)) %>%
    dplyr::mutate(group = factor(group, levels = group_levels(group)))
  dodge <- ggplot2::position_dodge(width = 0.9)
  ggplot2::ggplot(dat, ggplot2::aes(group, proportion, fill = farm)) +
    ggplot2::geom_col(position = dodge) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", proportion)),
                       position = dodge, vjust = -0.5, size = 3) +
    scale_fill_farm() +
    ggplot2::labs(x = "Haplotype group", y = "Proportion of farm reads (%)", fill = "Farm") +
    theme_pub()
}

plot_diversity <- function(gene) {
  label <- GENE_LABELS[[gene]]
  dat <- group_diversity %>%
    dplyr::filter(gene == label) %>%
    dplyr::mutate(group = factor(group, levels = group_levels(group)))
  # A group absent from one farm has no pi; its single bar is centred.
  ggplot2::ggplot(dat, ggplot2::aes(group, pi, fill = farm)) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single")) +
    scale_fill_farm() +
    ggplot2::labs(x = "Haplotype group", y = "Mean cluster diversity (pi)", fill = "Farm") +
    theme_pub()
}

save_fig(plot_distribution("B"), "S2_Fig", width = 6.5, height = 4.5)
save_fig(plot_diversity("B"),    "S3_Fig", width = 6.5, height = 4.5)
save_fig(plot_distribution("C"), "S4_Fig", width = 6.5, height = 4.5)
save_fig(plot_diversity("C"),    "S5_Fig", width = 6.5, height = 4.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
