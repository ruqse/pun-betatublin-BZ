#!/usr/bin/env Rscript
# ============================================================================
# Genomic context of the seven beta-tubulin genes in ASM225920v1
# ============================================================================
#
# Purpose:
#   Draws each beta-tubulin gene on its host scaffold (to scale) to show which
#   genes share a scaffold, and reports each host scaffold's length relative
#   to the assembly N50. Gene labels are joined from gene_map.tsv by WormBase
#   ParaSite gene ID.
#
# Inputs:
#   data/reference/btubulin_coords.tsv          gene coordinates (scripts/06_extract_btubulin_coords.sh)
#   data/reference/assembly_stats.csv  assembly N50 (same script)
#   data/metadata/gene_map.tsv                  gene labels
#
# Outputs:
#   figures/S8_Fig.tif (+ .png)   Figure S8: scaffold map of the beta-tubulin genes
#
# Usage (from the repository root):
#   bash scripts/06_extract_btubulin_coords.sh <annotations.gff3>   # optional, regenerates inputs
#   Rscript scripts/06_scaffold_map.R
# ============================================================================

source(file.path("scripts", "functions.R"))
check_packages(c("ggrepel"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Load coordinates                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

assembly <- readr::read_csv(file.path("data", "reference", "assembly_stats.csv"),
                            show_col_types = FALSE)

coords <- readr::read_tsv(file.path(DATA_DIR, "reference", "btubulin_coords.tsv"),
                          show_col_types = FALSE) %>%
  dplyr::left_join(load_gene_map() %>% dplyr::select(label, wbps_gene_id),
                   by = c("wormbase_gene_id" = "wbps_gene_id"))
# Guards against a coordinates table written with a different nomenclature.
stopifnot(!anyNA(coords$label), identical(coords$pun_name, coords$label))

# Scaffolds ordered by number of beta-tubulin genes, then by first gene start.
scaffold_order <- coords %>%
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(n_genes = dplyr::n(), min_start = min(start), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_genes), min_start)

coords <- coords %>%
  dplyr::group_by(scaffold) %>%
  dplyr::mutate(linkage = ifelse(dplyr::n() > 1, "Co-scaffolded", "Isolated")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(scaffold = factor(scaffold, levels = scaffold_order$scaffold),
                midpoint = (start + end) / 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.S8 - Scaffold map           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Thin grey line = full scaffold extent, so scaffold size differences
# (PgE153 ~20 kb vs PgR003 ~4.7 Mb) are visible next to each gene.
scaffold_extents <- dplyr::distinct(coords, scaffold, scaffold_length)

p_map <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = scaffold_extents,
                        ggplot2::aes(x = 0, xend = scaffold_length / 1e6,
                                     y = scaffold, yend = scaffold),
                        linewidth = 0.7, colour = "grey78", alpha = 0.9) +
  ggplot2::geom_segment(data = coords,
                        ggplot2::aes(x = start / 1e6, xend = end / 1e6,
                                     y = scaffold, yend = scaffold, colour = linkage),
                        linewidth = 3.4, alpha = 0.95) +
  ggrepel::geom_text_repel(data = coords,
                           ggplot2::aes(x = midpoint / 1e6, y = scaffold,
                                        label = label, colour = linkage),
                           size = 3.3, fontface = "bold", direction = "y",
                           nudge_y = 0.32, segment.size = 0.25, segment.alpha = 0.5,
                           min.segment.length = 0, seed = 123, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(`Co-scaffolded` = FARM_COLORS[["Resistant"]],
                                          Isolated = "#7A7A7A"), name = NULL) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.005, 0.02))) +
  ggplot2::labs(x = "Position on scaffold (Mb)", y = NULL) +
  theme_pub() +
  ggplot2::theme(legend.position = "bottom",
                 panel.grid.major.y = ggplot2::element_blank())

save_fig(p_map, "S8_Fig", width = 7.5, height = 0.6 + 0.6 * nrow(scaffold_order) + 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Linkage summary                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n50_bp <- assembly$n50_mb * 1e6
print(as.data.frame(assembly), row.names = FALSE)
coords %>%
  dplyr::arrange(label) %>%
  dplyr::transmute(label, wormbase_gene_id, scaffold, scaffold_length, start, end,
                   vs_n50 = ifelse(scaffold_length >= n50_bp, "above", "below")) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# Intergenic gap = next gene start - previous gene end on the same scaffold.
coords %>%
  dplyr::arrange(scaffold, start) %>%
  dplyr::group_by(scaffold) %>%
  dplyr::mutate(next_gene = dplyr::lead(label), gap_bp = dplyr::lead(start) - end) %>%
  dplyr::filter(!is.na(gap_bp)) %>%
  dplyr::select(scaffold, label, next_gene, gap_bp) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
