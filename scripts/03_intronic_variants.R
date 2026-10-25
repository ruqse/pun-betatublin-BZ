#!/usr/bin/env Rscript
# ============================================================================
# Farm-exclusive intronic variants (Figures 4 and 5)
# ============================================================================
#
# Purpose:
#   Merge the VEP-annotated read-based variant calls of the resistant and
#   susceptible farms per gene, classify each variant as shared or exclusive
#   to one farm, and summarise the farm-exclusive intronic and
#   splice-proximal variants by consequence (Fig 4) and by intron position
#   along Pun-bt-B and Pun-bt-C (Fig 5). Also exports the exonic variant
#   calls used by scripts/04_coding_variants.R to cross-check haplotype
#   amino-acid changes against VEP.
#
# Inputs:
#   data/variants/reads_annot/Pun-bt-{A..G}_{Resistant,Susceptible}.reads.txt
#                                   Ensembl VEP v111 tab output per gene and farm
#   data/metadata/gene_map.tsv      gene label <-> target transcript ID
#
# Outputs:
#   figures/Fig4.tif (+ .png)   Figure 4: farm-exclusive intronic and
#                               splice-proximal variant counts per gene
#   figures/Fig5.tif (+ .png)   Figure 5: density of farm-exclusive intronic
#                               variant positions, A: Pun-bt-B, B: Pun-bt-C
#   results/intronic_variants/farm_exclusive_intronic_counts.csv
#                               gene, farm, consequence, n (Fig 4 data)
#   results/read_variants/exonic_variants.csv
#                               merged exonic VEP calls, input to script 04
#
# Usage (from the repository root):
#   Rscript scripts/03_intronic_variants.R
# ============================================================================

source(file.path("scripts", "functions.R"))

VEP_DIR        <- file.path(DATA_DIR, "variants", "reads_annot")
INTRON_OUT_DIR <- file.path("results", "intronic_variants")
READ_VAR_DIR   <- file.path("results", "read_variants")
FIG5_GENES     <- c("B", "C")

dir.create(INTRON_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(READ_VAR_DIR,   recursive = TRUE, showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Load and merge VEP calls        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gene_map <- load_gene_map()

# Each VEP file annotates several β-tubulin loci; only rows on the target
# gene's transcript (gene_map.tsv) are kept.
read_vep <- function(gene, farm) {
  info <- gene_map[gene_map$gene == gene, ]
  path <- file.path(VEP_DIR, sprintf("%s_%s.reads.txt", info$label, farm))
  readr::read_delim(path, delim = "\t", col_types = readr::cols(.default = "c")) %>%
    dplyr::mutate(Gene_name = gene, FarmStatus = farm) %>%
    dplyr::filter(Feature == info$transcript_id)
}

all_variants <- do.call(rbind, Map(read_vep, rep(GENES, each = length(FARMS)), FARMS))

# A variant (same position and allele) seen in both farms is "Shared".
merged_variants <- dplyr::full_join(
  all_variants %>% dplyr::filter(FarmStatus == "Resistant"),
  all_variants %>% dplyr::filter(FarmStatus == "Susceptible"),
  by = c("Gene_name", "Location", "Allele"),
  suffix = c(".res", ".sus")
) %>%
  dplyr::mutate(VariantStatus = dplyr::case_when(
    !is.na(FarmStatus.res) & !is.na(FarmStatus.sus) ~ "Shared",
    !is.na(FarmStatus.res) ~ "Resistant-only",
    TRUE ~ "Susceptible-only"
  ))

# One row per VEP consequence term: a variant annotated
# "splice_polypyrimidine_tract_variant,intron_variant" is counted once under
# each term in Fig 4 and contributes two points to the Fig 5 densities.
classified_variants <- merged_variants %>%
  dplyr::mutate(Consequence_type = dplyr::case_when(
    stringr::str_detect(EXON.res, "[0-9]+/[0-9]+") |
      stringr::str_detect(EXON.sus, "[0-9]+/[0-9]+") ~ "EXON",
    stringr::str_detect(INTRON.res, "[0-9]+/[0-9]+") |
      stringr::str_detect(INTRON.sus, "[0-9]+/[0-9]+") ~ "INTRON",
    TRUE ~ "OTHER"
  )) %>%
  tidyr::separate_rows(Consequence.res, Consequence.sus, sep = ",") %>%
  dplyr::mutate(dplyr::across(c(Consequence.res, Consequence.sus), trimws))

classified_variants %>%
  dplyr::filter(Consequence_type == "EXON") %>%
  readr::write_csv(file.path(READ_VAR_DIR, "exonic_variants.csv"))

# Farm-exclusive variants carry their annotation in that farm's columns only.
exclusive_variants <- classified_variants %>%
  dplyr::filter(VariantStatus != "Shared") %>%
  dplyr::mutate(
    farm        = ifelse(VariantStatus == "Resistant-only", "Resistant", "Susceptible"),
    consequence = ifelse(farm == "Resistant", Consequence.res, Consequence.sus),
    intron      = ifelse(farm == "Resistant", INTRON.res, INTRON.sus)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.4 - intronic variants       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

intronic_counts <- exclusive_variants %>%
  dplyr::filter(Consequence_type == "INTRON", !is.na(consequence)) %>%
  dplyr::count(gene = Gene_name, farm, consequence)

readr::write_csv(intronic_counts,
                 file.path(INTRON_OUT_DIR, "farm_exclusive_intronic_counts.csv"))

fig4 <- intronic_counts %>%
  dplyr::mutate(gene = factor(gene, levels = GENES)) %>%
  ggplot2::ggplot(ggplot2::aes(x = gene, y = n, fill = consequence)) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(~ farm, labeller = ggplot2::as_labeller(
    c(Resistant = "Resistant-only", Susceptible = "Susceptible-only"))) +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::labs(x = "Gene", y = "Variant count", fill = "Consequence") +
  theme_pub() +
  ggplot2::theme(legend.position = "bottom",
                 strip.text = ggplot2::element_text(face = "bold")) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))

save_fig(fig4, "Fig4", width = 7.5, height = 4.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.5 - intron positions        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# VEP INTRON field is "intron_number/total_introns"; numbered 5' -> 3'.
intron_positions <- exclusive_variants %>%
  dplyr::filter(Gene_name %in% FIG5_GENES, !is.na(intron), intron != "-") %>%
  tidyr::separate(intron, into = c("intron_num", "total_introns"), sep = "/") %>%
  dplyr::mutate(intron_num = as.numeric(intron_num),
                total_introns = as.numeric(total_introns))

# Shared y-axis for both panels: highest stats::density() peak of any
# gene x farm group with more than one variant, plus 10%.
density_peaks <- vapply(
  split(intron_positions$intron_num, list(intron_positions$Gene_name, intron_positions$farm)),
  function(x) if (length(x) > 1) max(stats::density(x)$y) else 0,
  numeric(1)
)
y_limit <- 1.1 * max(density_peaks)

# x-axis spans the observed introns, capped at the gene's last intron.
plot_intron_density <- function(gene) {
  df <- intron_positions[intron_positions$Gene_name == gene, ]
  upper <- min(ceiling(max(df$intron_num)) + 0.5, max(df$total_introns) + 0.5)
  lower <- max(1, floor(min(df$intron_num)) - 0.5)
  ggplot2::ggplot(df, ggplot2::aes(x = intron_num, fill = farm)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_x_continuous(breaks = seq(2, floor(upper), by = 2),
                                limits = c(lower, upper)) +
    ggplot2::coord_cartesian(ylim = c(0, y_limit)) +
    scale_fill_farm() +
    ggplot2::labs(x = "Intron position", y = "Density", fill = "Farm") +
    theme_pub()
}

fig5 <- patchwork::wrap_plots(lapply(FIG5_GENES, plot_intron_density), nrow = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")

save_fig(fig5, "Fig5", width = 7.5, height = 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
