#!/usr/bin/env Rscript
# ============================================================================
# Beta-tubulin expression and response to thiabendazole in adult worms
# ============================================================================
#
# Purpose:
#   (A) Constitutive expression (TPM) of the seven beta-tubulin genes in the
#       anterior end and intestine of untreated adult P. univalens.
#   (B) DESeq2 log2 fold change of the beta-tubulin genes in anterior-end
#       samples after 24 h in vitro exposure to thiabendazole (TBZ) at
#       1e-5, 1e-7 and 1e-9 M, each versus control.
#
# Inputs:
#   data/rnaseq/quant/<sample>/quant.sf.gz           salmon v0.11.3, TBZ experiment (24 lane-level runs)
#   data/rnaseq/TBZ_samples.txt                      TBZ sample sheet
#   data/rnaseq/quant_tissues/<sample>/quant.sf.gz   salmon v1.9.0, untreated anterior end and intestine (20 runs)
#   data/rnaseq/tissue_samples.txt                   tissue sample sheet
#   data/rnaseq/tx2gene.tsv                          WBPS17 transcript -> gene (scripts/07a_make_tx2gene.R)
#   data/metadata/gene_map.tsv                       gene labels
#
# Outputs:
#   figures/Fig6.tif (+ .png)                    Figure 6: TPM by tissue (A) and TBZ log2 fold change (B)
#   results/expression/beta_tubulin_tpm.csv      wbps_gene_id, label, sample, tissue, condition, tpm
#   results/expression/beta_tubulin_deseq2.csv   wbps_gene_id, label, contrast, baseMean, log2FoldChange, lfcSE, pvalue, padj
#
# Usage (from the repository root):
#   Rscript scripts/07_expression_tbz.R
# ============================================================================

source(file.path("scripts", "functions.R"))
check_packages(c("tximport", "DESeq2"))

RNASEQ_DIR  <- file.path(DATA_DIR, "rnaseq")
RESULTS_DIR <- file.path("results", "expression")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_MAP       <- load_gene_map()
GENE_LABELS    <- with(GENE_MAP, stats::setNames(label, wbps_gene_id))
CONCENTRATIONS <- c(CD = "Control", TBZ5 = "10^-5 M", TBZ7 = "10^-7 M", TBZ9 = "10^-9 M")
TX2GENE        <- utils::read.delim(file.path(RNASEQ_DIR, "tx2gene.tsv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Functions                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Gene-level counts, TPM and average transcript lengths; columns = samples.
import_salmon <- function(quant_dir, samples) {
  files <- stats::setNames(file.path(quant_dir, samples, "quant.sf.gz"), samples)
  stopifnot(all(file.exists(files)))
  tximport::tximport(files, type = "salmon", tx2gene = TX2GENE)
}

# A png device measures the legend text; Rscript's default pdf device cannot
# render the beta glyph and would leave an Rplots.pdf behind.
legend_grob <- function(p) {
  grDevices::png(tempfile(fileext = ".png"))
  on.exit(grDevices::dev.off())
  g <- ggplot2::ggplotGrob(p + ggplot2::theme(legend.position = "bottom"))
  g$grobs[[which(g$layout$name == "guide-box-bottom")]]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.6A - TPM by tissue          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tissue_meta <- utils::read.delim(file.path(RNASEQ_DIR, "tissue_samples.txt"))
tissue_txi  <- import_salmon(file.path(RNASEQ_DIR, "quant_tissues"), tissue_meta$sample)

tpm_mat <- tissue_txi$abundance[names(GENE_LABELS), ]
tpm <- data.frame(wbps_gene_id = rownames(tpm_mat), tpm_mat, check.names = FALSE) %>%
  tidyr::pivot_longer(-wbps_gene_id, names_to = "sample", values_to = "tpm") %>%
  dplyr::left_join(tissue_meta, by = "sample") %>%
  dplyr::transmute(wbps_gene_id,
                   label     = unname(GENE_LABELS[wbps_gene_id]),
                   sample,
                   tissue    = ifelse(type == "anterior", "Anterior", "Intestine"),
                   condition = unname(CONCENTRATIONS[concentration]),
                   tpm)

readr::write_csv(tpm, file.path(RESULTS_DIR, "beta_tubulin_tpm.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           TBZ differential expression     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Second-batch controls (P25051_*, expbatch 2) are excluded so that controls
# and TBZ-exposed worms share experiment and sequencing batch.
tbz_meta <- utils::read.delim(file.path(RNASEQ_DIR, "TBZ_samples.txt")) %>%
  dplyr::filter(!grepl("P25051", sample)) %>%
  dplyr::mutate(concentration = factor(concentration, levels = names(CONCENTRATIONS)))
rownames(tbz_meta) <- tbz_meta$sample

tbz_txi <- import_salmon(file.path(RNASEQ_DIR, "quant"), tbz_meta$sample)

# Each sequencing lane is its own column (n = 6 lanes from 3 worms per group),
# as in the original analysis: its collapseReplicates() call grouped by lane
# (intrareps), so no lanes were merged.
dds <- DESeq2::DESeqDataSetFromTximport(tbz_txi, colData = tbz_meta,
                                        design = ~ concentration)
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
dds  <- DESeq2::DESeq(dds[keep, ])

deseq2 <- dplyr::bind_rows(lapply(names(CONCENTRATIONS)[-1], function(conc) {
  contrast <- paste0("concentration_", conc, "_vs_CD")
  res <- as.data.frame(DESeq2::results(dds, name = contrast))
  ids <- intersect(names(GENE_LABELS), rownames(res))
  data.frame(wbps_gene_id = ids,
             label        = unname(GENE_LABELS[ids]),
             contrast     = contrast,
             res[ids, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")],
             row.names    = NULL)
}))
message("Beta-tubulin genes removed by the count filter: ",
        paste(setdiff(GENE_LABELS, deseq2$label), collapse = ", "))

readr::write_csv(deseq2, file.path(RESULTS_DIR, "beta_tubulin_deseq2.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.6 - assembly                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(123)
jitter <- ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
gene_colours <- ggplot2::scale_colour_brewer(palette = "Set1", limits = unname(GENE_LABELS),
                                             name = "\u03b2-tubulin")
no_strips <- ggplot2::theme(strip.text.x = ggplot2::element_blank(),
                            legend.title = ggplot2::element_text(face = "bold"))

p_tpm <- ggplot2::ggplot(tpm, ggplot2::aes(tissue, log(tpm + 1), colour = label)) +
  # Outliers hidden: every replicate is already drawn as a jittered point.
  ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
  ggplot2::geom_point(position = jitter, size = 0.8) +
  ggplot2::facet_wrap(~ tissue, scales = "free_x", ncol = 2) +
  gene_colours +
  ggplot2::labs(x = "Tissue type", y = "Transcript per million (log scale)") +
  theme_pub(base_size = 9) + no_strips +
  ggplot2::theme(legend.position = "none")

lfc <- deseq2 %>%
  dplyr::mutate(condition   = CONCENTRATIONS[sub("concentration_(.*)_vs_CD", "\\1", contrast)],
                significant = padj < 0.05)

p_lfc <- ggplot2::ggplot(lfc, ggplot2::aes(condition, log2FoldChange, colour = label)) +
  ggplot2::geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(shape = significant), position = jitter, size = 2.2) +
  ggplot2::facet_wrap(~ condition, scales = "free_x", ncol = 3) +
  gene_colours +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                              name = expression(bold(italic(p.adj) ~ "< 0.05"))) +
  ggplot2::scale_x_discrete(labels = function(x) parse(text = sub(" ", " * ", x))) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::guides(shape = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 2.5))) +
  ggplot2::labs(x = "Thiabendazole concentration", y = expression(Log[2] ~ "Fold Change")) +
  theme_pub(base_size = 9) + no_strips +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

# One gene-colour legend centred under both panels (taken from panel A, which
# has all seven genes; Pun-bt-G is absent from B); the p.adj legend under B.
colour_legend <- legend_grob(p_tpm + ggplot2::guides(
  colour = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 2.5))))
shape_legend <- legend_grob(p_lfc + ggplot2::guides(colour = "none"))

fig6 <- patchwork::wrap_plots(
  p_tpm,
  p_lfc + ggplot2::theme(legend.position = "none"),
  patchwork::wrap_elements(full = colour_legend, ignore_tag = TRUE),
  patchwork::wrap_elements(full = shape_legend, ignore_tag = TRUE),
  design = "AB\nCC\n#D", heights = c(10, 0.8, 0.8)
) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))

save_fig(fig6, "Fig6", width = 7.5, height = 4.8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
