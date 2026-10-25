#!/usr/bin/env Rscript
# ============================================================================
# QA: regenerated results versus the numbers stated in the manuscript
# ============================================================================
#
# Purpose:
#   Reads the CSVs written by scripts 01-08 and compares each value quoted in
#   the manuscript text with the regenerated value, at the precision the
#   manuscript reports it. Prints PASS/FAIL per claim and exits 1 if any fail.
#   A FAIL means the code or the manuscript text needs fixing; never edit the
#   expected value here unless the manuscript text itself changes.
#
# Inputs:
#   results/**/*.csv, tables/*.csv
#
# Usage (from the repository root, after scripts/run_all.sh):
#   Rscript scripts/qa/check_manuscript_numbers.R
# ============================================================================

source(file.path("scripts", "functions.R"))

results <- list()
check <- function(claim, observed, expected, digits = NULL) {
  obs <- if (is.null(digits)) observed else round(observed, digits)
  if (is.numeric(obs)) obs <- as.vector(obs)
  ok  <- isTRUE(all.equal(obs, expected, check.attributes = FALSE))
  results[[length(results) + 1]] <<- data.frame(
    status = ifelse(ok, "PASS", "FAIL"), claim = claim,
    expected = paste(expected, collapse = "/"), observed = paste(obs, collapse = "/"))
}
read <- function(...) readr::read_csv(file.path(...), show_col_types = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.1 - n95 and dominance       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dom <- read("results", "haplotype_diversity", "dominance_summary.csv")
dom_val <- function(g, col) dom[[col]][match(paste0("Pun-bt-", g, c("Resistant", "Susceptible")),
                                             paste0(dom$gene, dom$farm))]
for (g in names(n95 <- list(A = c(9, 13), B = c(2, 19), C = c(12, 27), D = c(18, 21)))) {
  check(sprintf("n95 Pun-bt-%s Farm-R vs Farm-S", g), dom_val(g, "n95"), n95[[g]])
}
top1 <- list(B = c(95, 13), C = c(25, 8), F = c(52, 26), G = c(46, 37), A = c(39, 31))
for (g in names(top1)) {
  check(sprintf("Dominant haplotype %% Pun-bt-%s Farm-R vs Farm-S", g), 100 * dom_val(g, "top1"), top1[[g]], 0)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Fig.2/3, S2-S5 - meta-haplotypes     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pca <- read("results", "metahaplotype", "pca_variance.csv")
check("PCA Pun-bt-B PC1/PC2 %", unlist(pca[pca$gene == "Pun-bt-B", c("PC1_pct", "PC2_pct")]), c(96.0, 3.7), 1)
check("PCA Pun-bt-C PC1/PC2 %", unlist(pca[pca$gene == "Pun-bt-C", c("PC1_pct", "PC2_pct")]), c(74.3, 11.3), 1)

thr <- read("results", "metahaplotype", "threshold_selection.csv")
check("Pun-bt-B groups and silhouette threshold", unlist(thr[thr$gene == "Pun-bt-B", c("n_groups", "chosen_height")]), c(2, 1.4), 1)
check("Pun-bt-C groups and silhouette threshold", unlist(thr[thr$gene == "Pun-bt-C", c("n_groups", "chosen_height")]), c(4, 0.8), 1)

enr <- read("results", "metahaplotype", "group_enrichment.csv")
fc <- function(g, grp) enr$fold_change[enr$gene == paste0("Pun-bt-", g) & enr$group == grp]
check("Pun-bt-B Group-1 fold enrichment", fc("B", "Group-1"), 4.2, 1)
check("Pun-bt-C Group-1 / Group-3 fold enrichment", c(fc("C", "Group-1"), fc("C", "Group-3")), c(4.3, 2.1), 1)
check("Pun-bt-C Group-2 fold depletion", 1 / fc("C", "Group-2"), 3.3, 1)
check("Pun-bt-E Group-2 fold depletion", 1 / fc("E", "Group-2"), 30, 0)
check("Pun-bt-D Group-7 enrichment / Group-4 depletion", c(fc("D", "Group-7"), 1 / fc("D", "Group-4")), c(4.6, 3.6), 1)
check("Groups significant (FDR < 0.01, FC >= 2) or farm-exclusive", sum(enr$class != "Not significant"), 15)

asg <- read("results", "metahaplotype", "group_assignments.csv")
share <- function(g, farm) {
  x <- asg[asg$gene == paste0("Pun-bt-", g) & asg$farm == farm, ]
  100 * tapply(x$read_count, x$group, sum) / sum(x$read_count)
}
check("Pun-bt-B Farm-S Group-2 / Group-1 %", share("B", "Susceptible")[c("Group-2", "Group-1")], c(76, 24), 0)
check("Pun-bt-C Farm-R Group-3 / Group-1 %", share("C", "Resistant")[c("Group-3", "Group-1")], c(46.5, 24.8), 1)
check("Pun-bt-C Farm-S Group-2 / Group-4 %", share("C", "Susceptible")[c("Group-2", "Group-4")], c(41.4, 30.9), 1)

div <- read("results", "metahaplotype", "group_diversity.csv")
pi_of <- function(g, grp, farm) div$pi[div$gene == paste0("Pun-bt-", g) & div$group == grp & div$farm == farm]
check("pi Pun-bt-B Group-1 Farm-R", pi_of("B", "Group-1", "Resistant"), 0.01, 2)
check("pi Pun-bt-B Group-2 Farm-S", pi_of("B", "Group-2", "Susceptible"), 0.035, 3)
check("pi Pun-bt-B Group-1 Farm-S", pi_of("B", "Group-1", "Susceptible"), 0.005, 3)
check("Pun-bt-B Farm-S pi Group-2 / Group-1 (7-fold)", pi_of("B", "Group-2", "Susceptible") / pi_of("B", "Group-1", "Susceptible"), 7, 0)
check("pi Pun-bt-C Group-3 Farm-R / Farm-S", c(pi_of("C", "Group-3", "Resistant"), pi_of("C", "Group-3", "Susceptible")), c(0.046, 0.009), 3)
check("pi Pun-bt-C Group-1 Farm-R / Farm-S", c(pi_of("C", "Group-1", "Resistant"), pi_of("C", "Group-1", "Susceptible")), c(0.001, 0.033), 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Table S2 - coding variants      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

s2 <- read("tables", "S2_Table.csv")
freq <- function(g, aa) 100 * s2$freq[s2$gene_letter == g & s2$aa_key == aa][order(s2$farm[s2$gene_letter == g & s2$aa_key == aa])]
row1 <- function(g, aa, col) s2[[col]][s2$gene_letter == g & s2$aa_key == aa][1]
check("V315M % Farm-R / Farm-S", freq("B", "V315M"), c(94.96, 0), 2)
check("V315M OR lower 95% bound > 8400", row1("B", "V315M", "or_low") > 8400, TRUE)
check("F167S % Farm-R / Farm-S", freq("B", "F167S"), c(5.04, 0), 2)
check("F167S OR lower 95% bound > 2100", row1("B", "F167S", "or_low") > 2100, TRUE)
check("Q447P % Farm-R / Farm-S", freq("C", "Q447P"), c(28.02, 8.99), 2)
check("Q447P odds ratio", row1("C", "Q447P", "odds_ratio"), 3.94, 2)
check("D328T-Q329R-M330* % Farm-S", c(freq("C", "D328T")[2], freq("C", "Q329R")[2], freq("C", "M330*")[2]), rep(10.99, 3), 2)
check("No canonical F167Y/E198A/F200Y", any(s2$aa_key %in% c("F167Y", "E198A", "F200Y")), FALSE)
check("No differential substitutions at Pun-bt-A or Pun-bt-E", any(s2$gene_letter %in% c("A", "E")), FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Table S3 / Fig.S7 - G_ST             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

s3 <- read("tables", "S3_Table.csv")
gst <- stats::setNames(s3$Gst, s3$gene)
expected_gst <- c(A = 0.098, B = 0.322, C = 0.051, D = 0.037, E = 0.146, F = 0.136, G = 0.137)
check("G_ST per gene A-G", gst[paste0("Pun-bt-", names(expected_gst))], unname(expected_gst), 3)
check("G_ST Pun-bt-B / highest non-candidate", gst[["Pun-bt-B"]] / max(gst[paste0("Pun-bt-", c("A", "D", "E", "F", "G"))]), 2.2, 1)
check("H_S Pun-bt-C", s3$H_S[s3$gene == "Pun-bt-C"], 0.96, 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Fig.S8 - genomic context             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

asm <- read("data", "reference", "assembly_stats.csv")
check("Assembly scaffolds / total Mb / N50 Mb / L50",
      c(asm$n_scaffolds, round(asm$total_mb, 1), round(asm$n50_mb, 2), asm$l50), c(1274, 253.4, 1.83, 44))
coords <- readr::read_tsv(file.path("data", "reference", "btubulin_coords.tsv"), show_col_types = FALSE)
gene_map <- load_gene_map()
scaffold_of <- function(lbl) coords[[3]][match(gene_map$wbps_gene_id[gene_map$label == lbl], coords[[2]])]
check("Pun-bt-B scaffold (manuscript: PgB10)", scaffold_of("Pun-bt-B"), "PgB10")
check("Pun-bt-C scaffold", scaffold_of("Pun-bt-C"), "PgE153")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Fig.6 - expression              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tpm <- read("results", "expression", "beta_tubulin_tpm.csv")
med <- tpm %>% dplyr::group_by(tissue, label) %>% dplyr::summarise(m = stats::median(tpm), .groups = "drop")
rank_in <- function(tis) med$label[med$tissue == tis][order(-med$m[med$tissue == tis])]
check("Anterior: A highest, B second, G lowest", rank_in("Anterior")[c(1, 2, 7)], c("Pun-bt-A", "Pun-bt-B", "Pun-bt-G"))
check("Intestine: A highest, B second, C third, G lowest", rank_in("Intestine")[c(1, 2, 3, 7)],
      c("Pun-bt-A", "Pun-bt-B", "Pun-bt-C", "Pun-bt-G"))

de <- read("results", "expression", "beta_tubulin_deseq2.csv")
sig <- de[!is.na(de$padj) & de$padj < 0.05, ]
check("Significant TBZ responses (gene, contrast)", sort(paste(sig$label, sig$contrast)),
      sort(c("Pun-bt-C concentration_TBZ5_vs_CD", "Pun-bt-E concentration_TBZ9_vs_CD")))
check("Pun-bt-E fold change at 1e-9 M", 2^sig$log2FoldChange[sig$label == "Pun-bt-E"], 2, 0)
check("Pun-bt-C fold reduction at 1e-5 M", 2^-sig$log2FoldChange[sig$label == "Pun-bt-C"], 92, 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      Table S1 - sequencing depth          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

s1 <- read("tables", "S1_Table.csv")
rs <- read("results", "sequencing", "read_summary.csv")
check("Demultiplexed reads (14 barcodes)", rs$demultiplexed_reads, 2281170)
check("Per-barcode range", c(rs$per_barcode_min, rs$per_barcode_max), c(13295, 320220))
check("Per-barcode mean / SD", c(rs$per_barcode_mean, rs$per_barcode_sd), c(162941, 88190), 0)
check("Susceptible / Resistant mean depth (2.5-fold)", rs$susceptible_over_resistant, 2.5, 1)
check("Post-normalisation range", range(s1$reads_normalised_lower_farm), c(13295, 155928))
check("Post-normalisation median / mean / SD", c(rs$normalised_median, rs$normalised_mean, rs$normalised_sd),
      c(105378, 94271, 47582), 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Report                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

report <- do.call(rbind, results)
print(report, right = FALSE, row.names = FALSE)
n_fail <- sum(report$status == "FAIL")
cat(sprintf("\n%d checks, %d PASS, %d FAIL\n", nrow(report), nrow(report) - n_fail, n_fail))
readr::write_csv(report, file.path("results", "qa_manuscript_numbers.csv"))
if (n_fail > 0) quit(status = 1)
