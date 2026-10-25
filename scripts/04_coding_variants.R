#!/usr/bin/env Rscript
# ============================================================================
# Amino-acid variants in haplotype CDS (Table S2)
# ============================================================================
#
# Purpose:
#   Call amino-acid differences between each pbaa haplotype CDS and the
#   reference protein, weight them by the haplotype read count to obtain the
#   frequency of variant-carrying reads per farm, and test farm association
#   with Fisher's exact test (OR, 95% CI, p, BH FDR within gene). Variants are
#   flagged at canonical benzimidazole-resistance codons (167, 198, 200), in
#   the C-terminal region (codon >= 440), as premature stops, and as also
#   called by VEP from the reads of the same farm.
#
# Inputs:
#   data/cds_alignments/Pun-bt-{A..G}-align.fas
#                                   pbaa haplotype + reference protein
#                                   alignments; haplotype names carry farm,
#                                   cluster and read count
#   results/read_variants/exonic_variants.csv
#                                   VEP exonic calls from scripts/03_intronic_variants.R
#
# Outputs:
#   tables/S2_Table.csv         Table S2: amino-acid variants per gene and farm
#
# Usage (from the repository root):
#   Rscript scripts/03_intronic_variants.R   # writes the VEP exonic calls
#   Rscript scripts/04_coding_variants.R
# ============================================================================

source(file.path("scripts", "functions.R"))
check_packages(c("Biostrings"))

CANONICAL_POSITIONS <- c(167L, 198L, 200L)
CTAIL_START         <- 440L
CDS_ALIGN_DIR       <- file.path(DATA_DIR, "cds_alignments")
EXONIC_VARIANTS     <- file.path("results", "read_variants", "exonic_variants.csv")

if (!file.exists(EXONIC_VARIANTS)) {
  stop("Missing ", EXONIC_VARIANTS, "; run scripts/03_intronic_variants.R first.",
       call. = FALSE)
}

gene_map <- load_gene_map()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Helpers                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Residues after the first stop are not translated, so they are ignored.
trim_at_first_stop <- function(s) {
  pos <- regexpr("\\*", s)
  ifelse(pos > 0, substr(s, 1, pos), s)
}

# Alignment columns where the haplotype differs from the reference; gaps in
# either sequence are not counted as substitutions.
diffs_vs_ref <- function(seq_trim, ref_vec) {
  qv <- strsplit(seq_trim, "")[[1]]
  n  <- min(length(qv), length(ref_vec))
  pos <- which(qv[seq_len(n)] != ref_vec[seq_len(n)] &
                 qv[seq_len(n)] != "-" & ref_vec[seq_len(n)] != "-")
  dplyr::tibble(position = pos, ref = ref_vec[pos], alt = qv[pos],
                aa_key = paste0(ref_vec[pos], pos, qv[pos]))
}

# 2x2 table of variant vs other reads, Resistant vs Susceptible;
# OR > 1 = enriched in the resistant farm.
fisher_by_variant <- function(df) {
  df$farm <- factor(df$farm, levels = FARMS)
  df <- df[order(df$farm), ]
  r <- df[1, ]; s <- df[2, ]
  mat <- matrix(c(r$var_reads, r$total_reads - r$var_reads,
                  s$var_reads, s$total_reads - s$var_reads), nrow = 2)
  ft <- tryCatch(stats::fisher.test(mat, conf.int = TRUE), error = function(e) NULL)
  if (is.null(ft)) {
    return(dplyr::tibble(p_fisher = NA_real_, odds_ratio = NA_real_,
                         or_low = NA_real_, or_high = NA_real_))
  }
  dplyr::tibble(p_fisher = ft$p.value, odds_ratio = as.numeric(ft$estimate),
                or_low = ft$conf.int[1], or_high = ft$conf.int[2])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           VEP amino-acid keys             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

exonic_variants <- readr::read_csv(EXONIC_VARIANTS, col_types = readr::cols(.default = "c"))

# VEP Amino_acids "V/M" at Protein_position "315" -> "V315M". Multi-residue
# changes keep only the first residue on each side.
vep_keys_for_farm <- function(farm, status, suffix) {
  exonic_variants %>%
    dplyr::filter(VariantStatus %in% status) %>%
    dplyr::mutate(
      pos = suppressWarnings(as.integer(stringr::str_extract(.data[[paste0("Protein_position", suffix)]], "^\\d+"))),
      aa_pair = .data[[paste0("Amino_acids", suffix)]],
      ref = ifelse(stringr::str_detect(aa_pair, "/"), stringr::str_extract(aa_pair, "^[A-Z\\*]"), NA_character_),
      alt = ifelse(stringr::str_detect(aa_pair, "/"), stringr::str_extract(aa_pair, "(?<=/)[A-Z\\*]"), NA_character_),
      aa_key = ifelse(!is.na(ref) & !is.na(pos) & !is.na(alt), paste0(ref, pos, alt), NA_character_)
    ) %>%
    dplyr::filter(!is.na(aa_key)) %>%
    dplyr::distinct(gene_letter = Gene_name, aa_key) %>%
    dplyr::mutate(farm = farm, in_vep = TRUE)
}

vep_keys <- dplyr::bind_rows(
  vep_keys_for_farm("Resistant",   c("Resistant-only", "Shared"),   ".res"),
  vep_keys_for_farm("Susceptible", c("Susceptible-only", "Shared"), ".sus")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Table S2 - coding variants      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyze_gene <- function(gene) {
  label <- gene_map$label[gene_map$gene == gene]
  aa <- Biostrings::readAAStringSet(file.path(CDS_ALIGN_DIR, sprintf("%s-align.fas", label)))

  # The reference is the only record named "<label>:<transcript>".
  ref_idx <- grep(":", names(aa))
  stopifnot(length(ref_idx) == 1)
  ref_vec <- strsplit(trim_at_first_stop(as.character(aa[[ref_idx]])), "")[[1]]

  # Farm and read count come from the haplotype names. No chimera filter is
  # applied: the submitted Table S2 kept the pbaa-flagged chimeric clusters
  # (1 at Pun-bt-A, 4 at Pun-bt-C) although Methods state they were excluded.
  clusters <- dplyr::tibble(name = names(aa)[-ref_idx]) %>%
    dplyr::mutate(
      gene_full   = label,
      gene_letter = gene,
      farm        = ifelse(grepl("Resistant", name), "Resistant", "Susceptible"),
      cluster     = as.integer(stringr::str_match(name, "cluster-(\\d+)")[, 2]),
      read_count  = as.integer(stringr::str_match(name, "ReadCount-(\\d+)")[, 2]),
      aa_trim     = trim_at_first_stop(as.character(aa[-ref_idx]))
    )
  stopifnot(!anyNA(clusters$cluster), !anyNA(clusters$read_count))

  totals <- clusters %>%
    dplyr::group_by(gene_full, gene_letter, farm) %>%
    dplyr::summarise(total_reads = sum(read_count), .groups = "drop")

  cluster_diffs <- clusters %>%
    dplyr::mutate(diffs = lapply(aa_trim, diffs_vs_ref, ref_vec = ref_vec)) %>%
    dplyr::select(gene_full, gene_letter, farm, cluster, read_count, diffs) %>%
    tidyr::unnest(diffs)

  if (nrow(cluster_diffs) == 0) {
    message("Pun-bt-", gene, ": no amino-acid differences from the reference")
    return(NULL)
  }

  variant_reads <- cluster_diffs %>%
    dplyr::group_by(gene_full, gene_letter, aa_key, farm) %>%
    dplyr::summarise(var_reads = sum(read_count), .groups = "drop") %>%
    tidyr::complete(gene_full, gene_letter, aa_key, farm = FARMS,
                    fill = list(var_reads = 0L)) %>%
    dplyr::left_join(totals, by = c("gene_full", "gene_letter", "farm")) %>%
    dplyr::mutate(freq = var_reads / total_reads)

  # BH correction across the variants of this gene only.
  variant_stats <- variant_reads %>%
    dplyr::group_by(gene_full, gene_letter, aa_key) %>%
    dplyr::group_modify(~ fisher_by_variant(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_adj = stats::p.adjust(p_fisher, method = "BH"))

  variant_reads %>%
    dplyr::left_join(variant_stats, by = c("gene_full", "gene_letter", "aa_key")) %>%
    dplyr::mutate(
      freq_pct     = sprintf("%.2f%%", 100 * freq),
      pos          = as.integer(stringr::str_extract(aa_key, "\\d+")),
      is_canonical = pos %in% CANONICAL_POSITIONS,
      is_stop      = stringr::str_detect(aa_key, "\\*"),
      is_ctail     = !is.na(pos) & pos >= CTAIL_START
    ) %>%
    dplyr::left_join(vep_keys, by = c("gene_letter", "farm", "aa_key")) %>%
    dplyr::mutate(in_vep = tidyr::replace_na(in_vep, FALSE))
}

s2_table <- dplyr::bind_rows(lapply(GENES, analyze_gene))

write_table(s2_table, "S2_Table")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
