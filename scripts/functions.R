#!/usr/bin/env Rscript
# ============================================================================
# Shared helpers for the Parascaris univalens β-tubulin analyses
# ============================================================================
#
# Purpose:
#   Constants, input parsers, plotting theme and figure-saving helper used by
#   every script in scripts/. Sourced from the repository root:
#     source(file.path("scripts", "functions.R"))
#
# Inputs:
#   data/metadata/gene_map.tsv          gene label <-> WormBase ParaSite ID
#   data/haplotypes/{A..G}_cluster_info.txt  pbaa cluster headers
#   data/haplotypes/{A..G}_matrix.txt        Clustal Omega percent-identity
#
# Conventions:
#   - No library(); every non-base call is written pkg::fun.
#   - All paths are relative to the repository root.
#   - Farm colours: Resistant = #440154, Susceptible = #FDE725 (viridis ends).
# ============================================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Packages and pipe               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         "\nCreate the environment with: micromamba create -f envs/r-figures.yml",
         call. = FALSE)
  }
  invisible(TRUE)
}

check_packages(c("magrittr", "ggplot2", "dplyr", "tidyr", "readr", "stringr",
                 "patchwork", "scales"))

`%>%` <- magrittr::`%>%`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Constants                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (!file.exists(file.path("scripts", "functions.R"))) {
  stop("Run scripts from the repository root, e.g. Rscript scripts/01_haplotype_diversity.R",
       call. = FALSE)
}

DATA_DIR    <- "data"
HAPLO_DIR   <- file.path(DATA_DIR, "haplotypes")
FIGURES_DIR <- "figures"
TABLES_DIR  <- "tables"
GENES       <- LETTERS[1:7]
FARMS       <- c("Resistant", "Susceptible")
FARM_COLORS <- c(Resistant = "#440154", Susceptible = "#FDE725")

dir.create(FIGURES_DIR, showWarnings = FALSE)
dir.create(TABLES_DIR,  showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Input loaders                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Single source of truth for gene names (Roose et al. 2021 nomenclature).
load_gene_map <- function() {
  readr::read_tsv(file.path(DATA_DIR, "metadata", "gene_map.tsv"),
                  col_types = readr::cols(.default = "c", amplicon_kb = "d"))
}

# One row per pbaa cluster, parsed from a FASTA-style header such as
# >Pun-bt-B-Susceptible_cluster-0_ReadCount-20318 uchime_score:-1 ... cluster_freq:0.130699 ...
read_cluster_info <- function(gene, drop_chimeras = TRUE) {
  path <- file.path(HAPLO_DIR, sprintf("%s_cluster_info.txt", gene))
  headers <- grep("^>", readLines(path), value = TRUE)
  field <- function(key) {
    as.numeric(stringr::str_match(headers, paste0(key, ":([-0-9.eE]+)"))[, 2])
  }
  name <- stringr::str_match(headers, "^>(\\S+)")[, 2]
  info <- data.frame(
    name         = name,
    gene         = gene,
    farm         = ifelse(grepl("Resistant", name), "Resistant", "Susceptible"),
    cluster      = as.integer(stringr::str_match(name, "cluster-(\\d+)")[, 2]),
    read_count   = as.numeric(stringr::str_match(name, "ReadCount-(\\d+)")[, 2]),
    uchime_score = field("uchime_score"),
    cluster_freq = field("cluster_freq"),
    diversity    = field("diversity"),
    avg_quality  = field("avg_quality"),
    seq_length   = field("seq_length"),
    stringsAsFactors = FALSE
  )
  # pbaa flags putative chimeras with a positive uchime score; -1 = not flagged.
  info$is_chimeric <- !is.na(info$uchime_score) & info$uchime_score > 0
  if (drop_chimeras) info <- info[!info$is_chimeric, , drop = FALSE]
  stopifnot(!anyNA(info$read_count), !anyNA(info$cluster_freq))
  info
}

# Clustal Omega --percent-id --full matrix: first line = N, then one row per
# sequence (name followed by N identities). Returns a symmetric numeric matrix.
read_identity_matrix <- function(gene) {
  path  <- file.path(HAPLO_DIR, sprintf("%s_matrix.txt", gene))
  lines <- readLines(path)
  n     <- as.integer(lines[1])
  rows  <- strsplit(trimws(lines[-1][nzchar(lines[-1])]), "\\s+")
  stopifnot(length(rows) == n, all(lengths(rows) == n + 1))
  mat <- matrix(as.numeric(unlist(lapply(rows, `[`, -1))), nrow = n, byrow = TRUE,
                dimnames = list(vapply(rows, `[`, "", 1), vapply(rows, `[`, "", 1)))
  stopifnot(isSymmetric(unname(mat)))
  mat
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Theme and saving                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95"))
}

scale_fill_farm <- function(...) ggplot2::scale_fill_manual(values = FARM_COLORS, ...)
scale_color_farm <- function(...) ggplot2::scale_color_manual(values = FARM_COLORS, ...)

# PLOS ONE wants TIFF (LZW, 300-600 dpi); the PNG is a preview for the README.
save_fig <- function(p, name, width, height, dpi = 600) {
  tif <- file.path(FIGURES_DIR, paste0(name, ".tif"))
  png <- file.path(FIGURES_DIR, paste0(name, ".png"))
  ggplot2::ggsave(tif, p, width = width, height = height, dpi = dpi,
                  device = "tiff", compression = "lzw", bg = "white")
  ggplot2::ggsave(png, p, width = width, height = height, dpi = 150, bg = "white")
  message("Wrote ", tif)
  invisible(tif)
}

write_table <- function(df, name) {
  path <- file.path(TABLES_DIR, paste0(name, ".csv"))
  readr::write_csv(df, path)
  message("Wrote ", path)
  invisible(path)
}
