#!/usr/bin/env Rscript
# ============================================================================
# Transcript-to-gene table for tximport (one-off helper for script 07)
# ============================================================================
#
# Purpose:
#   Builds the transcript -> gene mapping used to summarise salmon transcript
#   quantifications to gene level, from the WormBase ParaSite WBPS17
#   P. univalens annotation. Run once; its output is committed so that
#   scripts/07_expression_tbz.R does not need the 130 MB GFF3.
#   The mapping is taken from a GenomicFeatures TxDb (TXNAME -> GENEID), as in
#   the original analysis.
#
# Inputs:
#   <gff3>   parascaris_univalens.PRJNA386823.WBPS17.annotations.gff3, from
#            https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/parascaris_univalens/PRJNA386823/parascaris_univalens.PRJNA386823.WBPS17.annotations.gff3.gz
#            (decompress before use)
#
# Outputs:
#   data/rnaseq/tx2gene.tsv   TXNAME, GENEID
#
# Usage (from the repository root):
#   Rscript scripts/07a_make_tx2gene.R path/to/parascaris_univalens.PRJNA386823.WBPS17.annotations.gff3
# ============================================================================

source(file.path("scripts", "functions.R"))
check_packages(c("GenomicFeatures", "AnnotationDbi"))

gff3 <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(gff3) || !file.exists(gff3)) {
  stop("Give the path to the WBPS17 GFF3 as the first argument.", call. = FALSE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           TxDb and tx2gene                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

txdb <- GenomicFeatures::makeTxDbFromGFF(gff3)
tx2gene <- AnnotationDbi::select(txdb,
                                 keys    = AnnotationDbi::keys(txdb, keytype = "TXNAME"),
                                 columns = "GENEID",
                                 keytype = "TXNAME")
stopifnot(!anyDuplicated(tx2gene$TXNAME))

out <- file.path(DATA_DIR, "rnaseq", "tx2gene.tsv")
readr::write_tsv(tx2gene, out)
message("Wrote ", out, " (", nrow(tx2gene), " transcripts, ",
        length(unique(tx2gene$GENEID)), " genes)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Session info                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(sessionInfo())
