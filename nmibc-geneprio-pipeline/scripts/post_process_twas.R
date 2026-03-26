#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script Name: post_process_twas.R
# Description:
#   This new script combines per-chromosome TWAS results and applies Transcriptome-wide significance threshold using Bonferroni correction 
#   correction based on the total number of genes actually tested across all chromosomes
#
# Usage:
#   Rscript post_process_twas.R <results_dir> <out_prefix>
#
# Arguments:
#   results_dir  : Directory containing per-chromosome TWAS result files
#   out_prefix   : Prefix for output files
# ------------------------------------------------------------------------------
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript post_process_twas.R <results_dir> <out_prefix>")
}
results_dir <- args[1]
out_prefix  <- args[2]

# Read all per-chromosome result files
chr_files <- list.files(results_dir, pattern="chr[0-9]+\\.txt$", full.names=TRUE)
cat(sprintf("Found %d chromosome result files\n", length(chr_files)))

# Combine all results
all_results <- rbindlist(lapply(chr_files, fread, sep="\t"), fill=TRUE)
cat(sprintf("Total genes attempted: %d\n", nrow(all_results)))

# Count only actually tested genes (TWAS.P is not NA)
tested  <- !is.na(all_results$TWAS.P)
N_eff   <- sum(tested)
cat(sprintf("Total genes tested (genome-wide): %d\n", N_eff))

# Apply genome-wide Bonferroni
Pcrit <- 0.05 / N_eff
cat(sprintf("Genome-wide Bonferroni threshold: P < %.2e (0.05 / %d)\n", Pcrit, N_eff))

# Identify significant genes
significant <- all_results[tested & TWAS.P < Pcrit]
cat(sprintf("Significant genes: %d\n", nrow(significant)))

# Write outputs
all_out  <- paste0(out_prefix, "_all_chromosomes.txt")
sig_out  <- paste0(out_prefix, "_significant.txt")

fwrite(all_results, all_out,  sep="\t", quote=FALSE)
fwrite(significant, sig_out,  sep="\t", quote=FALSE)

cat(sprintf("Written: %s\n", all_out))
cat(sprintf("Written: %s\n", sig_out))
