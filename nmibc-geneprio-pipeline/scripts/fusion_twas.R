# ------------------------------------------------------------------------------
# Script Name: fusion_twas.R
# Description: 
#   This script performs a robust FUSION Transcriptome-Wide Association Study (TWAS)
#   analysis on GWAS summary statistics using precomputed gene expression weights 
#   and a PLINK LD reference panel, all on the hg38 genome build. It outputs a 
#   comprehensive 18-column result table per gene, including TWAS statistics, 
#   heritability, best GWAS SNP, top eQTL, and model performance metrics.
#
# Usage:
#   Rscript fusion_twas.R --sumstats <GWAS_sumstats> --weights <weight_list> \
#     --weights_dir <weights_directory> --ref_ld_chr <plink_ld_prefix> \
#     --chr <chromosome> --out <output_file> [--hsq_p <min_hsq_pvalue>]
#
# Command-line Options:
#   --sumstats      : Path to GWAS summary statistics file (required, tab-delimited)
#   --weights       : Path to weight list file (required)
#   --weights_dir   : Directory containing weight files (required)
#   --ref_ld_chr    : Prefix for PLINK LD reference files (required)
#   --chr           : Chromosome to analyze (required, e.g., "1", "2", ..., "X")
#   --out           : Output file path (required)
#   --hsq_p         : Minimum heritability p-value threshold (default: 0.01)
#
# Main Steps:
#   1. Load and filter GWAS summary statistics for the specified chromosome.
#   2. Load PLINK LD reference panel for the chromosome.
#   3. Load gene expression weight list and filter to the chromosome.
#   4. For each gene:
#      - Load weight file and extract SNP weights.
#      - Subset LD and GWAS data to the gene region (±500kb).
#      - Match SNPs across weights, LD, and GWAS, performing allele QC.
#      - Compute TWAS Z-score and p-value using the best model.
#      - Extract heritability, best GWAS SNP, top eQTL, and model CV metrics.
#      - Record results in the output table.
#   5. Write the results table to the specified output file.
#
# Dependencies:
#   - plink2R
#   - optparse
#
# Notes:
#   - Assumes all input files are formatted correctly and on the hg38 build.
#   - Designed for use on HPC clusters with custom R library paths.
#   - Output table contains 18 columns per gene for downstream analysis.
# ------------------------------------------------------------------------------
#' @author Bashir Hussein 
#' @date 2025-01-07
#'
#!/usr/bin/env Rscript

# ------------------------------
# Robust FUSION TWAS Script on hg38 (full 18‐column output)
# ------------------------------
#!/usr/bin/env Rscript

# ------------------------------
# Robust FUSION TWAS Script on hg38 (with logging of skipped genes)
# ------------------------------
.libPaths("/gpfs/home2/hbashir1/Rlibs")
suppressMessages({
  library(plink2R)
  library(optparse)
})

# Command-line options
opt_list <- list(
  make_option("--sumstats",    type="character", help="GWAS sumstats [required]"),
  make_option("--weights",     type="character", help="Weight list file [required]"),
  make_option("--weights_dir", type="character", help="Directory for weight files"),
  make_option("--ref_ld_chr",  type="character", help="PLINK LD reference prefix [required]"),
  make_option("--chr",         type="character", help="Chromosome to analyze [required]"),
  make_option("--out",         type="character", help="Output file [required]")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# Allele‐QC helper
allele.qc <- function(a1,a2,ref1,ref2) {
  a1<-toupper(a1); a2<-toupper(a2)
  ref1<-toupper(ref1); ref2<-toupper(ref2)
  flip_ref <- function(x) { f<-x; f[x=="A"]<-"T"; f[x=="T"]<-"A"; f[x=="C"]<-"G"; f[x=="G"]<-"C"; f }
  flip1 <- flip_ref(ref1); flip2 <- flip_ref(ref2)
  keep <- a1%in%c("A","C","G","T") & a2%in%c("A","C","G","T") &
          !( (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|
             (a1=="C"&a2=="G")|(a1=="G"&a2=="C") )
  flip <- (a1==ref2 & a2==ref1) | (a1==flip2 & a2==flip1)
  list(keep=keep, flip=flip)
}

# 1) LOAD GWAS sumstats & filter to chromosome
raw <- read.table(opt$sumstats, header=TRUE, sep="\t", as.is=TRUE)
raw$MarkerName <- trimws(raw$MarkerName)
raw$chr <- sub(":.*","", raw$MarkerName)
raw_chr <- subset(raw, chr==paste0("chr",opt$chr))
sumstat <- data.frame(
  locus = trimws(raw_chr$MarkerName),
  A1    = toupper(raw_chr$Allele1),
  A2    = toupper(raw_chr$Allele2),
  Z     = raw_chr$Effect / raw_chr$StdErr,
  stringsAsFactors=FALSE
)

# 2) LOAD LD reference (hg38)
ld_pref <- paste0(opt$ref_ld_chr, opt$chr)
ld <- read_plink(ld_pref, impute="avg")
colnames(ld$bim) <- c("chr","locus_id","cM","pos","ref","alt")
ld$bim$locus <- with(ld$bim, paste0(chr, ":", pos, ":", ref, ":", alt))
ld_mat <- as.matrix(ld$bed)

# 3) LOAD weight list
wgtlist <- read.table(opt$weights, header=TRUE, stringsAsFactors=FALSE)
wgtlist <- subset(wgtlist, CHR==opt$chr)

# Prepare output with full 18 columns
out <- data.frame(
  FILE         = wgtlist$WGT,
  ID           = wgtlist$ID,
  CHR          = wgtlist$CHR,
  P0           = wgtlist$P0,
  P1           = wgtlist$P1,
  HSQ          = NA_real_,
  BEST.GWAS.ID = NA_character_,
  BEST.GWAS.Z  = NA_real_,
  EQTL.ID      = NA_character_,
  EQTL.R2      = NA_real_,
  EQTL.Z       = NA_real_,
  EQTL.GWAS.Z  = NA_real_,
  NSNP         = NA_integer_,
  MODEL        = NA_character_,
  MODELCV.R2   = NA_real_,
  MODELCV.PV   = NA_real_,
  TWAS.Z       = NA_real_,
  TWAS.P       = NA_real_,
  stringsAsFactors=FALSE
)

# Prepare skip log
skip_log <- data.frame(ID=character(), Reason=character(), stringsAsFactors=FALSE)

# PER-GENE LOOP
for(i in seq_len(nrow(wgtlist))) {
  gid <- wgtlist$ID[i]
  wfile <- file.path(opt$weights_dir, wgtlist$WGT[i])
  if (!file.exists(wfile)) {
    skip_log <- rbind(skip_log, data.frame(ID=gid, Reason="Missing weight file")); next
  }

  env <- new.env(); load(wfile, envir=env)
  if (!"snps" %in% ls(env)) {
    skip_log <- rbind(skip_log, data.frame(ID=gid, Reason="No SNPs in weight file")); next
  }
  snps <- env$snps; snps$id <- snps$id_hg38
  w_locus <- snps$id; wmat <- env$wgt.matrix; cvp <- env$cv.performance

  P0 <- max(1, wgtlist$P0[i] - 500000)
  P1 <- wgtlist$P1[i] + 500000
  sel <- which(ld$bim$pos >= P0 & ld$bim$pos <= P1)
  if (length(sel)==0) {
    skip_log <- rbind(skip_log, data.frame(ID=gid, Reason="No SNPs in LD window")); next
  }
  bim_sub <- ld$bim[sel,]; bed_sub <- ld_mat[,sel,drop=FALSE]

  sumstat_sub <- sumstat
  w_coord <- sub("(:[^:]+){2}$","", w_locus)
  ld_coord <- paste0(bim_sub$chr, ":", bim_sub$pos)
  m_pos <- match(w_coord, ld_coord)
  if (sum(!is.na(m_pos)) < 1) {
    skip_log <- rbind(skip_log, data.frame(ID=gid, Reason="No LD-GWAS overlap"))
    out[i, "NSNP"] <- 0
    next
  }
  valid <- which(!is.na(m_pos))
  w2 <- wmat[valid,,drop=FALSE]; curb <- bim_sub[m_pos[valid],]; geno2 <- bed_sub[,m_pos[valid],drop=FALSE]
  w_coord <- ld_coord[m_pos[valid]]

  g_coord <- sub("(:[^:]+){2}$","", sumstat_sub$locus)
  m2 <- match(w_coord, g_coord)
  if (sum(!is.na(m2)) < 1) {
    skip_log <- rbind(skip_log, data.frame(ID=gid, Reason="No GWAS SNPs matched"))
    out[i, "NSNP"] <- length(valid)
    next
  }
  z_all <- sumstat_sub$Z[m2]; a1_all <- sumstat_sub$A1[m2]; a2_all <- sumstat_sub$A2[m2]
  ref_all <- curb$ref; alt_all <- curb$alt
  qc <- allele.qc(a1_all,a2_all,ref_all,alt_all)
  flip <- qc$flip; flip[is.na(flip)] <- FALSE
  keep <- qc$keep; keep[is.na(keep)] <- FALSE
  z_all[flip] <- -z_all[flip]; a1_all[flip] <- ref_all[flip]; a2_all[flip] <- alt_all[flip]
  w2 <- w2[keep,,drop=FALSE]; geno2 <- geno2[,keep,drop=FALSE]; z_all <- z_all[keep]
  if (length(z_all) < 1) {
    skip_log <- rbind(skip_log, data.frame(ID=gid, Reason="No SNPs after allele QC"))
    out[i, "NSNP"] <- sum(keep)
    next
  }

  geno2_s <- scale(geno2)
  LDmat <- crossprod(geno2_s) / (nrow(geno2_s)-1)
  pidx <- grep("pval",rownames(cvp)); if (!length(pidx)) pidx <- grep("^p",rownames(cvp))
  best <- which.min(apply(cvp[pidx,,drop=FALSE],2,min,na.rm=TRUE))
  wv <- w2[,best]; denom <- sqrt(as.numeric(t(wv)%*%LDmat%*%wv))
  tZ <- as.numeric((wv %*% z_all)/denom); tP <- 2*pnorm(-abs(tZ))

  hsq <- env$hsq[1]
  loci_win <- sumstat_sub$locus[m2]
  best_g_idx <- which.max(abs(z_all)); best_gwas_id <- loci_win[best_g_idx]; best_gwas_z <- z_all[best_g_idx]
  top1_w <- w2[,"top1"]; eqtl_idx <- which.max(abs(top1_w)); eqtl_id <- w_coord[eqtl_idx]
  eqtl_r2 <- cvp["rsq","top1"]; eqtl_z <- top1_w[eqtl_idx]; eqtl_gwas_z <- z_all[eqtl_idx]
  model <- colnames(w2)[best]; model_cv_r2 <- cvp["rsq", best]; model_cv_pv <- cvp["pval",best]

  out[i,] <- list(
    wgtlist$WGT[i], wgtlist$ID[i], wgtlist$CHR[i], wgtlist$P0[i], wgtlist$P1[i],
    hsq, best_gwas_id, best_gwas_z, eqtl_id, eqtl_r2, eqtl_z, eqtl_gwas_z,
    nrow(LDmat), model, model_cv_r2, model_cv_pv, tZ, tP
  )
}

# Bonferroni correction and output

# Get tested genes
tested <- !is.na(out$TWAS.P)
N_eff <- sum(tested)
Pcrit <- 0.05 / N_eff
significant <- subset(out, tested & TWAS.P < Pcrit)

# Save significant results
sign_file <- sub("\\.txt$", ".significant.txt", opt$out)
write.table(significant, sign_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Save all tested results only
tested_file <- opt$out
write.table(out[tested,], tested_file, sep="\t", quote=FALSE, row.names=FALSE)

# Save skip log
skip_file <- sub("\\.txt$", ".skipped_genes.tsv", opt$out)
write.table(skip_log, skip_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Print summary
cat("\nNumber of genes attempted:", nrow(wgtlist), "\n")
cat("Number of genes tested:", N_eff, "\n")
cat("Bonferroni threshold:", format(Pcrit, scientific=TRUE), "\n")
cat("Number of significant genes:", nrow(significant), "\n")

cat("\n--- Gene Skipping Summary ---\n")
print(as.data.frame(table(skip_log$Reason)))
