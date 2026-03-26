#' coloc.R - Colocalization Analysis Script
#'
#' This script performs colocalization analysis between eQTL and GWAS summary statistics using the `coloc` R package.
#' It identifies loci with evidence for shared genetic signals between gene expression and disease/trait association.
#'
#' ## Usage
#' ```
#' eqtl[gz], gwas, outdir, gwas_p, window_kb
#' ```
#' - `eqtl[gz]`      : Path to eQTL summary statistics file (TSV or GZ).
#' - `gwas`          : Path to GWAS summary statistics file.
#' - `outdir`        : Output directory for results.
#' - `gwas_p`        : GWAS p-value threshold for lead SNP selection.
#' - `window_kb`     : Window size (in kilobases) around lead SNPs for locus definition.
#'
#' ## Workflow
#' 1. **Input Parsing**: Reads command-line arguments and checks input files.
#' 2. **Data Loading**: Loads eQTL and GWAS summary statistics using `data.table::fread`.
#' 3. **Preprocessing**:
#'    - Parses variant IDs to extract chromosome, position, reference, and alternate alleles.
#'    - Computes effect sizes, variances, and minor allele frequencies.
#'    - Estimates standard deviation of the eQTL trait.
#' 4. **Lead SNP Identification**: Selects lead SNPs based on p-value thresholds in GWAS dataset.
#' 5. **Locus Definition**: For each lead SNP, defines a locus window and extracts overlapping variants from both datasets.
#' 6. **Colocalization Analysis**: Runs `coloc.abf` for each locus to compute posterior probabilities for five hypotheses (H0-H4).
#' 7. **Results Output**:
#'    - Writes locus-level hypothesis summary (`coloc_hypotheses_summary.txt`).
#'    - Writes SNP-level colocalization results (`coloc_snp_level.txt`).
#' 8. **Candidate Gene Extraction**:
#'    - Identifies loci with strong colocalization evidence (PP.H4 ≥ 0.8).
#'    - Writes candidate gene list (`coloc_candidate_genes.txt`).
#'
#' ## Output Files
#' - `coloc_hypotheses_summary.txt` : Locus-level summary of posterior probabilities for each hypothesis.
#' - `coloc_snp_level.txt`          : SNP-level colocalization results.
#' - `coloc_candidate_genes.txt`    : Candidate genes and credible sets for loci with strong colocalization.
#' - `coloc.log`                    : Log file with progress and messages.
#'
#' ## Dependencies
#' - R packages: `data.table`, `coloc`
#'
#' ## Notes
#' - Assumes eQTL and GWAS files contain variant IDs in the format `chr:pos:ref:alt`.
#' - The eQTL sample size (`N_eqtl`) is hardcoded to 134 (bladder tissue).
#' - Designed for use with GTEx eQTL and GWAS summary statistics.
#'
#'

#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(coloc)
})

# === Parse command-line arguments ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: coloc.R eqtl[gz] gwas outdir eqtl_p gwas_p window_kb")
}
eqtl_path        <- args[1]
gwas_path        <- args[2]
output_dir       <- args[3]
eqtl_pval_thresh <- as.numeric(args[4])
gwas_pval_thresh <- as.numeric(args[5])
window_kb        <- as.numeric(args[6])

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# === Logging ===
logfile <- file.path(output_dir, "coloc.log")
con     <- file(logfile, open = "wt")
sink(con, type = "output"); sink(con, type = "message")
message("=== Starting coloc.R ==="); flush.console()

# === Read & prep eQTL data ===
eqtl <- fread(eqtl_path)
setnames(eqtl, "variant_id", "snp")
eqtl[, c("chr", "pos", "ref", "alt") := tstrsplit(snp, ":", fixed = TRUE)]
eqtl[, pos := as.integer(pos)]

N_eqtl <- 80
eqtl[, `:=`(
  beta_eqtl    = as.numeric(slope),
  varbeta_eqtl = (slope_se)^2,
  MAF_eqtl     = as.numeric(af)
)]

eqtl_sdY <- coloc:::sdY.est(
  vbeta = eqtl$varbeta_eqtl,
  maf   = eqtl$MAF_eqtl,
  n     = N_eqtl
)

# === Read & prep GWAS data ===
gwas <- fread(gwas_path)
setnames(gwas, "MarkerName", "snp")
gwas[, c("chr", "pos", "ref", "alt") := tstrsplit(snp, ":", fixed = TRUE)]
gwas[, pos := as.integer(pos)]

gwas[, `:=`(
  beta_gwas    = as.numeric(Effect),
  varbeta_gwas = (StdErr)^2,
  MAF_gwas     = as.numeric(MinFreq),
  N_gwas       = as.numeric(TotalSampleSize),
  cases_gwas   = as.numeric(TotalEvents)
)]
gwas[, controls_gwas := N_gwas - cases_gwas]

# === CHANGE 1: Only GWAS lead SNPs define loci ===
lead_snps <- gwas[`P-value` < gwas_pval_thresh, unique(snp)]
message(sprintf("Number of GWAS lead SNPs: %d", length(lead_snps)))

# Initialize containers
locus_list <- vector("list", length(lead_snps))
snp_list   <- vector("list", length(lead_snps))

# === Main loop: per GWAS lead SNP, run coloc ===
for (i in seq_along(lead_snps)) {
  lead <- lead_snps[i]
  cg <- gwas[snp == lead]
  if (nrow(cg) == 0) next

  chr    <- cg$chr[1]
  center <- cg$pos[1]
  start  <- center - window_kb * 1000
  end    <- center + window_kb * 1000

  re <- unique(eqtl[chr == chr & pos >= start & pos <= end], by = "snp")
  rg <- unique(gwas[chr == chr & pos >= start & pos <= end], by = "snp")
  m  <- merge(re, rg, by = "snp")
  m  <- m[complete.cases(m[, .(beta_eqtl, varbeta_eqtl, beta_gwas, varbeta_gwas)]), ]
  if (nrow(m) < 3) next

  gene <- unique(m$gene_name)[1]

  res <- coloc.abf(
    dataset1 = list(
      snp     = m$snp,
      beta    = m$beta_eqtl,
      varbeta = m$varbeta_eqtl,
      type    = "quant",
      N       = N_eqtl,
      sdY     = eqtl_sdY
    ),
    dataset2 = list(
      snp       = m$snp,
      beta      = m$beta_gwas,
      varbeta   = m$varbeta_gwas,
      type      = "cc",
      Ncases    = m$cases_gwas,
      Ncontrols = m$controls_gwas
    )
  )

  pp <- as.numeric(res$summary[paste0("PP.H", 0:4, ".abf")])
  locus_list[[i]] <- data.table(
    lead_snp  = lead,
    chr       = chr,
    center    = center,
    gene_name = gene,
    PP.H0     = pp[1], PP.H1 = pp[2], PP.H2 = pp[3],
    PP.H3     = pp[4], PP.H4 = pp[5]
  )

  tbl <- as.data.table(res$results)
  tbl[, `:=`(
    gene_name = gene,
    lead_snp  = lead,
    chr       = chr,
    center    = center
  )]
  snp_list[[i]] <- tbl

  message(sprintf("[coloc] %3d/%3d loci done", i, length(lead_snps)))
  flush.console()
}

# === Write summary outputs ===
all_loci <- rbindlist(locus_list)
fwrite(
  all_loci[PP.H4 >= 0.8],
  file.path(output_dir, "coloc_hypotheses_summary.txt"),
  sep = "\t"
)
fwrite(rbindlist(snp_list), file.path(output_dir, "coloc_snp_level.txt"), sep = "\t")
message("Wrote hypothesis and SNP-level summaries.")

# simply extract candidate genes from loci with PP.H4 >= 0.8
candidates <- all_loci[PP.H4 >= 0.8, .(
  lead_snp  = lead_snp,
  chr       = chr,
  center    = center,
  gene_name = gene_name,
  PP.H4     = PP.H4
)]

fwrite(candidates, file.path(output_dir, "coloc_candidate_genes.txt"), sep = "\t")
message("Wrote candidate gene list.")

message("=== Finished coloc.R ===")
sink(type="message"); sink(type="output"); close(con)
# EOF
