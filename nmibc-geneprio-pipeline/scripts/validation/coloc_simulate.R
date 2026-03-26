library(coloc)

# Load built-in test data
data(coloc_test_data)
attach(coloc_test_data)

# Mock eQTL (D1) — contains shared causal variant
eqtl_mock <- data.table(
  variant_id   = paste("1", D1$position, "A", "T", sep=":"),
  slope        = D1$beta,
  slope_se     = sqrt(D1$varbeta),
  pval_nominal = 2 * pnorm(-abs(D1$beta / sqrt(D1$varbeta))),
  af           = D1$MAF,
  gene_name    = "TESTGENE"
)

# Mock GWAS (D2) — same SNPs, shared causal variant
gwas_mock <- data.table(
  MarkerName      = paste("1", D2$position, "A", "T", sep=":"),
  Effect          = D2$beta,
  StdErr          = sqrt(D2$varbeta),
  `P-value`       = 2 * pnorm(-abs(D2$beta / sqrt(D2$varbeta))),
  MinFreq         = D2$MAF,
  TotalSampleSize = 1000,
  TotalEvents     = 500
)

outdir <- "C:/path/to/output"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fwrite(eqtl_mock, file.path(outdir, "mock_eqtl.tsv"), sep="\t")
fwrite(gwas_mock, file.path(outdir, "mock_gwas.tsv"), sep="\t")
