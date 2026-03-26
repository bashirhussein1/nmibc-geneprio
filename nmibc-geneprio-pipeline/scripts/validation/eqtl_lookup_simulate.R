library(data.table)

# Mock eQTL data with 10 variants
eqtl_mock <- data.table(
  variant_id = c("chr1:100:A:T", "chr1:200:C:G", "chr1:300:A:C",
                 "chr2:100:T:A", "chr2:200:G:C", "chr3:100:A:T",
                 "chr3:200:C:G", "chr4:100:A:C", "chr5:100:T:G", "chr6:100:A:G"),
  gene_id = paste0("GENE", 1:10),
  pval_nominal = c(1e-6, 1e-5, 1e-4, 1e-7, 1e-3, 1e-8, 1e-9, 1e-4, 1e-6, 1e-5),
  slope = rnorm(10)
)

# Mock GWAS data — 6 overlapping + 4 unique variants
gwas_mock <- data.table(
  MarkerName = c("chr1:100:A:T", "chr1:200:C:G", "chr2:100:T:A",
                 "chr3:100:A:T", "chr4:100:A:C", "chr5:100:T:G",
                 "chr7:999:A:T", "chr8:999:C:G", "chr9:999:A:C", "chr10:999:T:G"),
  Effect = rnorm(10),
  StdErr = abs(rnorm(10, 0.1))
)

fwrite(eqtl_mock, "mock_eqtl.txt", sep="\t")
fwrite(gwas_mock, "mock_gwas.txt", sep="\t")
cat("Expected overlapping variants: 6\n")
