library("monocle")
library("Matrix")
# Load:
path_data <- "."
metadata <- read.csv(file.path(path_data, "pdt_metadata.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
genes <- read.csv(file.path(path_data, "pdt_genes.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
hvg <- genes[genes$highly_variable == "True", ]$gene_short_name
c_sparse <- readMM(file.path(path_data, "pdt_counts.mtx"))
c_mat <- as.matrix(c_sparse)

# Convert to cds:
pd <- new("AnnotatedDataFrame", data = metadata)
fd <- new("AnnotatedDataFrame", data = genes)
cds <- newCellDataSet(
  c_mat,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size())
# Calculate:
cds <- estimateSizeFactors(cds)
pData(cds)$Pseudotime <- metadata$dpt_pseudotime

# Exclude the following section if you want to test all genes:
to_be_tested <- row.names(subset(fData(cds), gene_short_name %in% hvg))
cds_subset <- cds[to_be_tested, ]

diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")

saveRDS(diff_test_res, file = "diff_test_res.RDS")

sink("sessionInfo.txt")
sessionInfo()
sink()
