# Metadata
metadata <- seuratobject@meta.data[, selectedcolumns]
metadata <- cbind(rownames(metadata), metadata)
colnames(metadata) <- c(
  "Cell",
  "SampleID",
  "Sample",
  "Tissue_layer",
  "Flow_gate",
  "Cell_group",
  "Cell_type")
write.table(metadata, file = "arrayexpress_metadata.txt", row.names = F, quote = F, sep = "\t")

# Counts
exportcounts <- seuratobject@raw.data
exportcounts <- as.matrix(exportcounts)
exportdata <- cbind("Gene" = rownames(exportcounts), exportcounts)
write.table(exportdata, file = "arrayexpress_counts.txt", row.names = F, quote = F, sep = "\t")
