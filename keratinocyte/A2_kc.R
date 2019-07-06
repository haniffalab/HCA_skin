seuratobject <- readRDS(file = "keratinocyte.RDS")
# Add scanpy-calculated FDG to the Seurat object:
fdgvalues <- read.table("../../scanpy_scrubletremoved/keratinocyte/X_FA_FDG.txt", sep=' ')
scanpy.cell.names <- scan("../../scanpy_scrubletremoved/keratinocyte/cell_names.txt", what=character())
rownames(fdgvalues) <- scanpy.cell.names
seuratobject <- SetDimReduction(seuratobject, reduction.type="fdg_scanpy", slot="cell.embeddings", new.data=as.matrix(fdgvalues))
seuratobject <- SetDimReduction(seuratobject, reduction.type="fdg_scanpy", slot="key", new.data="fdg_scanpy")
kc.colours <- c("#E87D72", "#87AC34", "#0e6c8b", "#b8bae5")
DimPlot(object=seuratobject, reduction.use='fdg_scanpy', do.label=T, group.by = "harmony_annotation_11", cols.use = kc.colours)
# Prepend cluster numbers with text:
seuratobject <- SetAllIdent(seuratobject, id = 'Clustering_res.0.35_harmony_merged')
seuratobject@meta.data$Clustering_res.0.35_harmony_merged_txt <- paste0("cluster", seuratobject@meta.data$Clustering_res.0.35_harmony_merged)
saveRDS(seuratobject, file = "keratinocyte.RDS")
