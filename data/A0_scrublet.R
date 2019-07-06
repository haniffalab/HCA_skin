doublets_score <- read.csv("doublets_score.csv")
rownames(doublets_score) <- doublets_score$X
doublets_score <- doublets_score[, 2:4]
seuratobject <- AddMetaData(seuratobject, metadata = doublets_score)
png(file.path(plotDir, 'tSNE_scrublet_prd_hires.png'), width = 1500, height = 1500)
DimPlot(seuratobject, group.by = 'prd', reduction.use = "tsne_harmony")
dev.off()
sum(as.logical(seuratobject@meta.data$prd))
png(file.path(plotDir, 'tSNE_scrublet_ds_hires.png'), width = 1500, height = 1500)
FeaturePlot(seuratobject, features.plot = 'ds', no.legend = F, cols.use = c("grey", "red"))
dev.off()
write(x = rownames(seuratobject@meta.data)[seuratobject@meta.data$prd == "True"], file = "scrublet_doublet_cells_skin.txt")
all(seuratobject@meta.data$run == seuratobject@meta.data$Sanger_sample_ID)
saveRDS(seuratobject, file='seuratobject.RDS')
###############################################################################
# Subset singlets:
not.doublet <- rownames(seuratobject@meta.data)[seuratobject@meta.data$prd == "False"]
old.metadata <- seuratobject@meta.data[not.doublet, ]
subsetseurat <- SubsetData(seuratobject, cells.use = not.doublet, do.clean = T)
saveRDS(subsetseurat, 'seur.RDS')
