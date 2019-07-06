library("harmony")
dir.create(geneDir)
dir.create(plotDir)
dir.create(umapDir)
seuratobject <- readRDS('seur.RDS')
seuratobject <- NormalizeData(seuratobject)
seuratobject <- FindVariableGenes(
  object = seuratobject,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.0125,
  x.high.cutoff = 3,
  y.cutoff = 0.5)
write(seuratobject@var.genes, file = "vargenes.txt")
seuratobject <- ScaleData(seuratobject)
seuratobject <- RunPCA(object = seuratobject, pc.genes = seuratobject@var.genes, do.print = FALSE, pcs.print = 1:20, genes.print = 10)
seuratobject <- RunHarmony(object = seuratobject, group.by.vars = "Sample", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 20, max.iter.harmony = 5)
maxPC <- 20
seuratobject <- RunTSNE(seuratobject, reduction.use = "harmony", dims.use = 1:maxPC, do.fast = T, reduction.name="tsne", reduction.key="tsne", seed.use = 42)
seuratobject <- RunUMAP(seuratobject, dims.use=1:maxPC, seed.use = 42, reduction.use = "harmony")
###############################################################################
# Clustering
seuratobject <- FindClusters(seuratobject, reduction.type = "harmony", resolution=2, dims.use=1:maxPC, force.recalc=T, save.SNN=T)
seuratobject <- SetAllIdent(seuratobject, id = "res.2")
# Build node tree:
seuratobject <- BuildClusterTree(seuratobject, pcs.use = 1:maxPC)
seuratobject@cluster.tree
PlotClusterTree(seuratobject)
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, col.name = "All_clustering_res2")
###############################################################################
# Subset data:
# Myeloid
myeloid.subset <- SubsetData(seuratobject, ident.use = c(21, 18, 22, 7, 17, 41, 26), do.clean = T)
saveRDS(myeloid.subset, 'subset_myeloid_seurat.RDS')
table(myeloid.subset@meta.data$Labels_8_simple)
table(myeloid.subset@meta.data$Labels_7_backlabelled)
###############################################################################
# Nonimmune
nonimmune.subset <- SubsetData(seuratobject, ident.use = c(31, 6, 25, 23, 28, 39, 35, 32, 3, 14, 16, 29, 12, 36), do.clean = T)
saveRDS(nonimmune.subset, 'subset_nonimmune_seurat.RDS')
table(nonimmune.subset@meta.data$Labels_8_simple)
table(nonimmune.subset@meta.data$Labels_7_backlabelled)
###############################################################################
# Lymphoid
lymphoid.subset <- SubsetData(seuratobject, ident.use = c(27, 0, 15, 13, 5, 8), do.clean = T)
saveRDS(lymphoid.subset, 'subset_lymphoid_seurat.RDS')
table(lymphoid.subset@meta.data$Labels_8_simple)
table(lymphoid.subset@meta.data$Labels_7_backlabelled)
###############################################################################
# KC
kc.subset <- SubsetData(seuratobject, ident.use = c(1, 11, 9, 38, 37, 2, 4, 33, 10, 19, 30, 24, 20), do.clean = T)
saveRDS(kc.subset, 'subset_kc_seurat.RDS')
table(kc.subset@meta.data$Labels_8_simple)
table(kc.subset@meta.data$Labels_7_backlabelled)
###############################################################################
# Plasma-mast
plasmamast.subset <- SubsetData(seuratobject, ident.use = 34, do.clean = T)
saveRDS(plasmamast.subset, 'subset_plasmamast_seurat.RDS')
table(plasmamast.subset@meta.data$Labels_8_simple)
table(plasmamast.subset@meta.data$Labels_7_backlabelled)
###############################################################################
seuratobject <- SetAllIdent(seuratobject, id = "All_clustering_res2")
# Markers
markers.res2 <- FindAllMarkers(seuratobject,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.25,
  print.bar = TRUE,
  only.pos = FALSE,
  do.print = T,
  random.seed = 42)
markers.res2$FC <- exp(markers.res2$avg_logFC)
saveRDS(markers.res2, "markers.res2.RDS")
markers.res2 <- markerDescription(markers.res2, verbose = T)
write.table(markers.res2, file = "markers.res2.txt", sep = "\t", row.names = F)
saveRDS(markers.res2, 'markers.res2.RDS')
writeLines(capture.output(sessionInfo()), "sessioninfo.txt")
