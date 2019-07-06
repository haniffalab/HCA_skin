seuratobject <- readRDS("../data/subset_kc_seurat.RDS")
dir.create(plotDir)
dir.create(geneDir)
dir.create(umapDir)
seuratobject <- NormalizeData(seuratobject)
seuratobject <- ScaleData(seuratobject)
seuratobject <- FindVariableGenes(seuratobject, x.low.cutoff=0.1, x.high.cutoff=Inf, y.cutoff=0.5, y.high.cutoff=Inf)
length(seuratobject@var.genes)
# 1020 genes
write(seuratobject@var.genes, file = "hca.keratinocytes.vargenes.txt")
seuratobject <- RunPCA(seuratobject, pc.genes=seuratobject@var.genes, do.print=TRUE, pcs.print=1:5, genes.print=5)
library("harmony")
seuratobject <- RunHarmony(object = seuratobject, group.by.vars = "Sample", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 20, max.iter.harmony = 5)
maxPC <- 20
seuratobject <- RunTSNE(seuratobject, reduction.use="harmony", dims.use=1:maxPC, do.fast=T, reduction.name="tsne", reduction.key="tsne")
###############################################################################
seuratobject <- FindClusters(seuratobject, resolution=1.2, dims.use=1:maxPC, force.recalc=T, save.SNN=T, reduction.type="harmony")
resol <- 0.35
seuratobject <- FindClusters(seuratobject, resolution=resol, dims.use=1:maxPC, reuse.SNN=T, reduction.type="harmony")
selected.res <- 0.35
seuratobject <- SetAllIdent(seuratobject, id=paste0('res.', selected.res))
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, paste0('Clustering_res.', selected.res, '_harmony'))
# Build node tree:
seuratobject <- BuildClusterTree(seuratobject)
seuratobject@cluster.tree
PlotClusterTree(seuratobject)
###############################################################################
# Markers
seuratobject <- SetAllIdent(seuratobject, id = 'Clustering_res.0.35_harmony')
markers.res0.35 <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.25, min.diff.pct=0.2)
markers.res0.35$FC <- exp(markers.res0.35$avg_logFC)
saveRDS(markers.res0.35, "markers.res0.35.RDS")
markers.res0.35 <- markerDescription(markers.res0.35, verbose = T)
write.table(markers.res0.35, file = "markers.res0.35.txt", sep = "\t", row.names = F)
saveRDS(markers.res0.35, 'markers.res0.35.RDS')
###############################################################################
markers.cluster8 <- FindMarkers(seuratobject, ident.1 = 8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.cluster8$FC <- exp(markers.cluster8$avg_logFC)
saveRDS(markers.cluster8, 'markers.cluster8.RDS')
markers.cluster8$gene <- rownames(markers.cluster8)
markers.cluster8 <- markerDescription(markers.cluster8, verbose = T)
write.table(markers.cluster8, file = "markers.cluster8.txt", sep = "\t", row.names = F)
saveRDS(markers.cluster8, 'markers.cluster8.RDS')
# Based on cluster tree, gene expression and DE genes; cluster 8 was merged into cluster 0:
seuratobject <- SetIdent(seuratobject, cells.use = WhichCells(seuratobject, ident = 8), ident.use = 0)
levels(seuratobject@ident)
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, 'Clustering_res.0.35_harmony_merged')
# Identify clusters:
FeaturePlot(seuratobject, "CD83", reduction.use = "umap")
FeaturePlot(seuratobject, "MKI67", reduction.use = "umap")
FeaturePlot(seuratobject, "KRT1", reduction.use = "umap")
FeaturePlot(seuratobject, "KRT5", reduction.use = "umap")
seuratobject <- SetIdent(seuratobject, cells.use = WhichCells(seuratobject, ident = c(1, 3)), ident.use = "Premitotic_KC")
seuratobject <- SetIdent(seuratobject, cells.use = WhichCells(seuratobject, ident = c(6)), ident.use = "Mitotic_KC")
seuratobject <- SetIdent(seuratobject, cells.use = WhichCells(seuratobject, ident = c(5)), ident.use = "CD83_KC")
seuratobject <- SetIdent(seuratobject, cells.use = WhichCells(seuratobject, ident = c(0, 2, 4, 7)), ident.use = "Postmitotic_KC")
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, 'harmony_annotation_11')
saveRDS(seuratobject, file = "keratinocyte.RDS") # saved
#0e6c8b # postmitotic
#b8bae5 # premitotic
#87AC34 # mitotic (green)
#E87D72 # cd83 (red)
kc.colours <- c("#E87D72", "#87AC34", "#0e6c8b", "#b8bae5")

DimPlot(object=seuratobject, reduction.use='umap', do.label=T, group.by = "harmony_annotation_11", cols.use = kc.colours)

png(file.path(plotDir, paste0('umap_harmony_annotation_11__nolabel.png')), width=1000, height=1000)
DimPlot(object=seuratobject, reduction.use='umap', do.label=F, no.legend = T, no.axes = T, group.by='harmony_annotation_11', cols.use = kc.colours)
dev.off()

png(file.path(plotDir, paste0('umap_harmony_annotation_11.png')), width=1000, height=1000)
DimPlot(object=seuratobject, reduction.use='umap', do.label=F, no.legend = F, no.axes = F, group.by='harmony_annotation_11', cols.use = kc.colours)
dev.off()

###############################################################################
# Markers: clustering
seuratobject <- SetAllIdent(seuratobject, id = 'harmony_annotation_11')

markers.harmony_annotation_11 <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.25, min.diff.pct=0.2)

markers.harmony_annotation_11$FC <- exp(markers.harmony_annotation_11$avg_logFC)
saveRDS(markers.harmony_annotation_11, "markers.harmony_annotation_11.RDS")

markers.harmony_annotation_11 <- markerDescription(markers.harmony_annotation_11, verbose = T)
write.table(markers.harmony_annotation_11, file = "markers.harmony_annotation_11.txt", sep = "\t", row.names = F)
saveRDS(markers.harmony_annotation_11, 'markers.harmony_annotation_11.RDS')
###############################################################################
# Markers: clustering after merge of cluster #8:
seuratobject <- SetAllIdent(seuratobject, id = 'Clustering_res.0.35_harmony_merged')

markers.res0.35_merged <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)

markers.res0.35_merged$FC <- exp(markers.res0.35_merged$avg_logFC)
saveRDS(markers.res0.35_merged, "markers.res0.35_merged.RDS")

markers.res0.35_merged <- markerDescription(markers.res0.35_merged, verbose = T)
write.table(markers.res0.35_merged, file = "markers.res0.35_merged.txt", sep = "\t", row.names = F)
saveRDS(markers.res0.35_merged, 'markers.res0.35_merged.RDS')

###############################################################################
Convert(seuratobject, "anndata", filename = "kc.h5ad")
# Run a1_keratinocyte.py
