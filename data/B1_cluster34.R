seuratobject <- readRDS("subset_plasmamast_seurat.RDS")
dir.create(plotDir)
dir.create(geneDir)
dir.create(umapDir)
seuratobject <- NormalizeData(seuratobject)
seuratobject <- ScaleData(seuratobject)
seuratobject <- FindVariableGenes(seuratobject, x.low.cutoff=1, x.high.cutoff=Inf, y.cutoff=1, y.high.cutoff=Inf)
sort(seuratobject@var.genes)
length(seuratobject@var.genes)
seuratobject <- RunPCA(seuratobject, pc.genes=seuratobject@var.genes, do.print=TRUE, pcs.print=1:5, genes.print=5)
library("harmony")
seuratobject <- RunHarmony(object = seuratobject, group.by.vars = "Sample", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 20, max.iter.harmony = 5)
maxPC <- 20
seuratobject <- RunTSNE(seuratobject, reduction.use="harmony", dims.use=1:maxPC, do.fast=T, reduction.name="tsne", reduction.key="tsne")
FeaturePlot(seuratobject, features.plot = c("XBP1", "MS4A1"))
FeaturePlot(seuratobject, features.plot = c("IGHG1", "CD79A"))
FeaturePlot(seuratobject, features.plot = c("IGHA1", "AL928768.3"))
# See also :"CD79A", "CD79B", "IGKC", "JCHAIN", "IGHA1", "CD40LG"
PCElbowPlot(seuratobject, num.pc = 20)
PrintPCA(seuratobject, pcs.print = 1:20, genes.print = 5)
seuratobject <- RunUMAP(seuratobject, reduction.use = "harmony", dims.use=1:20, reduction.name="umap", reduction.key="umap", seed.use = 2)
###############################################################################
seuratobject <- FindClusters(seuratobject, resolution=0.5, dims.use=1:maxPC, force.recalc=T, save.SNN=T, reduction.type="harmony")
seuratobject <- SetAllIdent(seuratobject, id='res.0.5')
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, col.name = 'C34_clustering_res0.5')
###############################################################################
# Markers: clustering
seuratobject <- SetAllIdent(seuratobject, id = 'C34_clustering_res0.5')
markers.res0.5 <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.25, min.diff.pct=0.2)
markers.res0.5$FC <- exp(markers.res0.5$avg_logFC)
markers.res0.5 <- readRDS("markers.res0.5.RDS")
markers.res0.5 <- markerDescription(markers.res0.5, verbose = T)
write.table(markers.res0.5, file = "markers.res0.5.txt", sep = "\t", row.names = F)
saveRDS(markers.res0.5, 'markers.res0.5.RDS')
###############################################################################
# Add annotations:
celltypes <- seuratobject@ident
library(plyr)
celltypes <- mapvalues(celltypes,
  from=c(0:3),
  to=c("Mast", "Plasma_cell", "Plasma_cell", "Plasma_cell"))
seuratobject <- AddMetaData(seuratobject, celltypes, 'harmony_annotation_11')
seuratobject <- SetAllIdent(seuratobject, id='harmony_annotation_11')
saveRDS(seuratobject, 'subset_plasmamast.RDS')
