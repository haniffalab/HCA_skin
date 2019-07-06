seuratobject <- readRDS("../data/subset_lymphoid_seurat.RDS")
dir.create(plotDir)
dir.create(geneDir)
dir.create(umapDir)
seuratobject <- NormalizeData(seuratobject)
seuratobject <- ScaleData(seuratobject)
seuratobject <- FindVariableGenes(seuratobject, x.low.cutoff=0.1, x.high.cutoff=Inf, y.cutoff=0.5, y.high.cutoff=Inf)
seuratobject <- RunPCA(seuratobject, pc.genes=seuratobject@var.genes, do.print=TRUE, pcs.print=1:5, genes.print=5)
library("harmony")
seuratobject <- RunHarmony(object = seuratobject, group.by.vars = "Sample", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 20, max.iter.harmony = 5)
maxPC <- 20
seuratobject <- RunTSNE(seuratobject, reduction.use="harmony", dims.use=1:maxPC, do.fast=T, reduction.name="tsne", reduction.key="tsne")
PCElbowPlot(seuratobject, num.pc = 20)
PrintPCA(seuratobject, pcs.print = 1:20, genes.print = 5)
seuratobject <- RunUMAP(seuratobject, reduction.use = "harmony", dims.use=1:10, reduction.name="umap", reduction.key="umap", seed.use = 2) # 10 PCs used for visualization
###############################################################################
seuratobject <- FindClusters(seuratobject, resolution=1.2, dims.use=1:maxPC, force.recalc=T, save.SNN=T, reduction.type="harmony")
seuratobject <- FindClusters(seuratobject, resolution=1.5, dims.use=1:maxPC, reuse.SNN=T, reduction.type="harmony")
seuratobject <- SetAllIdent(seuratobject, id='res.1.5')
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, col.name = 'Lymphoid_clustering_res1.5')
###############################################################################
# Markers: clustering
seuratobject <- SetAllIdent(seuratobject, id = 'Lymphoid_clustering_res1.5')
markers.res1.5 <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.25, min.diff.pct=0.2)
markers.res1.5$FC <- exp(markers.res1.5$avg_logFC)
saveRDS(markers.res1.5, "markers.res1.5.RDS")
markers.res1.5 <- readRDS("markers.res1.5.RDS")
markers.res1.5 <- markerDescription(markers.res1.5, verbose = T)
write.table(markers.res1.5, file = "markers.res1.5.txt", sep = "\t", row.names = F)
saveRDS(markers.res1.5, 'markers.res1.5.RDS')
###############################################################################
# Add annotations:
seuratobject <- SetAllIdent(seuratobject, id='Lymphoid_clustering_res1.5')
celltypes <- seuratobject@ident
library(plyr)
celltypes <- mapvalues(celltypes,
  from=c(0:19),
  to=c("Th", "Th", "Th", "Th", "Th", "Treg", "NK1", "NK1", "Tc", "Tc", "Th", "NK2", "Tc", "Treg", "Tc", "ILC", "NK3", "Th", "NK1", "Tc"))
seuratobject <- AddMetaData(seuratobject, celltypes, 'harmony_annotation_11')
seuratobject <- SetAllIdent(seuratobject, id='harmony_annotation_11')
###############################################################################
# Markers: annotation
seuratobject <- SetAllIdent(seuratobject, id = 'Lymphoid_clustering_harmony_annotation_11')
markers.harmony_annotation_11 <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.25, min.diff.pct=0.2)
markers.harmony_annotation_11$FC <- exp(markers.harmony_annotation_11$avg_logFC)
markers.harmony_annotation_11 <- markerDescription(markers.harmony_annotation_11, verbose = T)
write.table(markers.harmony_annotation_11, file = "markers.harmony_annotation_11.txt", sep = "\t", row.names = F)
saveRDS(markers.harmony_annotation_11, 'markers.harmony_annotation_11.RDS')
###############################################################################
# Create new metadata:
seuratobject@meta.data$Tissue_harmony_annotation_11 <- paste(seuratobject@meta.data$Tissue, seuratobject@meta.data$harmony_annotation_11, sep="_")
unique(seuratobject@meta.data$Tissue_harmony_annotation_11)
saveRDS(seuratobject, file='lymphoid.RDS')
###############################################################################
# Compare epidermis and dermis Tc and Th populations:
seuratobject <- SetAllIdent(seuratobject, id='Tissue_harmony_annotation_11')
# Epi_Th vs Derm_Th
markers.Epi_ThvsDerm_Th <- FindMarkers(seuratobject, ident.1='Epi_Th', ident.2='Derm_Th', only.pos=F, min.pct=0.3, logfc.threshold=0.25)
markers.Epi_ThvsDerm_Th$FC <- exp(markers.Epi_ThvsDerm_Th$avg_logFC)
markers.Epi_ThvsDerm_Th$gene <- rownames(markers.Epi_ThvsDerm_Th)
markers.Epi_ThvsDerm_Th <- markerDescription(markers.Epi_ThvsDerm_Th, verbose=T)
saveRDS(markers.Epi_ThvsDerm_Th, 'markers.Epi_ThvsDerm_Th.RDS')
write.table(markers.Epi_ThvsDerm_Th, file='markers.Epi_ThvsDerm_Th.txt', sep="\t", row.names=F)
###############################################################################
# Epi_Tc vs Derm_Tc
markers.Epi_TcvsDerm_Tc <- FindMarkers(seuratobject, ident.1='Epi_Tc', ident.2='Derm_Tc', only.pos=F, min.pct=0.3, logfc.threshold=0.25)
markers.Epi_TcvsDerm_Tc$FC <- exp(markers.Epi_TcvsDerm_Tc$avg_logFC)
markers.Epi_TcvsDerm_Tc$gene <- rownames(markers.Epi_TcvsDerm_Tc)
markers.Epi_TcvsDerm_Tc <- markerDescription(markers.Epi_TcvsDerm_Tc, verbose=T)
saveRDS(markers.Epi_TcvsDerm_Tc, 'markers.Epi_TcvsDerm_Tc.RDS')
write.table(markers.Epi_TcvsDerm_Tc, file='markers.Epi_TcvsDerm_Tc.txt', sep="\t", row.names=F)
###############################################################################
saveRDS(seuratobject, 'lymphoid.RDS') # saved
