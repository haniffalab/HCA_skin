seuratobject <- readRDS("../data/subset_nonimmune_seurat.RDS")
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
seuratobject <- RunUMAP(seuratobject, reduction.use = "harmony", dims.use=1:20, reduction.name="umap", reduction.key="umap", seed.use=15)
###############################################################################
seuratobject <- FindClusters(seuratobject, resolution=1.2, dims.use=1:maxPC, force.recalc=T, save.SNN=T, reduction.type="harmony")
seuratobject <- FindClusters(seuratobject, resolution=1.5, dims.use=1:maxPC, reuse.SNN=T, reduction.type="harmony")
seuratobject <- FindClusters(seuratobject, resolution=2, dims.use=1:maxPC, reuse.SNN=T, reduction.type="harmony")
saveRDS(seuratobject, file = "nonimmune.RDS")
seuratobject <- SetAllIdent(seuratobject, id='res.1.2')
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, 'Clustering_res.1.2_harmony')
###############################################################################
# Markers: clustering
seuratobject <- SetAllIdent(seuratobject, id = 'Clustering_res.1.2_harmony')
markers.res1.2 <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.25, min.diff.pct=0.2)
markers.res1.2$FC <- exp(markers.res1.2$avg_logFC)
markers.res1.2 <- markerDescription(markers.res1.2, verbose = T)
write.table(markers.res1.2, file = "markers.res1.2.txt", sep = "\t", row.names = F)
saveRDS(markers.res1.2, 'markers.res1.2.RDS')
###############################################################################
# Add annotations:
celltypes <- seuratobject@ident
library(plyr)
celltypes <- mapvalues(celltypes,
  from=c(0:19),
  to=c("VE2", "VE2", "VE1", "VE1", "F3", "Pericyte", "F1", "VE1", "F1", "F1", "LE2", "VE1", "F2", "VE2", "F2", "LE1", "F2", "VE3", "Mel", "Stroma_Schwann"))
seuratobject <- AddMetaData(seuratobject, celltypes, 'harmony_annotation_10')
seuratobject <- SetAllIdent(seuratobject, id='harmony_annotation_10')
Schwann <- DimPlot(seuratobject, group.by='harmony_annotation_10', do.label=T, reduction.use = "tsne",do.identify = T)
seuratobject <- SetIdent(seuratobject, cells.use = Schwann, ident.use = "Schwann")
seuratobject <- AddMetaData(seuratobject, seuratobject@ident, 'harmony_annotation_11')
seuratobject <- SetAllIdent(seuratobject, id='harmony_annotation_11')
###############################################################################
# Markers: Schwann
seuratobject <- SetAllIdent(seuratobject, id = "harmony_annotation_11")
markers.Schwann <- FindMarkers(seuratobject, ident.1 = "Schwann", only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25, min.diff.pct = 0.2)
markers.Schwann$FC <- exp(markers.Schwann$avg_logFC)
markers.Schwann$gene <- rownames(markers.Schwann)
markers.Schwann <- markerDescription(markers.Schwann, verbose = T)
write.table(markers.Schwann, file = "markers.Schwann.txt", sep = "\t", row.names = F)
saveRDS(markers.Schwann, 'markers.Schwann.RDS')
schwann.markers <- c("PLP1", "SOX10", "ERBB3", "NGFR", "MPZ")
schwann.markers %in% markers.Schwann$gene
markers.Schwann[schwann.markers, ]
###############################################################################
celltypes <- seuratobject@ident
library(plyr)
celltypes <- mapvalues(celltypes,
                       from=c("F1", "F2", "F3", "Pericyte", "Schwann", "Stroma_Schwann", "Mel", "LE1", "LE2", "VE1", "VE2", "VE3"),
                       to=c("A_F1", "B_F2", "C_F3", "D_Pericyte", "E_Schwann", "F_Stroma Schwann", "G_Melanocyte", "H_LE1", "I_LE2", "J_VE1", "K_VE2", "L_VE3"))
seuratobject <- AddMetaData(seuratobject, celltypes, 'Alphabetical_11')
seuratobject <- SetAllIdent(seuratobject, id='Alphabetical_11')
###############################################################################
# Markers: Alphabetical_11
seuratobject <- readRDS(file = "nonimmune.RDS")
seuratobject <- SetAllIdent(seuratobject, id='Alphabetical_11')
markers.Alphabetical_11 <- FindAllMarkers(seuratobject, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.25, min.diff.pct=0.2)
markers.Alphabetical_11$FC <- exp(markers.Alphabetical_11$avg_logFC)
saveRDS(markers.Alphabetical_11, "markers.Alphabetical_11.RDS")
markers.Alphabetical_11 <- readRDS("markers.Alphabetical_11.RDS")
markers.Alphabetical_11 <- markerDescription(markers.Alphabetical_11, verbose = T)
write.table(markers.Alphabetical_11, file = "markers.Alphabetical_11.txt", sep = "\t", row.names = F)
saveRDS(markers.Alphabetical_11, 'markers.Alphabetical_11.RDS')
###############################################################################
# Markers: LE1 vs LE2
seuratobject <- SetAllIdent(seuratobject, id = "Alphabetical_11")

markers.LE1vsLE2 <- FindMarkers(seuratobject, ident.1 = "H_LE1", ident.2 = "I_LE2", only.pos = FALSE, min.pct = 0.3, logfc.threshold = 0.25)
markers.LE1vsLE2$FC <- exp(markers.LE1vsLE2$avg_logFC)
saveRDS(markers.LE1vsLE2, 'markers.LE1vsLE2.RDS')

markers.LE1vsLE2$gene <- rownames(markers.LE1vsLE2)
markers.LE1vsLE2 <- markerDescription(markers.LE1vsLE2, verbose = T)
write.table(markers.LE1vsLE2, file = "markers.LE1vsLE2.txt", sep = "\t", row.names = F)
saveRDS(markers.LE1vsLE2, 'markers.LE1vsLE2.RDS')

"SOD2" %in% markers.LE1vsLE2$gene
###############################################################################
