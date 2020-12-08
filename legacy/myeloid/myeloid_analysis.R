## ----setup, include=FALSE------------------------------------------------
library(Seurat)
library(reticulate)
library(dplyr)
library(magrittr)
library(harmony)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(BiocParallel)

setwd("~/Desktop/Skin DC analysis")
subset_myeloid_seurat <- readRDS("../data/subset_myeloid_seurat.RDS")

subset_myeloid_seurat <- NormalizeData(object = subset_myeloid_seurat)

subset_myeloid_seurat <- FindVariableGenes(object = subset_myeloid_seurat, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.125, x.high.cutoff = 5, y.cutoff = 0.5)

subset_myeloid_seurat <- ScaleData(object = subset_myeloid_seurat)

subset_myeloid_seurat <- RunPCA(object = subset_myeloid_seurat, pc.genes = subset_myeloid_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

subset_myeloid_seurat %<>% RunHarmony("Sample", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 20, max.iter.harmony = 5)

subset_myeloid_seurat %<>% RunUMAP(reduction.use = "harmony", dims.use = 1:12, do.fast = T)
subset_myeloid_seurat %<>% FindClusters(reduction.type = "harmony", resolution = 0.85, dims.use = 1:12, force.recalc = TRUE, print.output = FALSE)
DimPlot(subset_myeloid_seurat, reduction.use = "umap", do.label=T)

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)
new.cluster.ids <- c("moDC2", "Macrophages1", "moDC3", "LC3", "LC2", "Macrophage2", "moDC1", "LC1", "IL23_DC", "mig_cDC", "DC2", "LC4", "KLF10_LC", "DC1")
subset_myeloid_seurat@ident <- plyr::mapvalues(x = subset_myeloid_seurat@ident, from = current.cluster.ids, to = new.cluster.ids)

my_levels <- c("KLF10_LC", "LC1", "LC2", "LC3", "LC4", "IL23_DC", "DC1", "DC2", "mig_cDC", "Macrophages1", "Macrophages2", "moDC1", "moDC2", "moDC3")

subset_myeloid_seurat@ident <- factor(x = subset_myeloid_seurat@ident, levels = my_levels)

subset_myeloid_seurat <- StashIdent(object = subset_myeloid_seurat, save.name = "combined")



## ------------------------------------------------------------------------

a<-read.csv("mouse_CCR7_genes.cgi",header=T,sep='\t')
b<-((a[!(a$logFC<1),]))
c<-((b[!(b$adj.P.Val>0.05),]))
d<-as.character(c$Gene.symbol)

 
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "external_gene_name", values = d , mart = mouse, attributesL = c("external_gene_name"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])

subset_myeloid_seurat<-AddModuleScore(object=subset_myeloid_seurat,genes.list = list(humanx),ctrl.size=5,enrich.name='mouse_CCR7')

my_levels <- c("moDC3", "moDC2", "moDC1", "Monomac", "Macrophages", "DC2", "DC1", "IL23_DC", "LC4", "LC3", "LC2", "LC1", "KLF10_LC")

subset_myeloid_seurat@ident <- factor(x = subset_myeloid_seurat@ident, levels = my_levels)

DotPlot(object = subset_myeloid_seurat, genes.plot = c("mouse_CCR71"),plot.legend = TRUE, scale.by='radius', cols.use = c("black","red"), dot.min=0, do.return = T)


## ------------------------------------------------------------------------

a<-read.csv("XCR1_DC.txt",header=T,sep='\t')
b<-((a[(a$Geneset=="common_all"),]))
c<-as.character(b$Gene)

 
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "external_gene_name", values = c, mart = mouse, attributesL = c("external_gene_name"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])

subset_myeloid_seurat<-AddModuleScore(object=subset_myeloid_seurat,genes.list = list(humanx),ctrl.size=5,enrich.name='mouse_DC1')

DotPlot(object = subset_myeloid_seurat, genes.plot = c("mouse_DC11"),plot.legend = TRUE, scale.by='radius', cols.use = c("black","red"), dot.min=0, do.return = T)


## ------------------------------------------------------------------------

tonsil_CCR7_markers<-read.csv("./tonsil_CCR7_markers.csv", header=TRUE, row.names=1)
tonsil_CCR7_markers$gene<-row.names(tonsil_CCR7_markers)
tonsil_CCR7_markers<-tonsil_CCR7_markers[!(tonsil_CCR7_markers$p_val_adj>0.05),]
tonsil_CCR7_markers<-list(tonsil_CCR7_markers$gene)

ascites_CCR7_markers<-read.csv("./ascites_CCR7_markers.csv", header=TRUE, row.names=1)
ascites_CCR7_markers$gene<-row.names(ascites_CCR7_markers)
ascites_CCR7_markers<-ascites_CCR7_markers[!(ascites_CCR7_markers$p_val_adj>0.05),]
ascites_CCR7_markers<-list(ascites_CCR7_markers$gene)

SF_CCR7_markers<-read.csv("./SF_subset_markers.csv", header=TRUE, row.names=1)
SF_CCR7_markers<-subset(SF_CCR7_markers, cluster == "7")
SF_CCR7_markers$gene<-row.names(SF_CCR7_markers)
SF_CCR7_markers<-SF_CCR7_markers[!(SF_CCR7_markers$p_val_adj>0.05),]
SF_CCR7_markers<-list(SF_CCR7_markers$gene)

subset_myeloid_seurat<-AddModuleScore(object=subset_myeloid_seurat, genes.list = tonsil_CCR7_markers, ctrl.size=5, enrich.name='tonsil_CCR7')

subset_myeloid_seurat<-AddModuleScore(object=subset_myeloid_seurat, genes.list = ascites_CCR7_markers, ctrl.size=5, enrich.name='ascites_CCR7')

subset_myeloid_seurat<-AddModuleScore(object=subset_myeloid_seurat, genes.list = SF_CCR7_markers, ctrl.size=5, enrich.name='SF_CCR7')


Dotplot_other_tissue_labelled<-DotPlot(object = subset_myeloid_seurat, genes.plot = c("tonsil_CCR71", "ascites_CCR71", "SF_CCR71"), plot.legend = TRUE,x.lab.rot = T, cols.use = c("deepskyblue", "firebrick2"), dot.scale = 2.2, dot.min = 0.1, do.return = T)


## ------------------------------------------------------------------------

moDC<-FindMarkers(subset_myeloid_seurat, ident.1="moDC3",ident.2="moDC1", only.pos=T, min.pct=0.25)
LC<-FindMarkers(subset_myeloid_seurat, ident.1="LC4",ident.2="LC1", only.pos=T, min.pct=0.25)
cDC<-FindMarkers(subset_myeloid_seurat, ident.1="mig_cDC",ident.2=c("DC1","DC2"), only.pos=T, min.pct=0.25)

library(VennDiagram)

LC<-LC[!(LC$p_val_adj>0.05),]
moDC<-moDC[!(moDC$p_val_adj>0.05),]
cDC<-cDC[!(cDC$p_val_adj>0.05),]

venn.diagram(x<-list("moDC"=row.names(moDC),"LC"=row.names(LC),"mDC"=row.names(cDC)),"plot_venn",height = 3000, width = 4000, resolution =500, imagetype = "tiff",fill=rainbow(3),cat.col="white", cex=1.5)

common<-Reduce(intersect, list("moDC"=row.names(moDC),"LC"=row.names(LC),"mDC"=row.names(cDC)))
write.csv(common, "common_CCR7.csv")

