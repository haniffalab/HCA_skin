library(Seurat)
library(dplyr)
library(reticulate)
library(Matrix)
library(stringr)
library(tximport)
library(readr)
library(biomaRt)
library(harmony)
library(dplyr)
library(ggplot2)
library(reshape)

###Import adult skin immune cells as a Seurat object

setwd("~/hcaskin/blood")
adult_immune<-readRDS("adult_immune_Seurat.RDS")

###Import fetal skin immune cells as a Seurat object

setwd("~/hcaskin/SKIN_FETAL_mtx")
fetal_immune<-readRDS("fetal_immune_Seurat.RDS")

###Perform standard pre-processing

fetal_immune<-NormalizeData(fetal_immune)
fetal_immune<-FindVariableFeatures(fetal_immune, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)

adult_immune<- NormalizeData(adult_immune)
adult_immune<-FindVariableFeatures(adult_immune, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)

###Integrate both datasets

anchors <- FindTransferAnchors(reference = adult_immune, query = fetal_immune, reduction="pcaproject",
    dims = 1:30)

predictions <- TransferData(anchorset = anchors, refdata = adult_immune$final, 
    dims = 1:30)

fetal_immune <- AddMetaData(fetal_immune, metadata = predictions)

###Prepare data for plotting (Figure 1E)

a<-fetal_immune@meta.data
a$orig.ident<-NULL
a$nCount_RNA<-NULL
a$nFeature_RNA<-NULL
a$index<-NULL
a$sample_id<-NULL
a$sort_gate<-NULL
a$COHORT<-NULL
a$Status<-NULL
a$DNA.SOURCE<-NULL
a$Organ<-NULL
a$predicted.id<-NULL
a$prediction.score.max<-NULL
a$Site<-NULL
a$Tissue<-NULL
a$Enrichment<-NULL
a$Location<-NULL
a$Donor.ID<-NULL
a$Sex<-NULL
a$mad_prd<-NULL
a$stage<-NULL


data<-aggregate(.~anno_final,data=a,FUN=median)

colnames(data)<-c("clustering", "Th", "Tc", "LC","Treg", "NK", "Macro1", "ILC1_NK", "ILC2_3", "ILC1", "Macro2",
                 "DC2", "Infmono", "Mono", "moDC", "migDC", "DC1", "Mast_cell")

data<-melt(data, id = c("clustering"))

colnames(data)<-c("clustering","value","variable")

###Reorder for plotting
data$value <- factor(data$value, levels = rev(c("ILC1", "ILC1_NK", "NK", "ILC2_3", "Tc", "Th", "Treg", "Mast_cell",
                                             "Macro1", "Macro2", "Infmono", "DC1", "DC2", "LC", "Mono",
                                             "migDC", "moDC")))

data$clustering <- factor(data$clustering, levels = (c("fs_ILC","fs_NK", "fs_Mast cell", "fs_Macrophage",
                                                      "fs_DC1", "fs_DC2", "fs_LC", "fs_Monocyte")))

plot1<-ggplot(data = data, mapping = aes(x = clustering, y = value, fill = variable)) + 
geom_tile(width=0.975, height=0.975) + theme_bw() + coord_equal() + theme(axis.text.x = element_text(angle = 90)) +
scale_fill_gradient(low = "#568dba", high = "#f42222")
ggsave("fetal_adult_immune_comparison.pdf", plot = plot1, device = NULL, path = NULL,
  scale = 1, dpi = 300)
plot1

setwd("~/hcaskin/SKIN_FETAL_mtx")

fetal_nonimmune_counts<-readMM("fetal_nonimmune_counts.mtx")
fetal_nonimmune_genes<-read.csv("fetal_nonimmune_genes.csv")
fetal_nonimmune_metadata<-read.csv('fetal_nonimmune_metadata.csv')
fetal_nonimmune_cells<-read.csv('fetal_nonimmune_cell_names.csv')

row.names(fetal_nonimmune_counts)<-fetal_nonimmune_genes$gene_short_name
colnames(fetal_nonimmune_counts)<-fetal_nonimmune_cells$index

fetal_nonimmune <- CreateSeuratObject(counts = fetal_nonimmune_counts, project = "skin", min.cells = 3, min.features = 200)
row.names(fetal_nonimmune_metadata)<-fetal_nonimmune_metadata$index
fetal_nonimmune <- AddMetaData(fetal_nonimmune, fetal_nonimmune_metadata)
fetal_nonimmune@meta.data["stage"]<-"fetal"

saveRDS(fetal_nonimmune, "fetal_nonimmune_Seurat.RDS")

setwd("~/hcaskin/blood")

adult_nonimmune_counts<-readMM("adult_nonimmune_counts.mtx")
adult_nonimmune_genes<-read.csv("adult_nonimmune_genes.csv")
adult_nonimmune_metadata<-read.csv('adult_nonimmune_metadata.csv')
adult_nonimmune_cells<-read.csv('adult_nonimmune_cell_names.csv')

row.names(adult_nonimmune_counts)<-adult_nonimmune_genes$gene_short_name
colnames(adult_nonimmune_counts)<-adult_nonimmune_cells$index

adult_nonimmune <- CreateSeuratObject(counts = adult_nonimmune_counts, project = "skin", min.cells = 3, min.features = 200)
row.names(adult_nonimmune_metadata)<-adult_nonimmune_metadata$index
adult_nonimmune <- AddMetaData(adult_nonimmune, adult_nonimmune_metadata)
adult_nonimmune@meta.data["stage"]<-"adult"

saveRDS(adult_nonimmune, "adult_nonimmune_Seurat.RDS")

setwd("~/hcaskin/blood")
adult_nonimmune<-readRDS("adult_nonimmune_Seurat.RDS")

Idents(adult_nonimmune)<-"final"

adult_nonkc<-subset(adult_nonimmune, idents = c("Pre-proliferation_KC","Post_proliferation_KC"), invert=TRUE)
adult_kc<-subset(adult_nonimmune, idents = c("Pre-proliferation_KC","Post_proliferation_KC"))
kc.downsample = subset(adult_kc, cells = sample(Cells(adult_kc), 6000))
adult_nonimmune<-merge(adult_nonkc, y=kc.downsample)

setwd("~/hcaskin/SKIN_FETAL_mtx")
fetal_nonimmune<-readRDS("fetal_nonimmune_Seurat.RDS")

Idents(fetal_nonimmune) <- "anno_final"

fetal_nonimmune<-NormalizeData(fetal_nonimmune)
fetal_nonimmune<-FindVariableFeatures(fetal_nonimmune, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)

adult_nonimmune<- NormalizeData(adult_nonimmune)
adult_nonimmune<-FindVariableFeatures(adult_nonimmune, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)

anchors <- FindTransferAnchors(reference = adult_nonimmune, query = fetal_nonimmune, reduction="pcaproject",
    dims = 1:30)

predictions <- TransferData(anchorset = anchors, refdata = adult_nonimmune$final, 
    dims = 1:30)

fetal_nonimmune <- AddMetaData(fetal_nonimmune, metadata = predictions)

a<-fetal_nonimmune@meta.data
setwd("~/hcaskin/SKIN_FETAL_mtx")
write.csv(a, "fetal_adult_nonimmune_Seurat_comparison.csv")

a<-fetal_nonimmune@meta.data
a$orig.ident<-NULL
a$nCount_RNA<-NULL
a$nFeature_RNA<-NULL
a$index<-NULL
a$sample_id<-NULL
a$sort_gate<-NULL
a$COHORT<-NULL
a$Status<-NULL
a$DNA.SOURCE<-NULL
a$Organ<-NULL
a$predicted.id<-NULL
a$prediction.score.max<-NULL
a$Site<-NULL
a$Tissue<-NULL
a$Enrichment<-NULL
a$Location<-NULL
a$Donor.ID<-NULL
a$Sex<-NULL
a$mad_prd<-NULL
a$stage<-NULL


data<-aggregate(.~anno_final,data=a,FUN=median)

colnames(data)<-c("clustering", "Melanocyte", "LE2", "F2", "LE1", "F1", "F3", "VE1", "VE2", "Peri1",
                  "Peri2", "Schwann2","VE3", "Schwann1", "Post_prolif_KC", "Pre-prolif_KC")

data<-melt(data, id = c("clustering"))

colnames(data)<-c("clustering","value","variable")

levels(data$value)
data$value <- factor(data$value, levels = rev(c("Pre-prolif_KC","Post_prolif_KC", "Melanocyte", "Schwann1","Schwann2",
                                               "F1", "F2", "F3", "Peri1", "Peri2", "VE1", "VE2", "VE3",
                                               "LE1", "LE2")))
levels(data$value)

data$clustering<-as.factor(data$clustering)

levels(data$clustering)
data$clustering <- factor(data$clustering, levels = (c("fs_Keratinocyte", "fs_Melanocyte", "fs_Schwann cell",
                                                      "fs_Fibroblast", "fs_Smooth muscle", "fs_Vascular endothelium",
                                                      "fs_Lymphatic endothelium")))
levels(data$clustering)

setwd("~/hcaskin/blood")
plot1<-ggplot(data = data, mapping = aes(x = clustering, y = value, fill = variable)) + 
geom_tile(width=0.975, height=0.975) + theme_bw() + coord_equal() + theme(axis.text.x = element_text(angle = 90)) +
scale_fill_gradient(low = "#568dba", high = "#f42222")
ggsave("fetal_adult_nonimmune_comparison.pdf", plot = plot1, device = NULL, path = NULL,
  scale = 1, dpi = 300)
plot1
