if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu95av2.db")

library(qusage,rlist)#Category)
library(Category)
library(rlist)
library(pipeR)
library(biomaRt)
library(plyr)

################################
#Get list of all TFS
tf_list<-read.delim("/Users/Documents/Rscripts/List of all transcription factors/Homo_sapiens_transcription_factors_gene_list.txt", header = TRUE, sep = "\t", dec = ".")
tf_names <- unique(as.character(tf_list$Ensembl.ID))
  
###Check chromosomal location of all human TFS
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes<-listAttributes(ensembl)
all_tf_positions<-getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand','chromosomal_region',"hgnc_symbol","ensembl_gene_id"),
                   filters=c('ensembl_gene_id'),
                   values=tf_names,
                   mart=ensembl)
write.csv(all_tf_positions,"~/Downloads/all_human_tf_positions.csv")

## Read Gary's Rdata instead of csv::
#try = g

#Read DEGS
g = read.csv("~/Downloads/downs_DEGs.csv",row.names = 1)
g <- reshape2::melt(g,measure.vars = c(colnames(g)))
colnames(g)<-c("cell_label","gene")
concat<-g
concat<-na.omit(concat)

#for (i in colnames(g)){
#  print(i)
#  temp = g[i]
#  temp<-na.omit(temp)
#  colnames(temp)<-"gene"
#  temp$cell_label<-i
#  
#  if(!"concat"%in%ls()){
#    concat <- temp 
#  }else{
#    concat <- rbind(concat,temp)
#  }
#}

###Check chromosomal location of all genes
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
TF_position<-getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand',"hgnc_symbol"),
                   filters=c('hgnc_symbol'),
                   values=concat$gene,
                   mart=ensembl)
#Push locations back to gene list
concat$chromosome <- "NA"
concat$chromosome <- mapvalues(concat$gene,TF_position$hgnc_symbol,TF_position$chromosome_name,)

##Get only genes which are in chrom 19
#concat<-concat[(concat$chromosome==19),]
#concat

##TF predictions acquired from Encode database (https://amp.pharm.mssm.edu/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets)
#TF gene list from "https://amp.pharm.mssm.edu/Enrichr/#stats"
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X","~/Downloads/chea_encode_consesnus_tfs")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=ENCODE_TF_ChIP-seq_2015","~/Downloads/chip_seq_tfs")


# Choose option for database:
db = 2
{
#option1 -- chea-encode consensus 108 tf
if (db == 1){
x<- read.delim("~/Downloads/chea_encode_consesnus_tfs", header = FALSE, sep = "\t", dec = ".",stringsAsFactors = FALSE)
x$V2<-NULL
x$V1<-gsub("(.*) .*","\\1",x$V1)
x$V1<-gsub(" ", "",x$V1)
#x.list <- setNames(split(x, seq(nrow(x))), rownames(x$V1))
x.list <- as.list(as.data.frame(t(x)))
names(x.list) <- sapply(x.list, `[[`, 1)
x.list <- lapply(x.list, `[`, -1)
x.list<-lapply(x.list, function(z){ z[!is.na(z) & z != ""]})
database <- x.list
}
#option2 <800 tf-gene associations
if (db == 2){
x<- read.delim("~/Downloads/chip_seq_tfs", header = FALSE, sep = "\t", dec = ".",stringsAsFactors = FALSE)
x$V2<-NULL
x$V1<-sub(" .*", "", x$V1)
x$V1<-sub("[.].*", "", x$V1)
#x$V1<-gsub(" ", "",x$V1)
#x.list <- setNames(split(x, seq(nrow(x))), rownames(x$V1))
x <- data.frame(lapply(x, as.character), stringsAsFactors=FALSE)
x.list <- as.list(as.data.frame(t(x)))
names(x.list) <- sapply(x.list, `[[`, 1)
x.list <- lapply(x.list, `[`, -1)
x.list<-lapply(x.list, function(z){ z[!is.na(z) & z != ""]})
database <- x.list
}
#option3 -- 181 curated gene-tf associations
if (db == 3){
database <- read.gmt("/Users/Downloads/gene_set_library_crisp.gmt")
}
}

test = unlist(concat$gene)
test1 = unlist(database)
remove <- as.character(test[!test%in%test1])


#Remove genes not in universe
concat <- concat[(!concat$gene%in%remove),]

#Predict overall tfs regulating all DE genes
try<-list(concat$gene)
sigsets <- Map(function(x, y) x[x %in% y], database, MoreArgs=list(y=as.character(unlist(try))))
universe <- unlist(database, use.names=FALSE)
overall_tfs <- hyperg(database, sigsets, universe)
overall_tfs$tf <- rownames(overall_tfs)
rownames(overall_tfs)<-NULL
overall_tfs<-overall_tfs[(overall_tfs$p<0.05)&(overall_tfs$p!=0),]
overall_tfs$tf<-sub("[.].*", "", overall_tfs$tf)

###Check chromosomal location of overall TFS
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
TF_position<-getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand',"hgnc_symbol"),
                   filters=c('hgnc_symbol'),
                   values=overall_tfs$tf,
                   mart=ensembl)
#Push locations back to gene list
overall_tfs$chromosome <- "NA"
overall_tfs$chromosome <- mapvalues(overall_tfs$tf,TF_position$hgnc_symbol,TF_position$chromosome_name)


#Find TFS for each celltype
#Initialise empty frame
results_concat <- data.frame(module=character(),tf=character(),assayed=character(), significant=character(), universe=character(), representation=character(), p=character(), odds=character(), expected=character(),stringsAsFactors=FALSE)
for (z in unique(concat$cell_label)){
  print(z)
  temp_gene <- concat[(concat$cell_label==z),]
  universe <- unlist(database, use.names=FALSE)
  try<-list(temp_gene$gene)
  sigsets <- Map(function(x, y) x[x %in% y], database, MoreArgs=list(y=as.character(unlist(try))))
  temp_celltype_tf <- hyperg(database, sigsets, universe)
  temp_celltype_tf$module <- z
  temp_celltype_tf$tf <- rownames(temp_celltype_tf)
  temp_celltype_tf<-temp_celltype_tf[(temp_celltype_tf$p<0.05)&(temp_celltype_tf$p!=0),]
  row.names(temp_celltype_tf) <- NULL
  results_concat<- rbind(temp_celltype_tf,results_concat)
}

results_export <- results_concat[, c("module","tf","assayed", "p", "significant", "universe", "representation", "odds", "expected")]
results_export <- as.data.frame(t(apply(results_export, 1, unlist)))

###Check chromosomal location of all TFS
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
TF_position<-getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand',"hgnc_symbol"),
                   filters=c('hgnc_symbol'),
                   values=results_export$tf,
                   mart=ensembl)
#Push locations back to gene list
results_export$chromosome <- "NA"
results_export$chromosome <- mapvalues(results_export$tf,TF_position$hgnc_symbol,TF_position$chromosome_name)


write.csv(results_export,"/Users/Documents/Projects/Fetal Bone Marrow/network_gsea/HSC_adult_fetal/ABM_F_comp_merged_networks_summary/FBM_DOWNS_DE_TF_regs.csv")


#Check for expression
expressed_var <- read.csv("/Users/Downloads/genes_expressed.csv")
results_concat[(!results_concat$tf%in%expressed_var$X),]

#Get top 5 tfs for each module by pvalue
results_export<-as.data.frame(results_export)
filtered_export<-results_export %>%
  group_by(module) %>%
  top_n(-5, p)
unique(filtered_export$tf)

write.csv(filtered_export,"/Users/Documents/Projects/Fetal Bone Marrow/network_gsea/HSC_adult_fetal/ABM_F_comp_merged_networks_summary/filtered_Downs_FBM_TF_regs.csv")

