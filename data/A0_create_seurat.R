# Metadata
metadata <- read.delim('sample_metadata.txt', header = T, sep="\t", stringsAsFactors=FALSE)
metadata
megamatrix <- NA
for (i in 1:dim(metadata)[1])
{
  #the sanger study ID, which we used up to this point to identify samples, lives in the SECOND column
  minimatrix = readRDS(file.path('../processed', metadata[i, 2], 'final-count-matrix', 'cellranger-count-matrix.RDS'))
  #add on a tag to the cell barcode to easily identify the sample from which it came from
  #this tag lives in the FIFTH column of the metadata CSV
  for (j in 1:dim(minimatrix)[2])
  {
    colnames(minimatrix)[j] = paste(metadata[i, 5], colnames(minimatrix)[j], sep='_')
  }
  #create Seurat object, imposing the standard 200 gene minimum
  minimatrix = CreateSeuratObject(raw.data = minimatrix, min.cells = 3, min.genes = 200, project = "hca")
  #add all of the metadata from the CSV file
  #you can access it later at megamatrix@meta.data[,'Foetus'], for example
  for (j in 1:dim(metadata)[2])
  {
    #copy the value as many times as we have cells, and name everything so Seurat understands it
    holder = rep(metadata[i, j], dim(minimatrix@data)[2])
    minimatrix = AddMetaData(minimatrix, setNames(holder, colnames(minimatrix@data)), colnames(metadata)[j])
  }

  #add in the mitochondrial percentage, do Seurat's normalisation and filter to no more than 20% mitochondrial
  mito.genes = grep(pattern = "^MT-", x = rownames(x = minimatrix@data), value = TRUE)
  percent.mito = Matrix::colSums(minimatrix@raw.data[mito.genes, ]) / Matrix::colSums(minimatrix@raw.data)
  minimatrix = AddMetaData(minimatrix,percent.mito, 'percent.mito')
  minimatrix = FilterCells(object = minimatrix, subset.names = "percent.mito", high.thresholds=0.2)
  minimatrix = NormalizeData(object = minimatrix, normalization.method = "LogNormalize", scale.factor = 10000)

  if (is.na(megamatrix))
  {
    megamatrix = minimatrix
  }
  else
  {
    megamatrix = MergeSeurat(object1 = megamatrix, object2 = minimatrix)
  }
}
saveRDS(megamatrix, 'all_cells.RDS')

# Subsequently, Run SVM then export data for scrublet doublet detection:
# Check if rawdata and data dims are the same (we want the svm-processed cells only)
dim(seuratobject@data) == dim(seuratobject@raw.data)
# [1] TRUE TRUE
exportSeurat(seuratobject, fileprefix = "export_skin_all_", overwrite = T)
