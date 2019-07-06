exportSeurat <- function(seuratobject, fileprefix = "export_", overwrite = F) {
  library("Matrix")

  filenames <- c("rownames.txt", "colnames.txt", "metadata.txt", "data.txt")
  filenames <- paste0(fileprefix, filenames)

  if (!overwrite & any(file.exists(filenames))) {
    stop("Outfile already exists.")
  }

  if (!identical(dim(seuratobject@raw.data), dim(seuratobject@data))) {
    warning("Raw and processed data dimensions are not the same.")
  }

  write(rownames(seuratobject@raw.data), file = filenames[1])
  write(colnames(seuratobject@raw.data), file = filenames[2])
  write.table(seuratobject@meta.data,
    file = filenames[3], row.names = T, col.names = T, sep = "\t", quote = F)
  out.mtx <- t(seuratobject@raw.data)
  writeMM(out.mtx, file=filenames[4])

}
