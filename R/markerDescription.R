# Copyright 2018 Dorin-Mirel Popescu, Peter Vegh, Newcastle University
# Takes a dataframe (markers) with a 'gene' column containing HGNC gene IDs (i.e. Seurat::FindAllMarkers() output) and returns it with added Entrez ID, Gene name (description), Cluster occurrences, Summary columns.
# Requires biomaRt, rentrez, jsonlite packages
markerDescription <- function(markers, verbose = F) {
  # Define subfunction that converts HUGO IDs to Entrez IDs:
  hgncEntrez <- function(hgnc.genes) {
    library(biomaRt)
    while(TRUE) {
      mart <- try(useDataset("hsapiens_gene_ensembl", useMart("ensembl")))
      if (class(mart) == "try-error") {
        print("Retrieving mart error. Trying again")
      } else {
        break
      }
    }
    gene.ids <- getBM(
      filters = "hgnc_symbol",
      attributes = c("hgnc_symbol", "entrezgene"),
      values = unique(hgnc.genes),
      mart = mart)
    return(gene.ids)
  }


  library(jsonlite)
  library(rentrez)

  if (!("gene" %in% colnames(markers))) {
    stop("The input dataframe should have a \"gene\" column!")
  }


  # Retrieve entrez IDs for query with rentrez
  genelist <- hgncEntrez(markers$gene)
  if (verbose) {print('Entrez IDs retrieved.')}

  # query only 400 genes at a time because querying all at once produces errors
  queryMax <- 399
  startIndex <- 1
  highestIndex <- length(genelist$entrezgene)

  summaries <- c()
  descriptions <- c()
  while (TRUE) {
    stopIndex <- min(startIndex + queryMax, highestIndex)
    if (verbose) {print(paste('Querying', stopIndex, 'genes'))}
    genes.batch <- entrez_summary(db = "gene", id = genelist$entrezgene[startIndex:stopIndex], always_return_list = T)
    summaries.batch <- extract_from_esummary(genes.batch, "summary")
    descriptions.batch <- extract_from_esummary(genes.batch, "description")
    summaries <- c(summaries, summaries.batch)
    descriptions <- c(descriptions, descriptions.batch)

    startIndex <- stopIndex + 1
    if (startIndex > highestIndex) {
      break
    }
  }

  # Add entrez IDs:
  markers$entrez <- markers$gene
  markers$entrez <- with(genelist, entrezgene[match(markers$entrez, hgnc_symbol)])


  # Add description (i.e. name):
  markers$description <- as.character(markers$entrez)
  for (index in 1:nrow(markers)) {
    markers$description[index] <- descriptions[markers$description[index]]
  }


  # Calculate occurrences:
  countClusters <- function(gene.name, clusters, genes) {
    poss <- clusters[genes == gene.name]
    return(paste(poss, collapse = ","))
  }
  # identify clusters for each marker genes
  unsym <- as.vector(unique(markers$gene))
  unsym.count <- lapply(unsym, countClusters, clusters = markers$cluster, genes = markers$gene)
  unsym.count <- unlist(unsym.count)
  names(unsym.count) <- unsym
  poss <- match(markers$gene, names(unsym.count))
  # clusters per marker genes is written in the 'occurrences' column
  occurrences <- as.vector(unsym.count)[poss]

  # Add occurrences:
  markers$occurrences <- unname(occurrences)

  # Add summaries:
  markers$summary <- as.character(markers$entrez)
  for (index in 1:nrow(markers)) {
    markers$summary[index] <- summaries[markers$summary[index]]
  }


  return(markers)
}
