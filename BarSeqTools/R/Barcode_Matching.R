defineConsensus <- function(barcodeCluster) {
  
  barcodeConsensus <- as.character(DECIPHER::ConsensusSequence(Biostrings::DNAStringSet(barcodeCluster), threshold=0.9999))
  
  barcodeConsensus <- gsub("-", "", barcodeConsensus)
  barcodeConsensus <- strsplit(barcodeConsensus, "")[[1]]
  barcodeConsensus <- sapply(barcodeConsensus, function(i) {
    x <- Biostrings::IUPAC_CODE_MAP[names(Biostrings::IUPAC_CODE_MAP)==i]
    x <- strsplit(x, "")[[1]]
    x <- sample(x, 1)
    x
  })
  barcodeConsensus <- paste0(barcodeConsensus, collapse = "")
  
  data.frame(barcodeConsensus=barcodeConsensus, n=length(barcodeCluster))
  
}

#'@export
#'@title Create consensus barcodes from a vector of barcodes
#'@param barcodes A character vector containing the barcodes
#'@param maxMistmatch The maximum number of mismatches which will be allowed to call two barcodes the same
#'@param clusterPlot Should a dendrogram of the clustering be plotted?
#'In the case of \code{NULL} (the default), there will be no plot.
#'In order to save a dendrogram, please specify a file path here.
#'This should be a pdf.
#'@param plotWidth The plot width of the dendrogram (this is usuall very large).
#'@return A \code{data.frame} is returned showing the consensus barcodes and how many times a barcode has been matched to the consensus
#'Defauts to 3
createConsensus <- function(barcodes, maxMismatch=1, clusterPlot=NULL, plotWidth=1000) {
  
  barcodes <- barcodes[order(barcodes)]
  
  tagsMatrix <- stringdist::stringdistmatrix(unique(barcodes), useNames = T, method = "lv")
  
  freq <- table(as.character(barcodes))
  
  clusters <- hclust(tagsMatrix, method="average", members = freq)
  
  barcodeAssignments <- cutree(clusters, h=maxMismatch)
  
  if(is.null(clusterPlot)==F) {
    
    clusters <- as.dendrogram(clusters)
    
    colors <- colorspace::rainbow_hcl(3)
    
    clusters <- dendextend::color_branches(clusters, h = maxMismatch, col=colors)
    clusters <- dendextend::color_labels(clusters, h = maxMismatch, col=colors)
    clusters <- dendextend::assign_values_to_branches_edgePar(clusters, 10, "lwd")

    pdf(clusterPlot, height = 10, width = plotWidth)
    par(mar=c(15,4,4,2))
    plot(clusters)
    dev.off()
    
  }
  
  barcodeTable <- data.frame(barcodeConsensus=character(0), n=integer(0))
  
  for(i in 1:length(unique(barcodeAssignments))) {
    barcodeCluster <- unique(barcodes)[barcodeAssignments==i]
    barcodeCluster <- barcodes[barcodes %in% barcodeCluster]
    barcodeTable <- rbind(barcodeTable, defineConsensus(barcodeCluster))
  }
  
  return(barcodeTable)
  
}

#'@export
#'@title Matching Barcodes to a Barcode Library
#'@param barcodeTable A \code{data.frame} containing the barcodes to be matched
#'@param barcodeColumn Which column in \code{barcodeTable} contains the barcodes?
#'@param barcodeLibrary A vector containing the barcodes to be matched against
#'@return \code{matchBarcodes} returns the original \code{data.frame} supplied to the function, but with two additional columns containing the best matched barcode and the number of mismatches between the supplied barcode and the best match.
matchBarcodes <- function(barcodeTable, barcodeColumn, barcodeLibrary) {

  barcodeColumn <- deparse(substitute(barcodeColumn))

  barcodes <- unique(barcodeTable[[barcodeColumn]])

  matchMatrix <- as.matrix(stringdist::stringdistmatrix(barcodes, barcodeLibrary, method = "lv"))

  bestMatch <- apply(matchMatrix, 1, which.min)

  barcodes <- data.frame(barcodes, barcodeBestMatch=barcodeLibrary[bestMatch], nMismatches=apply(matchMatrix, 1, min))
  colnames(barcodes)[1] <- barcodeColumn
  colnames(barcodes)[2] <- paste0(barcodeColumn, "_BestMatch")

  matchedBarcodeTable <- merge(barcodeTable, barcodes, by = barcodeColumn)

  return(matchedBarcodeTable)

}

