#'@export
#'@title Calculation of Barcode Frequencies
#'@description \code{calculateFrequency} will calculate the total number of times each barcode-gene combination pair appears.
#'@param barcodeTable A \code{data.frame} containing the barcode and gene data
#'@param barcodeColumn The column containing the barcodes
#'@param geneColumn The column containing the genes
calculateFrequency <- function(barcodeTable, barcodeColumn, geneColumn) {
  
  barcodeColumn <- deparse(substitute(barcodeColumn))
  geneColumn <- deparse(substitute(geneColumn))
  
  plyr::count(barcodeTable, c(barcodeColumn, geneColumn))
  
}

#'@export
#'@rdname calculateFrequency
#'@description \code{aggregateProportions} will calculate the best match of genes to barcodes and provide information of the number of matches, mismatches, total matches and proportion of matches.
#'@param barcodeFrequencyTable A \code{data.frame} as produced by \code{calculateFrequency}.
aggregateProportions <- function(barcodeFrequencyTable, barcodeColumn, geneColumn) {
  
  barcodeColumn <- deparse(substitute(barcodeColumn))
  geneColumn <- deparse(substitute(geneColumn))
  
  output <- data.frame(barcode=character(0), gene=character(0), freq=integer(0), totalBarcode=integer(0), totalGene=integer(0))
  
  for(i in 1:nrow(barcodeFrequencyTable)) {
    barcode <- barcodeFrequencyTable[[barcodeColumn]][i]
    gene <- barcodeFrequencyTable[[geneColumn]][i]
    totBar <- sum(barcodeFrequencyTable$freq[which(barcodeFrequencyTable[[barcodeColumn]]==barcode)])
    totGene <- sum(barcodeFrequencyTable$freq[which(barcodeFrequencyTable[[geneColumn]]==gene)])
    output <- rbind(output, data.frame(barcode=barcode, gene=gene, freq=barcodeFrequencyTable$freq[i], totalBarcode=totBar, totalGene=totGene))
  }
  
  output$barcodeProportion <- output$freq/output$totalBarcode
  output$geneProportion <- output$freq/output$totalGene
  
  return(output)
  
}
