#'@export
#'@title Importing data from bedtools
#'@description \code{ExtractBarcodefromBEDTOOL} will read a file outputted from bedtools, extract the barcode information and remove unnecessary columns.
#'@param file The file path to the bedtools output.
#'@param outputFile The file path to the csv to be created.
#'@return  Nothing is returned from this function. Instead, a new csv file is written.
ExtractBarcodefromBEDTOOL<-function(file, outputfile){
  
  Barcodes<-read.csv(file, sep = "\t", header=F)
  Barcodes<-Barcodes[,c(1:4, 6, 10)]
  colnames(Barcodes)<-c("chr", "start", "end", "tag", "strand", "gene")
  
  Barcodes$tag <- stringr::str_match(Barcodes$tag, "2:N:0:\\[((?:A|T|C|G)*)\\]")[,2]
  
  write.csv(Barcodes, outputfile, row.names = F)
  
  return()
  
}