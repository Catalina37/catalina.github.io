library(devtools)
library(xlsx)

BioneerV2 <- read.xlsx("BarSeqTools/data-raw/Bioneer collection Version 2.xls", 1, rowIndex = 1:3007, colIndex = 1:4, stringsAsFactors = F)

BioneerV5 <- read.xlsx("BarSeqTools/data-raw/Bioneer_collection_version_5.xlsx", 1, rowIndex = 1:3421, colIndex = 1:6, stringsAsFactors = F)

pombeGenes <- suppressWarnings(tryCatch({
  read.table("ftp://ftp.pombase.org/pombe/names_and_identifiers/sysID2product.tsv", sep="\t", quote="", stringsAsFactors=F)
}, error=function(e) {return(NULL)}))
if(is.null(pombeGenes)==F) {
  colnames(pombeGenes) <- c("Systematic.Name", "Common.Name", "Aliases", "Product")
}

pombeRemovedGenes <- read.csv("BarSeqTools/data-raw/pombeRemovedGenes.csv", stringsAsFactors=F)[,-1]

use_data(BioneerV2, BioneerV5, pombeGenes, pombeRemovedGenes, pkg = "BarSeqTools", overwrite = T)
