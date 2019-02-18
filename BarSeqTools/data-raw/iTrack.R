library(Gviz)

options(ucscChromosomeNames=FALSE)

iTrack <- list()
for(i in 1:3) {
  
  ideogramData <- read.table(file = paste0("BarSeqTools/data-raw/cytoband_chr", i, ".txt"),
                             header = TRUE,
                             sep = " ",
                             stringsAsFactors = FALSE)
  
  itrack <- IdeogramTrack(chromosome=as.roman(i),
                               genome = "Schizosaccharomyces pombe",
                               name=NULL,
                               bands=ideogramData)
  
  iTrack[[as.character(as.roman(i))]] <- itrack
  
}

library(devtools)

use_data(dnTags, upTags, iTrack, internal = T, overwrite = T, pkg = "BarSeqTools")
