#'@export
deconstructConsensus <- function(tagsCsv, tag) {
  
  tags <- read.csv(tagsCsv, stringsAsFactors = F)
  
  tags <- tags[tags$tag_BestMatch == tag,]
  
  x <- plyr::ddply(tags, plyr::.(tag, gene), plyr::summarise, n=length(tag))
  
  x <- x[rev(order(x$n)),]
  
  x$tagID <- as.numeric(as.factor(x$tag))
  
  x
  
}
