#'@export
#'@title Conversion of common/alias names to systematic IDs
#'@param common A character vector of the names to be checked
#'@return A character vector the same length as the input is returned with updated IDs.
#'In the case that a \code{NULL} or empty string is provided, \code{"NA"} will  be returned.
#'In the case that no match is found in current names, the name will be checked to see if it has been removed, in which case \code{"Gene Removed"} will be returned.
#'If this also fails, \code{"Gene Not Found"} will be returned.
common2systematic <- Vectorize(function(common) {
  #Import data
  genes <- pombeGenes
  #If nothing has been entered, return NA
  if(is.null(common)) {return("NA")}
  if(common=="") {return("NA")}
  #Check if the user has entered a systematic name
  if(common %in% genes$Systematic.Name) {
    return(common)
  }
  if(toupper(common) %in% genes$Systematic.Name) {
    return(toupper(common))
  }
  #If the last letter is a character, set to lower case and check if systematic (i.e. for cases like SPCC1919.03c)
  if(paste(c(toupper(substring(common, 0, nchar(common)-1)), tolower(substring(common, nchar(common)))), collapse="") %in% genes$Systematic.Name) {
    name <- paste(c(toupper(substring(common, 0, nchar(common)-1)), tolower(substring(common, nchar(common)))), collapse="")
    return(name)
  }
  #Search if provided name is in common names
  systematic <- genes$Systematic.Name[match(common, genes$Common.Name)]
  #If no match, search in aliases
  if(is.na(systematic)) {
    systematic <- genes$Systematic.Name[match(common, genes$Aliases)]
  }
  #If no match, split the aliases by vertical bar comma (depending on species), and search the aliases again
  if(is.na(systematic)) {
    aliases <- genes$Aliases
    aliases <- strsplit(aliases, ",")
    systematic <- genes$Systematic.Name[sapply(1:length(aliases), function(i) {common %in% aliases[[i]]})]
  }
  #Return the retrieved systematic name
  if(length(systematic)!=0) {return(systematic)}
  #Convert to lower case and try again
  common <- tolower(common)
  #Search if provided name is in common names
  systematic <- genes$Systematic.Name[match(common, genes$Common.Name)]
  #If no match, search in aliases
  if(is.na(systematic)) {
    systematic <- genes$Systematic.Name[match(common, genes$Aliases)]
  }
  #If no match, split the aliases by vertical bar comma (depending on species), and search the aliases again
  if(is.na(systematic)) {
    aliases <- genes$Aliases
    aliases <- strsplit(aliases, ",")
    systematic <- genes$Systematic.Name[sapply(1:length(aliases), function(i) {common %in% aliases[[i]]})]
  }
  #Return the retrieved systematic name
  if(length(systematic)==1) {return(systematic)}
  #Check if in removed genes
  return(ifelse(common %in% pombeRemovedGenes$Gene, "Gene Removed", "Gene Not Found"))
}, USE.NAMES = F)