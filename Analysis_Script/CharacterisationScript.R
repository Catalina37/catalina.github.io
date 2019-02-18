---
Title: "Deletion Library Characterisation: Script "
---
# Load the custom built R package called BarSeqTools.
# All of the custom-built scripts are embeded within BarSeqTools.
library(BarSeqTools)

####### Sample Preparation #######
#1. Set working directory, one for uptag and one for dntag files
setwd("~/Desktop/Cat_BionnerV5_Recharc_May18-77598553/FASTQ_Generation_2018-05-22_07_28_46Z-96698246/uptag/")
setwd("~/Desktop/Cat_BionnerV5_Recharc_May18-77598553/FASTQ_Generation_2018-05-22_07_28_46Z-96698246/dntag/")

#2. Combine all uptag and dntag R1 and R2 fastq files
# uptag 
system("cat */*R1* > upTagsR1.fastq")
system("cat */*R2* > upTagsR2.fastq")
# dntag
system("cat */*R1* > dnTagsR1.fastq")
system("cat */*R2* > dnTagsR2.fastq")

#3. Load MAGIC: cutsom-built script linking barcodes & genomic sequences together
MAGIC("upTagsR1.fastq", "upTagsR2.fastq", "up")
MAGIC("dnTagsR1.fastq", "dnTagsR2.fastq", "dn")

####### Sample analysis is performed with Barcount #######
#1. Set working directory for Barcount analysis
setwd("~/Desktop/Cat_BionnerV5_Recharc_May18-77598553/FASTQ_Generation_2018-05-22_07_28_46Z-96698246/up_dn_UniqueMatches/")

#2. Extract barcodes from bed files & save as csv
ExtractBarcodefromBEDTOOL("upBedtools.Unique.bed", "upTags.Unq.csv")
ExtractBarcodefromBEDTOOL("dnBedtools.Unique.bed", "dnTags.Unq.csv")

#3. Name & load the cvs files 
upTags <- read.csv("upTags.Unq.csv")
dnTags <- read.csv("dnTags.Unq.csv")

#4. Load package & count barcode presence and sequence length
library(plyr)

upTagsCount <- ddply(upTags, .(tag), summarise, n=length(tag))
dnTagsCount <- ddply(dnTags, .(tag), summarise, n=length(tag))
# Files are too large for R to process, discard reads.

#5. Discard barcodes appearing <=2 to a unique gene
upTagsUncommon <- upTagsCount$tag[upTagsCount$n<=2]
dnTagsUncommon <- dnTagsCount$tag[dnTagsCount$n<=2]

# 6. Filter for the barcode reads mapped to a unique gene >2 
upTags <- upTags[-which(upTags$tag %in% upTagsUncommon),]
dnTags <- dnTags[-which(dnTags$tag %in% dnTagsUncommon),]

#7. Create barcode consensus by allowing 2 bp mismatches
upConsensus <- createConsensus(upTags$tag, 2, "upTags.pdf")
dnConsensus <- createConsensus(dnTags$tag, 2, "dnTags.pdf")

#8. Plot the distribution of the created barcode consensuses
hist(upConsensus$n, breaks=1500, xlim=c(0,10000), ylim=c(0, 100), main="Uptags", xlab="Frequency of Consensus Barcode")
hist(dnConsensus$n, breaks=1500, xlim=c(0,10000), ylim=c(0, 100), main="Dntags", xlab="Frequency of Consensus Barcode")

#9. Match barcodes to genes based on the barcode consenuses and actual barcodes
# Max difference allowed is 3 bp mismatches
upTags <- matchBarcodes(upTags, tag, upConsensus$barcodeConsensus)
dnTags <- matchBarcodes(dnTags, tag, dnConsensus$barcodeConsensus)

#10. Write out these tables as csv
write.csv(upTags, "upTags.csv", row.names = F)
write.csv(dnTags, "dnTags.csv", row.names = F)

#11. Calculate the number of mismatches between consensus & actual barcodes
table(upTags$nMismatches)
table(dnTags$nMismatches)

#12. Calculate the frequency for each gene-barcode consensus 
upTagsFreq <- calculateFrequency(upTags, tag_BestMatch, gene)
dnTagsFreq <- calculateFrequency(dnTags, tag_BestMatch, gene)

#13. Calculate aggregation for each gene-barcode consensus as follows:
# - each unique gene against every best barcode consensus 
# - each unique best barcode consensus against every gene
upTagsAgg <- aggregateProportions(upTagsFreq, tag_BestMatch, gene)
dnTagsAgg <- aggregateProportions(dnTagsFreq, tag_BestMatch, gene)

#14. Write out these tables as csv
write.csv(upTagsAgg, "upTagsAgg.csv", row.names = F)
write.csv(dnTagsAgg, "dnTagsAgg.csv", row.names = F)

#15. Create safe gene-barcodes with the highest frequency occurance based on:
# - min 10 occurences for each gene linked to that particular barcode
# - min 80% of the time for which that gene is linked to that barcode 
minFreq <- 10
minProportion <- 0.8

#16. Write a function to visualise these two cut offs 
  plotSafes <- function(tagsAgg, y, minFreq, minProportion, ylab, main) {
  plot(tagsAgg$freq, tagsAgg[[y]], cex=0.1, col=ifelse(tagsAgg$freq>minFreq & tagsAgg[[y]]>minProportion, "Blue", "Black"), xlab="Barcode-Gene Frequency", ylab=ylab, main=main)
  abline(h=minProportion, col="Red")
  abline(v=minFreq, col="Red")
}

### Plot to visualise the uptag safe list
par(mfrow=c(1,2)) 
plotSafes(upTagsAgg, "barcodeProportion", minFreq, minProportion, "Proportion of Barcode in Barcode-Gene Pair", main="Up Tags")
plotSafes(upTagsAgg, "geneProportion", minFreq, minProportion, "Proportion of Gene in Barcode-Gene Pair", main="Up Tags")

### Plot to visualise the dntag safe list
par(mfrow=c(1,2))
plotSafes(dnTagsAgg, "barcodeProportion", minFreq, minProportion, "Proportion of Barcode in Barcode-Gene Pair", main="Dn Tags")
plotSafes(dnTagsAgg, "geneProportion", minFreq, minProportion, "Proportion of Gene in Barcode-Gene Pair", main="Dn Tags")

#17. Subset these automatically detected gene-barcode pairs 
upTagsSafe <- subset(upTagsAgg, upTagsAgg$freq>=minFreq & upTagsAgg$barcodeProportion>=minProportion & upTagsAgg$geneProportion>=minProportion)
dnTagsSafe <- subset(dnTagsAgg, dnTagsAgg$freq>=minFreq & dnTagsAgg$barcodeProportion>=minProportion & dnTagsAgg$geneProportion>=minProportion)

#18. Save the sutomatically detected genes by writing these tables out as csv
upTagsSafe$Method <- "auto"
dnTagsSafe$Method <- "auto"

write.csv(dnTagsSafe, "dnTagsAutoSafe.csv", row.names = F)
write.csv(upTagsSafe, "upTagsAutoSafe.csv", row.names = F)

####### Analysis Quality Check #######
#1. Identify the number of safe genes, unique genes and unique barcodes
# uptag
nrow(upTagsSafe)
length(unique(upTagsSafe$barcode))
length(unique(upTagsSafe$gene))

#dntag
nrow(dnTagsSafe)
length(unique(dnTagsSafe$barcode))
length(unique(dnTagsSafe$gene))

####### Gene Browser #######
#1. Load the prerequisites for running the Gene Browser
myTxDb <- makeTxDbFromBiomart(biomart = "fungi_mart", dataset = "spombe_eg_gene", host="fungi.ensembl.org")
options(ucscChromosomeNames=FALSE)

#2. Set directory
setwd("~/Desktop/Cat_BionnerV5_Recharc_May18-77598553/FASTQ_Generation_2018-05-22_07_28_46Z-96698246/up_dn_UniqueMatches/")

#3. Load the dependency files to visualise the genes of interest
BahlerBarcodeBrowser("upTags.csv", "dnTags.csv", "upTagsAgg.csv", "dnTagsAgg.csv", "upTagsAutoSafe.csv", "dnTagsAutoSafe.csv", "upTagsManualSafe.csv", "dnTagsManualSafe.csv", myTxDb)

