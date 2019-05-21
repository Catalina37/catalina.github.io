#'@export
MAGIC <- function(R1inputFastq, R2inputFastq, tag, max_distance_flanks=1, max_distance_barcode=3, barcode_table="/usr/local/bin/Empty_ref.csv") {
  
  system <- function(...) {base::system('bash -l', input=c("shopt -s expand_aliases", ...))}
  
  dir.create("Analysis", showWarnings = F)
  
  if(tag=="up") {
    flanking_left <- "CAAGCTAAGATATC"
    flanking_right <- "TTTAAATGCGAAGTAA"
    primer <- "CGCTCCCGCCTTACTTCGCATTTAAA"
    bedtoolsArgs <- "-iu"
  }
  if(tag=="dn") {
    flanking_left <- "AGTGTCGAAAAGTATC"
    flanking_right <- "TTTAAAATCCCCCCTA"
    primer <- "TTGCGTTGCGTAGGGGGGATTTTAAA"
    bedtoolsArgs <- "-id"
  }
  
  system(paste0("
    barcount --fastq ", R1inputFastq, " --flanking_left ", flanking_left, " --flanking_right ", flanking_right, " --max_distance_flanks ", max_distance_flanks, " --max_distance_barcode ", max_distance_barcode, " --barcode_table ", barcode_table, " --debug --verbose --save_extracted_barcodes --out Analysis/R1.Filter.fastq
    
    python /usr/local/bin/fastqCombinedPairedEnd.py Analysis/R1.Filter.fastq_barcodes.fastq ", R2inputFastq, " 
    
    mv ", R2inputFastq, "_pairs_R2.fastq Analysis/", R2inputFastq, "_pairs_R2.fastq

    sed -i 's/\\(@.*2:N:0:\\).*/\\1/' Analysis/", R2inputFastq, "_pairs_R2.fastq 
    
    paste -d '~' Analysis/", R2inputFastq, "_pairs_R2.fastq Analysis/R1.Filter.fastq_barcodes.fastq | perl -F'~' -lane 'push(@buffer, $F[0]); if($line == 1){@buffer[0] .= \"[\".$F[1].\"]\"}; if(($line == 3) && @buffer){print join(\"\\n\",@buffer); @buffer = ()}; $line = ($line+1) % 4;' - > Analysis/R2_barcode.fastq
    
    fastx_trimmer -f 33 -l 75 -i Analysis/R2_barcode.fastq -Q33 -o Analysis/R2_Filtered_Trimmed.fastq
    
    fastq-grep -v '", primer, "' Analysis/R2_Filtered_Trimmed.fastq > Analysis/R2_1filtered.fastq
    
    bowtie2 --sam-no-qname-trunc -N 1 -x /usr/local/bin/pombe/pombe \\ -q -U Analysis/R2_1filtered.fastq -S Analysis/Aligned.bam
    
    grep -v 'XS:i:' Aligned.bam  > Aligned.bam

    samtools view -b Analysis/Aligned.bam | bedtools bamtobed > Analysis/Aligned.bed
    
    bedtools  annotate -names -counts -both -i Analysis/Aligned.bed -files /usr/local/bin/genes.gtf > Analysis/Annotated.bed
    
    sort -k1,1 -k2,2n Analysis/Annotated.bed > Analysis/Sorted.bed
    
    cat Analysis/Sorted.bed | sed -e 's/chromosome_1/I/g' | sed -e 's/chromosome_2/II/g' | sed -e 's/chromosome_3/III/g' > Analysis/Chr.bed
    
    awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$5,$6,$7}' Analysis/Chr.bed| sort-bed - > Analysis/Column.bed
  "))
  
  if(tag=="dn") {
    bedFile <- read.delim("Analysis/Column.bed", header=F, stringsAsFactors=F)
    strand <- bedFile[,6]
    strandPlus <- strand=="+"
    strandNeg <- strand=="-"
    strandNew <- rep(".", length(strand))
    strandNew[strandPlus] <- "-"
    strandNew[strandNeg] <- "+"
    bedFile[,6] <- strandNew
    readr::write_delim(bedFile, "Analysis/Column.bed", col_names = F, delim="\t")
  }
  
  system(paste0("
    bedtools closest ", bedtoolsArgs, " -a Analysis/Column.bed -b /usr/local/bin/BioneerV5.bed -s -D a > bedtools.bed
  "))
  
}
