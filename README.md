# catalina.github.io
 ---------------------------------------------------------------------------------------------------------------------------
#Barcount: written by Stephan Kamrad (stephan.kamrad.15@ucl.ac.uk).
Python package for finding barcodes within sequencing reads and matching to reference database.

#BarSeqTools: written by StJohn Townsend (stjohn.townsend.11@ucl.ac.uk).
R package containing utility functions for barcode matching and visualisation.

#Bioneer V5.0 characterisation: folder where all characterisation related files are found. 
 #Browser: folder containing all the required files to visualize the genes in the Browser.
 
 ---------------------------------------------------------------------------------------------------------------------------

# Instructions for installing & running BarSeqTools & Barcount
BarSeqTools is a custom-built R package containing utility function to characterise the 
Bioneer version 5.0 of the fission yeast deletion collection library by barcode matching 
and visualisation.

### Install BarSeqTools or Barcount:
------------------------------------------------------------------------------------------
1. Open R studio
2. Install devtools package
   -install.packages("devtools")
3. Load devtools package
   -library(devtools)
4. Install github package
   -install_github("Catalina37/catalina.github.io", subdir = "BarSeqTools/Barcount")
------------------------------------------------------------------------------------------


### Check the following R package dependencies are available:
------------------------------------------------------------------------------------------
-shiny
-Gviz
-S4Vectors
-tats4
-BiocGenerics
-parallel
-IRanges
-GenomicRanges
-GenomeInfoDb
-grid
-GenomicFeatures
-AnnotationDbi
-Biobase

------------------------------------------------------------------------------------------
BarSeqTools is ready to use!

------------------------------------------------------------------------------------------


### To perform the deletion collection characterisation, check the following dependencies 
### are available & installed on your machine:
------------------------------------------------------------------------------------------
1. Anaconda
For details on how to install it, check out this link: https://docs.anaconda.com/anaconda
/install/

2. Barcount
A custom-built Python script used to search for the barcodes within the read using flanks 
on either side of the barcode.


3. Bedtools
Install with conda install -c bioconda bedtools

4. Bowtie2
Install with conda: conda install -c bioconda bowtie2

5. fastqCombinedPairedEnd.py 
Python script used to filter read 2 based on the header information of read 1. 
Read 1 header contains the barcode, thus matching read 1 barcode containing reads to 
their respective read 2 reads.
The script can be downloaded from Github code page.

6. Fastq-grep
The package is part of fastx-toolkit. 
Install with conda: conda install -c bioconda fastq-tools

7. Fastx_trimmer
The package is part of fastq-tools. 
Install with conda: conda install -c biobuilds fastx-toolkit

8. Samtools
Install with conda: conda install -c bioconda samtools 

9. PEAR 
For details on how to install it, check out this link: 
https://cme.h-its.org/exelixis/web/software/pear/doc.html)

10. Python
For details on how to install it, check out this link: https://www.python.org/downloads/
------------------------------------------------------------------------------------------


### In addition to the above, the following file dependencies are also required:
------------------------------------------------------------------------------------------

1. All protein coding genes in gtf file format.

2. A bedfile format of the Bioneer library containing all the genomic coordinates and 
strands for all of the protein coding genes.

3. An empty csv file required to run barcount since barcount was made to run with the 
barcode_table argument. 
Note, the file is not being used thus, the overall result is not affected by this.

4. The indexed S.pobme genome folder.
------------------------------------------------------------------------------------------
